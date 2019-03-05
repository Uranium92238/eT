!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module diis_tool_class
!
!!
!!    DIIS solver class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!
!!    Objects of this class can be instantiated in order to solve a non-linear set of equations,
!!    O(X) = 0, where O and X are vectors of length N_O and N_X.
!!
!!    The squared norm of the averaged error vector
!!
!!       < O(X) > = sum_k w_k O(X)_k
!!
!!    is minimized with respect to the condition sum_k w_k = 1. The k index
!!    runs over previous iterations, usually limited to the last 8 iterations
!!    (determined by diis_dimension).
!!
!!    The standard update of the parameters is
!!
!!       X = sum_k w_k XdX_k
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
!
!
   type :: diis_tool
!
      character(len=40) :: name ! Solver name; determines the prefix of all DIIS files
!
      type(file), private :: dx   ! File containing the O vector from previous iterations
      type(file), private :: x_dx ! File containing the X vector from previous iterations
!
      type(file), private :: diis_matrix ! File containing the previously computed elements of the DIIS matrix
!
      integer, private :: iteration = 1 ! Variable keeping track of the current DIIS iteration
                                        ! Note: defined to increment by +1 each time 'update' is called.
!
      integer, private :: diis_dimension = 8   ! Standard is 8, though it might be useful to change this value
!
      integer, private :: n_parameters ! The length of the X vector
      integer, private :: n_equations  ! The length of the O vector
!
      logical :: accumulative = .true. ! Variable keeping track of wether we are cummulatively building the diis history or not
!                                      ! Note: default is to accumulate
!
   contains
!
      procedure :: init     => init_diis_tool     ! Called to set up a DIIS solver
      procedure :: update   => update_diis_tool   ! Called to perform a DIIS step
      procedure :: finalize => finalize_diis_tool
!
      procedure :: get_current_index => get_current_index_diis_tool
!
      procedure :: get_dx     => get_dx_diis_tool
      procedure :: get_x_dx   => get_x_dx_diis_tool
!
   end type diis_tool
!
!
contains
!
!
   subroutine finalize_diis_tool(solver)
!!
!!    Finalize DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
      implicit none
!
      class(diis_tool) :: solver
!
      solver%iteration = 1
      solver%diis_dimension = 8
!
      call disk%open_file(solver%x_dx, 'readwrite')
      call disk%open_file(solver%dx, 'readwrite')
      call disk%open_file(solver%diis_matrix, 'readwrite')
!
      call disk%close_file(solver%x_dx, 'delete')
      call disk%close_file(solver%dx, 'delete')
      call disk%close_file(solver%diis_matrix, 'delete')
!
   end subroutine finalize_diis_tool
!
!
   subroutine init_diis_tool(solver, name, n_parameters, n_equations, diis_dimension, accumulate)
!!
!!    Init DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
      class(diis_tool) :: solver
!
      character(len=*), intent(in) :: name
!
      integer, intent(in) :: n_parameters
      integer, intent(in) :: n_equations
!
      integer, intent(in), optional :: diis_dimension
!
      logical, intent(in), optional :: accumulate
!
      solver%name = trim(name)
      solver%n_parameters = n_parameters
      solver%n_equations = n_equations
!
      solver%accumulative = .true. 
      if (present(accumulate)) solver%accumulative = accumulate
!
      if (present(diis_dimension)) solver%diis_dimension = diis_dimension
!
      call solver%x_dx%init(trim(solver%name) // '_x_dx', 'sequential', 'unformatted')
      call solver%dx%init(trim(solver%name) // '_dx', 'sequential', 'unformatted')
      call solver%diis_matrix%init(trim(solver%name) // '_diis_matrix', 'sequential', 'unformatted')
!
   end subroutine init_diis_tool
!
!
    subroutine update_diis_tool(solver, dx, x_dx)
!!
!!    Update (DIIS_tool)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Based on the CCSD specific routine written by Sarai D. Folkestad
!!    and Eirik F. Kjønstad, May 2017.
!!
!!    Solves the DIIS problem, as described in the documentation of DIIS class module
!!    (see the top of this file).
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations, 1)  :: dx
      real(dp), dimension(solver%n_parameters, 1) :: x_dx
!
      real(dp), dimension(:,:), allocatable :: dx_i   ! To hold previous Δ x_i temporarily
      real(dp), dimension(:,:), allocatable :: x_dx_i ! To hold previous x_dx_i temporarily
!
      real(dp) :: ddot
!
      integer :: i = 0, j = 0
!
      integer      :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer :: current_index = 0 ! Progressing as follows: 1,2,...,7,8,1,2,...
!
      integer :: dummy = 0
!
      real(dp), dimension(:,:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
!     :: Open DIIS files
!
!     Ask disk manager to open the files
!
      call disk%open_file(solver%dx, 'readwrite')
      call disk%open_file(solver%x_dx, 'readwrite')
      call disk%open_file(solver%diis_matrix, 'readwrite')
!
!     :: Compute the current index
!     (1,2,...,7,8,1,2,...) for the standard diis_dimension = 8
!
      current_index = solver%iteration - &
               ((solver%diis_dimension)-1)*((solver%iteration-1)/((solver%diis_dimension)-1))
!
!     :: Save (Δ x_i) and (x_i + Δ x_i) to files
!
      if (current_index .eq. 1) then
!
         rewind(solver%dx%unit)
         rewind(solver%x_dx%unit)
!
      else
!
         do dummy = 1, current_index - 1 ! Get to the correct line in the file
!
            read(solver%dx%unit)
            read(solver%x_dx%unit)
!
         enddo
!
      endif
!
      write(solver%dx%unit)   (dx(i,1), i = 1, solver%n_equations)
      write(solver%x_dx%unit) (x_dx(i,1), i = 1, solver%n_parameters)
!
!     :: Solve the least squares problem, G * w = H
!
!        G : DIIS matrix, G_ij = Δ t_i Δ t_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., current_index. To enforce normality
!     of the solution, G is extended with a row & column of -1's
!     and H with a -1 at the end.
!
!     First set the DIIS vector to one
!
      call mem%alloc(diis_vector,current_index+1,1)
      diis_vector = zero
!
!     Allocate the DIIS matrix and read in previous matrix elements
!
      call mem%alloc(diis_matrix, current_index+1, current_index+1)
      diis_matrix = zero
!
      if (current_index .gt. 1) then
!
         rewind(solver%diis_matrix%unit)
!
         do j = 1, current_index - 1
            do i = 1, current_index - 1
!
               read(solver%diis_matrix%unit) diis_matrix(i,j)
!
            enddo
         enddo
!
      endif
!
!     Get the parts of the DIIS matrix G not constructed in
!     the previous iterations
!
      call mem%alloc(dx_i, solver%n_equations, 1) ! Allocate temporary holder of quasi-Newton estimates
      dx_i = zero
!
      rewind(solver%dx%unit)
!
      do i = 1, current_index
!
         read(solver%dx%unit) (dx_i(j,1), j = 1, solver%n_equations)
!
         diis_matrix(current_index,i) = ddot(solver%n_equations, dx, 1, dx_i, 1)
         diis_matrix(i,current_index) = diis_matrix(current_index,i)
!
         diis_matrix(current_index+1,i) = -one
         diis_matrix(i,current_index+1) = -one
!
      enddo
!
      diis_vector(current_index+1,1) = -one
!
!     Write the current DIIS matrix to file
!
      rewind(solver%diis_matrix%unit)
!
      do j = 1, current_index
         do i = 1, current_index
!
            write(solver%diis_matrix%unit) diis_matrix(i,j)
!
         enddo
      enddo
!
!     Solve the DIIS equation
!
      call mem%alloc(ipiv, solver%diis_dimension)
      ipiv = 0
!
!     Note: on exit, the solution is in the diis_vector,
!     provided info = 0 (see LAPACK documentation for more)
!
      call dgesv(current_index+1,  &
                  1,               &
                  diis_matrix,     &
                  current_index+1, &
                  ipiv,            &
                  diis_vector,     &
                  current_index+1, &
                  info)
!
      call mem%dealloc(ipiv, solver%diis_dimension)
!
!     :: Update the parameters (place in x_dx on exit)
!
      x_dx = zero 
!
      call mem%alloc(x_dx_i, solver%n_parameters, 1)
!
      rewind(solver%x_dx%unit)
!
      do i = 1, current_index
!
!        Read the x_i + Δ x_i vector
!
         x_dx_i = zero
         read(solver%x_dx%unit) (x_dx_i(j, 1), j = 1, solver%n_parameters)
!
!        Add w_i (x_i + Δ x_i) to the amplitudes
!
         call daxpy(solver%n_parameters, diis_vector(i, 1), x_dx_i, 1, x_dx, 1)
!
      enddo
!
      call mem%dealloc(x_dx_i, solver%n_parameters, 1)
!
!     Deallocations
!
      call mem%dealloc(dx_i, solver%n_equations, 1)
      call mem%dealloc(diis_vector, current_index + 1, 1)
      call mem%dealloc(diis_matrix, current_index + 1, current_index+1)
!
!     Close files
!
      call disk%close_file(solver%dx)
      call disk%close_file(solver%x_dx)
      call disk%close_file(solver%diis_matrix)
!
      solver%iteration = solver%iteration + 1
!
   end subroutine update_diis_tool
!
!
   function get_current_index_diis_tool(solver)
!
      implicit none
!
      class(diis_tool) :: solver
!
      integer :: get_current_index_diis_tool
!
      get_current_index_diis_tool = solver%iteration - &
               ((solver%diis_dimension)-1)*((solver%iteration-1)/((solver%diis_dimension)-1))
!
   end function get_current_index_diis_tool
!
!
   subroutine get_x_dx_diis_tool(solver, x_dx)
!!
!!    Get x_dx 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
!!    Assumes the position in the file is the correct possition
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_parameters), intent(out) :: x_dx
!
      read(solver%x_dx%unit) x_dx
!
   end subroutine get_x_dx_diis_tool
!
!
   subroutine get_dx_diis_tool(solver, dx)
!!
!!    Get dx 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
!!    Assumes the position in the file is the correct possition
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(out) :: dx
!
      read(solver%dx%unit) dx
!
   end subroutine get_dx_diis_tool
!
!
   subroutine prepare_to_set_x_dx_diis_tool(solver)
!!
!!    Prepare to set x_dx
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      if (solver%accumulative) call output%error_msg('can not set x_dx for acummulative diis.')
!
      call solver%x_dx%init_rewrite()
      call disk%open_rewrite_file(solver%x_dx, 'write')
!
      rewind(solver%x_dx%unit)
      rewind(solver%x_dx%rewrite_unit)
!
   end subroutine prepare_to_set_x_dx_diis_tool
!
!
   subroutine finnish_to_set_x_dx_diis_tool(solver)
!!
!!    Finnish to set x_dx
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      if (solver%accumulative) call output%error_msg('can not set x_dx for acummulative diis.')
!
!     1. Close file and rewrite file
!
      call disk%close_rewrite_file(solver%x_dx)
      call disk%close_file(solver%x_dx)
!
!     2. Overwrite file
!
      call disk%overwrite_file(solver%x_dx)
!
!     3. Finalize rewrite and open file
!
      call solver%x_dx%finalize_rewrite()
!
      call disk%open_file(solver%x_dx, 'readwrite')
!
   end subroutine finnish_to_set_x_dx_diis_tool
!
!
   subroutine set_x_dx_diis_tool(solver, x_dx)
!!
!!    Set x_dx
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
!!    Assumes the position in the file is the correct possition
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(in) :: x_dx
!
!     Sanity checks
!
      if (solver%accumulative) call output%error_msg('can not set dx for acummulative diis.')
      if (.not. solver%x_dx%rewrite_opened) call output%error_msg('rewrite file has not been initialized for x_dx.')
!
!     Write to rewrite file
!
      write(solver%x_dx%rewrite_unit) x_dx
!
   end subroutine set_x_dx_diis_tool
!
!
   subroutine prepare_to_set_dx_diis_tool(solver)
!!
!!    Prepare to set dx
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      if (solver%accumulative) call output%error_msg('can not set dx for acummulative diis.')
!
      call solver%dx%init_rewrite()
      call disk%open_rewrite_file(solver%dx, 'write')
!
      rewind(solver%dx%unit)
      rewind(solver%dx%rewrite_unit)
!
   end subroutine prepare_to_set_dx_diis_tool
!
!
   subroutine finnish_to_set_dx_diis_tool(solver)
!!
!!    Finnish to set dx
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      if (solver%accumulative) call output%error_msg('can not set dx for acummulative diis.')
!
!     1. Close file and rewrite file
!
      call disk%close_rewrite_file(solver%dx)
      call disk%close_file(solver%dx)
!
!     2. Overwrite file
!
      call disk%overwrite_file(solver%dx)
!
!     3. Finalize rewrite and open file
!
      call solver%dx%finalize_rewrite()
!
      call disk%open_file(solver%dx, 'readwrite')
!
   end subroutine finnish_to_set_dx_diis_tool
!
!
   subroutine set_dx_diis_tool(solver, dx)
!!
!!    Set dx
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
!!    Assumes the position in the file is the correct possition
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(in) :: dx
!
!     Sanity checks
!
      if (solver%accumulative) call output%error_msg('can not set dx for acummulative diis.')
      if (.not. solver%dx%rewrite_opened) call output%error_msg('rewrite file has not been initialized for dx.')
!
!     Write to rewrite file
!
      write(solver%dx%rewrite_unit) dx
!
   end subroutine set_dx_diis_tool
!
end module diis_tool_class
