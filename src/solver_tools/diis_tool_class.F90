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
      type(file), dimension(:), allocatable, private :: dx   ! File containing the O vector from previous iterations
      type(file), dimension(:), allocatable, private :: x_dx ! File containing the X vector from previous iterations
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
      logical :: accumulate = .true. ! Variable keeping track of wether we are cummulatively building the diis history or not
!                                    ! Note: default is to accumulate
!
      logical :: erase_history = .false.  ! (DIIS space) 1,2,3,4,...,8,1,2,3,... if true 
                                          ! (DIIS space) 1,2,3,4,...,8,8,8,8,... if false 
!
   contains
!
      procedure :: init     => init_diis_tool     ! Called to set up a DIIS solver
      procedure :: update   => update_diis_tool   ! Called to perform a DIIS step
      procedure :: finalize => finalize_diis_tool
!
      procedure :: get_current_dim => get_current_dim_diis_tool
!
      procedure :: get_x_dx   => get_x_dx_diis_tool
      procedure :: get_dx     => get_dx_diis_tool
!
      procedure :: set_x_dx   => set_x_dx_diis_tool
      procedure :: set_dx     => set_dx_diis_tool
!
      procedure, private :: construct_diis_matrix => construct_diis_matrix_diis_tool
      procedure, private :: write_current_vecs    => write_current_vecs_diis_tool
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
      integer :: i
!
      do i = 1, solver%diis_dimension
!
         call disk%open_file(solver%x_dx(i), 'readwrite')
         call disk%close_file(solver%x_dx(i), 'delete')
!
         call disk%open_file(solver%dx(i), 'readwrite')
         call disk%close_file(solver%dx(i), 'delete')
!
      enddo
!
      call disk%open_file(solver%diis_matrix, 'readwrite')
      call disk%close_file(solver%diis_matrix, 'delete')
!
      solver%iteration = 1
      solver%diis_dimension = 8

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
      integer :: I
!
      character(len=200) :: name_dx, name_x_dx
!
      solver%name = trim(name)
      solver%n_parameters = n_parameters
      solver%n_equations = n_equations
      solver%iteration = 1 
!
      solver%accumulate = .true. 
      if (present(accumulate)) then 
!
         solver%accumulate = accumulate
!
      endif
!
      if (solver%accumulate) then 
!
         solver%erase_history = .false.
!
      else
!
         solver%erase_history = .true.
!
      endif 
!
      if (present(diis_dimension)) solver%diis_dimension = diis_dimension
!
!     File handling
!
      allocate(solver%x_dx(solver%diis_dimension))
      allocate(solver%dx(solver%diis_dimension))
!
      do I = 1, solver%diis_dimension
!
         write(name_x_dx, '(a, a1, i3.3)') trim(solver%name) // '_x_dx', '_', I
         write(name_dx, '(a, a1, i3.3)') trim(solver%name) // '_dx', '_', I
!
         call solver%x_dx(i)%init(trim(name_x_dx), 'sequential', 'unformatted')
         call solver%dx(i)%init(trim(name_dx) // '_dx', 'sequential', 'unformatted')
!
      enddo
!
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
      real(dp), dimension(solver%n_equations)  :: dx
      real(dp), dimension(solver%n_parameters) :: x_dx
!
      real(dp), dimension(:), allocatable :: x_dx_i ! To hold previous x_dx_i temporarily
!
      integer :: i = 0
!
      integer :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer :: current_dim 
!
      real(dp), dimension(:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
!     :: Open DIIS files
!
!     Ask disk manager to open the files
!
     do i = 1, solver%diis_dimension
!
        call disk%open_file(solver%dx(i), 'readwrite')
        call disk%open_file(solver%x_dx(i), 'readwrite')
!
     enddo
!
!     Compute the current dimensionality of the problem 
!     (1, 2,..., 7, 8, 8, 8, 8,...) for the standard diis_dimension = 8
!
      current_dim = solver%get_current_dim()
!
!     Write current x_dx and dx to file
!
      call solver%write_current_vecs(dx, x_dx, current_dim)
!
!     :: Solve the least squares problem, G w = H
!
!        G : DIIS matrix, G_ij = Δx_i Δx_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., current_dim. To enforce normality
!     of the solution, G is extended with a row & column of -1's
!     and H with a -1 at the end.
!
!     Set the DIIS vector H
!
      call mem%alloc(diis_vector, current_dim + 1)
!
      diis_vector = zero
      diis_vector(current_dim + 1) = -one
!
!     Set the DIIS matrix G
!
      call mem%alloc(diis_matrix, current_dim + 1, current_dim + 1)
!
      call solver%construct_diis_matrix(current_dim, diis_matrix, dx)
!
!     Solve the DIIS equation G w = H for the DIIS weights w 
!
!     Note: on exit, the solution is in the diis_vector,
!     provided info = 0 (see LAPACK documentation for more)
!
      call mem%alloc(ipiv, current_dim + 1)
      ipiv = 0
!
      call dgesv(current_dim + 1,  &
                  1,               &
                  diis_matrix,     &
                  current_dim + 1, &
                  ipiv,            &
                  diis_vector,     &
                  current_dim + 1, &
                  info)
!
      call mem%dealloc(ipiv, current_dim + 1)
      call mem%dealloc(diis_matrix, current_dim + 1, current_dim + 1)
!
      if (info .ne. 0) call output%error_msg('could not solver DIIS matrix equation!')
!
!     Update the parameters (place in x_dx on exit)
!
      x_dx = zero 
!
      call mem%alloc(x_dx_i, solver%n_parameters)
!
      do i = 1, current_dim
!
!        Read the x_i + Δ x_i vector
!
         rewind(solver%x_dx(i)%unit)
         read(solver%x_dx(i)%unit) x_dx_i
!
!        Add w_i (x_i + Δ x_i) to the amplitudes
!
         call daxpy(solver%n_parameters, diis_vector(i), x_dx_i, 1, x_dx, 1)
!
      enddo
!
      call mem%dealloc(x_dx_i, solver%n_parameters)
!
      call mem%dealloc(diis_vector, current_dim + 1)
!
!     Close files
!
      do i = 1, solver%diis_dimension

         call disk%close_file(solver%dx(i))
         call disk%close_file(solver%x_dx(i))

      enddo
!
      solver%iteration = solver%iteration + 1
!
   end subroutine update_diis_tool
!
!
   function get_current_dim_diis_tool(solver)
!!
!!    Get current dimension 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      integer :: get_current_dim_diis_tool
!
      if (solver%iteration .gt. solver%diis_dimension .and. &
            .not. solver%erase_history) then 
!
         get_current_dim_diis_tool = solver%diis_dimension
!
      else
!
         get_current_dim_diis_tool = solver%iteration - &
                           solver%diis_dimension*((solver%iteration-1)/solver%diis_dimension)
!
      endif
!
   end function get_current_dim_diis_tool
!
!
   subroutine write_current_vecs_diis_tool(solver, dx, x_dx, current_dim)
!!
!!    Write current vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Writes the current dx and x_dx to file. Each index 1, 2, ..., diis_dim
!!    has its own file. If more iterations have passed than we keep on file
!!    (i.e., solver%iteration > solver%diis_dimension), then we do a first-in/
!!    last-out replacement unless erase_history is true. 
!!
      implicit none 
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(in)  :: dx
      real(dp), dimension(solver%n_parameters), intent(in) :: x_dx
!
      integer, intent(in) :: current_dim
!
      integer :: k
!
      type(file) :: tmp_file
!
      if (current_dim .le. solver%diis_dimension &
            .and. current_dim .eq. solver%iteration &
            .or. solver%erase_history) then 
!
!        Just write to the correct file 
!
         rewind(solver%dx(current_dim)%unit)
         write(solver%dx(current_dim)%unit) dx
!
         rewind(solver%x_dx(current_dim)%unit)
         write(solver%x_dx(current_dim)%unit) x_dx
!
      else 
!
!        Redefine files such that the current entry is the last, 
!        while all the others move one position up 
!
         tmp_file = solver%dx(1)
!
         do k = 1, current_dim - 1 
!
            solver%dx(k) = solver%dx(k+1)
!
         enddo 
!
         solver%dx(current_dim) = tmp_file
!
         rewind(solver%dx(current_dim)%unit)
         write(solver%dx(current_dim)%unit) dx
!
         tmp_file = solver%x_dx(1)
!
         do k = 1, current_dim - 1 
!
            solver%x_dx(k) = solver%x_dx(k+1)
!
         enddo 
!
         solver%x_dx(current_dim) = tmp_file
!
         rewind(solver%x_dx(current_dim)%unit)
         write(solver%x_dx(current_dim)%unit) x_dx
!
      endif
!
   end subroutine write_current_vecs_diis_tool
!
!
   subroutine construct_diis_matrix_diis_tool(solver, current_dim, diis_matrix, dx)
!!
!!    Construct DIIS matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      integer, intent(in) :: current_dim
!
      real(dp), dimension(current_dim + 1, current_dim + 1), intent(inout) :: diis_matrix
!
      real(dp), dimension(solver%n_equations), intent(in) :: dx
!
      real(dp), dimension(:), allocatable :: dx_i  ! To hold previous Δ x_i temporarily
      real(dp), dimension(:), allocatable :: dx_j  ! To hold previous Δ x_j temporarily
!
      integer :: i, j
!
      real(dp), dimension(:,:), allocatable :: diis_matrix_tmp 
!
      real(dp) :: ddot
!
      diis_matrix = zero
!
      if (solver%accumulate) then
!
         call disk%open_file(solver%diis_matrix, 'readwrite')
!
         if (solver%iteration .gt. solver%diis_dimension) then
!
!           First entry goes out, eliminating the (old) 1st row & column.
!
            call mem%alloc(diis_matrix_tmp, current_dim, current_dim)
!
            rewind(solver%diis_matrix%unit)
!
            do j = 1, current_dim
               do i = 1, current_dim
!
                  read(solver%diis_matrix%unit) diis_matrix_tmp(i,j)
!
               enddo
            enddo
!
            diis_matrix(1 : (current_dim - 1), 1 : (current_dim - 1)) = diis_matrix_tmp(2 : current_dim, 2 : current_dim)      
!
            call mem%dealloc(diis_matrix_tmp, current_dim, current_dim)
!
         elseif (current_dim .gt. 1) then
!
!           Still doing the first set of vectors: just read the previous matrix 
!
            rewind(solver%diis_matrix%unit)
!
            do j = 1, current_dim - 1
               do i = 1, current_dim - 1
!
                  read(solver%diis_matrix%unit) diis_matrix(i,j)
!
               enddo
            enddo
!
         endif
!
!        Get the parts of the DIIS matrix G not constructed in
!        the previous iterations
!
         call mem%alloc(dx_i, solver%n_equations)
!
         do i = 1, current_dim
!
            rewind(solver%dx(i)%unit)
            read(solver%dx(i)%unit) dx_i
!
            diis_matrix(current_dim, i) = ddot(solver%n_equations, dx, 1, dx_i, 1)
            diis_matrix(i, current_dim) = diis_matrix(current_dim, i)
!
            diis_matrix(current_dim + 1, i) = -one
            diis_matrix(i, current_dim + 1) = -one
!
         enddo
!
         call mem%dealloc(dx_i, solver%n_equations)
!
!        Write the current DIIS matrix to file
!
         rewind(solver%diis_matrix%unit)
!
         do j = 1, current_dim
            do i = 1, current_dim
!
               write(solver%diis_matrix%unit) diis_matrix(i,j)
!
            enddo
         enddo
!
         call disk%close_file(solver%diis_matrix)
!
      else ! Non-cumulative diis -> build diis matrix from scratch
!         
         call mem%alloc(dx_i, solver%n_equations)
         call mem%alloc(dx_j, solver%n_equations)
!
         do i = 1, current_dim
!
            rewind(solver%dx(i)%unit)
            read(solver%dx(i)%unit) dx_i
!
            do j = 1, i
!
               rewind(solver%dx(j)%unit)
               read(solver%dx(j)%unit) dx_j
!
               diis_matrix(i, j) = ddot(solver%n_equations, dx_i, 1, dx_j, 1)
               diis_matrix(j, i) = diis_matrix(i, j)
!
            enddo
!
            diis_matrix(current_dim + 1, i) = -one
            diis_matrix(i, current_dim + 1) = -one
!
         enddo 
!
         call mem%dealloc(dx_i, solver%n_equations)
         call mem%dealloc(dx_j, solver%n_equations)
!
      endif
!
   end subroutine construct_diis_matrix_diis_tool
!
!
   subroutine get_x_dx_diis_tool(solver, x_dx_i, i)
!!
!!    Get x_dx_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_parameters), intent(out) :: x_dx_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (i .gt. current_dim) call output%error_msg('Asked for an x_dx not in use.')
!
      call disk%open_file(solver%x_dx(i), 'read')
      rewind(solver%x_dx(i)%unit)
!
      read(solver%x_dx(i)%unit) x_dx_i
!
      call disk%close_file(solver%x_dx(i))
!
   end subroutine get_x_dx_diis_tool
!
!
   subroutine get_dx_diis_tool(solver, dx_i, i)
!!
!!    Get dx_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(out) :: dx_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (i .gt. current_dim) call output%error_msg('Asked for an dx not in use.')
!
      call disk%open_file(solver%dx(i), 'read')
      rewind(solver%dx(i)%unit)
!
      read(solver%dx(i)%unit) dx_i
!
      call disk%close_file(solver%dx(i))
!
   end subroutine get_dx_diis_tool
!
!
   subroutine set_x_dx_diis_tool(solver, x_dx_i, i)
!!
!!    Set x_dx_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_parameters), intent(in) :: x_dx_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (solver%accumulate) call output%error_msg('Can not set x_dx for cumulative diis.')
      if (i .gt. current_dim) call output%error_msg('Asked for an x_dx not in use.')
!
      call disk%open_file(solver%x_dx(i), 'write')
      rewind(solver%x_dx(i)%unit)
!
      write(solver%x_dx(i)%unit) x_dx_i
!
      call disk%close_file(solver%x_dx(i))
!
   end subroutine set_x_dx_diis_tool
!
!
   subroutine set_dx_diis_tool(solver, dx_i, i)
!!
!!    Set dx_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(in) :: dx_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (solver%accumulate) call output%error_msg('Can not set dx for cumulative diis.')
      if (i .gt. current_dim) call output%error_msg('Asked for an dx not in use.')
!
      call disk%open_file(solver%dx(i), 'write')
      rewind(solver%dx(i)%unit)
!
      write(solver%dx(i)%unit) dx_i
!
      call disk%close_file(solver%dx(i))
!
   end subroutine set_dx_diis_tool
!
!
end module diis_tool_class
