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
!!    e(x) = 0, where e and x are vectors of length n_equations and n_parameters.
!!
!!    The squared norm of the averaged error vector
!!
!!       < e(x) > = sum_k w_k e(x)_k
!!
!!    is minimized with respect to the condition sum_k w_k = 1. The k index
!!    runs over previous iterations, by default limited to the last 8 iterations
!!    (determined by diis_dimension).
!!
!!    The standard update of the parameters is
!!
!!       x = sum_k w_k (x + e)_k
!!
!
   use global_out, only : output
!
   use parameters
   use array_utilities, only : zero_array
!
   use sequential_file_class, only : sequential_file
!
   use memory_manager_class, only : mem
!
!
   type :: diis_tool
!
      character(len=40) :: name_ ! Solver name; determines the prefix of all DIIS files
!
      type(sequential_file), dimension(:), allocatable, private :: e_files ! Array of files with error vectors
      type(sequential_file), dimension(:), allocatable, private :: x_files ! Array of files with X vectors
!
      type(sequential_file), private :: diis_matrix ! File containing the previously computed elements of the DIIS matrix
!
      integer, private :: iteration = 1 ! Variable keeping track of the current DIIS iteration
!                                       ! Note: defined to increment by +1 each time 'update' is called.
!
      integer, private :: diis_dimension = 8 ! Standard is 8, though it might be useful to change this value
!
      integer, private :: n_parameters ! The length of the X vector
      integer, private :: n_equations  ! The length of the O vector
!
      logical :: accumulate = .true. ! Variable keeping track of wether we are cummulatively 
!                                    ! building the diis history or not
!                                    ! Note: default is to accumulate
!
      logical :: erase_history = .false.  ! (DIIS space) 1,2,3,4,...,8,1,2,3,... if true 
                                          ! (DIIS space) 1,2,3,4,...,8,8,8,8,... if false 
!
   contains
!
      procedure :: update     => update_diis_tool   ! Called to perform a DIIS step
      procedure :: cleanup    => cleanup_diis_tool
!
      procedure :: get_current_dim => get_current_dim_diis_tool
!
      procedure :: read_x  => read_x_diis_tool
      procedure :: read_e  => read_e_diis_tool
!
      procedure :: write_x => write_x_diis_tool
      procedure :: write_e => write_e_diis_tool
!
      procedure, private :: construct_diis_matrix => construct_diis_matrix_diis_tool
      procedure, private :: write_current_vecs    => write_current_vecs_diis_tool
!
   end type diis_tool
!
!
   interface diis_tool 
!
      procedure :: new_diis_tool
!
   end interface diis_tool 
!
!
contains
!
!
   function new_diis_tool(name_, n_parameters, n_equations, diis_dimension, accumulate, erase_history) &
                          result(solver)
!!
!!    New DIIS tool
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
      type(diis_tool) :: solver
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters
      integer, intent(in) :: n_equations
!
      integer, intent(in), optional :: diis_dimension
!
      logical, intent(in), optional :: accumulate
      logical, intent(in), optional :: erase_history
!
      integer :: i
!
      character(len=200) :: name_e, name_x
!
      solver%name_ = trim(name_)
      solver%n_parameters = n_parameters
      solver%n_equations  = n_equations
      solver%iteration    = 1 
!
      solver%accumulate = .true. 
      solver%erase_history = .false. 
!
      if (present(accumulate)) then 
!
         solver%accumulate = accumulate
!
      endif
!
      if (present(erase_history)) then
!
         solver%erase_history = erase_history
!
      endif
!
      if (present(diis_dimension)) solver%diis_dimension = diis_dimension
!
!     File handling
!
      allocate(solver%e_files(solver%diis_dimension))
      allocate(solver%x_files(solver%diis_dimension))
!
      do i = 1, solver%diis_dimension
!
         write(name_e, '(a, a1, i3.3)') trim(solver%name_) // '_e', '_', i
         write(name_x, '(a, a1, i3.3)') trim(solver%name_) // '_x', '_', I
!
         solver%e_files(i) = sequential_file(trim(name_e))
         solver%x_files(i) = sequential_file(trim(name_x))
!
      enddo
!
      solver%diis_matrix = sequential_file(trim(solver%name_) // '_matrix')
!
   end function new_diis_tool
!
!
   subroutine cleanup_diis_tool(solver)
!!
!!    Cleanup
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
         call solver%x_files(i)%delete_()
         call solver%e_files(i)%delete_()
!
      enddo
!
      call solver%diis_matrix%delete_()
!
      solver%iteration = 1
      solver%diis_dimension = 8

   end subroutine cleanup_diis_tool
!
!
    subroutine update_diis_tool(solver, e, x)
!!
!!    Update 
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
      real(dp), dimension(solver%n_equations)  :: e
      real(dp), dimension(solver%n_parameters) :: x
!
      real(dp), dimension(:), allocatable :: x_i ! To hold previous x_i temporarily
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
         call solver%e_files(i)%open_()
         call solver%x_files(i)%open_()
!
      enddo
!
!     Compute the current dimensionality of the problem 
!     (1, 2,..., 7, 8, 8, 8, 8,...) for the standard diis_dimension = 8 without erasing history
!
      current_dim = solver%get_current_dim()
!
!     Write current e and x vector to file
!
      call solver%write_current_vecs(e, x, current_dim)
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
      call solver%construct_diis_matrix(current_dim, diis_matrix, e)
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
      if (info .ne. 0) call output%error_msg('could not solve DIIS matrix equation!', info)
!
!     Update the parameters (place in x on exit)
!
      call zero_array(x, solver%n_parameters)
!
      call mem%alloc(x_i, solver%n_parameters)
!
      do i = 1, current_dim
!
!        Read the x_i + Δ x_i vector
!
         call solver%x_files(i)%rewind_()
         call solver%x_files(i)%read_(x_i, solver%n_parameters)
!
!        Add w_i (x_i + Δ x_i) to the amplitudes
!
         call daxpy(solver%n_parameters, diis_vector(i), x_i, 1, x, 1)
!
      enddo
!
      call mem%dealloc(x_i, solver%n_parameters)
!
      call mem%dealloc(diis_vector, current_dim + 1)
!
!     Close files
!
      do i = 1, solver%diis_dimension

         call solver%e_files(i)%close_()
         call solver%x_files(i)%close_()

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
   subroutine write_current_vecs_diis_tool(solver, e, x, current_dim)
!!
!!    Write current vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Writes the current e and x to file. Each index 1, 2, ..., diis_dim
!!    has its own file. If more iterations have passed than we keep on file
!!    (i.e., solver%iteration > solver%diis_dimension), then we do a first-in/
!!    last-out replacement unless erase_history is true. 
!!
      implicit none 
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(in)  :: e
      real(dp), dimension(solver%n_parameters), intent(in) :: x
!
      integer, intent(in) :: current_dim
!
      integer :: k
!
      type(sequential_file) :: tmp_file
!
      if (current_dim .le. solver%diis_dimension &
            .and. current_dim .eq. solver%iteration &
            .or. solver%erase_history) then 
!
!        Just write to the correct file 
!
         call solver%e_files(current_dim)%rewind_()
         call solver%e_files(current_dim)%write_(e, solver%n_equations)
!
         call solver%x_files(current_dim)%rewind_()
         call solver%x_files(current_dim)%write_(x, solver%n_parameters)
!
      else 
!
!        Redefine files such that the current entry is the last, 
!        while all the others move one position up 
!
         tmp_file = solver%e_files(1)
!
         do k = 1, current_dim - 1 
!
            solver%e_files(k) = solver%e_files(k+1)
!
         enddo 
!
         solver%e_files(current_dim) = tmp_file
!
         call solver%e_files(current_dim)%rewind_()
         call solver%e_files(current_dim)%write_(e, solver%n_equations)
!
         tmp_file = solver%x_files(1)
!
         do k = 1, current_dim - 1 
!
            solver%x_files(k) = solver%x_files(k+1)
!
         enddo 
!
         solver%x_files(current_dim) = tmp_file
!
         call solver%x_files(current_dim)%rewind_()
         call solver%x_files(current_dim)%write_(x, solver%n_parameters)
!
      endif
!
   end subroutine write_current_vecs_diis_tool
!
!
   subroutine construct_diis_matrix_diis_tool(solver, current_dim, diis_matrix, e)
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
      real(dp), dimension(solver%n_equations), intent(in) :: e
!
      real(dp), dimension(:), allocatable :: e_i  ! To hold previous e_i temporarily
      real(dp), dimension(:), allocatable :: e_j  ! To hold previous e_j temporarily
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
         call solver%diis_matrix%open_()
!
         if (solver%iteration .gt. solver%diis_dimension .and. .not. solver%erase_history) then
!
!           First entry goes out, eliminating the (old) 1st row & column.
!
               call mem%alloc(diis_matrix_tmp, current_dim, current_dim)
!
               call solver%diis_matrix%rewind_()
!
               do j = 1, current_dim
                  do i = 1, current_dim
!
                     call solver%diis_matrix%read_(diis_matrix_tmp(i,j))
!
                  enddo
               enddo
!
               diis_matrix(1 : (current_dim - 1), 1 : (current_dim - 1)) = &
               diis_matrix_tmp(2 : current_dim, 2 : current_dim)      
!
               call mem%dealloc(diis_matrix_tmp, current_dim, current_dim)
!
         elseif (current_dim .gt. 1) then
!
!           Still doing the first set of vectors: just read the previous matrix 
!
            call solver%diis_matrix%rewind_()
!
            do j = 1, current_dim - 1
               do i = 1, current_dim - 1
!
                  call solver%diis_matrix%read_(diis_matrix(i,j))
!
               enddo
            enddo
!
         endif
!
!        Get the parts of the DIIS matrix G not constructed in
!        the previous iterations
!
         call mem%alloc(e_i, solver%n_equations)
!
         do i = 1, current_dim
!
            call solver%e_files(i)%rewind_()
            call solver%e_files(i)%read_(e_i, solver%n_equations)
!
            diis_matrix(current_dim, i) = ddot(solver%n_equations, e, 1, e_i, 1)
            diis_matrix(i, current_dim) = diis_matrix(current_dim, i)
!
            diis_matrix(current_dim + 1, i) = -one
            diis_matrix(i, current_dim + 1) = -one
!
         enddo
!
         call mem%dealloc(e_i, solver%n_equations)
!
!        Write the current DIIS matrix to file
!
         call solver%diis_matrix%rewind_()
!
         do j = 1, current_dim
            do i = 1, current_dim
!
               call solver%diis_matrix%write_(diis_matrix(i,j))
!
            enddo
         enddo
!
         call solver%diis_matrix%close_() 
!
      else ! Non-cumulative diis -> build diis matrix from scratch
!         
         call mem%alloc(e_i, solver%n_equations)
         call mem%alloc(e_j, solver%n_equations)
!
         do i = 1, current_dim
!
            call solver%e_files(i)%rewind_()
            call solver%e_files(i)%read_(e_i, solver%n_equations)
!
            do j = 1, i
!
               call solver%e_files(j)%rewind_()
               call solver%e_files(j)%read_(e_j, solver%n_equations)
!
               diis_matrix(i, j) = ddot(solver%n_equations, e_i, 1, e_j, 1)
               diis_matrix(j, i) = diis_matrix(i, j)
!
            enddo
!
            diis_matrix(current_dim + 1, i) = -one
            diis_matrix(i, current_dim + 1) = -one
!
         enddo 
!
         call mem%dealloc(e_i, solver%n_equations)
         call mem%dealloc(e_j, solver%n_equations)
!
      endif
!
   end subroutine construct_diis_matrix_diis_tool
!
!
   subroutine read_x_diis_tool(solver, x_i, i)
!!
!!    Get x_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_parameters), intent(out) :: x_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (i .gt. current_dim) call output%error_msg('Asked for an x not in use.')
!
      call solver%x_files(i)%open_('read', 'rewind')
!
      call solver%x_files(i)%read_(x_i, solver%n_parameters)
!
      call solver%x_files(i)%close_()
!
   end subroutine read_x_diis_tool
!
!
   subroutine read_e_diis_tool(solver, e_i, i)
!!
!!    Get e_i
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(out) :: e_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (i .gt. current_dim) call output%error_msg('Asked for an e not in use.')
!
      call solver%e_files(i)%open_('read', 'rewind')
!
      call solver%e_files(i)%read_(e_i, solver%n_equations)
!
      call solver%e_files(i)%close_()
!
   end subroutine read_e_diis_tool
!
!
   subroutine write_x_diis_tool(solver, x_i, i)
!!
!!    Set x_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_parameters), intent(in) :: x_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (solver%accumulate) call output%error_msg('Can not set x for cumulative diis.')
      if (i .gt. current_dim) call output%error_msg('Asked for an x not in use.')
!
      call solver%x_files(i)%open_('write', 'rewind')
!
      call solver%x_files(i)%write_(x_i, solver%n_parameters)
!
      call solver%x_files(i)%close_()
!
   end subroutine write_x_diis_tool
!
!
   subroutine write_e_diis_tool(solver, e_i, i)
!!
!!    Set e_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: solver
!
      real(dp), dimension(solver%n_equations), intent(in) :: e_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = solver%get_current_dim()
!
      if (solver%accumulate) call output%error_msg('Can not set e for cumulative diis.')
      if (i .gt. current_dim) call output%error_msg('Asked for an e not in use.')
!
      call solver%e_files(i)%open_('write', 'rewind')
!
      call solver%e_files(i)%write_(e_i, solver%n_equations)
!
      call solver%e_files(i)%close_()
!
   end subroutine write_e_diis_tool
!
!
end module diis_tool_class
