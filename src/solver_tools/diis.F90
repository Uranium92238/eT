module diis_solver_class
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
   type :: diis
!
      character(len=40), private :: name ! Solver name; determines the prefix of all DIIS files
!
      type(file), private :: dx   ! File containing the O vector from previous iterations
      type(file), private :: x_dx ! File containing the X vector from previous iterations
!
      type(file), private :: diis_matrix ! File containing the previously computed elements of the DIIS matrix
!
      integer(i15), private :: iteration = 1 ! Variable keeping track of the current DIIS iteration
                                             ! Note: defined to increment by +1 each time 'update' is called.
!
      integer(i15), private :: diis_dimension = 8   ! Standard is 8, though it might be useful to change this value
!
      integer(i15), private :: n_parameters ! The length of the X vector
      integer(i15), private :: n_equations  ! The length of the O vector

!
   contains
!
      procedure :: init     => init_diis     ! Called to set up a DIIS solver
      procedure :: update   => update_diis   ! Called to perform a DIIS step
      procedure :: finalize => finalize_diis
!
      procedure :: get_current_index => get_current_index_diis
!
   end type diis
!
!
contains
!
!
   subroutine finalize_diis(solver)
!!
!!    Finalize DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
      implicit none
!
      class(diis) :: solver
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
   end subroutine finalize_diis
!
!
   subroutine init_diis(solver, name, n_parameters, n_equations, diis_dimension)
!!
!!    Init DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
      class(diis) :: solver
!
      character(len=*), intent(in) :: name
!
      integer(i15), intent(in) :: n_parameters
      integer(i15), intent(in) :: n_equations
!
      integer(i15), intent(in), optional :: diis_dimension
!
      solver%name = trim(name)
      solver%n_parameters = n_parameters
      solver%n_equations = n_equations
!
      if (present(diis_dimension)) solver%diis_dimension = diis_dimension
!
      call solver%x_dx%init(trim(solver%name) // '_x_dx', 'sequential', 'unformatted')
      call solver%dx%init(trim(solver%name) // '_dx', 'sequential', 'unformatted')
      call solver%diis_matrix%init(trim(solver%name) // '_diis_matrix', 'sequential', 'unformatted')
!
   end subroutine init_diis
!
!
    subroutine update_diis(solver, dx, x_dx)
!!
!!    Update (DIIS)
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
      class(diis) :: solver
!
      real(dp), dimension(solver%n_equations, 1)  :: dx
      real(dp), dimension(solver%n_parameters, 1) :: x_dx
!
      real(dp), dimension(:,:), allocatable :: dx_i   ! To hold previous Δ x_i temporarily
      real(dp), dimension(:,:), allocatable :: x_dx_i ! To hold previous x_dx_i temporarily
!
      real(dp) :: ddot
!
      integer(i15) :: i = 0, j = 0
!
      integer      :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer(i15) :: current_index = 0 ! Progressing as follows: 1,2,...,7,8,1,2,...
!
      integer(i15) :: dummy = 0
!
      real(dp), dimension(:,:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer(i15), dimension(:,:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
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
      call mem%alloc_int(ipiv, solver%diis_dimension, 1)
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
   end subroutine update_diis
!
!
!    subroutine update_diis(solver, o, x)
! !!
! !!    Update (DIIS)
! !!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
! !!
! !!    Based on the CCSD specific routine written by Sarai D. Folkestad
! !!    and Eirik F. Kjønstad, May 2017.
! !!
! !!    Solves the DIIS problem, as described in the documentation of DIIS class module
! !!    (see the top of this file).
! !!
!       implicit none
! !
!       class(diis) :: solver
! !
!       real(dp), dimension(solver%n_equations, 1)  :: o
!       real(dp), dimension(solver%n_parameters, 1) :: x
! !
!       real(dp), dimension(:,:), allocatable :: o_i ! To hold previous O_i temporarily
!       real(dp), dimension(:,:), allocatable :: x_i ! To hold previous X_i temporarily
! !
!       real(dp) :: ddot
! !
!       integer(i15) :: i = 0, j = 0, ioerror = 0
! !
!       integer(kind=4) :: info = -1 ! Error integer for dgesv
!       integer(i15) :: current_index = 0 ! Progressing as follows: 1,2,...,7,8,1,2,...
! !
!       integer(i15) :: dummy = 0
! !
!       real(dp), dimension(:,:), allocatable :: diis_vector
!       real(dp), dimension(:,:), allocatable :: diis_matrix
! !
!       integer(i15), dimension(:,:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
! !
! !     :: Open DIIS files
! !
! !     Initialize files
! !
!       call solver%o%init(trim(solver%name) // '_o', 'sequential', 'unformatted')
!       call solver%x%init(trim(solver%name) // '_x', 'sequential', 'unformatted')
!       call solver%diis_matrix%init(trim(solver%name) // '_diis_matrix', 'sequential', 'unformatted')
! !
! !     Ask disk manager to open the files
! !
!       call disk%open_file(solver%o, 'readwrite')
!       call disk%open_file(solver%x, 'readwrite')
!       call disk%open_file(solver%diis_matrix, 'readwrite')
! !
! !     :: Compute the current index
! !     (1,2,...,7,8,1,2,...) for the standard diis_dimension = 8
! !
!       current_index = solver%iteration - &
!          ((solver%diis_dimension)-1)*((solver%iteration-1)/((solver%diis_dimension)-1))
! !
! !     :: Save O_i and X_i to files
! !
!       if (current_index .eq. 1) then
! !
!          rewind(solver%o%unit)
!          rewind(solver%x%unit)
! !
!       else
! !
!          do dummy = 1, current_index - 1 ! Get to the correct line in the file
! !
!             read(solver%o%unit)
!             read(solver%x%unit)
! !
!          enddo
! !
!       endif
! !
!       write(solver%o%unit) (o(i,1), i = 1, solver%n_equations)
!       write(solver%x%unit) (x(i,1), i = 1, solver%n_parameters)
! !
! !     :: Solve the least squares problem, G * w = H
! !
! !        G : DIIS matrix, G_ij = O_i O_j,
! !        H : DIIS vector, H_i = 0,
! !
! !     where i, j = 1, 2, ..., current_index. To enforce normality
! !     of the solution, G is extended with a row & column of -1's
! !     and H with a -1 at the end.
! !
! !     First set the DIIS vector to one
! !
!       call mem%alloc(diis_vector, current_index+1, 1)
!       diis_vector = zero
! !
! !     Allocate the DIIS matrix and read in previous matrix elements
! !
!       call mem%alloc(diis_matrix, current_index+1, current_index+1)
!       diis_matrix = zero
! !
!       if (current_index .gt. 1) then
! !
!          rewind(solver%diis_matrix%unit)
! !
!          do j = 1, current_index - 1
!             do i = 1, current_index - 1
! !
!                read(solver%diis_matrix%unit) diis_matrix(i,j)
! !
!             enddo
!          enddo
! !
!       endif
! !
! !     Get the parts of the DIIS matrix G not constructed in
! !     the previous iterations
! !
!       call mem%alloc(o_i, solver%n_equations, 1)
!       o_i = zero
! !
!       rewind(solver%o%unit)
! !
!       do i = 1, current_index
! !
!          ioerror = 0
!          read(solver%o%unit, iostat=ioerror) (o_i(j,1), j = 1, solver%n_equations)
! !
!          diis_matrix(current_index,i) = ddot(solver%n_equations, o, 1, o_i, 1)
!          diis_matrix(i,current_index) = diis_matrix(current_index,i)
! !
!          diis_matrix(current_index+1,i) = -one
!          diis_matrix(i,current_index+1) = -one
! !
!       enddo
! !
!       diis_vector(current_index+1,1) = -one
! !
! !     Write the current DIIS matrix to file
! !
!       rewind(solver%diis_matrix%unit)
! !
!       do j = 1, current_index
!          do i = 1, current_index
! !
!             write(solver%diis_matrix%unit) diis_matrix(i,j)
! !
!          enddo
!       enddo
! !
! !     Solve the DIIS equation
! !
!       call mem%alloc_int(ipiv, solver%diis_dimension, 1)
!       ipiv = 0
! !
! !     Note: on exit, the solution is in the diis_vector,
! !     provided info = 0 (see LAPACK documentation for more)
! !
!      ! write(output%unit, *) 'diis mat: '
! !
!      ! do i = 1, current_index + 1
!      !    write(output%unit, *) 'row ', i
!      !    write(output%unit, *) diis_matrix(i, :)
!      ! enddo
! !
!       call dgesv(current_index+1,  &
!                   1,               &
!                   diis_matrix,     &
!                   current_index+1, &
!                   ipiv,            &
!                   diis_vector,     &
!                   current_index+1, &
!                   info)
! !
!       if (info .ne. 0) then
! !
!          write(output%unit, *) 'Something went wrong in DIIS (dgesv). Info: ', info
!          stop
! !
!       endif
! !
!       write(output%unit, *) 'diis vec sol: ', diis_vector(:,1)
! !
! !     :: Update the parameters (placed in x on exit)
! !
!       call mem%dealloc(o_i, solver%n_equations, 1)
! !
!       call mem%alloc(x_i, solver%n_parameters, 1)
! !
!       rewind(solver%x%unit)
! !
!       do i = 1, solver%n_parameters
! !
!          x(i,1) = zero
! !
!       enddo
! !
!       o = zero
!       do i = 1, current_index
! !
! !        Add w_i x_i to the parameters x
! !
!          x_i = zero
!          read(solver%x%unit) (x_i(j, 1), j = 1, solver%n_parameters)
! !
!       !   call daxpy(solver%n_parameters, diis_vector(i, 1), x_i, 1, x, 1)
!          call daxpy(solver%n_parameters, diis_vector(i, 1), x_i, 1, o, 1)
! !
!       enddo
! !
!       call mem%dealloc(x_i, solver%n_parameters, 1)
! !
! !     Deallocations
! !
!       call mem%dealloc(diis_vector, current_index + 1, 1)
!       call mem%dealloc(diis_matrix, current_index + 1, current_index+1)
! !
! !     Close files
! !
!       call disk%close_file(solver%o)
!       call disk%close_file(solver%x)
!       call disk%close_file(solver%diis_matrix)
! !
!       solver%iteration = solver%iteration + 1
! !
!    end subroutine update_diis
!
!
   function get_current_index_diis(solver)
!
      implicit none
!
      class(diis) :: solver
!
      integer(i15) :: get_current_index_diis
!
      get_current_index_diis = solver%iteration - &
               ((solver%diis_dimension)-1)*((solver%iteration-1)/((solver%diis_dimension)-1))
!
   end function get_current_index_diis
!
!
end module diis_solver_class
