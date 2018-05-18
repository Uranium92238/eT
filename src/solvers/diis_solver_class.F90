module diis_solver_class
!
!!
!!                             DIIS solver class module                                 
!!             Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018         
!!
!! 
!!       Objects of this class can be instantiated in order to solve a non-linear set of equations
!!       of the form O(X) = 0, where O and X are vectors of length N. It is solved by updating the X as
!!
!!          X^(n+1) = sum_k=1^n w_k (X^(k) + dX^(k)),
!!
!!       where w_k are a set of weights determined so as to minimize the norm of the 
!!       average error vector 
!!
!!          F(w) = sum_k w_k dX^(k),
!!
!!       with the constraint that sum_k w_k = 1. The "perturbation estimate" dX^(k)
!!       may be chosen in several ways, but dX^(n) = 0 must imply that O(X) = 0.
!!
!!       The number of previous steps to include in the k-sums is the DIIS dimensionality,
!!       for which the optimal value usually depends on the equation that is being solved.
!!       
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
   use input_output
!
   use file_class
   use disk_manager_class
   use memory_manager_class
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the file class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
   type :: diis  
!
      character(len=40) :: name ! Solver name; determines the prefix of all DIIS files 
!
      type(file) :: dx          ! File containing the dX vector from previous iterations 
      type(file) :: x_dx        ! File containing the X + dX vector from previous iterations 
      type(file) :: diis_matrix ! File containing the previously computed elements of the DIIS matrix 
!
      integer(i15) :: iteration = 1        ! Variable keeping track of the current DIIS iteration
                                           ! Note: defined to increment by +1 each time 'update' is called.
!
      integer(i15) :: diis_dimension = 8   ! Standard is 8, though it might be useful to change this value
      integer(i15) :: n_parameters         ! The length of the O, X and dX vectors 
!
   contains 
!
      procedure :: init   => init_diis   ! Called to set up a DIIS solver 
      procedure :: update => update_diis ! Called to perform a DIIS step 
!
   end type diis                                                                             
!
!
contains
!
!
   module subroutine init_diis(solver, name, n_parameters, diis_dimension)
!!
!!    Init DIIS 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
      class(diis) :: solver 
!
      character(len=*), intent(in) :: name 
      integer(i15), intent(in)     :: n_parameters
!
      integer(i15), intent(in), optional :: diis_dimension
!
      solver%name = trim(name) 
      solver%n_parameters = n_parameters
!
      if (present(diis_dimension)) solver%diis_dimension = diis_dimension
!
   end subroutine init_diis
!
!
    module subroutine update_diis(solver, dx, x_dx, disk, mem)
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
      real(dp), dimension(solver%n_parameters, 1) :: dx
      real(dp), dimension(solver%n_parameters, 1) :: x_dx 
!
      class(disk_manager)   :: disk 
      class(memory_manager) :: mem  
!
      real(dp), dimension(:,:), allocatable :: dx_i ! To hold previous Δ x_i temporarily
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
!     Set file names 
!
      solver%dx%name          = trim(solver%name) // '_dx'
      solver%x_dx%name        = trim(solver%name) // '_x_dx' 
      solver%diis_matrix%name = trim(solver%name) // '_diis_matrix'
!
!     Ask disk manager to open the files 
!
      call disk%open_file(solver%dx, 'unformatted', 'readwrite', 'sequential')
      call disk%open_file(solver%x_dx, 'unformatted', 'readwrite', 'sequential')
      call disk%open_file(solver%diis_matrix, 'unformatted', 'readwrite', 'sequential')
!
!     :: Compute the current index 
!     (1,2,...,7,8,1,2,...) for the standard diis_dimension = 8 
!
      current_index = solver%iteration - ((solver%diis_dimension)-1)*((solver%iteration-1)/((solver%diis_dimension)-1)) 
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
      write(solver%dx%unit)   (dx(i,1), i = 1, solver%n_parameters)
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
      call mem%alloc(dx_i, solver%n_parameters, 1) ! Allocate temporary holder of quasi-Newton estimates
      dx_i = zero 
!
      rewind(solver%dx%unit)
!
      do i = 1, current_index
!
         read(solver%dx%unit) (dx_i(j,1), j = 1, solver%n_parameters) 
!
         diis_matrix(current_index,i) = ddot(solver%n_parameters, dx, 1, dx_i, 1) 
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
!     :: Update the parameters (placed in dx on exit) 
!
      dx = zero
!
      rewind(solver%x_dx%unit)
!
      do i = 1, current_index
!
!        Read the x_i + Δ x_i vector 
!
         x_dx = zero
         read(solver%x_dx%unit) (x_dx(j, 1), j = 1, solver%n_parameters)
!
!        Add w_i (x_i + Δ x_i) to the amplitudes 
!
         call daxpy(solver%n_parameters, diis_vector(i, 1), x_dx, 1, dx, 1)
!
      enddo
!
!     Deallocations 
!
      call mem%dealloc(dx_i, solver%n_parameters, 1)
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
end module diis_solver_class
