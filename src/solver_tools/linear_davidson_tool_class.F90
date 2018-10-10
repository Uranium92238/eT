module linear_davidson_tool_class
!
!!
!!    Eigenvalue davidson tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    A tool to help solve an eigenvalue equation A X_n = F X_n
!!    for symmetric A using the Davidson algorithm. It is tailored to
!!    be usable also in cases where A cannot be stored, but where the 
!!    transformation X -> A X is implemented. A typical Davidson loop
!!    will use this tool as follows:
!!
!!       To-write.
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
   use davidson_tool_class
   use array_utilities
   use array_analysis
!
   type, extends(davidson_tool) :: linear_davidson_tool 
!
!
   contains 
!
   !  procedure :: prepare                  => prepare_linear_davidson_tool 
   !  procedure :: cleanup                  => cleanup_linear_davidson_tool
!
      procedure :: construct_next_trial_vec => construct_next_trial_vec_linear_davidson_tool
!
      procedure :: solve_reduced_problem    => solve_reduced_problem_linear_davidson_tool
!
      procedure :: construct_residual       => construct_residual_linear_davidson_tool
!
   end type linear_davidson_tool
!
contains
!
!
   subroutine construct_reduced_gradient(davidson)
!!
!!    Construct reduced gradient
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!! 
!!    Constructs the gradient vector in the reduced space
!!
!!       F_i = c_i^T * F
!!
      implicit none
!
      class(linear_davidson_tool) :: davidson
!
!     Open file
!
      call disk%open_file(davidson%trials, 'read')
!
!     Rewind file
!
      rewind(davidson%trials%unit)
!
!     Construct reduced vector
!
      if (allocated(davidson%F_red)) then 
!
         call mem%alloc(F_red_copy, davidson%dim_red - davidson%n_new_trials, 1)
         F_red_copy = davidson%F_red
!
         call mem%dealloc(davidson%F_red, davidson%dim_red - davidson%new_trials, 1)
         call mem%alloc(davidson%F_red, davidson%dim_red, 1)
!
         davidson%F_red(1:davidson%dim_red - davidson%new_trials, 1) = F_red_copy
!
      else
!
         call mem%alloc(davidson%F_red, davidson%dim_red, 1)
!
      endif
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      do i = davidson%dim_red - davidson%new_trials + 1, davidson%dim_red 
!
         call davidson%read_trial(c_i, i)
         davidson%F_red(i, 1) = ddot(davidson%n_parameters, c_i, 1, davidson%F, 1)
!
      enddo
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
!
!     Close File
!
      call disk%close_file(davidson%trials)
!
   end subroutine construct_reduced_gradient
!
!
!
   subroutine solve_reduced_problem_linear_davidson_tool(davidson)
!!
!!    Solve reduced problem 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!       A X = F,  
!!
!!    On exit, the eigenvalues are stored in the vectors omega_re and omega_im,
!!    where the real and imaginary parts of eigenvalue n is stored in the nth 
!!    row of these vectors.
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      integer(i15), dimension(:,:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      integer :: info = -1 ! Error integer for dgesv routine (LU factorization)
!
!     Solve the linear problem
!
      call wf%mem%alloc_int(ipiv, reduced_dim, 1)
      ipiv = 0
      info = 0
!
      solver%X_red = solver%F_red
!
      call dgesv(davidson%dim_red,  &
                  1,                & ! Number of RHS 
                  davidson%A_red,   &
                  davidson%dim_red, &
                  ipiv,             &
                  solver%X_red,     &
                  davidson%dim_red, &
                  info)
!
      if (info .ne. 0) then 
!
         write(unit_output,*) 'Error: could not solve reduced response equation.', info
         stop
!
      endif 
!
   end subroutine solve_reduced_problem_linear_davidson_tool
!
!
   subroutine construct_residual_linear_davidson_tool(davidson, R, X, norm_X)
!!
!!    Construct residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Constructs the residual
!!
!!       R = (A*X - F)/|X|
!!
      implicit none 
!
      class(linear_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1), intent(out)  :: R 
      real(dp), dimension(davidson%n_parameters, 1), intent(in)   :: X 
!
      real(dp), intent(in)     :: norm_X 
!
      call davidson%construct_AX(R, 1)                               ! R = AX 
      call daxpy(davidson%n_parameters, -one, davidson%F, 1, R, 1)   ! R = AX - F
      call dscal(davidson%n_parameters, one/norm_X, R, 1)            ! R = (AX - F)/|X|
!
   end subroutine construct_residual_linear_davidson_tool
!
!
   subroutine construct_next_trial_vec_eigen_davidson_tool(davidson, residual_norm, iteration)
!!
!!    Construct next trial vector  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    This routine decides whether to construct a new trial vector based on 
!!    the current residual norm. The new trial vector is a preconditioned 
!!    residual orthogonalized against the search space (the existing trial
!!    vectors). If it turns out that the root is not converged and the 
!!    new trial vector contains no new components (linear dependence),
!!    the routine terminates in an error. 
!!
!!    Note that the residual norm is set on exit and need not to be 
!!    calculated beforehand. 
!!    
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), intent(out) :: residual_norm 
!
      integer(i15) :: iteration
!
      real(dp) :: norm_X, norm_new_trial, norm_residual, norm_precond_residual
!
      real(dp), dimension(:,:), allocatable :: R, X 
!
!     Construct full space solution vector X, 
!     and the associated residual R 
!
      call mem%alloc(X, davidson%n_parameters, 1)
      call mem%alloc(R, davidson%n_parameters, 1)
!
      call davidson%construct_X(X, 1) 
      norm_X = get_l2_norm(X, davidson%n_parameters) 
!
      call davidson%construct_residual(R, X, norm_X)
!
      call disk%open_file(davidson%X, 'write', 'rewind')
!
      call dscal(davidson%n_parameters, one/norm_X, X, 1)
!
      write(davidson%X%unit) X
!
      call disk%close_file(davidson%X)
!
      call mem%dealloc(X, davidson%n_parameters, 1)
!
!     Calculate the norm of the residual and test for convergence. If not
!     converged, the residual is preconditioned & we remove components already 
!     in the search space by orthogonalizing it against existing trial vectors.
!
      norm_residual = get_l2_norm(R, davidson%n_parameters)
      residual_norm = norm_residual
!
      if (norm_residual .gt. davidson%residual_threshold) then 
!
         call davidson%projection(R)
         !call davidson%precondition(R)
!
         norm_precond_residual = get_l2_norm(R, davidson%n_parameters)
         call dscal(davidson%n_parameters, one/norm_precond_residual, R, 1)
!
         call davidson%orthogonalize_against_trial_vecs(R)
         norm_new_trial = get_l2_norm(R, davidson%n_parameters)
!
         if (norm_new_trial .gt. davidson%residual_threshold) then
!
            davidson%n_new_trials = davidson%n_new_trials + 1
            call dscal(davidson%n_parameters, one/norm_new_trial, R, 1)
!
            call disk%open_file(davidson%trials, 'write', 'append')
            write(davidson%trials%unit) R
            call disk%close_file(davidson%trials)
!
         endif 
!
      endif 
!
      call mem%dealloc(R, davidson%n_parameters, 1)
!
   end subroutine construct_next_trial_vec_eigen_davidson_tool
!
!
end module linear_davidson_tool_class
