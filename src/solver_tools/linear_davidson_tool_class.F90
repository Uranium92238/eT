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
      real(dp), dimension(:,:), allocatable :: F_red
      real(dp), dimension(:,:), allocatable :: F
!
   contains 
!
      procedure :: prepare                      => prepare_linear_davidson_tool 
!
      procedure :: construct_next_trial_vec     => construct_next_trial_vec_linear_davidson_tool
!
      procedure :: solve_reduced_problem        => solve_reduced_problem_linear_davidson_tool
!
      procedure :: construct_residual           => construct_residual_linear_davidson_tool
!
      procedure :: construct_reduced_gradient   => construct_reduced_gradient_linear_davidson_tool
!
!     Routines for linear response solver
!
      procedure :: prepare_response              => prepare_response_linear_davidson_tool
!
   end type linear_davidson_tool
!
contains
!
!
   subroutine prepare_linear_davidson_tool(davidson, name, n_parameters, residual_threshold, F)
!!
!!    Initialize 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name
!
      integer, intent(in) :: n_parameters 
!
      real(dp), dimension(n_parameters, 1) :: F 
!
      real(dp), intent(in)     :: residual_threshold
!
      call mem%alloc(davidson%F, n_parameters, 1)
      davidson%F = F
!
      davidson%n_parameters = n_parameters
      davidson%n_solutions  = 1
!
      davidson%residual_threshold   = residual_threshold
!
      davidson%name = trim(name)
!
      call davidson%X%init(trim(davidson%name) // '_X', 'sequential', 'unformatted')
      call davidson%trials%init(trim(davidson%name) // '_trials', 'sequential', 'unformatted')
      call davidson%transforms%init(trim(davidson%name) // '_transforms', 'sequential', 'unformatted')
      call davidson%preconditioner%init(trim(davidson%name) // '_preconditioner', 'sequential', 'unformatted')
!
!     For safety, delete old files if they are on disk
!
       call disk%delete(davidson%trials)
       call disk%delete(davidson%transforms)
       call disk%delete(davidson%X)
!
      davidson%do_precondition   = .false.         ! Switches to true if 'set_preconditioner' is called
!
      davidson%dim_red           = davidson%n_solutions     ! Initial dimension equal to number of solutions
      davidson%n_new_trials      = davidson%n_solutions 
!
      davidson%max_dim_red = 150
!
      davidson%current_n_trials = 0  
!
   end subroutine prepare_linear_davidson_tool
!
!
   subroutine construct_reduced_gradient_linear_davidson_tool(davidson)
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
      real(dp), dimension(:,:), allocatable :: c_i
!
      integer :: i
!
      real(dp) :: ddot
!
!     Construct reduced vector
!
      if (allocated(davidson%F_red)) call mem%dealloc(davidson%F_red, davidson%dim_red - 1, 1) 
!
      call mem%alloc(davidson%F_red, davidson%dim_red, 1)
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      call disk%open_file(davidson%trials, 'read')
      rewind(davidson%trials%unit)
!
      do i = 1, davidson%dim_red 
!
         read(davidson%trials%unit)c_i
!
         davidson%F_red(i, 1) = ddot(davidson%n_parameters, c_i, 1, davidson%F, 1)
!
      enddo
!
      call disk%close_file(davidson%trials)
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
!
   end subroutine construct_reduced_gradient_linear_davidson_tool
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
      integer, dimension(:,:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      real(dp), dimension(:,:), allocatable :: A_red_copy
!
      integer :: info = -1 ! Error integer for dgesv routine (LU factorization)
!
!     Solve the linear problem
!
      call mem%alloc(ipiv, davidson%dim_red, 1)
      ipiv = 0
      info = 0
!
      if (allocated(davidson%X_red)) call mem%dealloc(davidson%X_red, davidson%dim_red - 1, 1 )
!
      call mem%alloc(davidson%X_red, davidson%dim_red, 1 )
!
      davidson%X_red = davidson%F_red
!
      call mem%alloc(A_red_copy, davidson%dim_red, 1 )
      A_red_copy = davidson%A_red
!
      call dgesv(davidson%dim_red,  &
                  1,                & ! Number of RHS 
                  A_red_copy,       &
                  davidson%dim_red, &
                  ipiv,             &
                  davidson%X_red,   &
                  davidson%dim_red, &
                  info)
!
      call mem%dealloc(A_red_copy, davidson%dim_red, 1 )
      call mem%dealloc(ipiv, davidson%dim_red, 1)
!
      if (info .ne. 0) then 
!
         call output%error_msg('could not solve reduced response equation.')
!
      endif 
!
   end subroutine solve_reduced_problem_linear_davidson_tool
!
!
   subroutine construct_residual_linear_davidson_tool(davidson, R, n)
!!
!!    Construct residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Constructs the residual
!!
!!       R = A*X - F
!!
      implicit none 
!
      class(linear_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: R 
!
      integer, intent(in) :: n
!
      if (n .ne. 1) then
!
         call output%error_msg('for linear equations with davidson n = 1')
!
      endif
!
      call davidson%construct_AX(R, 1)                               ! R = AX 
      call daxpy(davidson%n_parameters, -one, davidson%F, 1, R, 1)   ! R = AX - F
!
   end subroutine construct_residual_linear_davidson_tool
!
!
   subroutine construct_next_trial_vec_linear_davidson_tool(davidson, residual_norm)
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
      class(linear_davidson_tool) :: davidson 
!
      real(dp), intent(out) :: residual_norm 
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
      call davidson%construct_residual(R, 1)
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
         call davidson%precondition(R)
!
         norm_precond_residual = get_l2_norm(R, davidson%n_parameters)
         call dscal(davidson%n_parameters, one/norm_precond_residual, R, 1)
!
         call davidson%orthogonalize_against_trial_vecs(R)
!
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
         else
!
            call output%error_msg('did not find any new trials.')
!
         endif 
!
      endif 
!
      call mem%dealloc(R, davidson%n_parameters, 1)
!
   end subroutine construct_next_trial_vec_linear_davidson_tool
!
!
!  Routines below for linear response solver
!
   subroutine prepare_response_linear_davidson_tool(davidson, name, n_parameters, residual_threshold, F, dim_f)
!!
!!    Initialize for response solver davidson object 
!!    Written by Josefine H. Andersen, March 2019
!!
      implicit none
!
      class(linear_davidson_tool) :: davidson
!
      character(len=*), intent(in) :: name
!
      integer, intent(in) :: n_parameters
!
      integer, intent(in) :: dim_f
!
      real(dp), dimension(n_parameters, dim_f) :: F
!
      real(dp), intent(in)     :: residual_threshold
!
      call mem%alloc(davidson%F, n_parameters, dim_f)
      davidson%F = F
!
      davidson%n_parameters = n_parameters
      davidson%n_solutions  = 1
!
      davidson%residual_threshold   = residual_threshold
!
      davidson%name = trim(name)
!
      call davidson%X%init(trim(davidson%name) // '_X', 'sequential', 'unformatted')
      call davidson%trials%init(trim(davidson%name) // '_trials', 'sequential', 'unformatted')
      call davidson%transforms%init(trim(davidson%name) // '_transforms', 'sequential', 'unformatted')
      call davidson%preconditioner%init(trim(davidson%name) // '_preconditioner', 'sequential', 'unformatted')
!
!     For safety, delete old files if they are on disk
!
       call disk%delete(davidson%trials)
       call disk%delete(davidson%transforms)
       call disk%delete(davidson%X)
!
      davidson%do_precondition   = .false.         ! Switches to true if 'set_preconditioner' is called
!
      davidson%dim_red           = davidson%n_solutions     ! Initial dimension equal to number of solutions
      davidson%n_new_trials      = davidson%n_solutions
!
      davidson%max_dim_red = 150
!
      davidson%current_n_trials = 0
!
   end subroutine prepare_response_linear_davidson_tool
!
!
end module linear_davidson_tool_class
