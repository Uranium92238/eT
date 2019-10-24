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
!!    Linear davidson tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    A tool to help solve an eigenvalue equation A X_n = F X_n
!!    for symmetric A using the Davidson algorithm. It is tailored to
!!    be usable also in cases where A cannot be stored, but where the 
!!    transformation X -> A X is implemented. 
!!
!!    Using the tool: after using the constructor, 
!!
!!       1. Set the initial set of trial vectors by constructing 
!!          them and calling davidson%set_trial(c). When you are 
!!          done writing trials, call davidson%orthonormalize_trial_vecs() 
!!          to make sure the trial vectors are orthonormalized.
!!
!!       2. (optional) Set preconditioner by call davidson%set_preconditioner(P),
!!          where P^-1 is a diagonal preconditioner for A. For the CC Jacobian, P
!!          is the vector of orbital differences.
!!
!!       3. Set up the iterative loop, which consists of four stages:
!!
!!          do while (.not. converged ...)
!!    
!!             I. Let the tool know a new iteration has begun 
!!
!!             call davidson%iterate()
!!
!!             II. Read new trial, transform it, & store the result 
!!
!!             call davidson%get_trial(c, davidson%dim_red)
!!             Transform: c <- A c 
!!             call davidson%set_transform(c, davidson%dim_red)
!!
!!             III. Solve reduced problem 
!!             
!!             call davidson%solve_reduced_problem()
!!
!!             IV. Construct residual and use it to generate new trial vectors 
!!
!!             call davidson%construct_residual(R)
!!                
!!             ...
!!
!!             if (residual_norm >= thr) call davidson%construct_next_trial(R)
!!
!!          enddo ! end of iterative loop 
!!
!
   use davidson_tool_class
!
   use memory_manager_class, only: mem 
!
   type, extends(davidson_tool) :: linear_davidson_tool 
!
      real(dp), dimension(:), allocatable :: F_red
      real(dp), dimension(:), allocatable :: F
!
   contains 
!
!     Prepare and Cleanup
!
      procedure :: cleanup                      => cleanup_linear_davidson_tool 
!
!     Linear Davidson specific routines
!
      procedure :: construct_next_trial         => construct_next_trial_linear_davidson_tool
      procedure :: construct_residual           => construct_residual_linear_davidson_tool
      procedure :: solve_reduced_problem        => solve_reduced_problem_linear_davidson_tool
!
      procedure :: construct_reduced_gradient   => construct_reduced_gradient_linear_davidson_tool
!
   end type linear_davidson_tool
!
!
   interface linear_davidson_tool
!
      procedure :: new_linear_davidson_tool
!
   end interface linear_davidson_tool
!
!
contains
!
!
   function new_linear_davidson_tool(name_, n_parameters, &
            lindep_threshold, max_dim_red, F, records_in_memory) result(davidson)
!!
!!    New linear Davidson tool 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    name_ :            Name of solver tool (used for temporary files)
!!
!!    n_parameters:      Dimensionality of the full vector space (A is n_parameters x n_parameters)
!!
!!    lindep_threshold:  Norm threshold for new trial vector after being added to the trial 
!!                       space. If the vector, after orthonormalization against the existing
!!                       trial space, has a norm below this threshold, then we assume that the 
!!                       trial basis has become linearly dependent. This causes an error stop. 
!!
!!    max_dim_red:       Maximum dimension of the reduced space. When exceeding this dimensionality,
!!                       the solutions are set as the basis for the new trial space. 
!!
!!    F:                 The linear equation is A X = F. Here F is n_parameters long. F is passed to 
!!                       the tool and kept in the tool until the tool is deallocated.
!!
!!    records_in_memory: If .true., the trials and transforms are kept in memory. Otherwise they 
!!                       are stored on file. 
!!
      implicit none 
!
      type(linear_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters, max_dim_red 
!
      real(dp), dimension(n_parameters) :: F 
!
      real(dp), intent(in) :: lindep_threshold
!
      logical, intent(in) :: records_in_memory 
!
      call mem%alloc(davidson%F, n_parameters)
      davidson%F = F
!
      davidson%n_parameters = n_parameters
      davidson%n_solutions  = 1
!
      davidson%lindep_threshold = lindep_threshold
      davidson%max_dim_red = max_dim_red
!
      davidson%name_ = trim(name_)
!
      davidson%do_precondition   = .false. ! Switches to true if 'set_preconditioner' is called
!
      davidson%dim_red      = 0    
      davidson%n_new_trials = davidson%n_solutions 
!
!     Set up array/or file array for trials and transforms    
!
      call davidson%prepare_trials_and_transforms(records_in_memory)
!
   end function new_linear_davidson_tool
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
      real(dp), dimension(:), allocatable :: c_i
!
      integer :: i
!
      real(dp) :: ddot
!
!     Construct reduced vector
!
      if (allocated(davidson%F_red)) call mem%dealloc(davidson%F_red, davidson%dim_red - 1) 
!
      call mem%alloc(davidson%F_red, davidson%dim_red)
      call mem%alloc(c_i, davidson%n_parameters)
!
      do i = 1, davidson%dim_red 
!
         call davidson%get_trial(c_i, i)
!
         davidson%F_red(i) = ddot(davidson%n_parameters, c_i, 1, davidson%F, 1)
!
      enddo
!
      call mem%dealloc(c_i, davidson%n_parameters)
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
!!    The routine performs three actions: 1) constructs the reduced space quantities,
!!    2) uses these to solve the reduced space problem, and 3) sets number of new trials 
!!    to zero (prepares the tool for the next task, which is the receiving of new residuals
!!    to make another new set of trial vectors).
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      real(dp), dimension(:,:), allocatable :: A_red_copy
!
      integer :: info = -1 ! Error integer for dgesv routine (LU factorization)
!
!     Construct reduced space quantities 
!
      call davidson%construct_reduced_matrix()
      call davidson%construct_reduced_gradient()
!
!     Solve the linear problem
!
      call mem%alloc(ipiv, davidson%dim_red)
      ipiv = 0
      info = 0
!
      if (allocated(davidson%X_red)) call mem%dealloc(davidson%X_red, davidson%dim_red - 1, 1)
!
      call mem%alloc(davidson%X_red, davidson%dim_red, 1)
!
      davidson%X_red(:, 1) = davidson%F_red
!
      call mem%alloc(A_red_copy, davidson%dim_red, davidson%dim_red)
!
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
      call mem%dealloc(A_red_copy, davidson%dim_red, davidson%dim_red)
!
      call mem%dealloc(ipiv, davidson%dim_red)
!
      if (info .ne. 0) then 
!
         call output%error_msg('could not solve reduced response equation.')
!
      endif 
!
      davidson%n_new_trials = 0
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
      class(linear_davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(davidson%n_parameters) :: R 
!
      integer, intent(in) :: n
!
      if (n /= 1) then
!
         call output%error_msg('for linear equations, n must be 1 in construct_residual_linear_davidson_tool')
!
      endif
!
      call davidson%construct_AX(R, 1)                               ! R = AX 
      call daxpy(davidson%n_parameters, -one, davidson%F, 1, R, 1)   ! R = AX - F
!
   end subroutine construct_residual_linear_davidson_tool
!
!
   subroutine construct_next_trial_linear_davidson_tool(davidson, R, alpha)
!!
!!    Construct next trial  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    R:     residual, AX - F
!!    alpha: shift in preconditioner (optional with deafult value 0)
!!
!!    The routine constructs a new trial vector from the residual R by 
!!    1) preconditioning it and 2) orthonormalizing the result with respect 
!!    to previous trial vectors. If the trial vector thus generated is 
!!    too close to zero, the routine exits with an error. Otherwise,
!!    a new trial vector has been added to the search space.
!!    
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: R 
!
      real(dp), optional :: alpha
!
      real(dp) :: norm_trial, alpha_local
      real(dp), dimension(:), allocatable :: trial 
!
!     Precondition 
!
      alpha_local = zero
      if (present(alpha)) alpha_local = alpha
!
      call mem%alloc(trial, davidson%n_parameters)
      call dcopy(davidson%n_parameters, R, 1, trial, 1)
!
      if (davidson%do_precondition) call davidson%preconditioner%do_(trial, shift=alpha_local)
!
!     Renormalize 
!
      norm_trial = get_l2_norm(trial, davidson%n_parameters)
      call dscal(davidson%n_parameters, one/norm_trial, trial, 1)
!
!     Add to trial space 
!
      davidson%n_new_trials = davidson%n_new_trials + 1
      call davidson%set_trial(trial, davidson%dim_red + davidson%n_new_trials)
!
      call mem%dealloc(trial, davidson%n_parameters)
!
   end subroutine construct_next_trial_linear_davidson_tool
!
!
   subroutine cleanup_linear_davidson_tool(davidson)
!!
!!    Finalize 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      if (allocated(davidson%A_red)) call mem%dealloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
      if (allocated(davidson%X_red)) call mem%dealloc(davidson%X_red, davidson%dim_red, davidson%n_solutions)
      if (allocated(davidson%F)) call mem%dealloc(davidson%F, davidson%n_parameters)
      if (allocated(davidson%F_red)) call mem%dealloc(davidson%F_red, davidson%dim_red)
!
   end subroutine cleanup_linear_davidson_tool
!
!
end module linear_davidson_tool_class
