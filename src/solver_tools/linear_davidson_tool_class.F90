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
!!          them and calling davidson%set_trial(c, n) for n = 1,2,... When you are 
!!          done writing trials, call davidson%orthonormalize_trial_vecs() 
!!          to make sure the trial vectors are orthonormalized.
!!
!!       2. (Optional) Set preconditioner by call davidson%set_preconditioner(P),
!!          where P^-1 is a diagonal preconditioner for A. For the CC Jacobian, P
!!          is the vector of orbital differences.
!!
!!       3. Set up the iterative loop, which consists of four stages:
!!          (Here shown for one linear equation.)
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
!!             call davidson%construct_residual(R, 1)
!!                
!!             ...
!!
!!             if (residual_norm >= thr) call davidson%construct_next_trial(R)
!!
!!          enddo ! end of iterative loop 
!!
!! 
!!    Modified by Eirik F. Kjønstad and Josefine H. Andersen, 2019. 
!!
!!    Generalization of tool routines to handle multiple frequencies 
!!    and multiple right-hand-sides. Adapted from routines written by 
!!    Josefine H. Andersen.
!!
!
   use davidson_tool_class
!
   use memory_manager_class, only: mem 
   use math_utilities, only: delta
!
   type, extends(davidson_tool) :: linear_davidson_tool 
!
      real(dp), dimension(:,:), allocatable :: F_red
      real(dp), dimension(:,:), pointer :: F 
!
      integer :: n_rhs
      real(dp), dimension(:), allocatable :: frequencies
!
   contains 
!
!     Linear Davidson specific routines
!
      procedure :: construct_next_trial               => construct_next_trial_linear_davidson_tool
      procedure :: construct_residual                 => construct_residual_linear_davidson_tool
      procedure :: solve_reduced_problem              => solve_reduced_problem_linear_davidson_tool
!
      procedure :: construct_reduced_gradient         => construct_reduced_gradient_linear_davidson_tool
!
      procedure, private :: general_preparations      => general_preparations_linear_davidson_tool
!
      procedure :: set_trials_to_preconditioner_guess => set_trials_to_preconditioner_guess_linear_davidson
!
      final :: destructor_linear_davidson_tool
!
   end type linear_davidson_tool
!
!
   interface linear_davidson_tool
!
      procedure :: new_linear_davidson_tool_one_rhs
      procedure :: new_linear_davidson_tool_multiple_rhs
!
   end interface linear_davidson_tool
!
!
contains
!
!
   function new_linear_davidson_tool_one_rhs(name_, n_parameters, lindep_threshold, &
            max_dim_red, F, n_equations, frequencies) result(davidson)
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
!!    F:                 If n_equations = 1, then F is the vector for the 
!!                       right-hand-side of the single linear equation 
!!
!!                         (A - freq) X = F.   (*)
!!
!!                       If n_equations > 1, then F is the vector for the 
!!                       right-hand-side of the set of equations
!!
!!                         (A - freq_k) X_k = F,    k = 1, 2, ..., n_frequencies.  (**)
!!    
!!    n_equations:       Number of equations / number of solutions. 
!!
!!    frequencies:       (Optional) Array of frequencies with length equal to n_equations. 
!!                       Default is an array of zeros. 
!!
      implicit none 
!
      type(linear_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters, max_dim_red, n_equations
!
      real(dp), dimension(n_equations), intent(in), optional :: frequencies
!
      real(dp), dimension(n_parameters), target, intent(in) :: F 
!
      real(dp), intent(in) :: lindep_threshold
!
!     Perform tasks common to constructor 
!
      call davidson%general_preparations(name_, n_parameters, n_equations, &
                           lindep_threshold, max_dim_red, 1)
!
!     Set F and frequencies 
!
      davidson%F(1:n_parameters, 1:davidson%n_rhs) => F(1:n_parameters)
!
      call mem%alloc(davidson%frequencies, davidson%n_solutions)
!
      if (present(frequencies)) then 
!
         davidson%frequencies = frequencies
!
      else 
!
         call zero_array(davidson%frequencies, davidson%n_solutions)
!
      endif
!
   end function new_linear_davidson_tool_one_rhs
!
!
   function new_linear_davidson_tool_multiple_rhs(name_, n_parameters, lindep_threshold, &
            max_dim_red, F, n_rhs, frequencies, n_frequencies) result(davidson)
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
!!    F:                 The columns of F contain the right-hand-side vectors of the set of equations
!!
!!                         (A - freq_k) X_k = F_k,   F_k = F(:,k),   k = 1, 2, ..., n_frequencies. (***)
!!                                                   
!!    n_equations:       Number of equations / number of solutions. 
!!
!!    frequencies:       Array of frequencies with length equal to n_frequencies.
!!
      implicit none 
!
      type(linear_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters, max_dim_red, n_rhs, n_frequencies
!
      real(dp), dimension(n_frequencies), intent(in) :: frequencies
!
      real(dp), dimension(n_parameters, n_rhs), target, intent(in) :: F 
!
      real(dp), intent(in) :: lindep_threshold
!
!     Perform tasks common to constructor 
!
      call davidson%general_preparations(name_, n_parameters, n_frequencies, &
                              lindep_threshold, max_dim_red, n_rhs)
!
!     Set F and frequencies 
!
      davidson%F(1:n_parameters, 1:n_rhs) => F
!
      call mem%alloc(davidson%frequencies, davidson%n_solutions)
!
      davidson%frequencies = frequencies
!
   end function new_linear_davidson_tool_multiple_rhs
!
!
   subroutine general_preparations_linear_davidson_tool(davidson, name_, &
                                             n_parameters, n_equations, lindep_threshold, &
                                             max_dim_red, n_rhs)
!!
!!    General preparations  
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Performs tasks that are in common for different constructors of the tool. 
!!    
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name_ 
!
      integer, intent(in) :: n_parameters, n_equations, max_dim_red, n_rhs
!
      real(dp), intent(in) :: lindep_threshold
!
!     Set tool parameters 
!
      davidson%name_ = trim(name_)
      davidson%n_parameters = n_parameters
      davidson%n_solutions = n_equations
      davidson%max_dim_red = max_dim_red
      davidson%lindep_threshold = lindep_threshold
      davidson%n_rhs = n_rhs
!
      davidson%do_precondition = .false. ! Switches to true if 'set_preconditioner' is called
!
!     Set some initial values 
!
      davidson%dim_red      = 0    
      davidson%n_new_trials = davidson%n_solutions 
!
   end subroutine general_preparations_linear_davidson_tool
!
!
   subroutine set_trials_to_preconditioner_guess_linear_davidson(davidson)
!!
!!    Set trials to precondtioner guess 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    If the preconditiner Q is set, this routine determines 
!!    initial guesses for the solution by assuming Q = A:
!!
!!       (A - freq) X = F => X = (Q - freq)^-1 F 
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      real(dp), dimension(:), allocatable :: c 
!
      integer :: i 
!
      call mem%alloc(c, davidson%n_parameters)
!
      do i = 1, davidson%n_solutions
!
         if (davidson%n_rhs == 1) then 
!
            call dcopy(davidson%n_parameters, davidson%F(:,1), 1, c, 1)
!
         else
!
            call dcopy(davidson%n_parameters, davidson%F(:,i), 1, c, 1)
!
         endif 
!
         call davidson%preconditioner%do_(c,                               &
                                          shift=davidson%frequencies(i),   &
                                          prefactor=-one)
!
         call davidson%set_trial(c, i)
!
      enddo
!
   end subroutine set_trials_to_preconditioner_guess_linear_davidson
!
!
   subroutine construct_reduced_gradient_linear_davidson_tool(davidson)
!!
!!    Construct reduced gradient
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!! 
!!    Constructs the gradient vector(s) in the reduced space:
!!
!!       F_ik = c_i^T * F_k 
!!
!!    Modified by Eirik F. Kjønstad and Josefine H. Andersen, 2019.
!!
!!    Generalization to multiple frequencies and multiple right-hand-sides. 
!!
      implicit none
!
      class(linear_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: c_i
!
      integer :: i, k
!
      real(dp) :: ddot
!
!     Construct reduced vector
!
      if (allocated(davidson%F_red)) call mem%dealloc(davidson%F_red, &
                                                      davidson%dim_red - davidson%n_new_trials, &
                                                      davidson%n_rhs) 
!
      call mem%alloc(davidson%F_red, davidson%dim_red, davidson%n_rhs)
!
      call mem%alloc(c_i, davidson%n_parameters)
!
      do i = 1, davidson%dim_red 
!
         call davidson%get_trial(c_i, i)
!
         do k = 1, davidson%n_rhs 
!
            davidson%F_red(i,k) = ddot(davidson%n_parameters, c_i, 1, davidson%F(:,k), 1)
!
         enddo
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
!!    Solves the reduced space equations:
!!
!!       (A - freq) X = F,  
!!
!!    where there can be many frequencies and right-hand-sides (see *, **, and *** 
!!    in the constructors).
!!
!!    In more detail, the routine performs three actions: 1) constructs the reduced space 
!!    quantities, 2) uses these to solve the reduced space problem, and 3) sets number of new trials 
!!    to zero (prepares the tool for the next task, which is the receiving of new residuals
!!    to make another new set of trial vectors).
!!
!!    Modified by Eirik F. Kjønstad and Josefine H. Andersen, 2019.
!!
!!    Generalization to multiple frequencies and multiple right-hand-sides. 
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      real(dp), dimension(:,:), allocatable :: M_k
!
      integer :: i, j, k
      integer :: info ! Error integer for dgesv routine (LU factorization)
!
!     Construct reduced space quantities 
!
      call davidson%construct_reduced_matrix()
      call davidson%construct_reduced_gradient()
!
!     (Re-)allocate the reduced space solution vector 
!
      if (allocated(davidson%X_red)) call mem%dealloc(davidson%X_red, &
                                                      davidson%dim_red - davidson%n_new_trials, &
                                                      davidson%n_solutions)
!
      call mem%alloc(davidson%X_red, davidson%dim_red, davidson%n_solutions)
!
!     Solve the eigenvalue equation(s)
!
      call mem%alloc(ipiv, davidson%dim_red)
      call mem%alloc(M_k, davidson%dim_red, davidson%dim_red)
!
      do k = 1, davidson%n_solutions
!
!        Set X_k equal to the right-hand-side for dgesv 
!
         if (davidson%n_rhs == 1) then ! X_k = F 
!
            davidson%X_red(:,k) = davidson%F_red(:,1)
!
         else ! X_k = F_k 
!
            davidson%X_red(:,k) = davidson%F_red(:,k)
!
         endif
!
!        Make M_k = A - freq_k * I 
!
         do i = 1, davidson%dim_red
            do j = 1, davidson%dim_red
!
               M_k(i,j) = davidson%A_red(i,j) - delta(i,j)*davidson%frequencies(k)
!
            enddo
         enddo
!
!        Solve M_k X_k = F, or 
!              M_k X_k = F_k         
!
         ipiv = 0
         info = 0
!
         call dgesv(davidson%dim_red,     &
                     1,                   & ! Number of RHS 
                     M_k,                 &
                     davidson%dim_red,    &
                     ipiv,                &
                     davidson%X_red(:,k), & ! Solution on exit 
                     davidson%dim_red,    &
                     info)
!
         if (info .ne. 0) then 
!
            call output%error_msg('could not solve linear equation in Davidson tool.')
!
         endif 
!
      enddo
!
      call mem%dealloc(ipiv, davidson%dim_red)
      call mem%dealloc(M_k, davidson%dim_red, davidson%dim_red)
!
!     Reset the number of new trials
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
!!    Constructs the nth residual. 
!!
!!    If one RHS, this is:
!!
!!       R = (A - freq_n)*X_n - F 
!!
!!    If multiple RHS, this is: 
!!   
!!       R = (A - freq_n)*X_n - F_n  
!!
!!    Modified by Eirik F. Kjønstad and Josefine H. Andersen, 2019.
!!
!!    Generalization to multiple frequencies and multiple right-hand-sides. 
!!
      implicit none 
!
      class(linear_davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(davidson%n_parameters) :: R 
!
      integer, intent(in) :: n
!
      real(dp), dimension(:), allocatable :: X 
!
      call mem%alloc(X, davidson%n_parameters)
!
      call davidson%construct_solution(X, n) ! X = X_n
      call davidson%construct_AX(R, n) ! R = A X_n 
!
      call daxpy(davidson%n_parameters, -davidson%frequencies(n), X, 1, R, 1) ! R = (A - freq_n) X_n
!
      call mem%dealloc(X, davidson%n_parameters)
!
      if (davidson%n_rhs == 1) then ! R = (A - freq_n) X_n - F
!
         call daxpy(davidson%n_parameters, -one, davidson%F, 1, R, 1) 
!
      else ! R = (A - freq_n) X_n - F_n
!
         call daxpy(davidson%n_parameters, -one, davidson%F(:,n), 1, R, 1) 
!
      endif 
!
   end subroutine construct_residual_linear_davidson_tool
!
!
   subroutine construct_next_trial_linear_davidson_tool(davidson, R, n)
!!
!!    Construct next trial  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    R:     nth residual, (A - freq) X - F
!!    n:     State number
!!
!!    Preconditions the residual and adds it to the trial space. 
!!    
      implicit none 
!
      class(linear_davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: R 
!
      integer, intent(in) :: n 
!
      real(dp) :: norm_trial
      real(dp), dimension(:), allocatable :: trial 
!
!     Precondition 
!
      call mem%alloc(trial, davidson%n_parameters)
      call dcopy(davidson%n_parameters, R, 1, trial, 1)
!
      if (davidson%do_precondition) call davidson%preconditioner%do_(trial, shift=davidson%frequencies(n))
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
   subroutine destructor_linear_davidson_tool(davidson)
!!
!!    Destructor  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      type(linear_davidson_tool) :: davidson 
!
      if (allocated(davidson%A_red)) call mem%dealloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
      if (allocated(davidson%X_red)) call mem%dealloc(davidson%X_red, davidson%dim_red, davidson%n_solutions)
      if (allocated(davidson%F_red)) call mem%dealloc(davidson%F_red, davidson%dim_red, davidson%n_rhs)
!
      davidson%F => null()
!
   end subroutine destructor_linear_davidson_tool
!
!
end module linear_davidson_tool_class
