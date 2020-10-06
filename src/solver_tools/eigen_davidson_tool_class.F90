!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module eigen_davidson_tool_class
!
!!
!!    Eigenvalue davidson tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    A tool to help solve an eigenvalue equation A X_n = omega_n X_n
!!    for general A using the Davidson algorithm. It is tailored to
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
!!             I. Preparations of reduced space for iteration 
!!
!!             if (davidson%red_dim_exceeds_max()) call davidson%set_trials_to_solutions()
!!
!!             call davidson%update_reduced_dim()
!!
!!             call davidson%orthonormalize_trial_vecs() 
!!
!!             II. Read new trials, transform them, & store the result 
!!
!!             do trial = davidson%first_new_trial(), davidson%last_new_trial()
!!
!!                call davidson%get_trial(c, trial)
!!                Transform: c <- A c 
!!                call davidson%set_transform(c)
!!
!!             enddo 
!!
!!             III. Solve reduced problem 
!!             
!!             call davidson%solve_reduced_problem()
!!
!!             IV. Construct residual and use it to generate new trial vectors 
!!
!!             do state = 1, n_states 
!!
!!                call davidson%construct_residual(R, state)
!!                
!!                ...
!!
!!                if (residual_norm >= thr) call davidson%construct_next_trial(R, state)
!!
!!             enddo 
!!
!!          enddo ! end of iterative loop 
!!             
!
   use davidson_tool_class
!
   use memory_manager_class, only: mem 
   use array_utilities, only: get_n_lowest
!
   type, extends(davidson_tool) :: eigen_davidson_tool 
!
      real(dp), dimension(:), allocatable :: omega_re 
      real(dp), dimension(:), allocatable :: omega_im
!
      logical :: non_unit_metric                     ! If trial vectors are not orthonormal,
                                                     ! then the metric is not the unit matrix 
                                                     ! Default is false, but can be set to true
                                                     ! in constructor
!
      real(dp), dimension(:,:), allocatable :: S_red ! Metric if non-orthonormal trial vectors 
!
   contains 
!
!     Procedures a user of the tool may need to use (see ancestor also)
!
      procedure :: solve_reduced_problem          => solve_reduced_problem_eigen_davidson_tool
      procedure :: construct_residual             => construct_residual_eigen_davidson_tool
      procedure :: construct_next_trial           => construct_next_trial_eigen_davidson_tool
!
!     Other routines 
!
      procedure, private :: construct_re_residual => construct_re_residual_eigen_davidson_tool
      procedure, private :: construct_im_residual => construct_im_residual_eigen_davidson_tool
!  
      procedure, private :: initialize_omega_re   => initialize_omega_re_eigen_davidson_tool
      procedure, private :: initialize_omega_im   => initialize_omega_im_eigen_davidson_tool
!  
      procedure, private :: destruct_omega_re     => destruct_omega_re_eigen_davidson_tool
      procedure, private :: destruct_omega_im     => destruct_omega_im_eigen_davidson_tool
!  
      procedure :: construct_reduced_metric           &
                => construct_reduced_metric_eigen_davidson_tool
!
      procedure :: construct_reduced_submetric        &
                => construct_reduced_submetric_eigen_davidson_tool
!
      procedure :: destruct_reduced_space_quantities  &
                => destruct_reduced_space_quantities_eigen_davidson_tool
!
      final :: destructor_eigen_davidson_tool 
!
   end type eigen_davidson_tool
!
!
   interface eigen_davidson_tool
!
      procedure :: new_eigen_davidson_tool
!
   end interface eigen_davidson_tool
!
!
contains
!
!
   function new_eigen_davidson_tool(name_, n_parameters, n_solutions, &
                        lindep_threshold, max_dim_red, non_unit_metric) result(davidson)
!!
!!    New eigen Davidson tool 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    name_ :            Name of solver tool (used for temporary files)
!!
!!    n_parameters:      Dimensionality of the full vector space (A is n_parameters x n_parameters)
!!
!!    n_solutions:       Number of solutions of the equation to solve for 
!!
!!    lindep_threshold:  Norm threshold for new trial vector after being added to the trial 
!!                       space. If the vector, after orthonormalization against the existing
!!                       trial space, has a norm below this threshold, then we assume that the 
!!                       trial basis has become linearly dependent. This causes an error stop. 
!!
!!    max_dim_red:       Maximum dimension of the reduced space. When exceeding this dimensionality,
!!                       the solutions are set as the basis for the new trial space. 
!!
!!    non_unit_metric:   (optional) Specifies that the metric can be non-unit (i.e., the trial 
!!                       can be non-orthonormal). Default: false. Use only when there is a reason 
!!                       to. Assumes that vectors are linearly independent. 
!!
      implicit none 
!
      type(eigen_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in)  :: n_parameters, n_solutions, max_dim_red 
      real(dp), intent(in) :: lindep_threshold  
!
      logical, optional, intent(in) :: non_unit_metric
!
      davidson%n_parameters = n_parameters
      davidson%n_solutions  = n_solutions
      davidson%max_dim_red  = max_dim_red
!
      davidson%lindep_threshold = lindep_threshold  
!
      davidson%name_ = trim(name_)
!
      davidson%non_unit_metric = .false.
      if (present(non_unit_metric)) davidson%non_unit_metric = non_unit_metric
!
      davidson%do_precondition   = .false. ! Switches to true if 'set_preconditioner' is called
!
      davidson%dim_red      = 0
      davidson%n_new_trials = n_solutions
!
   end function new_eigen_davidson_tool
!
!
   subroutine destructor_eigen_davidson_tool(davidson)
!!
!!    Destructor  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      type(eigen_davidson_tool), intent(inout) :: davidson 
!
      call davidson%destruct_reduced_space_quantities()
!
      call davidson%destruct_omega_re()
      call davidson%destruct_omega_im()
!
      if (davidson%do_precondition) call davidson%preconditioner%destruct_precondition_vector()
!
   end subroutine destructor_eigen_davidson_tool
!  
!
   subroutine initialize_omega_im_eigen_davidson_tool(davidson)
!!
!!    Initialize omega im 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      if (.not. allocated(davidson%omega_im)) call mem%alloc(davidson%omega_im, davidson%n_solutions)
!
   end subroutine initialize_omega_im_eigen_davidson_tool
!  
!
   subroutine initialize_omega_re_eigen_davidson_tool(davidson)
!!
!!    Initialize omega re 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      if (.not. allocated(davidson%omega_re)) call mem%alloc(davidson%omega_re, davidson%n_solutions)
!
   end subroutine initialize_omega_re_eigen_davidson_tool
!  
!
   subroutine destruct_omega_im_eigen_davidson_tool(davidson)
!!
!!    Destruct omega im 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      if (allocated(davidson%omega_im)) call mem%dealloc(davidson%omega_im, davidson%n_solutions)
!
   end subroutine destruct_omega_im_eigen_davidson_tool
!  
!
   subroutine destruct_omega_re_eigen_davidson_tool(davidson)
!!
!!    Destruct omega re 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      if (allocated(davidson%omega_re)) call mem%dealloc(davidson%omega_re, davidson%n_solutions)
!
   end subroutine destruct_omega_re_eigen_davidson_tool
!
!
   subroutine solve_reduced_problem_eigen_davidson_tool(davidson)
!!
!!    Solve reduced problem 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    For the current reduced matrix A, this routine solves the eigenvalue problem 
!!
!!       A X_n = omega_n X_n,    n = 1, 2, ..., dim_red.
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
      use array_utilities, only: invert
!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), dimension(:), allocatable :: work
!
      real(dp), dimension(:), allocatable :: omega_re
      real(dp), dimension(:), allocatable :: omega_im
!
      integer, dimension(:), allocatable :: index_list
!
      real(dp), dimension(:,:), allocatable :: X_red
      real(dp), dimension(:,:), allocatable :: A_red ! Safe copy to avoid BLAS overwrite
      real(dp), dimension(:,:), allocatable :: S_red_inv 
!
      integer :: info = 0, j = 0, i = 0, worksize
      real(dp) :: dummy = 0.0
      real(dp), dimension(1) :: optwork
!
!     Construct reduced space matrix A 
!
      call davidson%construct_reduced_matrix()
!
!     Solve reduced eigenvalue problem
!
      info = 0
!
      call mem%alloc(X_red, davidson%dim_red, davidson%dim_red)
      call mem%alloc(A_red, davidson%dim_red, davidson%dim_red)
!
      call zero_array(X_red, davidson%dim_red**2)
      call dcopy(davidson%dim_red**2, davidson%A_red, 1, A_red, 1)
!
      if (davidson%non_unit_metric) then 
!
!        Construct reduced space metric S and invert it 
!
         call davidson%construct_reduced_metric()
!
         call mem%alloc(S_red_inv, davidson%dim_red, davidson%dim_red)
         call invert(S_red_inv, davidson%S_red, davidson%dim_red)
!
!        To solve A X = w S X, we transform the system to (S^-1 A) X = w X
!        by redefining the A matrix as S^-1 A.  
!
         call dgemm('N', 'N',             &
                     davidson%dim_red,    &
                     davidson%dim_red,    &
                     davidson%dim_red,    &
                     one,                 &
                     S_red_inv,           &
                     davidson%dim_red,    &
                     davidson%A_red,      &
                     davidson%dim_red,    &
                     zero,                &
                     A_red,               &
                     davidson%dim_red)
!
         call mem%dealloc(S_red_inv, davidson%dim_red, davidson%dim_red)
!
      endif 
!
      call mem%alloc(omega_re, davidson%dim_red)
      call mem%alloc(omega_im, davidson%dim_red)
!
!     Find optimal work size, lwork = -1
!
      call dgeev('N','V',             &
                  davidson%dim_red,   &
                  A_red,              &
                  davidson%dim_red,   &
                  omega_re,           &
                  omega_im,           &
                  dummy,              &              
                  1,                  &
                  X_red,              &
                  davidson%dim_red,   &
                  optwork,            &   
                  -1,                 &
                  info)
!
      worksize = ceiling( optwork(1) )
!
      call mem%alloc(work, worksize)
!      
      call dgeev('N','V',             &
                  davidson%dim_red,   &
                  A_red,              &
                  davidson%dim_red,   &
                  omega_re,           &
                  omega_im,           &
                  dummy,              &              
                  1,                  &
                  X_red,              &
                  davidson%dim_red,   &
                  work,               &   
                  worksize,           &
                  info)
!
      call mem%dealloc(A_red, davidson%dim_red, davidson%dim_red)
!
      if (info .ne. 0) call output%error_msg('could not solve reduced equation in Davidson davidson.')
!
      call mem%dealloc(work, worksize)
!
!     Find lowest n_solutions eigenvalues and sort them (the corresponding indices
!     are placed in the integer array index_list)
!
      call davidson%initialize_omega_im()
      call davidson%initialize_omega_re()
!
      call mem%alloc(index_list, davidson%n_solutions)
!
      call get_n_lowest(davidson%n_solutions, davidson%dim_red, &
                        omega_re, davidson%omega_re, index_list)
!
!     Pick out solution vectors and imaginary parts of eigenvalues according to index_list
!
      if (allocated(davidson%X_red)) &
         call mem%dealloc(davidson%X_red, davidson%dim_red - davidson%n_new_trials, davidson%n_solutions)
!
      call mem%alloc(davidson%X_red, davidson%dim_red, davidson%n_solutions)
!
!$omp parallel do private(i, j) collapse(2)
      do j = 1, davidson%n_solutions
         do i = 1, davidson%dim_red
!
            davidson%X_red(i, j) = X_red(i, index_list(j))
            davidson%omega_im(j) = omega_im(index_list(j))
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_red, davidson%dim_red, davidson%dim_red)
      call mem%dealloc(omega_re, davidson%dim_red)
      call mem%dealloc(omega_im, davidson%dim_red)
      call mem%dealloc(index_list, davidson%n_solutions)
!
      davidson%n_new_trials = 0
!
   end subroutine solve_reduced_problem_eigen_davidson_tool
!
!
   subroutine construct_residual_eigen_davidson_tool(davidson, R, X, n)
!!
!!    Construct residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Constructs the nth full space residual 
!!
!!       R = A X - omega X,
!!
!!    without any preconditioning (this is applied after a call
!!    to construct residual from the construct next trial vector
!!    routine). Here, X is the nth eigenvector and omega is the 
!!    eigenvalue. The residual is normalized by the norm of the 
!!    solution X.
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(davidson%n_parameters), intent(out)  :: R 
      real(dp), dimension(davidson%n_parameters), intent(in)   :: X 
!
      integer, intent(in), optional :: n
!
      integer :: k
!
      real(dp) :: norm_X 
!
      if (present(n)) then
! 
         k = n
!
      else
!
         k = 1
!
      endif   
!
      norm_X = get_l2_norm(X, davidson%n_parameters)
!
      if (davidson%omega_im(k) == zero) then  ! standard case: the nth root is not part of a complex pair
!
         call davidson%construct_AX(R, k) ! set R = AX 
!
        call daxpy(davidson%n_parameters, - davidson%omega_re(k), X, 1,  R, 1)
        call dscal(davidson%n_parameters, one/norm_X, R, 1)
!
      else
!
!        If it's the first root of the complex pair, construct the real residual; 
!        if it's the second, construct the imaginary residual
!
         if (k == 1) then
!
            if (davidson%n_solutions < 2) call output%error_msg('add one more root to treat the complex pair.')
!
            call davidson%construct_re_residual(R, X, norm_X, k)
!
         elseif (n == davidson%n_solutions) then 
!
            if (davidson%omega_re(k-1) /= davidson%omega_re(k)) &
                        call output%error_msg('add one more root to treat the complex pair.')
!
            call davidson%construct_im_residual(R, X, norm_X, k)
!
         else ! neither first or last, so it's safe to look at n + 1 and n - 1 
!
            if (davidson%omega_re(k) == davidson%omega_re(k - 1)) then
!
               call davidson%construct_im_residual(R, X, norm_X, k)
!
            elseif (davidson%omega_re(k) == davidson%omega_re(k + 1)) then 
!
               call davidson%construct_re_residual(R, X, norm_X, k)
!
            else ! should never happen, but just in case, let the user know there's a bug
!
               call output%error_msg('something went very wrong when trying to construct imaginary residual.')
!
            endif 
!
         endif
!
      endif
!
   end subroutine construct_residual_eigen_davidson_tool
!
!
   subroutine construct_re_residual_eigen_davidson_tool(davidson, R, X_re, norm_X_re, n)
!!
!!    Construct real residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    For a complex pair, 
!!
!!       X+ = X_re + i X_im 
!!       X- = X_re - i X_im,
!!
!!    we will find the real part X_re as the first of the two roots,
!!    and the imaginary part X_im as the second of the two roots 
!!    (in X_red). This routine constructs the residual associated
!!    with the real part. As such, the construct residual routine 
!!    makes sure the n-th root is X_re and (n + 1)-th is X_im (and
!!    exists, which might not be the case if too few roots are 
!!    requested). On exit, the real residual is placed in R: 
!!    
!!       R = (A X_re - omega_re X_re) + omega_im X_im 
!!
!!    The residual is divided by the norm of X+ (or, 
!!    equivalently, X-).
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(davidson%n_parameters)             :: R 
      real(dp), dimension(davidson%n_parameters), intent(in) :: X_re  
!
      integer, intent(in)  :: n 
      real(dp), intent(in) :: norm_X_re
!
      real(dp), dimension(:), allocatable :: X_im 
!
      real(dp) :: norm_X_im
!
      call mem%alloc(X_im, davidson%n_parameters)  
!
      call davidson%construct_solution(X_im, n + 1)   ! set X_im
      call davidson%construct_AX(R, n)                ! set R = A X_re 
!
      call daxpy(davidson%n_parameters, - davidson%omega_re(n), X_re, 1, R, 1) 
      call daxpy(davidson%n_parameters, davidson%omega_im(n), X_im, 1, R, 1) 
!
      norm_X_im = get_l2_norm(X_im, davidson%n_parameters)
!
      call dscal(davidson%n_parameters, one/(sqrt(norm_X_re**2 + norm_X_im**2)), R, 1)
!
      call mem%dealloc(X_im, davidson%n_parameters) 
!
   end subroutine construct_re_residual_eigen_davidson_tool
!
!
   subroutine construct_im_residual_eigen_davidson_tool(davidson, R, X_im, norm_X_im, n)
!!
!!    Construct imaginary residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    For a complex pair, 
!!
!!       X+ = X_re + i X_im 
!!       X- = X_re - i X_im,
!!
!!    we will find the real part X_re as the first of the two roots,
!!    and the imaginary part X_im as the second of the two roots 
!!    (in X_red). This routine constructs the residual associated
!!    with the imaginary part. As such, the construct residual 
!!    routine makes sure the n-th root is X_im and (n - 1)-th 
!!    is X_re. On exit, the imaginary residual is placed in R: 
!!    
!!       R = (A X_im - omega_re X_im) - omega_im X_re 
!!
!!    The residual is divided by the norm of X+ (or, equivalently, 
!!    X-). Note that omega_im of root n is the negative of omega_im 
!!    of root n - 1. To avoid inconsistency, we use omega_im from 
!!    the previous root (otherwise, we would change the sign of the
!!    last constribution).
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(davidson%n_parameters)             :: R 
      real(dp), dimension(davidson%n_parameters), intent(in) :: X_im  
!
      integer, intent(in)  :: n 
      real(dp), intent(in) :: norm_X_im
!
      real(dp), dimension(:), allocatable :: X_re  
!
      real(dp) :: norm_X_re
!
      call mem%alloc(X_re, davidson%n_parameters)  
!
      call davidson%construct_solution(X_re, n - 1) ! set X_re 
      call davidson%construct_AX(R, n)       ! set R = A X_im 
!
      call daxpy(davidson%n_parameters, - davidson%omega_re(n), X_im, 1, R, 1) 
      call daxpy(davidson%n_parameters, - davidson%omega_im(n - 1), X_re, 1, R, 1) 
!
      norm_X_re = get_l2_norm(X_re, davidson%n_parameters)
!
      call dscal(davidson%n_parameters, one/(sqrt(norm_X_re**2 + norm_X_im**2)), R, 1)
!
      call mem%dealloc(X_re, davidson%n_parameters)  
!
   end subroutine construct_im_residual_eigen_davidson_tool
!
!
   subroutine construct_next_trial_eigen_davidson_tool(davidson, R, n)
!!
!!    Construct next trial  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    R: nth residual, AX - omega*X 
!!    n: state number 
!!
!!    Preconditions the residual and adds it to the trial space.
!!    
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: R
!
      integer, intent(in) :: n 
!
      real(dp) :: norm_trial
!
      real(dp), dimension(:), allocatable :: trial 
!
!     Precondition 
!
      call mem%alloc(trial, davidson%n_parameters)
      call dcopy(davidson%n_parameters, R, 1, trial, 1)
!
      if (davidson%do_precondition) &
         call davidson%preconditioner%do_(trial,                        &
                                          shift=davidson%omega_re(n),   &
                                          prefactor=-one)
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
   end subroutine construct_next_trial_eigen_davidson_tool
!
!
   subroutine construct_reduced_metric_eigen_davidson_tool(davidson)
!!
!!    Construct reduced metric 
!!    Written by Eirik F. Kjønstad, Nov 2019 and Mar 2020
!!
!!    Constructs 
!!
!!       S_ij = c_i^T c_j     (davidson%S_red)
!!
!!    where {c_k} is the set of current trials. To avoid re-calculating 
!!    contributions, the routine computes the new elements in S in every 
!!    iteration and pads the existing S_red matrix.
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), dimension(:,:), allocatable :: S_red_copy
!
      integer :: req_0, req_1_i, req_1_j, req_2
!
      type(batching_index), allocatable :: batch_i, batch_j 
!
      integer :: prev_dim_red
!
      type(timings), allocatable :: timer 
!
      timer = timings('Davidson: time to construct reduced metric', 'v')
      call timer%turn_on()
!
!     Will batch over reduced space indices, i and j for trials
!
      req_0   = 0
      req_1_i = davidson%trials%required_to_load_record()
      req_1_j = req_1_i
      req_2   = 0
!
      if (davidson%dim_red .eq. davidson%n_solutions) then 
!
!        First iteration or reset of space: calculate all elements
!
         call mem%alloc(davidson%S_red, davidson%dim_red, davidson%dim_red)
!
         batch_i = batching_index(davidson%dim_red)
         batch_j = batching_index(davidson%dim_red)
!
         call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2)
!
         call davidson%construct_reduced_submetric(batch_i, batch_j)
!
      else
!
!        Not first iteration: calculate new elements only
!
!        Transfer previously calculated elements and 
!        extend the dimensionality of the reduced matrix 
!
         prev_dim_red = davidson%dim_red - davidson%n_new_trials
!
         call mem%alloc(S_red_copy, prev_dim_red, prev_dim_red)
!
         call dcopy(prev_dim_red**2, davidson%S_red, 1, S_red_copy, 1)
!
         call mem%dealloc(davidson%S_red, prev_dim_red, prev_dim_red)
!
         call mem%alloc(davidson%S_red, davidson%dim_red, davidson%dim_red)
!
         davidson%S_red(1:prev_dim_red, 1:prev_dim_red) = S_red_copy
!
         call mem%dealloc(S_red_copy, prev_dim_red, prev_dim_red)
!
!        Compute the new blocks in the reduced matrix 
!
!        i index new, j index all 
!
         batch_i = batching_index(dimension_=davidson%n_new_trials, &
                                  offset=prev_dim_red)
!
         batch_j = batching_index(dimension_=davidson%dim_red, &
                                  offset=0)
!
         call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2)
!
         call davidson%construct_reduced_submetric(batch_i, batch_j)
!
      endif
!
      call timer%turn_off()
!
   end subroutine construct_reduced_metric_eigen_davidson_tool
!
!
   subroutine construct_reduced_submetric_eigen_davidson_tool(davidson, batch_i, batch_j)
!!
!!    Construct reduced submetric 
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Constructs parts of the reduced metric S_ij based on the prepared 
!!    batching indices batch_i and batch_j.
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(inout) :: davidson 
!
      type(batching_index), intent(inout) :: batch_i, batch_j
!
      integer :: current_i_batch, current_j_batch, i, j
!
      real(dp), dimension(:,:), pointer, contiguous :: c_i, c_j 
!
      type(interval), allocatable :: j_interval
!
      real(dp) :: ddot
!
      call davidson%trials%prepare_records([batch_i, batch_j])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call davidson%trials%load(c_i, batch_i, 1)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            if (batch_j%first .gt. batch_i%last) cycle ! Nothing to calculate; 
                                                       ! go to next batch of j 
!
            j_interval = interval(first=batch_j%first, &
                                  last=min(batch_i%last, batch_j%last))
!
            call davidson%trials%load(c_j, j_interval, 2)
!
            do i = batch_i%first, batch_i%last
               do j = j_interval%first, min(batch_j%last, i) 
!
                  davidson%S_red(i, j) = ddot(davidson%n_parameters,             &
                                              c_i(:, i - batch_i%first + 1),     &
                                              1,                                 &
                                              c_j(:, j - j_interval%first + 1),  &
                                              1)
!
                  davidson%S_red(j, i) = davidson%S_red(i, j)
!
               enddo
            enddo
!
         enddo ! end of j batches 
      enddo ! end of i batches 
!
      call davidson%trials%free_records()
!
   end subroutine construct_reduced_submetric_eigen_davidson_tool
!
!
   subroutine destruct_reduced_space_quantities_eigen_davidson_tool(davidson)
!!
!!    Destruct reduced space quantities 
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
!!    Deallocates reduced space quantities, e.g. when re-setting the reduced space,
!!    or upon destruction of the Davidson tool.
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson
!
      if (allocated(davidson%A_red)) call mem%dealloc(davidson%A_red, &
         davidson%dim_red, davidson%dim_red)
      if (allocated(davidson%X_red)) call mem%dealloc(davidson%X_red, &
         davidson%dim_red, davidson%n_solutions)
      if (allocated(davidson%S_red)) &
         call mem%dealloc(davidson%S_red, davidson%dim_red, davidson%dim_red)
!
   end subroutine destruct_reduced_space_quantities_eigen_davidson_tool
!
!
end module eigen_davidson_tool_class
