!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module nonlinear_davidson_cc_es_class
!
!!
!! Non-linear Davidson coupled cluster excited state solver class module
!! Written by Eirik F. Kjønstad, 2019-2020
!!
!! Based in part on the original Davidson solver written by Sarai D. Folkestad
!! and Eirik F. Kjønstad, 2018-2019.
!!
!! The algorithm implemented here is described by C. Hättig and F. Weigend,
!! J. Chem. Phys. 113, 5154 (2000). It can be used on its own or for preconvergence
!! to produce sufficiently good starting guesses for the standard DIIS solver.
!!
!! Solves the CC excited state eigenvalue equation
!!
!!    A R = omega R   or  L^T A = omega L^T
!!
!! for a set of right or left states (R, L) and the associated
!! excitation energies omega. Here, A is the coupled cluster
!! Jacobian matrix:
!!
!!       A_mu,nu = < mu |[H-bar, tau_nu] | HF >,    H-bar = e-T H eT.
!!
!! Since A depends on the excitation energy of the state in the CCn
!! models (if equations are folded), this solver consists of a set of
!! micro iterations in which the Davidson algorithm is applied for
!! frozen energies. In each macro iteration, the energies are updated
!! with the solution obtained via the Davidson algorithm in the micro
!! iterations.
!!
!! See E. R. Davidson, J. Comput. Phys. 17, 87 (1975) for more details
!! regarding the Davidson algorithm.
!!
!! Davidson's algorithm expands a subspace using the residuals from
!! previous (micro) iterations.*
!!
!! * More precisely, using preconditioned residuals based on the
!!   orbital differences approximation of A,
!!
!!    A_mu,nu = < mu |[F, tau_nu] | HF > = epsilon_mu delta_mu,nu
!!
!!   Here,
!!
!!               epsilon_ai = orbital_energy(a) - orbital_energy(i)
!!
!!             epsilon_aibj = orbital_energy(a) + orbital_energy(b)
!!                          - orbital_energy(i) - orbital_energy(j),
!!
!!   and so on.
!!
!
   use parameters
   use ccs_class, only: ccs
   use timings_class, only: timings
   use global_out, only: output
   use memory_manager_class, only: mem
   use eigen_davidson_tool_class, only: eigen_davidson_tool
   use abstract_cc_es_class, only: abstract_cc_es!
!
   use start_vector_tool_class, only: start_vector_tool
   use abstract_projection_tool_class, only : abstract_projection_tool
   use convergence_tool_class, only: convergence_tool
!
   implicit none
!
!  Definition of solver class
!
   type, extends(abstract_cc_es) :: nonlinear_davidson_cc_es
!
!     Settings for micro iteration (Davidson subspace)
!
      integer  :: max_dim_red                         ! Maximum dimension of Davidson subspace
      integer  :: max_micro_iterations
      real(dp) :: relative_micro_residual_threshold   ! Micro iteration thresholds is this number
                                                      ! times the highest residual norm
!
!     Other settings
!
      logical :: prepare_wf                           ! In Davidson preconvergence, the wavefunction
                                                      ! is already prepared.
!
      class(eigen_davidson_tool), allocatable :: davidson ! tool to handle reduced space
!
   contains
!
      procedure :: run => run_nonlinear_davidson_cc_es
!
      procedure, private :: read_settings
      procedure, private :: print_settings
      procedure, private :: do_micro_iterations
!
   end type nonlinear_davidson_cc_es
!
!
   interface nonlinear_davidson_cc_es
!
      procedure :: new_nonlinear_davidson_cc_es
      procedure :: new_nonlinear_davidson_cc_es_preconvergence
!
   end interface nonlinear_davidson_cc_es
!
!
contains
!
!
   function new_nonlinear_davidson_cc_es(transformation, wf, start_vector, &
                                         projector, convergence_checker, n_states, &
                                         max_iterations, records_in_memory) result(this)
!!
!!    New non-linear Davidson CC ES
!!    Written by Eirik F. Kjønstad, Nov 2019 and Jan 2020
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      type(nonlinear_davidson_cc_es) :: this
      class(ccs), intent(inout), target :: wf
!
      character(len=*), intent(in) :: transformation
!
      class(start_vector_tool), intent(in) :: start_vector
      class(abstract_projection_tool), intent(in) :: projector
      type(convergence_tool), intent(in) :: convergence_checker
!
      integer, intent(in) :: n_states, max_iterations
      logical, intent(in) :: records_in_memory
!
      this%wf => wf
!
      this%total_timer = timings('Non-linear Davidson ES solver (' // &
                                 trim(transformation) //')', pl='normal')
      call this%total_timer%turn_on()
!
      this%name_ = 'Non-linear Davidson coupled cluster excited state solver'
      this%tag   = 'Davidson'
!
      this%description1 = 'A non-linear Davidson solver that calculates the lowest&
               & eigenvalues and the right or left eigenvectors of the Jacobian matrix A.&
               & The eigenvalue problem is solved in a reduced space for fixed energies,&
               & before the energies are updated and a new reduced space built.'
!
      this%description2 = 'See C. Hättig and F. Weigend, J. Chem. Phys. 113, 5154 (2000).'
!
      call this%print_banner()
!
      this%transformation   = trim(transformation)
      this%n_singlet_states = n_states
      this%max_iterations   = max_iterations
!
      this%projector = projector
      this%start_vectors = start_vector
      this%convergence_checker = convergence_checker
!
      this%max_dim_red                        = 0
      this%max_micro_iterations               = 100
      this%relative_micro_residual_threshold  = 1.0d-1
      this%prepare_wf                         = .true.
!
!
      call this%read_settings()
!
      call this%print_settings()
!
      if (this%n_singlet_states == 0) &
         call output%error_msg('number of excitations must be specified.')
!
      this%wf%n_singlet_states = this%n_singlet_states
!
      this%davidson = eigen_davidson_tool(name_ = 'nonlin_cc_es_davidson',        &
                                          n_parameters = this%wf%n_es_amplitudes, &
                                          n_solutions = this%n_singlet_states,    &
                                          lindep_threshold = 1.0d-11,             &
                                          max_dim_red = this%max_dim_red,         &
                                          records_in_memory = records_in_memory,  &
                                          non_unit_metric = .true.)
!
   end function new_nonlinear_davidson_cc_es
!
!
   function new_nonlinear_davidson_cc_es_preconvergence(transformation, wf,                  &
                                                        start_vector,                        &
                                                        projector,                           &
                                                        convergence_checker,                 &
                                                        n_singlet_states,                    &
                                                        max_iterations,                      &
                                                        records_in_memory,                   &
                                                        relative_micro_residual_threshold,   &
                                                        max_micro_iterations,                &
                                                        max_dim_red,                         &
                                                        prepare_wf) result(this)
!!
!!    New non-linear Davidson for DIIS preconvergence
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
!!    Constructor that can be called from another solver to perform e.g.
!!    preconvergence using non-linear Davidson.
!!
!!    Assumes that the wavefunction has been prepared for the excited state calculation
!!    (see "prepare_wf_for_excited_state" routine). Several tasks are currently re-done,
!!    but these are all cheap operations and should not affect performance.
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      type(nonlinear_davidson_cc_es) :: this
!
      class(ccs), target :: wf
!
      class(start_vector_tool), intent(in) :: start_vector
      class(abstract_projection_tool), intent(in) :: projector
      type(convergence_tool), intent(in) :: convergence_checker
!
      real(dp), intent(in) :: relative_micro_residual_threshold
      integer, intent(in) :: max_iterations, max_dim_red, n_singlet_states
      integer, intent(in) :: max_micro_iterations
      logical, intent(in) :: prepare_wf, records_in_memory
!
      character(len=*), intent(in) :: transformation
!
      this%wf => wf
!
      this%total_timer = timings('Non-linear Davidson ES solver (' // &
                                 trim(transformation) //')', pl='normal')
      call this%total_timer%turn_on()
!
      this%name_ = 'Non-linear Davidson coupled cluster excited state solver'
      this%tag   = 'Davidson'
!
      this%description1 = 'A Davidson solver that calculates the lowest eigenvalues and &
               & the right or left eigenvectors of the Jacobian matrix, A. The eigenvalue &
               & problem is solved in a reduced space, the dimension of which is &
               & expanded until the convergence criteria are met.'
!
      this%description2 = 'A complete description of the algorithm can be found in &
                                 & E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call this%print_banner()
!
      this%transformation   = trim(transformation)
      this%n_singlet_states = n_singlet_states
      this%max_iterations   = max_iterations
!
      this%start_vectors = start_vector
      this%projector = projector
      this%convergence_checker = convergence_checker
!
      this%max_micro_iterations               = max_micro_iterations
      this%relative_micro_residual_threshold  = relative_micro_residual_threshold
      this%max_dim_red                        = max_dim_red
      this%prepare_wf                         = prepare_wf
!
      call this%print_settings()
!
      this%davidson = eigen_davidson_tool(name_ = 'nonlin_cc_es_davidson',        &
                                          n_parameters = this%wf%n_es_amplitudes, &
                                          n_solutions = this%n_singlet_states,    &
                                          lindep_threshold = 1.0d-11,             &
                                          max_dim_red = this%max_dim_red,         &
                                          records_in_memory = records_in_memory,  &
                                          non_unit_metric = .true.)
!
   end function new_nonlinear_davidson_cc_es_preconvergence
!
!
   subroutine print_settings(this)
!!
!!    Print settings
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
      implicit none
!
      class(nonlinear_davidson_cc_es) :: this
!
      call this%print_es_settings()
!
      call output%printf('m', 'Max micro iterations:                (i4)',  &
                         ints=[this%max_micro_iterations], fs='(t6,a)')
!
      call output%printf('m', 'Relative micro threshold:     (e11.2)', &
                         reals=[this%relative_micro_residual_threshold], fs='(t6,a/)')
!
   end subroutine print_settings
!
!
   subroutine run_nonlinear_davidson_cc_es(this)
!!
!!    Run
!!    Written by Eirik F. Kjønstad, Nov 2019 and Jan 2020
!!
!!    Can be called, after constructing the solver object, to converge the
!!    excited state equations for the right or left coupled cluster states.
!!
!!    For a description of the algorithm employed, see the documentation in
!!    the top of the module. A more detailed description of the micro iterations
!!    is given in the associated routine - see "do_micro_iterations".
!!
!!    This driver routine simply calls the micro iterations routine so long as
!!    the states have not yet converged. Macro iterations consist in 1) updating
!!    the energies using the new guesses for the solutions and 2) checking whether
!!    they satisfy convergence criteria. If not, micro iterations are called again.
!!
!!    The solver is most suited for folded non-linear coupled cluster models, such
!!    as low-memory CC2 and CC3, for which the standard Davidson algorithm cannot be used.
!!    The solver is used in combination with DIIS as a preconvergence step to increase the
!!    robustness of DIIS. See "do_davidson_preconvergence" in the DIIS excited state solver.
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(nonlinear_davidson_cc_es), intent(inout) :: this
!
!     Convergence logical arrays, one component for each state
!
      logical, dimension(:), allocatable :: converged
!
!     Previous energies, current residuals, one component for each state
!
      real(dp), dimension(:), allocatable :: prev_energies
      real(dp), dimension(:), allocatable :: residual_norms
!
      integer :: iteration, state
!
      integer :: n_micro_iterations, n_transformations ! Information regarding micro iterations
!
      real(dp) :: ddot
!
!     Current guess of solutions
!     (columns stores the different states, i.e. X(:,k) is the guess for X_k)
!
      real(dp), dimension(:,:), allocatable :: X
!
!     Associated residuals, R(:,k) = A(omega_k) X(:,k) - omega_k X(:,k)
!     The current guess for the energy omega_k is saved in this%energies(k).
!
      real(dp), dimension(:,:), allocatable :: R
!
      type(timings) :: macro_iteration_time
      character(len=40) :: timer_name
!
      if (this%prepare_wf) call this%prepare_wf_for_excited_state()
!
!     Initialize energies, residual norms, and convergence arrays
!
      call this%initialize_energies()
      this%energies = zero
!
      call mem%alloc(prev_energies, this%n_singlet_states, set_zero=.true.)
      call mem%alloc(residual_norms, this%n_singlet_states, set_zero=.true.)
!
      call mem%alloc(converged, this%n_singlet_states)
      converged = .false.
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call mem%alloc(R, this%wf%n_es_amplitudes, this%n_singlet_states)
      call mem%alloc(X, this%wf%n_es_amplitudes, this%n_singlet_states)
!
      do state = 1, this%n_singlet_states
         call this%start_vectors%get(X(:,state), state, this%energies(state))
      end do
!
!     Enter iterative loop
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. this%max_iterations))
!
         iteration = iteration + 1
!
         write(timer_name, '(a,i0)') 'Non-linear Davidson: macro iteration ', iteration
         macro_iteration_time = timings(trim(timer_name), pl="n")
         call macro_iteration_time%turn_on()
!
!        Update energies and residuals for the micro-iterated X
!
         call output%printf('n', 'Macro iteration: (i0)', ints=[iteration], fs='(/t3,a)')
!
         call output%printf('n',&
                            'Root     Eigenvalue (Re)     Residual norm     Delta E', fs='(/t3,a)')
!
         call output%print_separator('n', 60,'-')
!
         do state = 1, this%n_singlet_states
!
!           Construct residual and compute residual norm for the given state,
!           and update the energy estimate
!
            call this%wf%construct_Jacobian_transform(this%transformation, &
                                                 X(:,state),            &
                                                 R(:,state),            &
                                                 this%energies(state))
!
            call this%projector%project(R(:,state))
!
            this%energies(state) = ddot(this%wf%n_es_amplitudes, X(:,state), 1, R(:,state), 1)
!
            call daxpy(this%wf%n_es_amplitudes, -this%energies(state), X(:,state), 1, R(:,state), 1)
!
            residual_norms(state) = get_l2_norm(R(:,state), this%wf%n_es_amplitudes)
!
!           Test the convergence
!
            converged(state) = this%convergence_checker%has_converged(residual_norms(state),&
                                               this%energies(state) - prev_energies(state), &
                                               iteration)
!
!           Print current residual and energy for the state
!
            call output%printf('n', '(i2)     (f16.12)     (e11.4)       (e11.4)',   &
                               reals=[this%energies(state), residual_norms(state), &
                                      abs(this%energies(state) - prev_energies(state))], &
                               ints=[state])
!
         enddo
!
         call output%print_separator('n', 60,'-')
!
         if (.not. all(converged)) then
!
!           Keep energies fixed and solve the eigenvalue problem
!           approximately by using the Davidson algorithm
!
            call this%do_micro_iterations(X, R,           &
                                            n_micro_iterations, &
                                            n_transformations,  &
                                            residual_norms)
!
!           Give the user some information from the micro iterations if non-verbose output
!
            call output%printf('n', 'Number of micro iterations: (i0)', &
                                 ints=[n_micro_iterations], fs='(/t6,a)')
!
            call output%printf('n', 'Number of transformations:  (i0)', &
                                 ints=[n_transformations], fs='(t6,a/)')
!
         else ! all converged!
!
            call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                               ints=[iteration], fs='(t3,a)')
!
            call this%print_summary(X)
!
         endif
!
         call this%wf%save_excited_state(X, 1, this%n_singlet_states, &
                                    this%transformation, this%energies)
!
         call dcopy(this%n_singlet_states, this%energies, 1, prev_energies, 1)
!
         call macro_iteration_time%turn_off()
         call macro_iteration_time%reset()
!
      enddo
!
      if (.not. all(converged)) &
         call output%error_msg('Unable to converge equations in the given maximum number&
                              & of macro iterations. Try to increase the value.')
!
      call this%wf%set_excitation_energies(this%energies, this%transformation)
!
      call mem%dealloc(R, this%wf%n_es_amplitudes, this%n_singlet_states)
      call mem%dealloc(X, this%wf%n_es_amplitudes, this%n_singlet_states)
      call mem%dealloc(prev_energies, this%n_singlet_states)
      call mem%dealloc(residual_norms, this%n_singlet_states)
!
      call mem%dealloc(converged, this%n_singlet_states)
!
   end subroutine run_nonlinear_davidson_cc_es
!
!
   subroutine do_micro_iterations(this, X, R,   &
                                                           n_micro_iterations, &
                                                           n_transformations,  &
                                                           residual_norms)
!!
!!    Do micro-iterations
!!    Written by Eirik F. Kjønstad, Nov 2019 and Jan 2020
!!
!!    Arguments:
!!
!!       X:  Column k contains the current estimates for solutions to (*)
!!           The corresponding excitation energy omega_k is assumed to be given by
!!           this%energies(k). On exit, the columns contain the updated guesses
!!           for the solutions to (*).
!!
!!       R:  Column k contains the residual corrensponding to X_k,
!!           namely A(omega_k) X_k - omega_k X_k. Here, omega_k is
!!           the previous estimate of the kth excitation energy.
!!
!!       n_micro_iterations: On exit, the number of micro iterations needed
!!       n_transformations:  On exit, the number of Jacobian transformations needed
!!
!!       residual_norms: Vector containing the norms of the columns of R.
!!
!!    Description:
!!
!!       Approximately solves the eigenvalue problem
!!
!!          A(omega_k) X_k = omega_k X_k     (*)
!!
!!       for the current solver energies omega_k.
!!
!!       The algorithm sets the first trial vectors equal to the current guesses for X_k
!!       and does not orthonormalize them (instead this is handled by a non-unit metric in
!!       the Davidson tool). The first set of preconditioned residuals are then given by
!!       the standard quasi-Newton estimates: (eps - omega_k)^{-1} R_k.
!!
!!       When performing the transformation of trial vectors, A(alpha) c,
!!       where c is some trial and alpha the corresponding energy,
!!       then alpha = omega_k when c is X_k or a trial
!!       generated from a residual corresponding to a guess for X_k.
!!
      use array_utilities, only: invert, get_l2_norm
!
      implicit none
!
      class(nonlinear_davidson_cc_es) :: this
!
!     Current solution guesses and associated residuals and residual norms
!
      real(dp), dimension(this%wf%n_es_amplitudes, this%n_singlet_states), intent(inout) :: X
      real(dp), dimension(this%wf%n_es_amplitudes, this%n_singlet_states), intent(in)    :: R
!
      real(dp), dimension(this%n_singlet_states), intent(in) :: residual_norms
!
!     On exit, number of micro iterations needed and number of transformations needed
!
      integer, intent(out) :: n_micro_iterations, n_transformations
!
      real(dp), dimension(:), allocatable :: c         ! stores trial vector temporarily
      real(dp), dimension(:), allocatable :: residual  ! stores residual temporarily
!
      integer :: state, trial, iteration
!
!     Vectors containing solution information (state norms and residual norms)
!     and whether it has converged yet
!
      real(dp), dimension(:), allocatable :: norm_X
      real(dp), dimension(:), allocatable :: micro_residual_norms
!
      logical,  dimension(:), allocatable :: converged
!
      real(dp) :: lindep_threshold                          ! linear dependence threshold
                                                            ! for Davidson
!
      integer :: corresponding_state                        ! state that corresponds to trial
      integer,  dimension(:), allocatable :: trial_to_state ! map from trial to corresponding state
!
      real(dp), dimension(:), allocatable :: eps            ! orbital differences
!
      real(dp) :: micro_residual_threshold
!
      real(dp), parameter :: default_lindep_threshold = 1.0d-11
!
      type(timings) :: micro_iteration_time
      character(len=40) :: timer_name
!
!     The residual threshold in the microiterations is taken to be proportional to the maximum
!     residual norm in the current macroiteration or the final residual threshold
!
      micro_residual_threshold = max(half*this%convergence_checker%residual_threshold, &
                                 this%relative_micro_residual_threshold*maxval(residual_norms))
      lindep_threshold         = min(default_lindep_threshold, micro_residual_threshold)
!
      if (lindep_threshold .lt. default_lindep_threshold) then
!
         call output%warning_msg('Linear dependence threshold was set to (e9.2), which is below &
                                 &the default value of (e9.2). May lead to instabilities. &
                                 &Consider loosening the relative micro threshold.', &
                                 reals=[lindep_threshold, default_lindep_threshold])
!
      endif
!
      call this%davidson%set_lindep_threshold(lindep_threshold)
!
      call output%printf('n', 'Starting on microiterations', fs='(/t6,a)')
!
      call this%davidson%initialize()
!
      call mem%alloc(eps, this%wf%n_es_amplitudes)
!
      call this%wf%get_orbital_differences(eps, this%wf%n_es_amplitudes)
      call this%davidson%set_preconditioner(eps)
!
      call mem%dealloc(eps, this%wf%n_es_amplitudes)
!
!     Set up initial trial space and save transforms
!
      call mem%alloc(trial_to_state, this%n_singlet_states)
      call mem%alloc(c, this%wf%n_es_amplitudes)
!
      do state = 1, this%n_singlet_states
!
         call this%davidson%set_trial(X(:,state), state)
         trial_to_state(state) = state
!
         call dcopy(this%wf%n_es_amplitudes, R(:,state), 1, c, 1)
         call daxpy(this%wf%n_es_amplitudes, this%energies(state), X(:,state), 1, c, 1)
!
         call this%davidson%set_transform(c, state)
!
      enddo
!
      call this%davidson%update_reduced_dim()
!
!     Initialize quantities used in the loop
!     and enter the iterative micro iterations Davidson loop
!
      call mem%alloc(converged, this%n_singlet_states)
      converged = .false.
!
      iteration         = 0
      n_transformations = 0
!
      call mem%alloc(residual, this%wf%n_es_amplitudes)
!
      call mem%alloc(norm_X, this%n_singlet_states, set_zero=.true.)
      call mem%alloc(micro_residual_norms, this%n_singlet_states, set_zero=.true.)
!
      do while (.not. all(converged) .and. iteration .le. this%max_micro_iterations)
!
         iteration = iteration + 1
!
         write(timer_name, '(a,i0)') 'Non-linear Davidson: micro iteration ', iteration
         micro_iteration_time = timings(trim(timer_name), pl="n")
         call micro_iteration_time%turn_on()
!
!        Solve reduced eigenvalue problem
!
         call this%davidson%solve_reduced_problem()
!
!        Loop over states and construct residuals, next trial vectors,
!        as well as check for convergence
!
         call output%printf('v', 'Micro iteration: (i0)', ints=[iteration], fs='(/t6,a)')
!
         call output%printf('v', 'Reduced space dimension: (i0)', &
                              ints=[this%davidson%dim_red], fs='(t6,a)')
!
         call output%printf('v', &
                            'Root     Eigenvalue (Re)        Eigenvalue (Im)      Residual norm', &
                            fs='(/t6,a)')
!
         call output%print_separator('v', 68,'-', fs='(t6,a)')
!
         trial = 0
!
         do state = 1, this%n_singlet_states
!
            call this%davidson%construct_solution(X(:,state), state)
!
            call this%davidson%construct_residual(residual, X(:,state), state)
!
            norm_X(state) = get_l2_norm(X(:,state), this%wf%n_es_amplitudes)
!
            micro_residual_norms(state) = get_l2_norm(residual, this%wf%n_es_amplitudes)/norm_X(state)
!
            converged(state) = micro_residual_norms(state) .le. micro_residual_threshold
!
            if (.not. converged(state)) then
!
               if (micro_residual_norms(state) .le. lindep_threshold) then
!
                  call output%warning_msg('Residual norm for root (i0) smaller than linear ' // &
                                          'dependence threshold, but energy and residual '   // &
                                          'thresholds have not yet been met. No new trial '  // &
                                          'added for this root.', ints=[state], fs='(/t6,a)')
!
               else
!
                  trial = trial + 1
!
                  trial_to_state(trial) = state
                  call this%davidson%add_new_trial(residual, state)
!
               endif
!
            endif
!
            call output%printf('v', '(i2)     (f16.12)       (f16.12)           (e11.4)', &
            ints=[state], reals=[this%davidson%omega_re(state), this%davidson%omega_im(state), &
                                 micro_residual_norms(state)], fs='(t6,a)')
!
         enddo
!
         call output%print_separator('v', 68,'-', fs='(t6,a)')
!
         if (.not. all(converged)) then
!
!           Reduced space preparations
!
            if (this%davidson%red_dim_exceeds_max()) then
               call this%davidson%set_trials_to_solutions()
            end if
!
            call this%davidson%update_reduced_dim()
!
            call this%davidson%orthonormalize_trial_vecs()
!
!           Transform new orthonormalized trial vectors
!
            do trial = this%davidson%first_new_trial(), this%davidson%last_new_trial()
!
               n_transformations = n_transformations + 1
!
               call this%davidson%get_trial(c, trial)
!
!              From which state residual did this trial vector originate?
!              We use that state's energy in the Jacobian transformation of the trial
!
               corresponding_state = trial_to_state(trial - this%davidson%first_new_trial() + 1)
!
               call this%wf%construct_Jacobian_transform(this%transformation,             &
                                                    c,                                    &
                                                    residual,                             &
                                                    this%energies(corresponding_state))
!
               call this%projector%project(residual)
!
               call this%davidson%set_transform(residual, trial)
!
            enddo
!
         endif
!
         call micro_iteration_time%turn_off()
         call micro_iteration_time%reset()
!
      enddo  ! End of iterative loop!
!
!     Note: the full space solutions are not orthonormal in general because the metric
!           is not equal to the identity matrix:
!
!           x^T x = sum_ij x_i b_i^T b_j x_j = x_red^T S_red x_red,   x_red^T x_red = 1.
!
!     Hence, we have to normalize the converged states. (This is also why we divide by
!     the norm of X when testing convergence above.)
!
      do state = 1, this%n_singlet_states
!
         call dscal(this%wf%n_es_amplitudes, one/norm_X(state), X(:, state), 1)
!
      enddo
!
      n_micro_iterations = iteration
!
      if (.not. all(converged)) &
         call output%error_msg('Unable to converge equations in the given maximum number&
                              & of micro iterations. Try to increase the value.')
!
!     Final deallocations before quitting micro iterations
!
      call mem%dealloc(c, this%wf%n_es_amplitudes)
      call mem%dealloc(residual, this%wf%n_es_amplitudes)
!
      call mem%dealloc(norm_X, this%n_singlet_states)
      call mem%dealloc(micro_residual_norms, this%n_singlet_states)
      call mem%dealloc(trial_to_state, this%n_singlet_states)
!
      call this%davidson%cleanup()
!
      call mem%dealloc(converged, this%n_singlet_states)
!
   end subroutine do_micro_iterations
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
      use global_in, only: input
!
      implicit none
!
      class(nonlinear_davidson_cc_es) :: this
!
      this%max_dim_red = max(100, 10*this%n_singlet_states)
!
      call input%get_keyword('max reduced dimension', &
                             'solver cc es',          &
                             this%max_dim_red)
!
      call input%get_keyword('max micro iterations', &
                             'solver cc es',         &
                             this%max_micro_iterations)
!
      call input%get_keyword('rel micro threshold', &
                             'solver cc es',        &
                             this%relative_micro_residual_threshold)
!
   end subroutine read_settings
!
!
end module nonlinear_davidson_cc_es_class
