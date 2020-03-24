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
!!       A_mu,nu = < mu | [H-bar, tau_nu] | HF >,    H-bar = e-T H eT. 
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
!!    A_mu,nu = < mu | [F, tau_nu] | HF > = epsilon_mu delta_mu,nu 
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
   use timings_class, only: timings
   use global_in, only: input
   use global_out, only: output
   use ccs_class, only: ccs
   use string_utilities, only: convert_to_uppercase
   use array_utilities, only: get_l2_norm, zero_array
   use memory_manager_class, only: mem
   use eigen_davidson_tool_class, only: eigen_davidson_tool
   use abstract_cc_es_class, only: abstract_cc_es
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
   contains
!     
      procedure :: run                               &
                => run_nonlinear_davidson_cc_es
!
      procedure :: read_settings                     &
                => read_settings_nonlinear_davidson_cc_es
!
      procedure :: read_davidson_settings            &
                => read_davidson_settings_nonlinear_davidson_cc_es
!
      procedure :: print_settings                    &
                => print_settings_nonlinear_davidson_cc_es
!
      procedure :: do_micro_iterations               &
                => do_micro_iterations_nonlinear_davidson_cc_es
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
   function new_nonlinear_davidson_cc_es(transformation, wf, restart) result(solver)
!!
!!    New non-linear Davidson CC ES 
!!    Written by Eirik F. Kjønstad, Nov 2019 and Jan 2020
!!
!!    Mostly identical to the constructor for the standard Davidson solver 
!!    (Sarai D. Folkestad and Eirik F. Kjønstad, 2018).
!!
      implicit none
!
      type(nonlinear_davidson_cc_es) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      logical, intent(in) :: restart
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) &
                        // ' excited state (' // trim(transformation) //')', pl='normal')
      call solver%timer%turn_on()
!
      solver%name_ = 'Non-linear Davidson coupled cluster excited state solver'
      solver%tag   = 'Davidson'
!
      solver%description1 = 'A non-linear Davidson solver that calculates the lowest&
               & eigenvalues and the right or left eigenvectors of the Jacobian matrix A.&
               & The eigenvalue problem is solved in a reduced space for fixed energies,&
               & before the energies are updated and a new reduced space built.'
!
      solver%description2 = 'See C. Hättig and F. Weigend, J. Chem. Phys. 113, 5154 (2000).'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%n_singlet_states                   = 0
      solver%max_iterations                     = 100
      solver%max_micro_iterations               = 100
      solver%eigenvalue_threshold               = 1.0d-6
      solver%residual_threshold                 = 1.0d-6
      solver%relative_micro_residual_threshold  = 1.0d-1
      solver%restart                            = restart
      solver%max_dim_red                        = 100 
      solver%transformation                     = trim(transformation)
      solver%es_type                            = 'valence'
      solver%records_in_memory                  = .false.  
      solver%storage                            = 'disk'
      solver%prepare_wf                         = .true.
!
      call solver%read_settings()
      call solver%print_settings()
!
      if (solver%n_singlet_states == 0) &
         call output%error_msg('number of excitations must be specified.')
!
      wf%n_singlet_states = solver%n_singlet_states
!
!     Determine whether to store records in memory or on file
!
      if (trim(solver%storage) == 'memory') then 
!
         solver%records_in_memory = .true.
!
      elseif (trim(solver%storage) == 'disk') then 
!
         solver%records_in_memory = .false.
!
      else 
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                 trim(solver%storage))
!
      endif 
!
   end function new_nonlinear_davidson_cc_es
!
!
   function new_nonlinear_davidson_cc_es_preconvergence(transformation,                      &
                                                        wf,                                  &
                                                        restart,                             &
                                                        energy_threshold,                    &
                                                        residual_threshold,                  &
                                                        max_iterations,                      &
                                                        relative_micro_residual_threshold,   &
                                                        max_micro_iterations,                &
                                                        max_dim_red,                         &
                                                        es_type,                             &
                                                        storage,                             &
                                                        n_singlet_states,                    &
                                                        prepare_wf) result(solver)
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
      implicit none 
!
      type(nonlinear_davidson_cc_es) :: solver 
!
      class(ccs) :: wf 
!
      integer, intent(in) :: max_iterations, max_micro_iterations, max_dim_red, n_singlet_states
!
      real(dp), intent(in) :: energy_threshold, residual_threshold
      real(dp), intent(in) :: relative_micro_residual_threshold
!
      character(len=*), intent(in) :: es_type, storage, transformation  
!
      logical, intent(in) :: restart, prepare_wf 
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) &
                        // ' excited state (' // trim(transformation) //')', pl='normal')
      call solver%timer%turn_on()
!
      solver%name_ = 'Non-linear Davidson coupled cluster excited state solver'
      solver%tag   = 'Davidson'
!
      solver%description1 = 'A Davidson solver that calculates the lowest eigenvalues and &
               & the right or left eigenvectors of the Jacobian matrix, A. The eigenvalue &
               & problem is solved in a reduced space, the dimension of which is &
               & expanded until the convergence criteria are met.'
!
      solver%description2 = 'A complete description of the algorithm can be found in &
                                 & E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%n_singlet_states                   = n_singlet_states
      solver%max_iterations                     = max_iterations
      solver%max_micro_iterations               = max_micro_iterations
      solver%eigenvalue_threshold               = energy_threshold
      solver%residual_threshold                 = residual_threshold
      solver%relative_micro_residual_threshold  = relative_micro_residual_threshold
      solver%restart                            = restart
      solver%max_dim_red                        = max_dim_red 
      solver%transformation                     = trim(transformation)
      solver%es_type                            = trim(es_type)
      solver%storage                            = trim(storage)
      solver%prepare_wf                         = prepare_wf
!
      call solver%print_settings()
!
!     Determine whether to store records in memory or on file
!
      if (trim(solver%storage) == 'memory') then 
!
         solver%records_in_memory = .true.
!
      elseif (trim(solver%storage) == 'disk') then 
!
         solver%records_in_memory = .false.
!
      else 
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                 trim(solver%storage))
!
      endif 
!
   end function new_nonlinear_davidson_cc_es_preconvergence
!
!
   subroutine print_settings_nonlinear_davidson_cc_es(solver)
!!
!!    Print settings    
!!    Written by Eirik F. Kjønstad, Jan 2020 
!!
      implicit none 
!
      class(nonlinear_davidson_cc_es) :: solver 
!
      call solver%print_es_settings()
!
      call output%printf('m', 'Max reduced space dimension:         (i4)',  &
                         ints=[solver%max_dim_red], fs='(/t6,a)')
!
      call output%printf('m', 'Max micro iterations:                (i4)',  &
                         ints=[solver%max_micro_iterations], fs='(t6,a)')
!
      call output%printf('m', 'Relative micro threshold:     (e11.2)', &
                         reals=[solver%relative_micro_residual_threshold], fs='(t6,a)')
!
   end subroutine print_settings_nonlinear_davidson_cc_es
!
!
   subroutine run_nonlinear_davidson_cc_es(solver, wf)
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
      implicit none
!
      class(nonlinear_davidson_cc_es) :: solver
!
      class(ccs) :: wf
!
!     Convergence logical arrays, one component for each state 
!
      logical, dimension(:), allocatable :: converged
      logical, dimension(:), allocatable :: converged_eigenvalue
      logical, dimension(:), allocatable :: converged_residual
!
!     Previous energies, current residuals, one component for each state 
!
      real(dp), dimension(:), allocatable :: prev_energies
      real(dp), dimension(:), allocatable :: residual_norms 
!
      integer :: n_solutions_on_file
      integer :: iteration, state
!
      integer :: n_micro_iterations, n_transformations ! Information regarding micro iterations 
!
      real(dp) :: ddot
!
!     Current guess of solutions and temporary norm variable
!     (columns stores the different states, i.e. X(:,k) is the guess for X_k)
!
      real(dp) :: norm_X
      real(dp), dimension(:,:), allocatable :: X
!
!     Associated residuals, R(:,k) = A(omega_k) X(:,k) - omega_k X(:,k)
!     The current guess for the energy omega_k is saved in solver%energies(k).
!
      real(dp), dimension(:,:), allocatable :: R
!
!     Initialize solver tools and prepare wavefunction
!
      call solver%initialize_start_vector_tool(wf)
      call solver%initialize_projection_tool(wf)
!
      if (solver%prepare_wf) &
         call solver%prepare_wf_for_excited_state(wf)
!
!     Initialize energies, residual norms, and convergence arrays 
!
      call solver%initialize_energies()
      solver%energies = zero
!
      call mem%alloc(prev_energies, solver%n_singlet_states)
      call mem%alloc(residual_norms, solver%n_singlet_states)
!
      call zero_array(prev_energies, solver%n_singlet_states)
      call zero_array(residual_norms, solver%n_singlet_states)
!
      call mem%alloc(converged, solver%n_singlet_states)
      call mem%alloc(converged_residual, solver%n_singlet_states)
      call mem%alloc(converged_eigenvalue, solver%n_singlet_states)
!
      converged            = .false.
      converged_residual   = .false.
      converged_eigenvalue = .false.
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call mem%alloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
      call mem%alloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
!
      do state = 1, solver%n_singlet_states
!
         call solver%start_vectors%get(X(:,state), state)
!
      enddo 
!
!     Overwrite some or all of the initial guesses if restart is requested 
!
      if (solver%restart) then 
!
         call solver%determine_restart_transformation(wf) ! Read right or left?
!
         n_solutions_on_file = wf%get_n_excited_states_on_file(solver%restart_transformation)
!
         call output%printf('n', 'Requested restart - there are (i0) (a0) eigenvectors on file.', &
                              ints=[n_solutions_on_file], chars=[solver%restart_transformation])
!
         do state = 1, n_solutions_on_file
!
            call wf%read_excited_state(X(:,state), state, solver%restart_transformation)
!
!           Normalize X in case it is not normalized
!
            norm_X = get_l2_norm(X(:,state), wf%n_es_amplitudes)
            call dscal(wf%n_es_amplitudes, one/norm_X, X(:, state), 1)
!
         enddo
!
         solver%energies = zero
         call wf%read_excitation_energies(n_solutions_on_file, &
                                          solver%energies(1:n_solutions_on_file))
!
      endif
!
!     Enter iterative loop
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1   
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
         do state = 1, solver%n_singlet_states
!
!           Construct residual and compute residual norm for the given state,
!           and update the energy estimate 
!
            call dcopy(wf%n_es_amplitudes, X(:,state), 1, R(:,state), 1)
!
            call wf%construct_Jacobian_transform(solver%transformation, R(:,state), &
                                                 solver%energies(state))
!
            if (solver%projector%active) call solver%projector%do_(R(:,state))            
!
            solver%energies(state) = ddot(wf%n_es_amplitudes, X(:,state), 1, R(:,state), 1)
!
            call daxpy(wf%n_es_amplitudes, -solver%energies(state), X(:,state), 1, R(:,state), 1)
!
            residual_norms(state) = get_l2_norm(R(:,state), wf%n_es_amplitudes)
!
!           Test the convergence of residual and energy 
!
            converged_residual(state) = residual_norms(state) .lt. solver%residual_threshold
!
            converged_eigenvalue(state) = abs(solver%energies(state) - prev_energies(state)) &
                                     .lt. solver%eigenvalue_threshold

!
            converged(state) = (converged_eigenvalue(state) .and. converged_residual(state)) &
                          .or. (converged_residual(state) .and. iteration .eq. 1)

!
!           Print current residual and energy for the state 
!
            call output%printf('n', '(i2)     (f16.12)     (e11.4)       (e11.4)',     &
                                       ints=[state],                                   &
                                       reals=[solver%energies(state),                  &
                                              residual_norms(state),                   &
                                              abs(solver%energies(state) - prev_energies(state))])
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
            call solver%do_micro_iterations(wf,                            &
                                            X,                             &
                                            R,                             &
                                            n_micro_iterations,            &
                                            n_transformations,             &
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
            call solver%print_summary(wf, X) 
!
         endif 
!
!        Save excited states and excitation energies
!
         do state = 1, solver%n_singlet_states
!
            call wf%save_excited_state(X(:,state), state, solver%transformation)
!
         enddo 
!
         call wf%save_excitation_energies(solver%n_singlet_states,   &
                                          solver%energies,           &
                                          solver%transformation)
!
         call dcopy(solver%n_singlet_states, solver%energies, 1, prev_energies, 1)
!
      enddo 
!
      if (.not. all(converged)) &
         call output%error_msg('Unable to converge equations in the given maximum number&
                              & of macro iterations. Try to increase the value.')
!
      call mem%dealloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
      call mem%dealloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
      call mem%dealloc(prev_energies, solver%n_singlet_states)
      call mem%dealloc(residual_norms, solver%n_singlet_states)
!
      call mem%dealloc(converged, solver%n_singlet_states)
      call mem%dealloc(converged_residual, solver%n_singlet_states)
      call mem%dealloc(converged_eigenvalue, solver%n_singlet_states)
!
   end subroutine run_nonlinear_davidson_cc_es
!
!
   subroutine do_micro_iterations_nonlinear_davidson_cc_es(solver, wf, X, R,     &
                                                         n_micro_iterations,  &
                                                         n_transformations,   &
                                                         residual_norms)
!!
!!    Do micro-iterations
!!    Written by Eirik F. Kjønstad, Nov 2019 and Jan 2020 
!!
!!    Arguments:
!!
!!       X:  Column k contains the current estimates for solutions to (*) 
!!           The corresponding excitation energy omega_k is assumed to be given by 
!!           solver%energies(k). On exit, the columns contain the updated guesses 
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
      use array_utilities, only: invert
!
      implicit none 
!
      class(ccs) :: wf 
!
      class(nonlinear_davidson_cc_es) :: solver 
!
!     Current solution guesses and associated residuals and residual norms
!
      real(dp), dimension(wf%n_es_amplitudes, solver%n_singlet_states), intent(inout) :: X 
      real(dp), dimension(wf%n_es_amplitudes, solver%n_singlet_states), intent(in)    :: R 
!
      real(dp), dimension(solver%n_singlet_states), intent(in) :: residual_norms
!
!     On exit, number of micro iterations needed and number of transformations needed 
!
      integer, intent(out) :: n_micro_iterations, n_transformations
!
      class(eigen_davidson_tool), allocatable :: davidson ! tool to handle reduced space 
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
      call output%printf('n', 'Starting on microiterations', fs='(/t6,a)')
!
!     The residual threshold in the microiterations is taken to be proportional to the maximum 
!     residual norm in the current macroiteration
!
      micro_residual_threshold = solver%relative_micro_residual_threshold*maxval(residual_norms) 
!
      lindep_threshold         = 1.0d-11
!
!     Initialize Davidson tool 
!
      davidson = eigen_davidson_tool('nonlin_cc_es_davidson',  &
                                     wf%n_es_amplitudes,       &
                                     solver%n_singlet_states,  &
                                     lindep_threshold,         &
                                     solver%max_dim_red,       &
                                     non_unit_metric=.true.)
!
      call davidson%initialize_trials_and_transforms(solver%records_in_memory)      
!
      call mem%alloc(eps, wf%n_es_amplitudes)
!
      call wf%get_es_orbital_differences(eps, wf%n_es_amplitudes)
      call davidson%set_preconditioner(eps)
!
      call mem%dealloc(eps, wf%n_es_amplitudes)
!
!     Set up initial trial space and save transforms 
!
      call mem%alloc(trial_to_state, solver%n_singlet_states)
      call mem%alloc(c, wf%n_es_amplitudes)
!
      do state = 1, solver%n_singlet_states
!
         call davidson%set_trial(X(:,state), state)
         trial_to_state(state) = state 
!
         call dcopy(wf%n_es_amplitudes, R(:,state), 1, c, 1)
         call daxpy(wf%n_es_amplitudes, solver%energies(state), X(:,state), 1, c, 1)
!
         call davidson%set_transform(c, state)
!
      enddo 
!
      call davidson%update_reduced_dim()
!
!     Initialize quantities used in the loop 
!     and enter the iterative micro iterations Davidson loop 
!
      call mem%alloc(converged, solver%n_singlet_states)
!
      converged = .false.
!
      iteration         = 0
      n_transformations = 0
!
      call mem%alloc(residual, wf%n_es_amplitudes)
!
      call mem%alloc(norm_X, solver%n_singlet_states)
      call mem%alloc(micro_residual_norms, solver%n_singlet_states)
!
      call zero_array(norm_X, solver%n_singlet_states)
      call zero_array(micro_residual_norms, solver%n_singlet_states)
!
      do while (.not. all(converged) .and. iteration .le. solver%max_micro_iterations)
!
         iteration = iteration + 1
!
!        Solve reduced eigenvalue problem 
!
         call davidson%solve_reduced_problem()
!
!        Loop over states and construct residuals, next trial vectors,
!        as well as check for convergence 
!
         call output%printf('v', 'Micro iteration: (i0)', ints=[iteration], fs='(/t6,a)')
!
         call output%printf('v', 'Reduced space dimension: (i0)', &
                              ints=[davidson%dim_red], fs='(t6,a)')
!
         call output%printf('v', &
                            'Root     Eigenvalue (Re)        Eigenvalue (Im)      Residual norm', &
                            fs='(/t6,a)')
!
         call output%print_separator('v', 68,'-', fs='(t6,a)')
!
         trial = 0
!
         do state = 1, solver%n_singlet_states
!
            call davidson%construct_solution(X(:,state), state)
!
            call davidson%construct_residual(residual, X(:,state), state)
!
            norm_X(state) = get_l2_norm(X(:,state), wf%n_es_amplitudes)
!
            micro_residual_norms(state) = get_l2_norm(residual, wf%n_es_amplitudes)/norm_X(state)
!
            converged(state) = micro_residual_norms(state) .le. micro_residual_threshold
!
            if (.not. converged(state)) then 
!
               if (micro_residual_norms(state) .le. lindep_threshold) then 
!
                  call output%warning_msg('Residual norm for root (i0) smaller than linear ' // &
                                          'dependence threshold, but energy and residual  '  // &
                                          'thresholds have not yet been met. No new trial ' // &
                                          'added for this root.', ints=[state], fs='(/t6,a)')
!
               else
!
                  trial = trial + 1
!
                  trial_to_state(trial) = state 
                  call davidson%construct_next_trial(residual, state)
!
               endif
!
            endif
!
            call output%printf('v', '(i2)     (f16.12)       (f16.12)           (e11.4)', &
            ints=[state], reals=[davidson%omega_re(state), davidson%omega_im(state), &
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
            if (davidson%red_dim_exceeds_max()) call davidson%set_trials_to_solutions()
!
            call davidson%update_reduced_dim()
!
            call davidson%orthonormalize_trial_vecs() 
!
!           Transform new orthonormalized trial vectors 
!
            do trial = davidson%first_new_trial(), davidson%last_new_trial()
!
               n_transformations = n_transformations + 1
!
               call davidson%get_trial(c, trial)
!
!              From which state residual did this trial vector originate? 
!              We use that state's energy in the Jacobian transformation of the trial
!
               corresponding_state = trial_to_state(trial - davidson%first_new_trial() + 1)
!
               call wf%construct_Jacobian_transform(solver%transformation,                &
                                                    c,                                    &
                                                    solver%energies(corresponding_state))
!
               if (solver%projector%active) call solver%projector%do_(c)
! 
               call davidson%set_transform(c, trial)
!
            enddo 
!
         endif 
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
      do state = 1, solver%n_singlet_states
!
         call dscal(wf%n_es_amplitudes, one/norm_X(state), X(:, state), 1)
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
      call mem%dealloc(c, wf%n_es_amplitudes)
      call mem%dealloc(residual, wf%n_es_amplitudes)
!
      call mem%dealloc(norm_X, solver%n_singlet_states)
      call mem%dealloc(micro_residual_norms, solver%n_singlet_states)
      call mem%dealloc(trial_to_state, solver%n_singlet_states)
!
      call davidson%finalize_trials_and_transforms()
!
      call mem%dealloc(converged, solver%n_singlet_states)
!
   end subroutine do_micro_iterations_nonlinear_davidson_cc_es
!
!
   subroutine read_settings_nonlinear_davidson_cc_es(solver)
!!
!!    Read settings 
!!    Written by Eirik F. Kjønstad, Jan 2020 
!!
      implicit none 
!
      class(nonlinear_davidson_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_davidson_settings()
!
   end subroutine read_settings_nonlinear_davidson_cc_es
!
!
   subroutine read_davidson_settings_nonlinear_davidson_cc_es(solver)
!!
!!    Read settings 
!!    Written by Eirik F. Kjønstad, Jan 2020 
!!
      implicit none 
!
      class(nonlinear_davidson_cc_es) :: solver 
!
      call input%get_keyword_in_section('max reduced dimension',  &
                                        'solver cc es',           &
                                        solver%max_dim_red)
!
      call input%get_keyword_in_section('max micro iterations',  &
                                        'solver cc es',           &
                                        solver%max_micro_iterations)
!
      call input%get_keyword_in_section('rel micro threshold',    &
                                        'solver cc es',           &
                                        solver%relative_micro_residual_threshold)
!
   end subroutine read_davidson_settings_nonlinear_davidson_cc_es
!
!
end module nonlinear_davidson_cc_es_class
