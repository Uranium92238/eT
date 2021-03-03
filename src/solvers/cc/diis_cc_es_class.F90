!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module diis_cc_es_class
!
!!
!! DIIS coupled cluster excited state solver class module
!! Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
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
!! The equation is solved using the direct inversion of the iterative 
!! subspace (DIIS) algorithm. See Pulay, P. Convergence acceleration 
!! of iterative sequences. The case of SCF iteration. Chem. Phys. Lett. 
!! 1980, 73, 393−398. This algorithm combines estimates for the parameters 
!! and the errors and finds a least-squares solution to the error being 
!! zero (see diis_tool solver tool for more details). 
!!
!! The update estimates used in DIIS are the following. Assuming we 
!! are solving for a right excited state R, the update estimate is
!!
!!    R_mu <- R_mu - Residual_mu/(epsilon_mu - omega),
!!
!! where Residual = A R - omega R and omega = R^T A R.
!! (R is normalized.) 
!!
!! This solver can handle the folded CC2 and CC3 equations, where 
!! the Jacobian transformation becomes non-linear 
!!
!!    A(omega) R = omega R   or    L^T A(omega) = omega L^T.
!!
!! In this case, Davidson's algorithm is not as well-suited. 
!! Literature for DIIS/CC2: C. Hättig and F. Weigend, 
!! J. Chem. Phys. 113, 5154 (2000).
!!
!! Note/warning for usage: if poor starting guesses for L and R,  
!! are used this solver often finds duplicate roots and may in some
!! cases not converge. 
!!
!
   use parameters
   use global_out, only: output
   use global_in, only: input
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use ccs_class, only: ccs
   use diis_tool_class, only: diis_tool
   use abstract_cc_es_class, only: abstract_cc_es
   use es_valence_start_vector_tool_class, only: es_valence_start_vector_tool
   use es_valence_projection_tool_class, only: es_valence_projection_tool
   use array_utilities, only: quicksort_with_index_ascending, get_l2_norm
   use string_utilities, only: convert_to_uppercase
   use precondition_tool_class, only: precondition_tool
   use convergence_tool_class, only: convergence_tool 
!
   implicit none
!
   type, extends(abstract_cc_es) :: diis_cc_es
!
      integer :: diis_dimension
!
      logical :: crop ! Standard DIIS if false; CROP variant of DIIS if true
!
      logical :: davidson_preconvergence     ! Perform non-linear Davidson first, then go over 
                                             ! to DIIS to converge below the residual threshold 
!
      real(dp) :: preconvergence_threshold   ! Non-linear Davidson threshold
!
   contains
!     
      procedure, non_overridable :: run            => run_diis_cc_es
!
      procedure :: read_settings                   => read_settings_diis_cc_es
      procedure :: read_diis_settings              => read_diis_settings_diis_cc_es
!
      procedure :: print_settings                  => print_settings_diis_cc_es
!
      procedure :: do_davidson_preconvergence      => do_davidson_preconvergence_diis_cc_es
!
   end type diis_cc_es
!
!
   interface diis_cc_es 
!
      procedure :: new_diis_cc_es
!
   end interface diis_cc_es
!
!
contains
!
!
   function new_diis_cc_es(transformation, wf, restart) result(solver)
!!
!!    New DIIS CC ES 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(diis_cc_es) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      logical, intent(in) :: restart
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' excited state (' // trim(transformation) //')')
      call solver%timer%turn_on()
!
!     Set printables
!
      solver%name_  = 'DIIS coupled cluster excited state solver'
      solver%tag    = 'DIIS'
!
      solver%description1 = 'A DIIS solver that solves for the lowest eigenvalues and &
                           & the right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                           & problem is solved by DIIS extrapolation of residuals for each &
                           & eigenvector until the convergence criteria are met.'
!
      solver%description2 = 'More on the DIIS algorithm can be found in &
                           &P. Pulay, Chemical Physics Letters, 73(2), 393-398 (1980).'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%n_singlet_states          = 0
      solver%max_iterations            = 100
      solver%diis_dimension            = 20
      solver%restart                   = restart
      solver%transformation            = trim(transformation)
      solver%es_type                   = 'valence'
      solver%records_in_memory         = .false.
      solver%storage                   = 'disk'
      solver%crop                      = .false.
      solver%davidson_preconvergence   = .false.
      solver%preconvergence_threshold  = 1.0d-2
!
!     Initialize convergence checker with default threshols
!
      solver%convergence_checker = convergence_tool(energy_threshold   = 1.0d-3,   &
                                                    residual_threshold = 1.0d-3,   &
                                                    energy_convergence = .false.)
!
      call solver%read_settings()
!
      call solver%print_settings()
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      wf%n_singlet_states = solver%n_singlet_states
! 
   end function new_diis_cc_es
!
!
   subroutine print_settings_diis_cc_es(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_es) :: solver 
!
      call solver%print_es_settings()
!
      call output%printf('m', 'DIIS dimension:               (i11)', &
                         ints=[solver%diis_dimension], fs='(/t6,a/)')
!
      if (solver%crop) then 
!
         call output%printf('m', 'Enabled CROP in the DIIS algorithm.', fs='(t6,a/)')
!
      endif
!
      if (solver%davidson_preconvergence) then 
!
         call output%printf('m', 'Enabled preconvergence using&
                                 & the non-linear Davidson solver.', fs='(t6,a)')
!
         call output%printf('m', 'Preconvergence threshold:   (e11.2)', &
                                 reals=[solver%preconvergence_threshold], fs='(/t6,a/)')
!
      endif
!
   end subroutine print_settings_diis_cc_es
!
!
   subroutine read_settings_diis_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_diis_settings()
!
   end subroutine read_settings_diis_cc_es
!
!
   subroutine read_diis_settings_diis_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_es) :: solver 
!
      call input%get_keyword('diis dimension', 'solver cc es', solver%diis_dimension)
!
      if (input%is_keyword_present('crop', 'solver cc es')) then 
!
         solver%crop = .true.
!
      endif
!
      if (input%is_keyword_present('davidson preconvergence', 'solver cc es')) then 
!
         solver%davidson_preconvergence = .true.
!
      endif
!
      call input%get_keyword('preconvergence threshold',  &
                                        'solver cc es',              &
                                        solver%preconvergence_threshold)
!
   end subroutine read_diis_settings_diis_cc_es
!
!
   subroutine run_diis_cc_es(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es) :: solver
!
      class(ccs) :: wf
!
      logical, dimension(:), allocatable :: converged
!
      real(dp), dimension(:), allocatable :: prev_energies
      real(dp), dimension(:), allocatable :: residual_norms 
!
      integer, dimension(:), allocatable :: X_order
!
      type(diis_tool), dimension(:), allocatable :: diis 
!
      integer :: iteration, state
!
      character(len=3) :: string_state
!
      real(dp) :: norm_X
!
      real(dp), dimension(:), allocatable   :: eps
      real(dp), dimension(:,:), allocatable :: X, R
!
      real(dp) :: ddot
!
!     Prepare wavefunction for excited state, and possibly do Davidson preconvergence 
!
      call solver%prepare_wf_for_excited_state(wf)
!
      if (solver%davidson_preconvergence) then 
!
         call solver%do_davidson_preconvergence(wf)   ! Use Davidson to get first guesses 
         solver%restart = .true.                      ! Restart from these guesses
!
      endif
!
!     Initialize solver tools 
!
      call mem%alloc(eps, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(eps, wf%n_es_amplitudes)
!
      solver%preconditioner = precondition_tool(eps, wf%n_es_amplitudes)
!
      call mem%dealloc(eps, wf%n_es_amplitudes)
!
      call solver%initialize_start_vector_tool(wf)
      call solver%initialize_projection_tool(wf)
!
!     Initialize energies, residual norms, and convergence arrays 
!
      call mem%alloc(prev_energies, solver%n_singlet_states)
      call mem%alloc(residual_norms, solver%n_singlet_states)
!
      prev_energies     = zero 
      residual_norms    = zero 
!
      call mem%alloc(converged, (solver%n_singlet_states))
!
      converged = .false.
!
!     Make DIIS tools array & initialize the individual DIIS tools 
!
      allocate(diis(solver%n_singlet_states))
!
      do state = 1, solver%n_singlet_states
!
         write(string_state, '(i3.3)') state
         diis(state) = diis_tool('diis_cc_es_' // string_state,      &
                                 wf%n_es_amplitudes,                 &
                                 wf%n_es_amplitudes,                 &
                                 dimension_=solver%diis_dimension,   &
                                 crop=solver%crop)
!
         call diis(state)%initialize_storers(solver%records_in_memory)
!
      enddo 
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call solver%initialize_energies()
!
      call mem%alloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
!
      call solver%set_initial_guesses(wf, X, 1, wf%n_singlet_states)
!
!     Enter iterative loop
!
      call mem%alloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1   
!
         call output%printf('n', 'Iteration: (i18)', ints=[iteration], fs='(/t3,a)')
!
         call output%printf('n', ' Root     Eigenvalue (Re)     Residual norm', fs='(/t3,a)')
         call output%print_separator('n', 47, '-')
!
         do state = 1, solver%n_singlet_states
!
            if (.not. converged(state)) then 
!
!              Construct R = AX 
!
               call wf%construct_Jacobian_transform(solver%transformation, &
                                                      X(:,state),          &
                                                      R(:,state),          &
                                                      solver%energies(state))
!
!              Project if relevant (CVS, IP)
!
               if (solver%projector%active) call solver%projector%do_(R(:,state))
!
!              Calculate energy E = X^T R = X^T A X
!
               solver%energies(state) = ddot(wf%n_es_amplitudes, X(:,state), 1, R(:,state), 1)
!
!              Construct residual, R = A X - energy X, and calculate its norm           
!
               call daxpy(wf%n_es_amplitudes, -solver%energies(state), X(:,state), 1, R(:,state), 1)
!
               residual_norms(state) = get_l2_norm(R(:, state), wf%n_es_amplitudes)
!
!              Update convergence logicals 
!
               converged(state) = solver%convergence_checker%has_converged(residual_norms(state), &
                                          solver%energies(state)-prev_energies(state), iteration)
!
!              Perform DIIS extrapolation to the optimal next guess for X,
!              then normalize it to avoid accumulating norm in X
!
               if (.not. converged(state)) then
!
!                 Form quasi-Newton estimate for X 
!
                  call solver%preconditioner%do_(R(:,state),                     &
                                                 shift=solver%energies(state),   &
                                                 prefactor=-one)
!
                  call daxpy(wf%n_es_amplitudes, one, R(:, state), 1, X(:, state), 1)
!
!                 DIIS extrapolate using previous quasi-Newton estimates
!
                  call diis(state)%update(R(:,state), X(:,state))
!
!                 Renormalize X 
!
                  norm_X = get_l2_norm(X(:,state), wf%n_es_amplitudes)
                  call dscal(wf%n_es_amplitudes, one/norm_X, X(:, state), 1)
!
               endif 
!
            endif 
!
            call output%printf('n', '(i4)  (f18.12)      (e11.4)', &
                               ints=[state], reals=[solver%energies(state), &
                               residual_norms(state)])
!
         enddo
!
!        Save excited states and excitation energies
!
         call wf%save_excited_state(X,                       &
                                    1,                       &
                                    solver%n_singlet_states, &
                                    solver%transformation,   &
                                    solver%energies)
!
         prev_energies = solver%energies 
!
         call output%print_separator('n', 47, '-')
!
      enddo 
!
      if (all(converged)) then 
!
         if (iteration .eq. 1) then 
!
            call output%printf('n', 'Note: residual of all states converged  in &
                               &first iteration.', fs='(/t3,a)')
            call output%printf('n', 'Energy convergence has not been tested.', fs='(t3,a/)')
!
         endif
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
         call output%printf('n', '- Resorting roots according to excitation energy.', &
                            fs='(/t3,a)')
!
!        Sort roots and store the original state number in X_order
!
         call mem%alloc(X_order, solver%n_singlet_states)
!
         call quicksort_with_index_ascending(solver%energies, X_order, solver%n_singlet_states)
!
         do state = 1, solver%n_singlet_states
!
            if(X_order(state) .ne. state) then
!
               call output%printf('v', 'Root number (i0) renamed to state (i0)' &
                                  &, ints=[X_order(state), state], fs='(t5,a)')
!
            end if
!
            call wf%save_excited_state(X(:, X_order(state)), &
                                       state,                &
                                       state,                &
                                       solver%transformation,&
                                       solver%energies(state))
!
         enddo
!
         call output%printf('n', '- Stored converged states to file.', fs='(/t3,a)')
!
         call wf%set_excitation_energies(solver%energies, solver%transformation)
!
         call solver%print_summary(wf, X, X_order)
!
         call mem%dealloc(X_order, solver%n_singlet_states)
!
         call wf%set_excitation_energies(solver%energies, solver%transformation)
!
      else 
!
         call output%error_msg("Did not converge in the max number of iterations.")
!
      endif
!
      call mem%dealloc(prev_energies, solver%n_singlet_states)
      call mem%dealloc(residual_norms, solver%n_singlet_states)
!
      call mem%dealloc(converged, (solver%n_singlet_states))
!
      call mem%dealloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
      call mem%dealloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
!
      do state = 1, solver%n_singlet_states
!
         call diis(state)%finalize_storers()
!
      enddo 
!
      call solver%preconditioner%destruct_precondition_vector()
!
   end subroutine run_diis_cc_es
!
!
   subroutine do_davidson_preconvergence_diis_cc_es(solver, wf)
!!
!!    Do Davidson preconvergence 
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
!!    Sets up non-linear Davidson solver and runs it to get first guesses
!!    for the eigenstates. 
!!
      use nonlinear_davidson_cc_es_class, only: nonlinear_davidson_cc_es
!
      implicit none 
!
      class(diis_cc_es), intent(in) :: solver 
!
      class(ccs) :: wf 
!
      class(nonlinear_davidson_cc_es), allocatable :: davidson_solver 
!
      real(dp) :: relative_micro_residual_threshold
      integer  :: max_dim_red, max_micro_iterations
!
      call output%printf('m', 'Running the non-linear Davidson solver to produce&
                              & first guesses for the DIIS solver. When finished,&
                              & the DIIS solver will restart from the preconverged&
                              & solutions.', ffs='(/t3,a)')
!
!     Set some defaults 
!
      max_micro_iterations              = 100
      max_dim_red                       = 100
      relative_micro_residual_threshold = 1.0d-1
!
!     Read non-default values, if provided by user 
!
      call input%get_keyword('max reduced dimension',  &
                                        'solver cc es',           &
                                        max_dim_red)
!
      call input%get_keyword('max micro iterations',  &
                                        'solver cc es',           &
                                        max_micro_iterations)
!
      call input%get_keyword('rel micro threshold',    &
                                        'solver cc es',           &
                                        relative_micro_residual_threshold)
!
!     Allocate non-linear Davidson solver & run it
!
      davidson_solver = nonlinear_davidson_cc_es(solver%transformation,             &
                                                 wf,                                &
                                                 solver%restart,                    &  
                                                 solver%preconvergence_threshold,   &
                                                 solver%preconvergence_threshold,   &
                                                 solver%max_iterations,             &
                                                 relative_micro_residual_threshold, & 
                                                 max_micro_iterations,              & 
                                                 max_dim_red,                       & 
                                                 solver%es_type,                    &
                                                 solver%storage,                    &
                                                 solver%n_singlet_states,           &
                                                 prepare_wf=.false.,                &
                                                 energy_convergence=solver%convergence_checker%energy_convergence)
!
      call davidson_solver%run(wf)
!
      call davidson_solver%cleanup(wf)
!
      call output%printf('m', 'Finished preconvergence! The DIIS solver will now restart&
                              & from the preconverged solutions.', ffs='(/t3,a)')
!
   end subroutine do_davidson_preconvergence_diis_cc_es
!
!
end module diis_cc_es_class
