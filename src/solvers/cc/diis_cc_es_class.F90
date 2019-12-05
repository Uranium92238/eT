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
   use precondition_tool_class, only: precondition_tool
   use es_valence_start_vector_tool_class, only: es_valence_start_vector_tool
   use es_valence_projection_tool_class, only: es_valence_projection_tool
   use array_utilities, only: quicksort_with_index_ascending, get_l2_norm
   use string_utilities, only: convert_to_uppercase
!
   implicit none
!
   type, extends(abstract_cc_es) :: diis_cc_es
!
      integer :: diis_dimension
!
      class(precondition_tool), allocatable :: preconditioner 
!
      logical :: crop ! Standard DIIS if false; CROP variant of DIIS if true
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
      real(dp), dimension(:), allocatable :: eps 
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
      solver%n_singlet_states     = 0
      solver%max_iterations       = 100
      solver%eigenvalue_threshold = 1.0d-6
      solver%residual_threshold   = 1.0d-6
      solver%transformation       = 'right'
      solver%diis_dimension       = 20
      solver%restart              = restart
      solver%transformation       = trim(transformation)
      solver%es_type              = 'valence'
      solver%records_in_memory    = .false.
      solver%storage              = 'disk'
      solver%crop                 = .false.
!
      call solver%read_settings()
      call solver%print_settings()
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      call solver%initialize_energies()
      solver%energies = zero
!
      wf%n_singlet_states = solver%n_singlet_states
!
      call solver%initialize_start_vector_tool(wf)
      call solver%initialize_projection_tool(wf)
!
      call solver%prepare_wf_for_excited_state(wf)

      if (wf%frozen_core .and. solver%es_type=='core') call output%error_msg('No support for frozen core with CVS yet.')
!
!     Initialize preconditioner 
!
      call mem%alloc(eps, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(eps, wf%n_es_amplitudes)
!
      solver%preconditioner = precondition_tool(eps, wf%n_es_amplitudes)
!
      call mem%dealloc(eps, wf%n_es_amplitudes)
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
      call output%printf('DIIS dimension:               (i11)', pl='minimal', &
                         ints=[solver%diis_dimension], fs='(/t6,a)')
!
      if (solver%crop) then 
!
         call output%printf('Enabled CROP in the DIIS algorithm.', pl='minimal', fs='(/t6,a)')
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
      call input%get_keyword_in_section('diis dimension', 'solver cc es', solver%diis_dimension)
!
      if (input%requested_keyword_in_section('crop', 'solver cc es')) then 
!
         solver%crop = .true.
!
      endif
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
      logical, dimension(:), allocatable :: converged_eigenvalue
      logical, dimension(:), allocatable :: converged_residual
!
      real(dp), dimension(:), allocatable :: prev_energies
      real(dp), dimension(:), allocatable :: residual_norms 
!
      integer, dimension(:), allocatable :: prev_state_numbers
!
      type(diis_tool), dimension(:), allocatable :: diis 
!
      integer :: iteration, state, n_solutions_on_file
!
      character(len=3) :: string_state
!
      real(dp) :: norm_X
!
      real(dp), dimension(:), allocatable   :: eps
      real(dp), dimension(:,:), allocatable :: X, R
      real(dp), dimension(:,:), allocatable :: X_sorted
!
      real(dp) :: ddot
!
!     Initialize energies, residual norms, and convergence arrays 
!
      call mem%alloc(prev_energies, solver%n_singlet_states)
      call mem%alloc(residual_norms, solver%n_singlet_states)
!
      prev_energies     = zero 
      residual_norms    = zero 
!
      allocate(converged(solver%n_singlet_states))
      allocate(converged_residual(solver%n_singlet_states))
      allocate(converged_eigenvalue(solver%n_singlet_states))
!
      converged            = .false.
      converged_residual   = .false.
      converged_eigenvalue = .false.
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
      call mem%alloc(eps, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(eps, wf%n_es_amplitudes)
!
      call mem%alloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
!
      do state = 1, solver%n_singlet_states
!
         call solver%start_vectors%get(X(:,state), state)
!
      enddo 
!
      if (solver%restart) then ! Overwrite all or some of the orbital differences 
!
         call solver%determine_restart_transformation(wf) ! Read right or left?
!
         n_solutions_on_file = wf%get_n_excited_states_on_file(solver%restart_transformation)
!
         call output%printf('Requested restart - there are (i0) (a0) eigenvectors on file.', &
                              ints=[n_solutions_on_file], chars=[solver%restart_transformation])
!
         do state = 1, n_solutions_on_file
!
            call wf%read_excited_state(X(:,state), state, solver%restart_transformation)
!
!           Normalize X in case it has been changed by FOP
!
            norm_X = get_l2_norm(X(:,state), wf%n_es_amplitudes)
            call dscal(wf%n_es_amplitudes, one/norm_X, X(:, state), 1)
!
         enddo
!
         solver%energies = zero
         call wf%read_excitation_energies(n_solutions_on_file, solver%energies(1:n_solutions_on_file))
!
      endif
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
         call output%printf('Iteration: (i18)', &
                     ints=[iteration], fs='(/t3,a)', pl='n')
!
         call output%printf('Root     Eigenvalue (Re)     Residual norm', &
                            fs='(/t3,a)', pl='n')
         call output%print_separator('n', 46, '-')
!
         do state = 1, solver%n_singlet_states
!
            if (.not. converged(state)) then 
!
!              Construct R = AX 
!
               call dcopy(wf%n_es_amplitudes, X(:,state), 1, R(:,state), 1)
!
               call wf%construct_Jacobian_transform(solver%transformation, &
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
               converged_eigenvalue(state) = abs(solver%energies(state)-prev_energies(state)) &
                                                      .lt. solver%eigenvalue_threshold
!
               converged_residual(state)   = residual_norms(state) .lt. solver%residual_threshold
!
               converged(state) = converged_eigenvalue(state) .and. converged_residual(state)
!
!              Perform DIIS extrapolation to the optimal next guess for X,
!              then normalize it to avoid accumulating norm in X
!
               if (converged_residual(state) .and. iteration .eq. 1) then 
!
                  converged(state) = .true.
!
               endif
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
            call output%printf('(i0)   (f19.12)      (e11.4)', ints=[state], &
                 reals=[solver%energies(state), residual_norms(state)], pl='n')
!
         enddo
!
!        Save excited states and excitation energies
!
         do state = 1, solver%n_singlet_states
!
            call wf%save_excited_state(X(:,state), state, solver%transformation)
!
         enddo 
!
         call wf%save_excitation_energies(solver%n_singlet_states, &
                           solver%energies, solver%transformation)
         prev_energies = solver%energies 
!
         call output%print_separator('n', 46, '-')
!
      enddo 
!
      if (all(converged)) then 
!
         if (iteration .eq. 1) then 
!
            call output%printf('Note: residual of all states converged &
                              & in first iteration.', fs='(/t3,a)',pl='n')
            call output%printf('Energy convergence has not been tested.', &
                                fs='(t3,a/)',pl='n')
!
         endif
!
         call output%printf('Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)',pl='m')
         call output%printf('- Resorting roots according to excitation energy.', &
                            fs='(/t3,a)', pl='n')
!
!        Sort roots and store the original state number in prev_state_numbers
!
         call mem%alloc(prev_state_numbers, solver%n_singlet_states)
!
         call quicksort_with_index_ascending(solver%energies, prev_state_numbers, solver%n_singlet_states)
!
         call mem%alloc(X_sorted, wf%n_es_amplitudes, solver%n_singlet_states)
!
         do state = 1, solver%n_singlet_states
!
!           Need X_sorted because print_summary needs the vectors in the correct order
            call dcopy(wf%n_es_amplitudes, X(:, prev_state_numbers(state)), 1, X_sorted(:, state), 1)
!
            if(prev_state_numbers(state) .ne. state) then
!
               call output%printf('Root number (i3) renamed to state (i3)',      &
                     ints=[prev_state_numbers(state), state], fs='(t5,a)', pl='v')
!
            end if
!
            call wf%save_excited_state(X_sorted(:, state), state, solver%transformation)
!
         enddo
!
         call mem%dealloc(prev_state_numbers, solver%n_singlet_states)
!
         call output%printf('- Stored converged states to file.', fs='(/t3,a)', pl='n')
!
         call solver%print_summary(wf, X_sorted)
!
         call mem%dealloc(X_sorted, wf%n_es_amplitudes, solver%n_singlet_states)
!
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies, solver%transformation)
!
         call wf%check_for_degeneracies(solver%transformation, solver%residual_threshold)
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
      deallocate(converged)
      deallocate(converged_residual)
      deallocate(converged_eigenvalue)
!
      call mem%dealloc(eps, wf%n_es_amplitudes)
      call mem%dealloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
      call mem%dealloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
!
      do state = 1, solver%n_singlet_states
!
         call diis(state)%finalize_storers()
!
      enddo 
!
   end subroutine run_diis_cc_es
!
!
end module diis_cc_es_class
