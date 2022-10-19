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
!!       A_mu,nu = < mu |[H-bar, tau_nu] | HF >,    H-bar = e-T H eT.
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
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use ccs_class, only: ccs
   use diis_tool_class, only: diis_tool
   use abstract_cc_es_class, only: abstract_cc_es
   use precondition_tool_class, only: precondition_tool
   use convergence_tool_class, only: convergence_tool
!
   use start_vector_tool_class, only: start_vector_tool
   use abstract_projection_tool_class, only : abstract_projection_tool
   use convergence_tool_class, only: convergence_tool
!
   implicit none
!
   type, extends(abstract_cc_es) :: diis_cc_es
!
      logical :: davidson_preconvergence     ! Perform non-linear Davidson first, then go over
                                             ! to DIIS to converge below the residual threshold
!
      real(dp) :: preconvergence_threshold   ! Non-linear Davidson threshold
!
      type(diis_tool), dimension(:), allocatable :: diis
!
      character(len=:), allocatable :: es_type
!
   contains
!
      procedure, non_overridable :: run => run_diis_cc_es
!
      procedure, private :: read_settings
      procedure, private :: print_settings
      procedure, private :: do_davidson_preconvergence
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
   function new_diis_cc_es(transformation, wf, start_vector, projector, &
                           convergence_checker, n_states, max_iterations, &
                           records_in_memory, es_type) result(this)
!!
!!    New DIIS CC ES
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      type(diis_cc_es) :: this
      class(ccs), intent(inout), target :: wf
!
      character(len=*), intent(in) :: transformation, es_type
!
      class(start_vector_tool), intent(in) :: start_vector
      class(abstract_projection_tool), intent(in) :: projector
      type(convergence_tool), intent(in) :: convergence_checker
!
      integer, intent(in) :: n_states, max_iterations
      logical, intent(in) :: records_in_memory
!
      character(len=3) :: string_state
      logical :: crop
      integer :: diis_dimension, state
!
      this%wf => wf
!
      this%iteration_timer = timings('DIIS CC ES solver iteration')
      this%total_timer = timings('DIIS CC ES solver (' // &
                                 trim(transformation) //')')
      call this%total_timer%turn_on()
!
!     Set printables
!
      this%name_  = 'DIIS coupled cluster excited state solver'
      this%tag    = 'DIIS'
!
      this%description1 = 'A DIIS solver that solves for the lowest eigenvalues and &
                           &the eigenvectors of the Jacobian matrix, A. The eigenvalue &
                           &problem is solved by DIIS extrapolation of residuals for each &
                           &eigenvector until the convergence criteria are met.'
!
      this%description2 = 'More on the DIIS algorithm can be found in &
                           &P. Pulay, Chemical Physics Letters, 73(2), 393-398 (1980).'
!
      call this%print_banner()
!
      this%transformation   = trim(transformation)
      this%es_type          = trim(es_type)
      this%n_singlet_states = n_states
      this%max_iterations   = max_iterations
!
      this%start_vectors = start_vector
      this%projector = projector
      this%convergence_checker = convergence_checker
!
      this%davidson_preconvergence  = .false.
      this%preconvergence_threshold = 1.0d-2
!
      crop           = .false.
      diis_dimension = 20
!
      call this%read_settings(crop, diis_dimension)
!
      call this%print_settings()
!
      if (this%n_singlet_states == 0) then
         call output%error_msg('number of excitations must be specified.')
      end if
!
      this%wf%n_singlet_states = this%n_singlet_states
!
!     Make DIIS tools array & initialize the individual DIIS tools
!
      allocate(this%diis(this%n_singlet_states))
!
      do state = 1, this%n_singlet_states
!
         write(string_state, '(i3.3)') state
         this%diis(state) = diis_tool('diis_cc_es_' // string_state,    &
                                       this%wf%n_es_amplitudes,         &
                                       this%wf%n_es_amplitudes,         &
                                       dimension_=diis_dimension,       &
                                       crop=crop,                       &
                                       records_in_memory=records_in_memory)
!
      enddo
!
!
   end function new_diis_cc_es
!
!
   subroutine print_settings(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(diis_cc_es) :: this
!
      call this%print_es_settings()
!
      if (this%davidson_preconvergence) then
!
         call output%printf('m', 'Enabled preconvergence using&
                                 & the non-linear Davidson solver.', fs='(t6,a)')
!
         call output%printf('m', 'Preconvergence threshold:   (e11.2)', &
                                 reals=[this%preconvergence_threshold], fs='(/t6,a/)')
!
      endif
!
   end subroutine print_settings
!
!
   subroutine read_settings(this, crop, diis_dimension)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      use global_in, only: input
!
      implicit none
!
      class(diis_cc_es)           :: this
      logical, intent(inout)      :: crop
      integer, intent(inout)      :: diis_dimension
!
      call input%get_keyword('diis dimension', 'solver cc es', diis_dimension)
      crop = input%is_keyword_present('crop', 'solver cc es')
!
      if (input%is_keyword_present('davidson preconvergence', 'solver cc es')) then
!
         this%davidson_preconvergence = .true.
!
      endif
!
      call input%get_keyword('preconvergence threshold',  &
                                        'solver cc es',              &
                                        this%preconvergence_threshold)
!
   end subroutine read_settings
!
!
   subroutine run_diis_cc_es(this)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use array_utilities, only: quicksort_with_index_ascending, get_l2_norm
!
      implicit none
!
      class(diis_cc_es), intent(inout) :: this
!
      logical, dimension(:), allocatable :: converged
!
      real(dp), dimension(:), allocatable :: prev_energies
      real(dp), dimension(:), allocatable :: residual_norms
!
      integer, dimension(:), allocatable :: X_order
!
      integer :: iteration, state
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
      call this%prepare_wf_for_excited_state()
!
      if (this%davidson_preconvergence) call this%do_davidson_preconvergence()   ! Use Davidson to get first guesses
!
!     Initialize solver tools
!
      do state = 1, this%n_singlet_states
!
         call this%diis(state)%initialize()
!
      enddo
!
      call mem%alloc(eps, this%wf%n_es_amplitudes)
      call this%wf%get_orbital_differences(eps, this%wf%n_es_amplitudes)
!
      this%preconditioner = precondition_tool(this%wf%n_es_amplitudes)
      call this%preconditioner%initialize_and_set_precondition_vector(eps)
!
      call mem%dealloc(eps, this%wf%n_es_amplitudes)
!
!     Initialize energies, residual norms, and convergence arrays
!
      call mem%alloc(prev_energies, this%n_singlet_states)
      call mem%alloc(residual_norms, this%n_singlet_states)
!
      prev_energies     = zero
      residual_norms    = zero
!
      call mem%alloc(converged, this%n_singlet_states, set_to=.false.)
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call this%initialize_energies()
!
      call mem%alloc(X, this%wf%n_es_amplitudes, this%n_singlet_states)
!
      do state = 1, this%n_singlet_states
         call this%start_vectors%get(X(:,state), state, this%energies(state))
      end do
!
!     Enter iterative loop
!
      call mem%alloc(R, this%wf%n_es_amplitudes, this%n_singlet_states)
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. this%max_iterations))
!
         call this%iteration_timer%turn_on()
!
         iteration = iteration + 1
!
         call output%printf('n', 'Iteration: (i18)', ints=[iteration], fs='(/t3,a)')
!
         call output%printf('n', ' Root     Eigenvalue (Re)     Residual norm', fs='(/t3,a)')
         call output%print_separator('n', 47, '-')
!
         do state = 1, this%n_singlet_states
!
            if (.not. converged(state)) then
!
!              Construct R = AX
!
               call this%wf%construct_Jacobian_transform(this%transformation, &
                                                      X(:,state),          &
                                                      R(:,state),          &
                                                      this%energies(state))
!
!              Project if relevant (CVS, IP)
!
               call this%projector%project(R(:,state))
!
!              Calculate energy E = X^T R = X^T A X
!
               this%energies(state) = ddot(this%wf%n_es_amplitudes, X(:,state), 1, R(:,state), 1)
!
!              Construct residual, R = A X - energy X, and calculate its norm
!
               call daxpy(this%wf%n_es_amplitudes, -this%energies(state), &
                          X(:,state), 1, R(:,state), 1)
!
               residual_norms(state) = get_l2_norm(R(:, state), this%wf%n_es_amplitudes)
!
!              Update convergence logicals
!
               converged(state) = this%convergence_checker%has_converged(residual_norms(state), &
                                          this%energies(state)-prev_energies(state), iteration)
!
!              Perform DIIS extrapolation to the optimal next guess for X,
!              then normalize it to avoid accumulating norm in X
!
               if (.not. converged(state)) then
!
!                 Form quasi-Newton estimate for X
!
                  call this%preconditioner%do_(R(:,state),                     &
                                                 shift=this%energies(state),   &
                                                 prefactor=-one)
!
                  call daxpy(this%wf%n_es_amplitudes, one, R(:, state), 1, X(:, state), 1)
!
!                 DIIS extrapolate using previous quasi-Newton estimates
!
                  call this%diis(state)%update(R(:,state), X(:,state))
!
!                 Renormalize X
!
                  norm_X = get_l2_norm(X(:,state), this%wf%n_es_amplitudes)
                  call dscal(this%wf%n_es_amplitudes, one/norm_X, X(:, state), 1)
!
               endif
!
            endif
!
            call output%printf('n', '(i4)  (f18.12)      (e11.4)', &
                               ints=[state], reals=[this%energies(state), &
                               residual_norms(state)])
!
         enddo
!
!        Save excited states and excitation energies
!
         call this%wf%save_excited_state(X,                       &
                                    1,                       &
                                    this%n_singlet_states, &
                                    this%transformation,   &
                                    this%energies)
!
         prev_energies = this%energies
!
         call output%print_separator('n', 47, '-')
!
         call this%iteration_timer%turn_off()
         call this%iteration_timer%reset()
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
         call mem%alloc(X_order, this%n_singlet_states)
!
         call quicksort_with_index_ascending(this%energies, X_order, this%n_singlet_states)
!
         do state = 1, this%n_singlet_states
!
            if(X_order(state) .ne. state) then
!
               call output%printf('v', 'Root number (i0) renamed to state (i0)' &
                                  &, ints=[X_order(state), state], fs='(t5,a)')
!
            end if
!
            call this%wf%save_excited_state(X(:, X_order(state)), &
                                       state,                &
                                       state,                &
                                       this%transformation,&
                                       this%energies(state))
!
         enddo
!
         call this%wf%set_excitation_energies(this%energies, this%transformation)
!
         call mem%dealloc(X_order, this%n_singlet_states)
!
         call this%wf%set_excitation_energies(this%energies, this%transformation)
!
      else
!
         call output%error_msg("Did not converge in the max number of iterations.")
!
      endif
!
      call mem%dealloc(prev_energies, this%n_singlet_states)
      call mem%dealloc(residual_norms, this%n_singlet_states)
!
      call mem%dealloc(converged, (this%n_singlet_states))
!
      call mem%dealloc(X, this%wf%n_es_amplitudes, this%n_singlet_states)
      call mem%dealloc(R, this%wf%n_es_amplitudes, this%n_singlet_states)
!
      do state = 1, this%n_singlet_states
!
         call this%diis(state)%finalize()
!
      enddo
!
      call this%preconditioner%destruct_precondition_vector()
!
   end subroutine run_diis_cc_es
!
!
   subroutine do_davidson_preconvergence(this)
!!
!!    Do Davidson preconvergence
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
!!    Sets up non-linear Davidson solver and runs it to get first guesses
!!    for the eigenstates.
!!
      use global_in, only: input
      use nonlinear_davidson_cc_es_class, only: nonlinear_davidson_cc_es
      use cc_es_start_vector_factory_class, only: cc_es_start_vector_factory
!
      implicit none
!
      class(diis_cc_es), intent(inout) :: this
!
      class(nonlinear_davidson_cc_es), allocatable :: davidson_solver
!
      type(convergence_tool) :: preconvergence_checker
      type(cc_es_start_vector_factory) :: start_vector_factory
!
      real(dp) :: relative_micro_residual_threshold
      integer  :: max_dim_red, max_micro_iterations
!
      logical :: records_in_memory, energy_convergence
!
      call output%printf('m', 'Running the non-linear Davidson solver to produce&
                              & first guesses for the DIIS solver. When finished,&
                              & the DIIS solver will restart from the preconverged&
                              & solutions.', ffs='(/t3,a)')
!
      max_micro_iterations              = 100
      max_dim_red                       = 100
      relative_micro_residual_threshold = 1.0d-1
!
      call input%get_keyword('max reduced dimension',  &
                             'solver cc es',&
                             max_dim_red)
!
      call input%get_keyword('max micro iterations',   &
                             'solver cc es',&
                             max_micro_iterations)
!
      call input%get_keyword('rel micro threshold',    &
                             'solver cc es',&
                             relative_micro_residual_threshold)
!
      energy_convergence = input%is_keyword_present('energy threshold', 'solver cc es')
!
      records_in_memory = .false.
      call input%place_records_in_memory('solver cc es', records_in_memory)
!
      preconvergence_checker = convergence_tool(this%preconvergence_threshold, &
                                                this%preconvergence_threshold, &
                                                energy_convergence)
!
      davidson_solver = nonlinear_davidson_cc_es(this%transformation,               &
                                                 this%wf,                           &
                                                 this%start_vectors,                &
                                                 this%projector,                    &
                                                 preconvergence_checker,            &
                                                 this%n_singlet_states,             &
                                                 this%max_iterations,               &
                                                 records_in_memory,                 &
                                                 relative_micro_residual_threshold, &
                                                 max_micro_iterations,              &
                                                 max_dim_red,                       &
                                                 prepare_wf=.false.)
!
      call davidson_solver%run()
      call davidson_solver%cleanup()
!
      call this%wf%print_es_summary(this%transformation, 'singlet')
!
      call output%printf('m', 'Finished preconvergence! The DIIS solver will now restart&
                              & from the preconverged solutions.', ffs='(/t3,a)')
!
!     Overwrite start vector tool to restart
!
      deallocate(this%start_vectors)
      start_vector_factory = cc_es_start_vector_factory(this%es_type, &
                                                        this%transformation, &
                                                        restart=.true.)
      call start_vector_factory%create(this%wf, this%start_vectors)
!
   end subroutine do_davidson_preconvergence
!
!
end module diis_cc_es_class
