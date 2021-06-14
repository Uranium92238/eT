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
module response_engine_class
!!
!!    Response coupled cluster engine class module
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Josefine H. Andersen, Apr 2019
!!
   use parameters
   use global_in,             only: input
   use global_out,            only: output
   use timings_class,         only: timings
   use memory_manager_class,  only: mem
   use sequential_file_class, only: sequential_file
   use task_list_class,       only: task_list
!
   use es_engine_class, only: es_engine
   use ccs_class,       only: ccs
!
   use davidson_cc_linear_equations_class, only: davidson_cc_linear_equations
!
   use array_utilities, only: copy_and_scale
!
   type, extends(es_engine) :: response_engine
!
!     Run equation-of-motion or linear response?
!
      logical :: eom, lr
!
!     Restart left or right vectors?
!
      logical :: es_restart_left
      logical :: es_restart_right
!
!     Do transition moments, permanent moments and/or polarizabilities?
!
      logical :: transition_moments, permanent_moments
      logical :: polarizabilities
!
      logical, dimension(3,3) :: compute_polarizability  ! which components to compute?
                                                         ! xx, xy, xz; yx, yy, yz; ...
!
      integer :: n_frequencies
      real(dp), dimension(:), allocatable :: frequencies ! array of frequencies for which
                                                         ! to evaluate the polarizability
!
      logical, dimension(3) :: t_response_components     ! compute amplitude response?
                                                         ! x, y, z; determined by true/false
                                                         ! in compute_polarizability
!
!     etaX_mu = < HF | e-T X eT | mu >
!     xiX_mu = < mu | e-T X eT | HF >
!
      real(dp), dimension(:,:), allocatable :: etaX ! columns correspond to components of X
      real(dp), dimension(:,:), allocatable :: xiX  ! columns correspond to components of X
!
!     File arrays to store vectors needed for linear response
!
      type(sequential_file), dimension(:), allocatable     :: M_vectors   ! used for transition moments
      type(sequential_file), dimension(:,:,:), allocatable :: t_responses ! used for polarizabilities
!
!     Operators
!
      logical :: dipole_length
!
!     Plotting
!
      logical :: plot_mn_densities, plot_es_densities
      integer :: n_states_to_plot
      integer, dimension(:), allocatable :: states_to_plot
!
      integer :: n_initial_states ! Number of initial states used in transition strength
!                                 ! calculations.
                                  ! Default is n=1 (calculate transitions from the GS)
!
      integer, dimension(:), allocatable :: initial_states ! The specified initial states for
!                                                          ! transition strength calculations.
!                                                          ! Default is initial_states = [0],
!                                                          ! meaning transition moments
!                                                          ! from the ground state
!
   contains
!
      procedure :: run => run_response_engine
!
!     Read settings
!
      procedure :: read_settings &
                => read_settings_response_engine
!
      procedure :: read_response_settings &
                => read_response_settings_response_engine
!
      procedure :: set_initial_states &
                => set_initial_states_response_engine
!
      procedure :: set_states_to_plot &
                => set_states_to_plot_response_engine
!
      procedure :: set_polarizability_components &
                => set_polarizability_components_response_engine
!
      procedure :: set_t_response_components &
                => set_t_response_components_response_engine
!
      procedure :: set_printables &
                => set_printables_response_engine
!
!     Compute properties
!
      procedure :: do_eom_transition_moments &
                => do_eom_transition_moments_response_engine
      procedure :: do_lr_transition_moments &
                => do_lr_transition_moments_response_engine
!
      procedure :: construct_etaX_and_xiX &
                => construct_etaX_and_xiX_response_engine
      procedure :: determine_M_vectors &
                => determine_M_vectors_response_engine
      procedure :: determine_amplitude_response &
                => determine_amplitude_response_response_engine
      procedure :: calculate_polarizabilities &
                => calculate_polarizabilities_response_engine
!
!     Summaries
!
      procedure :: print_lr_transition_moment_summary  &
                => print_lr_transition_moment_summary_response_engine
!
      procedure :: print_eom_transition_moment_summary &
                => print_eom_transition_moment_summary_response_engine
!
      procedure :: print_permanent_moments_summary &
                => print_permanent_moments_summary_response_engine
!
      procedure :: print_operator_per_state &
                => print_operator_per_state_response_engine
!
!     Visualization
!
      procedure :: visualize_cc_densities => visualize_cc_densities_response_engine
!
   end type response_engine
!
!
   interface response_engine
!
      procedure :: new_response_engine
!
   end interface response_engine
!
!
contains
!
!
   function new_response_engine(wf) result(engine)
!!
!!    New response engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
!     Needed for defaults and sanity checks
      class(ccs), intent(in)       :: wf
!
      type(response_engine) :: engine
!
!     Set standards and then read if nonstandard
!
      engine%gs_algorithm = 'diis'
!
      if (wf%name_ .eq. 'ccsd(t)'   .or. &
          wf%name_ .eq. 'mp2'       .or. &
          wf%name_ .eq. 'mlcc2'     .or. &
          wf%name_ .eq. 'mlccsd'    .or. &
          wf%name_ .eq. 'low memory cc2') then
!
         call output%error_msg("Response properties not implemented for (a0)", &
                               chars=[wf%name_])
!
      end if
!
      if (wf%name_ .eq. 'cc3' .or. &
          wf%name_ .eq. 'low memory cc2') then
!
         engine%multipliers_algorithm = 'diis'
         engine%es_algorithm          = 'non-linear davidson'
!
      else if (wf%name_ .eq. 'ccs' .or. &
               wf%name_ .eq. 'cc2' .or. &
               wf%name_ .eq. 'mlcc2') then
!
         engine%multipliers_algorithm = 'diis'
         engine%es_algorithm          = 'davidson'
!
      else
!
         engine%multipliers_algorithm = 'davidson'
         engine%es_algorithm          = 'davidson'
!
      end if
!
      engine%es_type                = 'valence'
      engine%lr                     = .false.
      engine%eom                    = .false.
      engine%polarizabilities       = .false.
      engine%transition_moments     = .false.
      engine%permanent_moments      = .false.
      engine%compute_polarizability = .false.
      engine%t_response_components  = .false.
!
      engine%gs_restart             = .false.
      engine%multipliers_restart    = .false.
      engine%es_restart             = .false.
      engine%es_restart_left        = .false.
      engine%es_restart_right       = .false.
!
      engine%dipole_length          = .false.
!
      engine%plot_density           = .false.
      engine%plot_mn_densities      = .false.
      engine%plot_es_densities      = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_response_engine
!
!
   subroutine run_response_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(response_engine) :: engine
      class(ccs) :: wf
!
      real(dp) :: energy_threshold, residual_threshold
!
      if (trim(wf%name_) == 'low memory cc2') then
!
         call output%printf('m', 'Requested response calculation with lowmem-cc2.', &
                            fs='(/t3,a)')
         call output%error_msg('Response not implemented for low memory cc2.')
!
      end if
!
!     Determine ground state | CC >
!
      call engine%gs_engine%run(wf)
!
!     Determine multipliers < Λ |
!
      call engine%do_multipliers(wf)
!
!     Construct etaX and xiX if they are going to be needed
!
      if (engine%polarizabilities .or. engine%lr) then
!
         call mem%alloc(engine%xiX, wf%n_es_amplitudes, 3)
         call mem%alloc(engine%etaX, wf%n_es_amplitudes, 3)
!
         call wf%prepare_for_properties
         call engine%construct_etaX_and_xiX(wf) ! etaX_mu = < Lambda | [X-bar, τ_mu] | HF >
                                                ! xiX_mu = < mu | X-bar | HF >
!
      endif
!
!     Calculate transition moments or permanent moments
!
      if (engine%transition_moments .or. engine%permanent_moments) then
!
!        Save intermediates from solving the mutlipliers
!        needed for the transition densities
!
         call wf%save_tbar_intermediates()
!
!        Excited state solutions, left & right
!
         call engine%do_excited_state(wf,                                                       &
                                      transformation='right',                                   &
                                      restart=(engine%es_restart .or. engine%es_restart_right))
!
         call engine%do_excited_state(wf,                                                       &
                                      transformation='left',                                    &
                                      restart = .true.)
!
         call engine%get_thresholds(energy_threshold, residual_threshold)
!
         if (engine%es_algorithm .eq. 'diis') then
            call wf%remove_parallel_states(residual_threshold, 'both')
         end if
!
!        Biorthornormalize and look for duplicate states
!
         call wf%biorthonormalize_L_and_R(energy_threshold, residual_threshold)
!
!        Compute transition moments
!
         if (engine%eom) then
!
            call engine%do_eom_transition_moments(wf)
!
         elseif (engine%lr) then
!
            if (engine%permanent_moments .and. .not. engine%transition_moments) then
               call output%error_msg('Only transition moments implemented for LR.')
            end if
!
            call engine%do_lr_transition_moments(wf)
!
         endif
!
      endif
!
!     Calculate polarizabilities
!
      if (engine%polarizabilities) then
!
         if (.not. engine%transition_moments) call wf%prepare_for_jacobian()
!
         call engine%determine_amplitude_response(wf)
         call engine%calculate_polarizabilities(wf)
!
      endif
!
!     Deallocate xi and eta if they were constructed
      if (engine%polarizabilities .or. &
            (engine%lr .and. engine%transition_moments)) then
!
         call mem%dealloc(engine%xiX, wf%n_es_amplitudes, 3)
         call mem%dealloc(engine%etaX, wf%n_es_amplitudes, 3)
!
      endif
!
      if (engine%polarizabilities) then
         call mem%dealloc(engine%frequencies, engine%n_frequencies)
      end if
!
   end subroutine run_response_engine
!
!
   subroutine construct_etaX_and_xiX_response_engine(engine, wf)
!!
!!    Construct etaX and xiX
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
!!    Allocates and constructs the vectors
!!
!!       etaX_mu = < Lambda | [X-bar, τ_mu] | HF > (+ EOM correction)
!!       xiX_mu = < mu | X-bar | HF >,
!!
!!    where
!!
!!       X-bar = e-T X eT.
!!
!!    These are kept as engine variables afterwards.
!!
      implicit none
!
      class(response_engine) :: engine
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:), allocatable :: X
!
      integer :: k
!
!     Compute the operator X = mu, the dipole operator
!
      call mem%alloc(X, wf%n_mo, wf%n_mo, 3)
      call wf%get_t1_oei('dipole', X)
!
!     Construct the right-hand-side vector for operator X,
!     i.e. xiX_mu = < mu | X-bar | HF >.
!
      do k = 1, 3
!
         call wf%construct_xiX(X(:,:,k), engine%xiX(:,k))
!
      enddo
!
!     Construct the left-hand-side vector for operator X,
!     i.e. etaX_nu = < Lambda | [X-bar, τ_nu] | HF >
!
      do k = 1, 3
!
         if (engine%eom) then
!
            call wf%construct_eom_etaX(X(:,:,k), engine%xiX(:,k), engine%etaX(:,k))
!
         else ! LR
!
            call wf%construct_etaX(X(:,:,k), engine%etaX(:,k))
!
         endif
!
      enddo
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo, 3)
!
   end subroutine construct_etaX_and_xiX_response_engine
!
!
   subroutine calculate_polarizabilities_response_engine(engine, wf)
!!
!!    Calculate polarizabilities
!!    Written by Josefine H. Andersen, spring 2019
!!
!!    Calculates requested polarizabilities,
!!
!!       << X, Y >>(omega) = 1/2*[ eta^X (t^Y(omega) + t^Y(-omega))
!!                               + eta^Y (t^X(omega) + t^X(-omega))
!!                               + t^X(-omega) F t^Y(omega)
!!                               + t^X(omega) F t^Y(-omega) ], (*)
!!
!!    for the set of requested frequencies omega. Note that only components
!!    requested by the user are computed; see keyword "polarizabilities".
!!
!!    This routine is called after the amplitude response
!!    vectors (t^X, t^Y)  needed for the polarizability components
!!    have been converged. These are solutions to the equations
!!
!!       (A - omega_k I) t^X(omega_k) = -xi^X
!!
!!    Currently, X and Y are components of the dipole operator.
!!
!!    The terms involving the F matrix (terms 3 and 4 in (*)) are only
!!    added in LR theory. These are thus ignored if EOM polarizabilities
!!    are requested.
!!
!!    New storage of t response and adapted to late-2019
!!    program structure by Eirik F. Kjønstad, Nov 2019.
!!
!!    Revision of cross terms June 2021 by Anna Kristina Schnack-Petersen
!!
      implicit none
!
      class(response_engine) :: engine
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: tk, tl ! holds amplitude response, components k,l
                                                      ! columns: [+, -] (positive/negative frequency)
!
      real(dp), dimension(:,:), allocatable :: Ftk ! holds F tk
                                                   ! columns: [+ -]
!
      integer :: k, l, freq, sign_
!
      real(dp) :: ddot, polarizability
!
      character(len=4), dimension(3) :: operator_ = ['mu_x', 'mu_y', 'mu_z']
!
      call engine%tasks%print_('polarizabilities')
      call output%printf('m', 'The convention applied here defines the polarizabilities as &
                              &the response functions, without negative sign.', fs='(t6,a)')
!
!     Allocate arrays to hold amplitude response vectors as well as
!     F-transformed response vectors (if LR)
!
      call mem%alloc(tk, wf%n_es_amplitudes, 2)
      call mem%alloc(tl, wf%n_es_amplitudes, 2)
!
      if (engine%lr) then
!
         call mem%alloc(Ftk, wf%n_es_amplitudes, 2)
!
      endif
!
!     Compute polarizabilities
!
      do freq = 1, engine%n_frequencies
!
         do k = 1, 3
!
!           Check if we need k
!
            if (any(engine%compute_polarizability(k,:))) then
!
               do sign_ = 1, 2 ! +, -
!
                  call engine%t_responses(freq, k, sign_)%open_('read', 'rewind')
                  call engine%t_responses(freq, k, sign_)%read_(tk(:,sign_), wf%n_es_amplitudes)
                  call engine%t_responses(freq, k, sign_)%close_()
!
               enddo
!
            endif
!
!
!           Compute & print diagonal polarizabilities
!
            if (engine%compute_polarizability(k,k)) then
!
               polarizability = zero
!
               do sign_ = 1, 2 ! +, -
!
                  polarizability = polarizability + &
                         ddot(wf%n_es_amplitudes, engine%etaX(:,k), 1, tk(:,sign_), 1)
!
               enddo
!
!
               if (engine%lr) then
!
                  call dcopy(wf%n_es_amplitudes, tk(:,1), 1, Ftk(:,1), 1)
                  call dcopy(wf%n_es_amplitudes, tk(:,2), 1, Ftk(:,2), 1)
!
                  call wf%F_transformation(Ftk(:,1))
                  call wf%F_transformation(Ftk(:,2))
!
                  polarizability = polarizability + &
                                        ddot(wf%n_es_amplitudes, Ftk(:,2), 1, tk(:,1), 1)
!
               endif
!
               call output%printf('m', '<< ' // operator_(k) // ', ' //  &
                                  operator_(k) // ' >>' // '((e8.2)): (f19.12)', &
                                  reals=[engine%frequencies(freq), &
                                  polarizability], fs='(t6,a)')
!
            endif
!
!           Compute & print non-diagonal polarizabilities
!
            do l = 1, k - 1
!
               if (engine%compute_polarizability(k,l)) then
!
                  polarizability = zero
!
                  do sign_ = 1, 2 ! +, -
!
                     call engine%t_responses(freq, l, sign_)%open_('read', 'rewind')
                     call engine%t_responses(freq, l, sign_)%read_(tl(:,sign_), wf%n_es_amplitudes)
                     call engine%t_responses(freq, l, sign_)%close_()
!
                     polarizability = polarizability +                            &
                         half*ddot(wf%n_es_amplitudes, engine%etaX(:,k), 1, tl(:,sign_), 1) + &
                         half*ddot(wf%n_es_amplitudes, engine%etaX(:,l), 1, tk(:,sign_), 1)
!
                  enddo
!
                  if (engine%lr) then
!  
                     polarizability = polarizability +                  &
                           half*ddot(wf%n_es_amplitudes, Ftk(:,1), 1, tl(:,2), 1) + &
                           half*ddot(wf%n_es_amplitudes, Ftk(:,2), 1, tl(:,1), 1)
!
                  endif
!
                  call output%printf('m', '<< ' // operator_(k) // ', ' //  &
                                     operator_(l) // ' >>' // '((e8.2)): (f19.12)', &
                                     reals=[engine%frequencies(freq), &
                                     polarizability], fs='(t6,a)')
!
               endif
!
            enddo
!
         enddo
!
      enddo
!
!     Deallocate arrays to hold amplitude response vectors as well as
!     F-transformed response vectors (if LR)
!
      call mem%dealloc(tk, wf%n_es_amplitudes, 2)
      call mem%dealloc(tl, wf%n_es_amplitudes, 2)
!
      if (engine%lr) then
!
         call mem%dealloc(Ftk, wf%n_es_amplitudes, 2)
!
      endif
!
   end subroutine calculate_polarizabilities_response_engine
!
!
   subroutine determine_amplitude_response_response_engine(engine, wf)
!!
!!    Determine amplitude response
!!    Written by Eirik F. Kjønstad and Josefine H. Andersen, 2019
!!
!!    Makes preparations for determining the amplitude response t^X(omega_k) and t^X(-omega_k),
!!    where X is the dipole components (x,y,z) and omega_k the set of frequencies for which
!!    the polarizability is to be evaluated.
!!
!!    These responses are solutions of the equations
!!
!!       (A - omega_k I) t^X(omega_k) = -xi^X
!!
!!    The amplitude response t^X(+- omega_k) are stored to file. In particular,
!!    these are stored in a sequential file array:
!!
!!       engine%t_responses(freq, component, sign_), where:
!!
!!          freq:       1,2,3,...; denotes the frequency number as specified in input
!!                      e.g., if freq = 2, the file contains an amplitude response vector
!!                      for the second frequency specified on input
!!
!!          component:  1,2,3 (x,y,z); the component of the dipole moment vector
!!
!!          sign_:      1,2 (+, -); whether it is t^X(+omega_k) or t^X(-omega_k)
!!
!!    This wrapper routine, and the file storage for amplitude response vectors, was made by
!!    Eirik F. Kjønstad, Nov 2019. It is adapted/based on the general structure set up to solve
!!    amplitude response originally written by Josefine H. Andersen, spring 2019.
!!
      implicit none
!
      class(response_engine) :: engine
!
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable     :: rhs ! right-hand-side
!
      integer :: k, freq, sign_
!
      real(dp) :: prefactor
!
      character(len=200) :: file_name
      character(len=1)   :: sign_character
!
      class(davidson_cc_linear_equations), allocatable :: t_response_solver
!
!     Initialize amplitude response files for storing solutions
!
      allocate(engine%t_responses(engine%n_frequencies, 3, 2))
!
      do k = 1, 3
         do freq = 1, engine%n_frequencies
            do sign_ = 1, 2
!
               write(file_name, '(a, i3.3, a, i3.3, a, i3.3)') 'dipole_t_response_component_', &
                                                            k, '_frequency_', freq, '_sign_', sign_
!
               engine%t_responses(freq,k,sign_) = sequential_file(file_name)
!
            enddo
         enddo
      enddo
!
!     Solve response equations for one component and all frequencies
!
      call mem%alloc(rhs, wf%n_es_amplitudes)
!
      do k = 1, 3
!
!        kth component not needed for requested polarizability calculations
         if (.not. engine%t_response_components(k)) cycle
!
         call copy_and_scale(-one, engine%xiX(:,k), rhs, wf%n_es_amplitudes)
!
         do sign_ = 1, 2
!
            sign_character = '-'
            if (sign_ == 2) sign_character = '+'
!
            call output%printf('v', 'Determining amplitude response...', fs='(/t3,a)')
!
            call output%printf('v', 'Asking solver to get response for ' //  &
                               'component (i0) and ((a1))-frequencies.', &
                               ints=[k], chars=[sign_character], fs='(t3,a)')
!
            t_response_solver = davidson_cc_linear_equations(wf,                                   &
                                                             section='cc response',                &
                                                             eq_description='Solving for the       &
                                                             &amplitude response vectors in CC     &
                                                             &response theory.',                   &
                                                             n_frequencies=engine%n_frequencies,   &
                                                             n_rhs=1)
!
            prefactor = real((-1)**sign_, kind=dp) ! for frequencies
!
            call t_response_solver%run(wf, rhs, prefactor*engine%frequencies, &
                              engine%t_responses(:,k,sign_), 'right')
!
            call t_response_solver%cleanup(wf)
!
         enddo
      enddo
!
      call mem%dealloc(rhs, wf%n_es_amplitudes)
!
   end subroutine determine_amplitude_response_response_engine
!
!
   subroutine determine_M_vectors_response_engine(engine, wf)
!!
!!    Determine M vectors
!!    Written by Eirik F. Kjønstad and Josefine H. Andersen, 2019
!!
!!    Determines the solutions to the equations
!!
!!       (A^T + omega_k I) M_k = -F R_k,
!!
!!    where F is the so-called F transformation,
!!    (F_mu,nu = < Lambda | [[H-bar,tau_mu],tau_nu] | HF > ), and
!!    the omega_k are the excitation energies.
!!
!!    The M vectors are required to compute the transition moments
!!    using CC response theory.
!!
!!    On exit, the converged M vectors are stored in the file array
!!    engine%M_vectors.
!!
!!    This solver wrapper routine, and the file storage for M vectors, was made by
!!    Eirik F. Kjønstad, Nov 2019. It is adapted/based on the general structure set up
!!    to determine M vectors originally written by Josefine H. Andersen, spring 2019.
!!
      implicit none
!
      class(response_engine) :: engine
!
      class(ccs) :: wf
!
      integer :: k
!
      real(dp), dimension(:,:), allocatable :: minus_FR ! [-F R_k], k = 1, 2, ..., n_singlet_states
!
      character(len=200) :: file_name
!
      class(davidson_cc_linear_equations), allocatable :: M_vectors_solver
!
!     Build the matrix FR of right-hand-sides, -FR = [-F R_k], k = 1, 2, ..., n_singlet_states
!
      call mem%alloc(minus_FR, wf%n_es_amplitudes, wf%n_singlet_states)
!
      call wf%read_excited_state(minus_FR,            &
                                 1,                   &
                                 wf%n_singlet_states, &
                                 'right')
!
      do k = 1, wf%n_singlet_states
!
         call wf%F_transformation(minus_FR(:,k))
!
      enddo
!
      call dscal(wf%n_es_amplitudes*wf%n_singlet_states, -one, minus_FR, 1)
!
!     Initialize file array to store the converged M_vectors
!
      allocate(engine%M_vectors(wf%n_singlet_states))
!
      do k = 1, wf%n_singlet_states
!
         write(file_name, '(a, i3.3)') 'M_vector_state_', k
!
         engine%M_vectors(k) = sequential_file(trim(file_name))
!
      enddo
!
!     Call solver to converge the M vectors
!
      M_vectors_solver = davidson_cc_linear_equations(wf,                                       &
                                                      section='cc response',                    &
                                                      eq_description='Solving for the M vectors &
                                                      &in CC response theory.',                 &
                                                      n_frequencies=wf%n_singlet_states,        &
                                                      n_rhs=wf%n_singlet_states)
!
      call M_vectors_solver%run(wf, minus_FR, -wf%right_excitation_energies, &
                              engine%M_vectors, 'left')
!
      call M_vectors_solver%cleanup(wf)
!
      call mem%dealloc(minus_FR, wf%n_es_amplitudes, wf%n_singlet_states)
!
   end subroutine determine_M_vectors_response_engine
!
!
   subroutine read_settings_response_engine(engine)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(response_engine) :: engine
!
      call engine%read_gs_settings()
      call engine%read_es_settings()
      call engine%read_response_settings()
!
   end subroutine read_settings_response_engine
!
!
   subroutine read_response_settings_response_engine(engine)
!!
!!    Read response settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(response_engine) :: engine
!
      integer :: n_states
!
      character(len=200) :: restart_string
!
      engine%eom = input%is_keyword_present('eom','cc response')
!
      engine%lr = input%is_keyword_present('lr','cc response')
!
      engine%transition_moments = input%is_keyword_present('transition moments','cc response')
!
      engine%permanent_moments = input%is_keyword_present('permanent moments','cc response')
!
      if (input%is_keyword_present('polarizabilities','cc response')) then
!
         engine%polarizabilities = .true.
!
         call engine%set_polarizability_components()
!
         call engine%set_t_response_components()
!
!        Read in for which frequencies to compute the polarizability
!
         if (input%is_keyword_present('frequencies', 'cc response')) then
!
            engine%n_frequencies = input%get_n_elements_for_keyword('frequencies', 'cc response')
!
            call mem%alloc(engine%frequencies, engine%n_frequencies)
!
            call input%get_array_for_keyword('frequencies', 'cc response', &
                                                   engine%n_frequencies, engine%frequencies)
!
         else
!
            call output%error_msg('Requested polarizabilities, but no &
                                  &frequencies were provided!')
!
         endif
      endif
!
      engine%dipole_length = input%is_keyword_present('dipole length','cc response')
!
!     Sanity checks
!
      if (engine%eom .and. engine%lr) &
         call output%error_msg('can not run lr and eom in same calculation.')
!
      if (.not. engine%eom .and. .not. engine%lr) &
         call output%error_msg('specify either eom or lr for response.')
!
!     Since dipole is the only operator we currently have, we will check that it is indeed true
!
      if (.not. engine%dipole_length) &
            call output%error_msg('no operator selected in response calculation')
!
!     Check if the user wants to restart the right or left states in particular,
!     and not do restart on both
!
      if (input%is_keyword_present('restart', 'solver cc es')) then
!
         call input%get_keyword('restart',      &
                                'solver cc es', &
                                restart_string)
!
         if (trim(restart_string) == 'right') then
!
            engine%es_restart = .false. ! Turn off right & left restart
!
            engine%es_restart_right = .true. ! Turn on right restart
!
         elseif (trim(restart_string) == 'left') then
!
            engine%es_restart = .false. ! Turn off right & left restart
!
            engine%es_restart_left = .true. ! Turn on left restart
!
         endif
!
      endif
!
      if (engine%eom .and. (engine%transition_moments .or. engine%permanent_moments)) then
!
         call input%get_required_keyword('singlet states', 'solver cc es', n_states)
!
         call engine%set_initial_states(n_states)
!
!        Plotting of transition densities
!        plot_density      - logical to enable plotting of the gs density
!        plot_mn_densities - logical to enable plotting of tranistion densities
!        plot_es_densities - logical to enable plotting of excited state densities
!        states_to_plot    - list of integers containing the state numbers of the states
!                            for which densities shall be plotted (default all initial states)
!
         engine%plot_density = input%is_keyword_present('plot cc density', 'visualization')
!
         engine%plot_mn_densities = input%is_keyword_present('plot transition densities', &
                                                             'visualization')
         engine%plot_es_densities = input%is_keyword_present('plot es densities', &
                                                             'visualization')
!
         if (engine%plot_mn_densities .or. engine%plot_es_densities) then
            call engine%set_states_to_plot(n_states, engine%n_initial_states, &
                                           engine%initial_states)
         end if
!
      end if
!
   end subroutine read_response_settings_response_engine
!
!
   subroutine set_polarizability_components_response_engine(engine)
!!
!!    Set polarizability components
!!    Written by Alexander C. Paul, May 2021
!!
      implicit none
!
      class(response_engine), intent(inout) :: engine
!
      integer, dimension(:), allocatable :: polarizabilities
      integer :: k, n_polarizabilities
!
      n_polarizabilities = input%get_n_elements_for_keyword('polarizabilities', 'cc response')
!
      if (n_polarizabilities == 0) then
!
!        Not specified means you get all components (xx, xy, xz, ...)
!        of the polarizability
!
         engine%compute_polarizability(:,:) = .true.
!
      else
!
!        Figure out which polarizabilities-components have been requested
!        to set true/false in the compute_polarizability array
!
         call mem%alloc(polarizabilities, n_polarizabilities)
!
         call input%get_array_for_keyword('polarizabilities', 'cc response', &
                                          n_polarizabilities, polarizabilities)
!
         do k = 1, n_polarizabilities
!
            if (polarizabilities(k) == 11) then
!
               engine%compute_polarizability(1,1) = .true.
!
            elseif (polarizabilities(k) == 12 .or. polarizabilities(k) == 21) then
!
               engine%compute_polarizability(1,2) = .true.
               engine%compute_polarizability(2,1) = .true.
!
            elseif (polarizabilities(k) == 13 .or. polarizabilities(k) == 31) then
!
               engine%compute_polarizability(1,3) = .true.
               engine%compute_polarizability(3,1) = .true.
!
            elseif (polarizabilities(k) == 22) then
!
               engine%compute_polarizability(2,2) = .true.
!
            elseif (polarizabilities(k) == 23 .or. polarizabilities(k) == 32) then
!
               engine%compute_polarizability(2,3) = .true.
               engine%compute_polarizability(3,2) = .true.
!
            elseif (polarizabilities(k) == 33) then
!
               engine%compute_polarizability(3,3) = .true.
!
            endif
!
         enddo
!
         call mem%dealloc(polarizabilities, n_polarizabilities)
!
      endif
!
   end subroutine set_polarizability_components_response_engine
!
!
   subroutine set_t_response_components_response_engine(engine)
!!
!!    Set t response components
!!    Written by Alexander C. Paul, May 2021
!!
!!    Determine which components of the response amplitudes will be
!!    necessary to compute the requested polarizabilities,
!!
      implicit none
!
      class(response_engine), intent(inout) :: engine
!
      integer :: k, l
!
         do k = 1, 3
            do l = 1, k
!
               if (engine%compute_polarizability(k,l)) then
!
                  engine%t_response_components(k) = .true.
                  engine%t_response_components(l) = .true.
!
               endif
!
            enddo
         enddo
!
   end subroutine set_t_response_components_response_engine
!
!
   subroutine set_initial_states_response_engine(engine, n_states)
!!
!!    Set initial states
!!    Written Alexander C. Paul and Sarai D. Folkestad, Apr 2021
!!
!!    Determine the states for which properties shall be calculated
!!    if we are doing transition moments or permanent moments:
!!
!!       Example of usage:
!!       initial states: {0, 1, 2}
!!       or: [0,2]
!!
!!       gives all transitions from the ground state
!!          0 -> 1, 2, 3, ..
!!       and excitations between excited state 1 and 2 (ordered according to energy)
!!       to all other excited states
!!          1 -> 2, 3, ...
!!          2 -> 3, 4, ...
!!       or properties of the excited states 1 and 2
!!
      implicit none
!
      class(response_engine), intent(inout) :: engine
      integer, intent(in) :: n_states
!
      integer, dimension(:), allocatable :: initial_states
!
      if (input%is_keyword_present('initial states','cc response')) then
!
         engine%n_initial_states = input%get_n_elements_for_keyword('initial states', &
                                                                    'cc response')
!
         call mem%alloc(initial_states, engine%n_initial_states)
!
         call input%get_array_for_keyword('initial states', 'cc response', &
                                          engine%n_initial_states, initial_states)
!
!        Transitions from the ground state shall always be considered
         if (any(initial_states == 0)) then
!
            call mem%alloc(engine%initial_states, engine%n_initial_states)
!
            engine%initial_states = initial_states
!
            call mem%dealloc(initial_states, engine%n_initial_states)
!
         else
!
!           Copy and add 0 as initial state
!
            call mem%alloc(engine%initial_states, engine%n_initial_states + 1)
!
            engine%initial_states(1)  = 0
            engine%initial_states(2:) = initial_states
!
            call mem%dealloc(initial_states, engine%n_initial_states)
!
            engine%n_initial_states = engine%n_initial_states + 1
!
         endif
!
         if (engine%n_initial_states > n_states + 1) then
            call output%error_msg('Requested properties for more states &
                                  &than requested in the calculation.')
         end if
!
         if (any(engine%initial_states > n_states)) then
            call output%error_msg('Requested properties for state with a number &
                                  &larger than the number of states')
         end if
!
         if (any(engine%initial_states < 0)) then
            call output%error_msg('Requested properties of state with a negative number.')
         end if
!
      else
!
!        Only calculate transitions from ground state (state 0)
!
         engine%n_initial_states = 1
!
         call mem%alloc(engine%initial_states, engine%n_initial_states)
!
         engine%initial_states(1) = 0
!
      endif
!
   end subroutine set_initial_states_response_engine
!
!
   subroutine set_states_to_plot_response_engine(engine, n_states, n_initial, initial_states)
!!
!!    Set states to plot
!!    Written Alexander C. Paul, Apr 2021
!!
!!    Determine the states for which densities shall be plotted.
!!    If not present the densities for all initial states are requested
!!
      use array_utilities, only: copy_integer
!
      implicit none
!
      class(response_engine), intent(inout) :: engine
!
      integer, intent(in) :: n_initial, n_states
      integer, dimension(n_initial), intent(in) :: initial_states
!
      integer :: k
!
      if (input%is_keyword_present('states to plot', 'visualization')) then
!
         engine%n_states_to_plot = input%get_n_elements_for_keyword('states to plot', &
                                                                    'visualization')
!
         if (engine%n_states_to_plot > n_initial) then
            call output%error_msg('Requested more states to plot than initial states.')
         end if
!
         call mem%alloc(engine%states_to_plot, engine%n_states_to_plot)
!
         call input%get_array_for_keyword('states to plot', 'visualization', &
                                           engine%n_states_to_plot, engine%states_to_plot)
!
         do k = 1, engine%n_states_to_plot
            if (.not. any(initial_states == engine%states_to_plot(k))) then
               call output%error_msg('Requested plotting for state that is not an initial state.')
            end if
         end do
!
         if (any(engine%states_to_plot > n_states)) then
            call output%error_msg('Requested plotting of state with a number &
                                  &larger than the number of states')
         end if
!
         if (any(engine%states_to_plot < 0)) then
            call output%error_msg('Requested plotting of state with a negative number')
         end if
!
      else
!
         engine%n_states_to_plot = n_initial
!
         call mem%alloc(engine%states_to_plot, engine%n_states_to_plot)
!
         call copy_integer(initial_states, engine%states_to_plot, n_initial)
!
      end if
!
   end subroutine set_states_to_plot_response_engine
!
!
   subroutine do_lr_transition_moments_response_engine(engine, wf)
!!
!!    Do LR transition moments
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
!!    Adapted from routines written by Josefine H. Andersen, spring 2019.
!!
!!    Routine that performs three tasks:
!!
!!       - Calls routine to determine the M vectors, which are needed to compute the
!!         transition moments in response theory
!!
!!       - Calls wavefunction routine that computes the moments and strengths from
!!         the etaX, xiX and M vectors.
!!
!!       - Calls routine to print the transition moments / strengths to output
!!
!!    This routine was based on a deprecated routine for EOM moments
!!    (Josefine H. Andersen, Sarai D. Folkestad, June 2019). Reintroduced
!!    with M vector contributions and making use of skip_states array
!!    (Eirik F. Kjønstad, Nov 2019).
!!    Removed skip_states array as the parallel states are removed by the engine
!!    (Alexander C. Paul Sep 2020)
!!
      implicit none
!
      class(response_engine) :: engine
      class(ccs)        :: wf
!
      real(dp), dimension(3) :: transition_strength, transition_moment_left, transition_moment_right
!
      integer :: state, component
!
      real(dp), dimension(:), allocatable :: M ! M vector associated with state
!
      type(timings) :: timer
!
      call engine%tasks%print_('properties')
!
      timer = timings('Time to converge M vectors and calculate LR transition moments.', pl='n')
!
      call timer%turn_on()
!
      call engine%determine_M_vectors(wf)
!
      call mem%alloc(M, wf%n_es_amplitudes)
!
      do state = 1, wf%n_singlet_states
!
         call engine%M_vectors(state)%open_('read', 'rewind')
         call engine%M_vectors(state)%read_(M, wf%n_es_amplitudes)
         call engine%M_vectors(state)%close_()
!
         do component = 1, 3
!
            call wf%calculate_lr_transition_strength(transition_strength(component),      &
                                                     engine%etaX(:,component),            &
                                                     engine%xiX(:,component),            &
                                                     state,                               &
                                                     transition_moment_left(component),   &
                                                     transition_moment_right(component),  &
                                                     M)
!
         enddo
!
         call engine%print_lr_transition_moment_summary(transition_strength,                 &
                                                     transition_moment_left,              &
                                                     transition_moment_right, state,      &
                                                     wf%right_excitation_energies(state))
!
      enddo
!
      call mem%dealloc(M, wf%n_es_amplitudes)
!
      call timer%turn_off()
!
   end subroutine do_lr_transition_moments_response_engine
!
!
   subroutine do_eom_transition_moments_response_engine(engine, wf)
!!
!!    Do EOM transition moments
!!    Written by Josefine H. Andersen, Sarai D. Folkestad
!!    and Alexander C. Paul, June 2019
!!
!!    Computes the EOM dipole transition moments using transition densities
!!    and dipole moment integrals.
!!
!!    Restructured by Alexander C. Paul and Sarai D. Folkestad, Apr 2020
!!
!!    Restructured for general (GS and ES) transition momnets and state
!!    properties and to allow for the calculation of transition moments
!!    of states generated in another calculation.
!!
      implicit none
!
      class(response_engine) :: engine
      class(ccs) :: wf
!
      real(dp), dimension(:,:,:), allocatable :: operator_
      real(dp), dimension(:,:,:), allocatable :: transition_moments
!
      type(timings) :: EOM_timer
!
      logical  :: visualize
!
      visualize = (engine%plot_density .or. engine%plot_mn_densities &
                .or. engine%plot_es_densities)
!
      call engine%tasks%print_('properties')
!
      EOM_timer = timings('Time to calculate EOM properties')
      call EOM_timer%turn_on()
!
      call wf%prepare_for_properties
      call wf%initialize_gs_density()
      call wf%construct_gs_density()
!
      call wf%initialize_density_intermediates
!
      call output%printf('m', ':: EOM properties calculation', fs='(/t3,a)')
!
      if (engine%dipole_length) then
!
         call mem%alloc(operator_, wf%n_mo, wf%n_mo, 3)
!
!        Constructs dipole operator in t1-transformed basis.
!
         call wf%get_t1_oei('dipole', operator_)
!
      else
!
         call output%error_msg('EOM transition moments requested, but operator not specified')
!
      endif
!
!     Construct transition density and calculate transition moments
!
      call mem%alloc(transition_moments, wf%n_singlet_states+1, wf%n_singlet_states+1, 3)
!
      call wf%compute_eom_transition_moments(operator_,                 &
                                             transition_moments,        &
                                             engine%n_initial_states,   &
                                             engine%initial_states,     &
                                             engine%transition_moments, &
                                             engine%permanent_moments,  &
                                             visualize) ! Print densities to file
                                                        ! only if visualization is requested
!
!     Print summary
!
      call engine%print_eom_transition_moment_summary(wf, transition_moments)
!
      if (engine%permanent_moments) then
!
         call engine%print_permanent_moments_summary(wf, transition_moments, 3)
!
      end if
!
      if (visualize) call engine%visualize_cc_densities(wf)
!
!     Cleanup
!
      call mem%dealloc(transition_moments, wf%n_singlet_states+1, wf%n_singlet_states+1, 3)
!
      call mem%dealloc(operator_, wf%n_mo, wf%n_mo, 3)
!
      call wf%destruct_density_intermediates
      call wf%destruct_gs_density()
!
      call mem%dealloc(engine%initial_states, engine%n_initial_states)
!
      if (visualize) call mem%dealloc(engine%states_to_plot, engine%n_states_to_plot)
!
      call EOM_timer%turn_off()
!
   end subroutine do_eom_transition_moments_response_engine
!
!
   subroutine print_lr_transition_moment_summary_response_engine(engine, transition_strength, &
                  transition_moment_left, transition_moment_right, state, excitation_energy)
!!
!!    Print transition moment summary
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(response_engine) :: engine
!
      real(dp), dimension(3), intent(in) :: transition_strength
      real(dp), dimension(3), intent(in) :: transition_moment_left
      real(dp), dimension(3), intent(in) :: transition_moment_right
      real(dp), intent(in) :: excitation_energy
      integer, intent(in)  :: state
!
      character(len=3) :: calculation_type
      character(len=1), dimension(3) :: components = ['X', 'Y', 'Z']
!
      integer :: component
!
      real(dp) :: sum_strength
!
      if (engine%lr) calculation_type  = 'LR'
      if (engine%eom) calculation_type = 'EOM'
!
      call output%printf('m', 'State (i0):', fs='(/t6,a)', ints=[state])
      call output%print_separator('m', 10, '-', fs='(t6,a)')
!
      call output%printf('m', 'Calculation type:             (a19)', &
                         chars=[calculation_type], fs='(t6,a)')
      call output%printf('m', 'Excitation energy [E_h]:      (f19.12)', &
                         reals=[excitation_energy], fs='(t6,a)')
      call output%printf('m', 'Excitation energy [eV]:       (f19.12)', &
                         reals=[excitation_energy*Hartree_to_eV], fs='(t6,a)')
!
      call output%printf('m', 'Hartree-to-eV (CODATA 2014):  (f19.8)', &
                         reals=[Hartree_to_eV], fs='(t6,a)')
!
      call output%printf('m', 'Transition moments [a.u.]         Transition &
                         &strength [a.u.]', ll=74, fs='(/t6,14X,a)')
      call output%print_separator('m', 74, '-', fs='(t6,a)')
            call output%printf('m', 'Comp. q     < n |q| m >       < m |q| n >    &
                              &    < n |q| m > < m |q| n >', ll=79, fs='(t6,a)')
      call output%print_separator('m', 74, '-', fs='(t6,a)')
!
      sum_strength = zero
!
      do component = 1, 3
!
         call output%printf('m', '(a0)      (f17.10) (f17.10)       (f17.10)', &
                            reals=[transition_moment_left(component), &
                            transition_moment_right(component), &
                            transition_strength(component)], &
                            chars=[components(component)], fs='(t6,a)')
!
         sum_strength = sum_strength + transition_strength(component)
!
      enddo
!
      call output%print_separator('m', 74, '-', fs='(t6,a)')
!
      call output%printf('m', 'Oscillator strength: (f19.12)', &
                         reals=[(two/three)*excitation_energy*sum_strength], fs='(t6,a)')
!
   end subroutine print_lr_transition_moment_summary_response_engine
!
!
   subroutine print_eom_transition_moment_summary_response_engine(engine, wf, &
                                                                  transition_moments)
!!
!!    Print transition moment summary
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(response_engine) :: engine
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_singlet_states+1, wf%n_singlet_states+1, 3), &
                                                intent(in) :: transition_moments
!
      integer :: state_i, state_f, component, i
!
      real(dp) :: energy_i, energy_f
!
      character(len=3) :: calculation_type
      character(len=1), dimension(3) :: components = ['X', 'Y', 'Z']
!
      real(dp) :: sum_strength
!
      call output%printf('m', '- Summary of EOM transition properties calculation:', fs='(/t3,a)')
!
      calculation_type = 'EOM'
!
      do i = 1, engine%n_initial_states
!
         state_i = engine%initial_states(i)
!
         energy_i = 0.0d0
         if (state_i > 0) then
            energy_i = wf%right_excitation_energies(state_i)
!
            if (.not. engine%transition_moments .and. engine%permanent_moments) exit
!
         end if
!
         do state_f = 1, wf%n_singlet_states
!
            if (state_i == state_f) cycle
!
            if (any(engine%initial_states == state_f) .and. state_i > state_f) cycle
!
            energy_f = wf%right_excitation_energies(state_f)
!
            call output%printf('m', 'States m = (i0) and n = (i0):', fs='(/t6,a)', ints=[state_i, state_f])
            call output%print_separator('m', 25, '-', fs='(t6,a)')
!
            call output%printf('m', 'Calculation type:             (a19)', &
                         chars=[calculation_type], fs='(t6,a)')
!
            call output%printf('m', 'Excitation energy [E_h]:      (f19.12)', &
                         reals=[energy_f - energy_i], fs='(t6,a)')
!
            call output%printf('m', 'Excitation energy [eV]:       (f19.12)', &
                         reals=[(energy_f - energy_i)&
                                 *Hartree_to_eV], fs='(t6,a)')
!
            call output%printf('m', 'Hartree-to-eV (CODATA 2014):  (f19.8)', &
                         reals=[Hartree_to_eV], fs='(t6,a)')
!
            call output%printf('m', 'Transition moments [a.u.]         Transition &
                         &strength [a.u.]', ll=74, fs='(/t6,14X,a)')
!
            call output%print_separator('m', 74, '-', fs='(t6,a)')
!
            call output%printf('m', 'Comp. q     < n |q| m >       < m |q| n >    &
                              &    < n |q| m > < m |q| n >', ll=79, fs='(t6,a)')
!
            call output%print_separator('m', 74, '-', fs='(t6,a)')
!
            sum_strength = zero
!
            do component = 1, 3
!
!              Index = 1 indicates the ground state in the transition moments array
!
               call output%printf('m', '(a0)      (f17.10) (f17.10)       (f17.10)',      &
                                  reals=[transition_moments(state_f+1, state_i+1, component), &
                                         transition_moments(state_i+1, state_f+1, component), &
                                         transition_moments(state_i+1, state_f+1, component)* &
                                         transition_moments(state_f+1, state_i+1, component)],&
                                  chars=[components(component)], fs='(t6,a)')
!
               sum_strength = sum_strength + transition_moments(state_i+1, state_f+1, component)* &
                                             transition_moments(state_f+1, state_i+1, component)
!
            enddo
!
            call output%print_separator('m', 74, '-', fs='(t6,a)')
!
            call output%printf('m', 'Oscillator strength: (f19.12)', fs='(t6,a)', &
                                reals=[(two/three)*(energy_f - energy_i)*sum_strength])
!
         enddo
      enddo
!
   end subroutine print_eom_transition_moment_summary_response_engine
!
!
   subroutine print_operator_per_state_response_engine(engine, electronic, nuclear, &
                                                       components, n_components)
!!
!!    Print operator per excited states
!!    Written by Alexander C. Paul, June 2020
!!
!!    Prints operator/property for every initial state
!!    and the corresponding nuclear contribution.
!!
!!    Information about the operator, conversion factors between units, ...
!!    shall be printed in the routine calling this routine.
!!
      use array_utilities, only: get_l2_norm
      use range_class
!
      implicit none
!
      class(response_engine), intent(in) :: engine
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(engine%n_initial_states, n_components), intent(in) :: electronic
      real(dp), dimension(n_components), intent(in) :: nuclear
!
      character(len=4), dimension(n_components), intent(in) :: components
!
      real(dp), dimension(:,:), allocatable :: operator_
      real(dp), dimension(:), allocatable :: norms
!
      integer :: n_printed, n_to_print, ll, n_columns
      integer :: p, c, table, n_tables
!
      character(len=15), dimension(:), allocatable :: labels
      character(len=100) :: line
!
      type(range_), allocatable :: subset
!
      n_columns = engine%n_initial_states + 1
!
      call mem%alloc(operator_, n_components, n_columns)
      call mem%alloc(norms, n_columns)
!
!     Prepare arrays to print: Nuclear, state 0, state 1, ...
!
      do c = 1, n_components
         operator_(c, 1) = nuclear(c)
      end do
!
      do p = 2, n_columns
         do c = 1, n_components
!
            operator_(c,p) = electronic(p-1,c) + nuclear(c)
!
         end do
      end do
!
      do p = 1, n_columns
         norms(p) = get_l2_norm(operator_(:,p), n_components)
      end do
!
      allocate(labels(n_columns))
!
      labels(1) = '        Nuclear'
      do c = 1, engine%n_initial_states
         write(labels(c+1), '(i15)') engine%initial_states(c)
      end do
!
!     Print permanent moments in tables for 4 states at a time
      n_tables = (n_columns - 1)/4 + 1
!
      n_printed = 0
!
      do table = 1, n_tables
!
         n_to_print = min(n_columns-n_printed, 4)
         ll = 7 + 15*n_to_print
!
         subset = range_(n_printed + 1, n_to_print)
!
!        Print header of the table (Comp.  Label1  Label2 ...)
!
         write(line,'(a,i0,a)') ' Comp.', subset%length, '(a15)'
         call output%printf('m', trim(line), fs='(//t6,a)', ll=80, &
                            chars=[labels(subset%first:subset%get_last())])
!
         call output%print_separator('m', ll,'-', fs='(t6,a)')
!
!        Print table/components of the operator for every state in the subset
!
         do c = 1, n_components
!
            write(line,'(a5,1x,i0,a)') components(c), subset%length, '(f15.6)'
            call output%printf('m', trim(line), fs='(t6,a)', &
                               reals=[operator_(c, subset%first:subset%get_last())])
!
         end do
!
         call output%print_separator('m', ll,'-', fs='(t6,a)')
!
!        Print Norm of the operator for every state in the subset
         write(line,'(a,i0,a)') 'Norm  ', subset%length, '(f15.6)'
         call output%printf('m', trim(line), fs='(t6,a)', &
                           reals=[norms(subset%first:subset%get_last())])
!
         call output%print_separator('m', ll,'-', fs='(t6,a)')
!
         n_printed = n_printed + n_to_print
!
      end do
!
      call mem%dealloc(operator_, n_components, n_columns)
      call mem%dealloc(norms, n_columns)
      deallocate(labels)
!
   end subroutine print_operator_per_state_response_engine
!
!
   subroutine print_permanent_moments_summary_response_engine(engine, wf, operator_, &
                                                              n_components)
!!
!!    Print permanent moments summary
!!    Written by Alexander C. Paul, May 2020
!!
      use parameters
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(response_engine) :: engine
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(wf%n_singlet_states+1, wf%n_singlet_states+1, n_components), &
                                                                           intent(in) :: operator_
!
      real(dp), dimension(:),   allocatable :: nuclear
      real(dp), dimension(:,:), allocatable :: electronic
!
      character(len=4), dimension(n_components) :: components
      integer :: state_i, i, c
!
      call output%printf('m', '- Summary of EOM permanent moments calculation:', fs='(/t3,a)')
!
      call mem%alloc(nuclear, n_components)
      call mem%alloc(electronic, engine%n_initial_states, n_components)
!
      if (engine%dipole_length) then
!
         call output%printf('m', 'Total permanent dipole moments in [a.u.]:', fs='(/t6,a)')
         call output%print_separator('m', 41, '=', fs='(t6,a)')
!
         call output%printf('m', 'Conversion factor from au to Debye: (f11.9)', &
                             reals=[au_to_debye], fs='(/t6,a)')
!
         nuclear = wf%get_nuclear_dipole()
!
         components = [ 'X', 'Y', 'Z']
!
      end if
!
!     Select only the states that shall be printed
!
      do i = 1, engine%n_initial_states
!
         state_i = engine%initial_states(i)
!
         do c = 1, 3
            electronic(i, c) = operator_(state_i+1, state_i+1, c)
         end do
!
      end do
!
      call engine%print_operator_per_state(electronic, nuclear,  &
                                           components, n_components)
!
      call mem%dealloc(nuclear, n_components)
      call mem%dealloc(electronic, engine%n_initial_states, n_components)
!
   end subroutine print_permanent_moments_summary_response_engine
!
!
   subroutine set_printables_response_engine(engine)
!!
!!    Set printables
!!    Written by sarai D. Folkestad, May 2019
!!
!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(response_engine) :: engine
!
      character(len=5) :: response_type
!
      engine%name_ = 'Coupled cluster response engine'
!
      engine%tag = 'response'
!
      if (engine%eom) then
!
         response_type = 'EOM'
!
      else
!
         response_type = 'LR'
!
      endif
!
!     Prepare the list of tasks
!
      engine%tasks = task_list()
!
      call engine%tasks%add(label='cholesky', &
                            description='Cholesky decomposition of the electron &
                                         &repulsion integrals')
!
      call engine%tasks%add(label='mo preparations',                             &
                            description='Preparation of MO basis and integrals')
!
      call engine%tasks%add(label='gs solver',                                &
                            description='Calculation of the ground state ('// &
                           trim((engine%gs_algorithm))//' algorithm)')
!
      call engine%tasks%add(label='multipliers solver',                       &
                            description='Calculation of the multipliers ('    &
                            //trim((engine%multipliers_algorithm))&
                            //' algorithm)')
!
      call engine%tasks%add(label='es solver',                                &
                           description='Calculation of the excited states ('// &
                           trim((engine%es_algorithm))//' algorithm)')
!
      if (engine%transition_moments .or. engine%permanent_moments) then
!
         call engine%tasks%add(label='properties', &
                           description='Calculation of excited state properties (' &
                           //trim(response_type)//')')
!
      endif
!
      if (engine%polarizabilities) then
!
         call engine%tasks%add(label='polarizabilities',          &
                               description='Calculation of the '  &
                               //trim(response_type)// ' polarizabilities')
!
      endif
!
      if (engine%plot_density .or. engine%plot_mn_densities &
         .or. engine%plot_es_densities) then
!
         call engine%tasks%add(label='plotting', description='Visualization of CC densities')
!
      end if
!
      engine%description = 'Calculates dipole transition moments and oscillator strengths between &
                           &the ground state and the excited states.'
!
   end subroutine set_printables_response_engine
!
!
   subroutine visualize_cc_densities_response_engine(engine, wf)
!!
!!    Visualize cc densities
!!    Written by Alexander C. Paul, Dec 2020
!!
      use visualization_class, only: visualization
!
      implicit none
!
      class(response_engine) :: engine
      class(ccs) :: wf
!
      type(visualization), allocatable :: visualizer
!
      real(dp), dimension(:,:), allocatable :: c_D_ct, density
!
      integer :: p, state_p, state_q
      character(len=10) :: tag
!
      call engine%tasks%print_('plotting')
!
      visualizer = visualization(wf%ao)
      call visualizer%initialize(wf%ao)
!
      call mem%alloc(c_D_ct, wf%ao%n, wf%ao%n)
!
!     GS density
      if (engine%plot_density .or. any(engine%states_to_plot == 0)) then
         call wf%add_t1_terms_and_transform(wf%density, c_D_ct)
         call visualizer%plot_density(wf%ao, c_D_ct, 'dm_000_000')
      end if
!
!     Transition densities
!
      call mem%alloc(density, wf%n_mo, wf%n_mo)
!
      if (engine%plot_mn_densities .and. engine%transition_moments) then
!
         do p = 1, engine%n_initial_states
!
            state_p = engine%initial_states(p)
!
            do state_q = 1, wf%n_singlet_states
!
               if (state_p == state_q) cycle
!
!              Have the density matrices been computed
               if ( .not. any(engine%states_to_plot == state_p) &
               .or. .not. any(engine%states_to_plot == state_q)) cycle
!
               c_D_ct = wf%get_density_for_plotting(density, state_p, state_q)
!
               write(tag, '(a, i3.3, a, i3.3)') 'dm_', state_p, '_', state_q
               call visualizer%plot_density(wf%ao, c_D_ct, tag)
!
               c_D_ct = wf%get_density_for_plotting(density, state_q, state_p)
!
               write(tag, '(a, i3.3, a, i3.3)') 'dm_', state_q, '_', state_p
               call visualizer%plot_density(wf%ao, c_D_ct, tag)
!
            end do
         end do
!
      end if
!
!     Excited state densities
!
      if (engine%plot_es_densities .and. engine%permanent_moments) then
!
         do p = 2, engine%n_initial_states ! first initial state is the GS
!
            state_p = engine%initial_states(p)
!
!           Has the density matrix been computed
            if (.not. any(engine%states_to_plot == state_p)) cycle
!
            c_D_ct = wf%get_density_for_plotting(density, state_p, state_p)
!
            write(tag, '(a, i3.3, a, i3.3)') 'dm_', state_p, '_', state_p
            call visualizer%plot_density(wf%ao, c_D_ct, tag)
!
         end do
!
      end if
!
      call mem%dealloc(density, wf%n_mo, wf%n_mo)
!
      call mem%dealloc(c_D_ct, wf%ao%n, wf%ao%n)
      call visualizer%cleanup()
!
   end subroutine visualize_cc_densities_response_engine
!
!
end module response_engine_class
