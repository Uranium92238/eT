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
!     Do transition moments and polarizabilities?
!
      logical :: transition_moments
      logical :: polarizabilities
!
      logical, dimension(3,3) :: compute_polarizability  ! which components to compute?
                                                         ! xx, xy, xz; yx, yy, yz; ... 
!
      integer :: n_frequencies
      real(dp), dimension(:), allocatable :: frequencies ! array of frequencies for which 
                                                         ! to evaluate the polarizability 
!
      logical, dimension(3) :: compute_t_response        ! compute amplitude response?
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
      logical :: plot_tdm
      integer :: n_states_to_plot
      integer, dimension(:), allocatable :: states_to_plot
!
   contains
!
      procedure :: run                              => run_response_engine
!
      procedure :: read_settings                    => read_settings_response_engine
      procedure :: read_response_settings           => read_response_settings_response_engine
!
      procedure, nopass :: get_thresholds           => get_thresholds_response_engine
!
      procedure :: print_transition_moment_summary  => print_transition_moment_summary_response_engine
!
      procedure :: set_printables                   => set_printables_response_engine
!
      procedure :: do_eom_transition_moments        => do_eom_transition_moments_response_engine
      procedure :: do_lr_transition_moments         => do_lr_transition_moments_response_engine
!
      procedure :: construct_etaX_and_xiX           => construct_etaX_and_xiX_response_engine
      procedure :: determine_M_vectors              => determine_M_vectors_response_engine
      procedure :: determine_amplitude_response     => determine_amplitude_response_response_engine 
      procedure :: calculate_polarizabilities       => calculate_polarizabilities_response_engine
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
      engine%gs_algorithm           = 'diis'
!
      if (wf%name_ .eq. 'ccsd(t)' .or. &
          wf%name_ .eq. 'mp2' .or.     &
          wf%name_ .eq. 'mlcc2' .or.   &
          wf%name_ .eq. 'mlccsd' .or.   &
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
         engine%es_algorithm          = 'diis'
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
      engine%compute_polarizability = .false.
      engine%compute_t_response     = .false.
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
      engine%plot_tdm               = .false.
!
      call engine%read_settings()
!
!
      engine%restart =  engine%gs_restart .or. &
                        engine%multipliers_restart .or. &
                        engine%es_restart
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
      class(ccs)         :: wf
!
      real(dp) :: energy_threshold, residual_threshold
!
      logical, dimension(:), allocatable :: skip_states
!
      if (trim(wf%name_) == 'low memory cc2') then
!
         call output%printf('m', 'Requested response calculation with lowmem-cc2.', &
                            fs='(/t3,a)')
         call output%error_msg('Response not implemented for low memory cc2.')
!
      end if
!
      call engine%tasks%print_('cholesky')
!
      call engine%do_cholesky(wf)
!
      call engine%tasks%print_('mo preparations')
!
      call wf%mo_preparations()
!
      call engine%restart_handling(wf)
!
!     Ground state solution
!
      call engine%do_ground_state(wf)
!
!     Determine multipliers
!
      call engine%do_multipliers(wf)
!
!     Construct etaX and xiX if they are going to be needed 
!
      if (engine%polarizabilities .or. &
            (engine%lr .and. engine%transition_moments)) then 
!
         call mem%alloc(engine%xiX, wf%n_es_amplitudes, 3)
         call mem%alloc(engine%etaX, wf%n_es_amplitudes, 3)
!
         call engine%construct_etaX_and_xiX(wf) ! etaX_mu = < Lambda | [X-bar, τ_mu] | HF >
                                                ! xiX_mu = < mu | X-bar | HF >         
!
      endif 
!
!     Calculate transition moments 
!
      if (engine%transition_moments) then 
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
                                      restart=(engine%es_restart .or. engine%es_restart_left))
!
!        Biorthornormalize and look for duplicate states 
!
         call mem%alloc(skip_states, wf%n_singlet_states)
!
         call engine%get_thresholds(energy_threshold, residual_threshold)
!
         call wf%biorthonormalize_L_and_R(energy_threshold, residual_threshold, skip_states)
!  
!        Compute transition moments 
!
         if (engine%eom) then 
!
            call engine%do_eom_transition_moments(wf, skip_states)
!
         elseif (engine%lr) then 
!
            call engine%do_lr_transition_moments(wf, skip_states)
!
         endif 
!
         call mem%dealloc(skip_states, wf%n_singlet_states)
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
      call wf%construct_mu(X)
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
      integer :: k, l, freq, sign 
!
      real(dp) :: ddot, polarizability
!
      character(len=4), dimension(3) :: operator = ['mu_x', 'mu_y', 'mu_z']
!
      call engine%tasks%print_('polarizabilities')
!
!     Allocate arrays to hold amplitude response vectors as well as 
!     F-transformed response vectors (if LR)
!
      call mem%alloc(tk, wf%n_gs_amplitudes, 2)
      call mem%alloc(tl, wf%n_gs_amplitudes, 2)
!
      if (engine%lr) then 
!
         call mem%alloc(Ftk, wf%n_gs_amplitudes, 2)
!
      endif 
!
!     Compute polarizabilities
!
      do freq = 1, engine%n_frequencies
!
         do k = 1, 3
!
!           Compute & print diagonal polarizabilities
!
            if (engine%compute_polarizability(k,k)) then 
!
               polarizability = zero 
!
               do sign = 1, 2 ! +, - 
!
                  call engine%t_responses(freq, k, sign)%open_('read', 'rewind')
                  call engine%t_responses(freq, k, sign)%read_(tk(:,sign), wf%n_gs_amplitudes)
                  call engine%t_responses(freq, k, sign)%close_()
!
                  polarizability = polarizability + &
                         ddot(wf%n_gs_amplitudes, engine%etaX(:,k), 1, tk(:,sign), 1)
!
               enddo
!
               if (engine%lr) then 
!
                  call dcopy(wf%n_gs_amplitudes, tk(:,1), 1, Ftk(:,1), 1)
                  call dcopy(wf%n_gs_amplitudes, tk(:,2), 1, Ftk(:,2), 1) 
!
                  call wf%F_transformation(Ftk(:,1))                           
                  call wf%F_transformation(Ftk(:,2))                           
!
                  polarizability = polarizability + &
                        ddot(wf%n_gs_amplitudes, Ftk(:,2), 1, tk(:,1), 1)
!
               endif 
!
               call output%printf('m', '<< ' // operator(k) // ', ' //  &
                                  operator(k) // ' >>' // '((e8.2)): (f19.12)', &
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
                  do sign = 1, 2 ! +, - 
!
                     call engine%t_responses(freq, l, sign)%open_('read', 'rewind')
                     call engine%t_responses(freq, l, sign)%read_(tl(:,sign), wf%n_gs_amplitudes)
                     call engine%t_responses(freq, l, sign)%close_()
!
                     polarizability = polarizability +                                       &
                         half*ddot(wf%n_gs_amplitudes, engine%etaX(:,k), 1, tl(:,sign), 1) + &
                         half*ddot(wf%n_gs_amplitudes, engine%etaX(:,l), 1, tk(:,sign), 1)
!
                  enddo    
!
                  if (engine%lr) then 
!
                     polarizability = polarizability +                              &
                           half*ddot(wf%n_gs_amplitudes, Ftk(:,1), 1, tl(:,2), 1) + &
                           half*ddot(wf%n_gs_amplitudes, Ftk(:,2), 1, tl(:,1), 1)
!
                  endif
!
                  call output%printf('m', '<< ' // operator(k) // ', ' //  &
                                     operator(l) // ' >>' // '((e8.2)): (f19.12)', &
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
      call mem%dealloc(tk, wf%n_gs_amplitudes, 2)
      call mem%dealloc(tl, wf%n_gs_amplitudes, 2)
!
      if (engine%lr) then 
!
         call mem%dealloc(Ftk, wf%n_gs_amplitudes, 2)
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
!!       engine%t_responses(freq, component, sign), where:
!!
!!          freq:       1,2,3,...; denotes the frequency number as specified in input
!!                      e.g., if freq = 2, the file contains an amplitude response vector 
!!                      for the second frequency specified on input 
!!
!!          component:  1,2,3 (x,y,z); the component of the dipole moment vector 
!!
!!          sign:       1,2 (+, -); whether it is t^X(+omega_k) or t^X(-omega_k) 
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
      integer :: k, freq, sign
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
            do sign = 1, 2 
!
               write(file_name, '(a, i3.3, a, i3.3, a, i3.3)') 'dipole_t_response_component_', &
                                                            k, '_frequency_', freq, '_sign_', sign
!
               engine%t_responses(freq,k,sign) = sequential_file(file_name)
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
         if (.not. engine%compute_t_response(k)) cycle ! skip; kth component not needed
                                                       ! for requested polarizability calculations! 
!
         call copy_and_scale(-one, engine%xiX(:,k), rhs, wf%n_es_amplitudes)
!
         do sign = 1, 2
!
            sign_character = '-'
            if (sign == 2) sign_character = '+'  
!
            call output%printf('m', 'Determining amplitude response...', fs='(/t3,a)')
!
            call output%printf('m', 'Asking solver to get response for ' //  &
                               'component (i0) and ((a1))-frequencies.', &
                               ints=[k], chars=[sign_character], fs='(t3,a)')
!
            t_response_solver = davidson_cc_linear_equations(wf,                             &
                                                             section='cc response',          &
                                                             eq_description='Solving for the &
                                                             &amplitude response vectors in CC &
                                                             &response theory.')
!
            prefactor = real((-1)**sign, kind=dp) ! for frequencies 
!
            call t_response_solver%run(wf, rhs, 1, prefactor*engine%frequencies, &
                              engine%n_frequencies, engine%t_responses(:,k,sign), 'right')
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
                                                      &in CC response theory.')
!
      call M_vectors_solver%run(wf, minus_FR, wf%n_singlet_states, -wf%right_excitation_energies, &
                           wf%n_singlet_states, engine%M_vectors, 'left')
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
      integer :: n_polarizabilities, k, l
      integer, dimension(:), allocatable :: polarizabilities
!
      character(len=200) :: restart_string
!
      if (input%requested_keyword_in_section('eom','cc response')) then 
!
         engine%eom = .true.
!
      endif 
!
      if (input%requested_keyword_in_section('lr','cc response')) then 
!
         engine%lr = .true.
!
      endif 
!
      if (input%requested_keyword_in_section('transition moments','cc response')) then 
!
         engine%transition_moments = .true.
!
      endif 
!
      if (input%requested_keyword_in_section('polarizabilities','cc response')) then 
!
         engine%polarizabilities = .true.
!
         n_polarizabilities = input%get_n_elements_for_keyword_in_section('polarizabilities', &
                                                                           'cc response')
!
         if (n_polarizabilities == 0) then 
!
!           Not specified means you get all components (xx, xy, xz, ...)
!           of thee polarizability
!
            engine%compute_polarizability(:,:) = .true.
!
         else 
!
!           Figure out which polarizabilities-components have been requested 
!           to set true/false in the compute_polarizability array
!
            call mem%alloc(polarizabilities, n_polarizabilities)
!
            call input%get_array_for_keyword_in_section('polarizabilities', 'cc response', &
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
!        Based on the requested polarizabilities, 
!        determine which amplitude response components 
!        will be necessary to determine 
!
         do k = 1, 3
            do l = 1, k
!
               if (engine%compute_polarizability(k,l)) then 
!
                  engine%compute_t_response(k) = .true.
                  engine%compute_t_response(l) = .true. 
!
               endif
!
            enddo
         enddo
!
!        Read in for which frequencies to compute the polarizability
!
         if (input%requested_keyword_in_section('frequencies', 'cc response')) then 
!
            engine%n_frequencies = input%get_n_elements_for_keyword_in_section('frequencies', 'cc response')
!  
            call mem%alloc(engine%frequencies, engine%n_frequencies)
!
            call input%get_array_for_keyword_in_section('frequencies', 'cc response', &
                                                   engine%n_frequencies, engine%frequencies) 
!
         else 
!
            call output%error_msg('Asked for polarizabilities, but no frequencies was ' // &
                                    'provided!')
!
         endif
      endif
!
!     Set operator
!
      engine%dipole_length = input%requested_keyword_in_section('dipole length','cc response')
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
      if (input%requested_keyword_in_section('restart', 'solver cc es')) then 
!
         call input%get_keyword_in_section('restart',       &
                                           'solver cc es',  &
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
!     Plotting of transition densities
!     plot_density   - logical to enable plotting of the gs density
!     plot_tdm       - logical to enable plotting of tranistion densities
!     states_to_plot - list of integers containing states of which 
!                      the transition densities shall be plotted (default all states)
!
      engine%plot_density = &
               input%requested_keyword_in_section('plot cc density', 'visualization')
!
      engine%plot_tdm = &
               input%requested_keyword_in_section('plot transition densities', 'visualization')
!
      if (engine%plot_tdm) then
!
         if (input%requested_keyword_in_section('states to plot', 'visualization')) then
!
            engine%n_states_to_plot = input%get_n_elements_for_keyword_in_section(&
                                      'states to plot', 'visualization')
!
            call mem%alloc(engine%states_to_plot, engine%n_states_to_plot)
!
            call input%get_array_for_keyword_in_section('states to plot',         &
                                                        'visualization',          &
                                                         engine%n_states_to_plot, &
                                                         engine%states_to_plot)
!
         else
!
            call input%get_required_keyword_in_section('singlet states',   &
                                                       'solver cc es',     &
                                                        engine%n_states_to_plot)
!
            call mem%alloc(engine%states_to_plot, engine%n_states_to_plot)
!
            do k = 1, engine%n_states_to_plot
               engine%states_to_plot(k) = k
            end do
!
         end if
!
      end if
!
   end subroutine read_response_settings_response_engine
!
!
   subroutine do_lr_transition_moments_response_engine(engine, wf, skip_states)
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
!!
      implicit none
!
      class(response_engine) :: engine
      class(ccs)        :: wf
!
      logical, dimension(wf%n_singlet_states), intent(in) :: skip_states
!
      real(dp), dimension(3) :: transition_strength, transition_moment_left, transition_moment_right
!
      integer :: state, component
!
      real(dp), dimension(:), allocatable :: M ! M vector associated with state 
!
      type(timings) :: timer
!     
      call engine%tasks%print_('transition moments')
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
         if (skip_states(state)) then
!
            call output%warning_msg('Skipped state (i0) because it is parallel to the &
                                    &previous state', ints=[state], fs='(/t3,a)')
!
            cycle
!
         endif
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
         call engine%print_transition_moment_summary(transition_strength,                 &
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
   subroutine do_eom_transition_moments_response_engine(engine, wf, skip_states)
!!
!!    Do EOM transition moments 
!!    Written by Josefine H. Andersen, Sarai D. Folkestad 
!!    and Alexander C. Paul, June 2019
!!
!!    Computes the EOM dipole transition moments using transition densities 
!!    and dipole moment integrals. 
!!
      use visualization_class, only: visualization
!
      implicit none
!
      class(response_engine) :: engine
      class(ccs)             :: wf
!
      logical, dimension(wf%n_singlet_states), intent(in) :: skip_states
!
      real(dp), dimension(:,:,:), allocatable :: operator
      real(dp), dimension(:,:),   allocatable :: c_D_ct
!
      real(dp), dimension(3) :: transition_strength, transition_moment_left, transition_moment_right
!
      integer :: state, component
!
      real(dp) :: trace_r_tdm, trace_l_tdm
      integer :: p
!
      type(sequential_file) :: density_file
      character(len=15)     :: file_name
!
      type(timings) :: EOM_timer
!
      type(visualization), allocatable :: visualizer
!
      call engine%tasks%print_('transition moments')
!
      EOM_timer = timings('Total time to calculate EOM transition moments.')
!
      call EOM_timer%turn_on()
!
      call wf%prepare_for_density()
      call wf%initialize_gs_density()
      call wf%construct_gs_density()
!
      call wf%initialize_transition_densities()
!
      call output%printf('m', ':: EOM properties calculation', fs='(/t3,a)')
!
      if (engine%dipole_length) then
!
         call mem%alloc(operator, wf%n_mo, wf%n_mo, 3)
!
!        Constructs dipole operator in t1-transformed basis.
!
         call wf%construct_mu(operator)
!
!        Loop over excited states, construct transition density
!        and calculate transition strength
!
         call output%printf('m', '- Summary of EOM properties calculation:', &
                            fs='(/t3,a)')
!
         do state = 1, wf%n_singlet_states
!
            if (skip_states(state)) then
!
               call output%warning_msg('Skipped state (i0) because it is parallel to the &
                                       &previous state', ints=[state], fs='(/t3,a)')
               cycle
!
            end if
!
!           Construct right tdm and write to file
!
            call wf%construct_right_transition_density(state)
!
            write(file_name, '(a, i3.3)') 'right_tdm_', state
            density_file = sequential_file(trim(file_name))
            call density_file%open_('write')
            call density_file%write_(wf%right_transition_density, wf%n_mo**2)
            call density_file%close_()
!
!           Construct left tdm and write to file
!
            call wf%construct_left_transition_density(state)
!
            write(file_name, '(a, i3.3)') 'left_tdm_', state
            density_file  = sequential_file(trim(file_name))
            call density_file%open_('write')
            call density_file%write_(wf%left_transition_density, wf%n_mo**2)
            call density_file%close_()
!
            do component = 1, 3
!
               transition_moment_left(component) = wf%calculate_expectation_value(operator(:,:,component),   &
                                                                              wf%left_transition_density)
!
               transition_moment_right(component) = wf%calculate_expectation_value(operator(:,:,component),  &
                                                                              wf%right_transition_density)
!
               transition_strength(component) = transition_moment_left(component)* &
                                                transition_moment_right(component)
!
            enddo
!
!           Print results
!
            call engine%print_transition_moment_summary(transition_strength, transition_moment_left, &
                                          transition_moment_right, state,              &
                                          wf%right_excitation_energies(state))
!
!           Print traces of the density matrices for Debugging
!
            trace_r_tdm = zero
            trace_l_tdm = zero
!
            do p = 1, wf%n_mo
!
               trace_r_tdm = trace_r_tdm + wf%right_transition_density(p,p)
               trace_l_tdm = trace_l_tdm + wf%left_transition_density(p,p)
!
            end do
!
            call output%printf('debug', 'Trace left transition density:  (f15.12)', &
                               reals=[trace_l_tdm], fs='(t6,a)')
            call output%printf('debug', 'Trace right transition density: (f15.12)', &
                               reals=[trace_r_tdm], fs='(t6,a/)')
!
         enddo
!
         call mem%dealloc(operator, wf%n_mo, wf%n_mo, 3)
!
      endif
!
      if (engine%plot_density .or. engine%plot_tdm) then
!
         call engine%tasks%print_('plotting')
!
         visualizer = visualization(wf%system, wf%n_ao)
!
         call visualizer%initialize(wf%system)
         call mem%alloc(c_D_ct, wf%n_ao, wf%n_ao)
!
         if (engine%plot_density) then
!
            call wf%add_t1_terms_and_transform(wf%density, c_D_ct)
            call visualizer%plot_density(wf%system, c_D_ct, 'cc_gs_density')
!
         end if
!
         if (engine%plot_tdm) then
!
            do p = 1, engine%n_states_to_plot
!
               state = engine%states_to_plot(p)
!
               write(file_name, '(a, i3.3)') 'right_tdm_', state
               density_file  = sequential_file(trim(file_name))
               call density_file%open_('read')
               call density_file%read_(wf%right_transition_density, wf%n_mo**2)
!
               call wf%add_t1_terms_and_transform(wf%right_transition_density, c_D_ct)
               call visualizer%plot_density(wf%system, c_D_ct, file_name)
!
               call density_file%close_()
!
               write(file_name, '(a, i3.3)') 'left_tdm_', state
               density_file  = sequential_file(trim(file_name))
               call density_file%open_('read')
               call density_file%read_(wf%left_transition_density, wf%n_mo**2)
!
               call wf%add_t1_terms_and_transform(wf%left_transition_density, c_D_ct)
               call visualizer%plot_density(wf%system, c_D_ct, file_name)
!
               call density_file%close_
!            
            end do
!
         end if
!         
         call mem%dealloc(c_D_ct, wf%n_ao, wf%n_ao)
         call visualizer%cleanup()
!
      end if
!
      call wf%destruct_transition_densities()
      call wf%destruct_gs_density()
!
      call EOM_timer%turn_off()
!
      if (engine%plot_tdm) call mem%dealloc(engine%states_to_plot, engine%n_states_to_plot)
!
   end subroutine do_eom_transition_moments_response_engine
!
!
   subroutine get_thresholds_response_engine(energy_threshold, residual_threshold)
!!
!!    Get thresholds from input
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    Get thresholds from input to perform checks for parallel states and 
!!    to check that the left and right states are consistent
!!
      implicit none
!
      real(dp), intent(out) :: energy_threshold, residual_threshold
!
!     Set default values and overwrite if specified in input
!
      energy_threshold = 1.0d-6
      residual_threshold = energy_threshold
!
!     Set thresholds for the sanity checks if the roots are ordered correctly
!
      if (input%requested_keyword_in_section('energy threshold', 'solver cc es') .and. &
          input%requested_keyword_in_section('residual threshold', 'solver cc es')) then 
!
        call input%get_keyword_in_section('energy threshold', 'solver cc es', energy_threshold)
        call input%get_keyword_in_section('residual threshold', 'solver cc es', residual_threshold)
!
      else if (input%requested_keyword_in_section('residual threshold', 'solver cc es')) then 
!
        call input%get_keyword_in_section('residual threshold', 'solver cc es', residual_threshold)
        energy_threshold = residual_threshold
!
      else if (input%requested_keyword_in_section('energy threshold', 'solver cc es')) then 
!
         call input%get_keyword_in_section('energy threshold', 'solver cc es', energy_threshold)
         residual_threshold = energy_threshold
!
      endif
!
   end subroutine get_thresholds_response_engine
!
!
   subroutine print_transition_moment_summary_response_engine(engine, transition_strength, &
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
      call output%printf('m', 'Comp. q     < k |q| 0 >       < 0 |q| k >        &
                         &< 0 |q| k > < k |q| 0 >  ', ll=79, fs='(t6,a)')
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
   end subroutine print_transition_moment_summary_response_engine
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
                           description='Calculation of the excited state ('// &
                           trim((engine%es_algorithm))//' algorithm)')
!
      if (engine%transition_moments) then
!
         call engine%tasks%add(label='transition moments',                       &
                           description='Calculation of the transition moments (' &
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
      if (engine%plot_density .or. engine%plot_tdm) then
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
end module response_engine_class
