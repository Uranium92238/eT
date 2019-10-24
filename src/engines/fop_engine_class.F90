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
module fop_engine_class
!!
!!    First order coupled cluster engine class module
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Josefine H. Andersen, Apr 2019
!!
   use parameters
   use global_in,             only: input
   use global_out,            only: output
   use timings_class,         only: timings
   use memory_manager_class,  only: mem
   use sequential_file_class, only: sequential_file
!
   use es_engine_class, only: es_engine
   use ccs_class,       only: ccs
!
   type, extends(es_engine) :: fop_engine
!
!     Property solver types
!
      logical :: eom
      logical :: lr
!
!     Operators
!
      logical :: dipole_length
!
   contains
!
      procedure :: run                       => run_fop_engine
!
      procedure :: read_settings             => read_settings_fop_engine
      procedure :: read_fop_settings         => read_fop_settings_fop_engine
!
      procedure :: do_eom                    => do_eom_fop_engine
!
      procedure, nopass :: get_thresholds    => get_thresholds_fop_engine
!
      procedure, nopass :: print_summary_eom => print_summary_eom_fop_engine
!
      procedure :: set_printables            => set_printables_fop_engine
!
   end type fop_engine
!
!
   interface fop_engine
!
      procedure :: new_fop_engine 
!
   end interface fop_engine
!
!
contains
!
!
   function new_fop_engine() result(engine)
!!
!!    New FOP engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(fop_engine) :: engine
!
!     Set standards and then read if nonstandard
!
      engine%es_algorithm           = 'davidson'
      engine%gs_algorithm           = 'diis'
      engine%multipliers_algorithm  = 'davidson'
      engine%es_type                = 'valence'
      engine%lr                     = .false.
      engine%eom                    = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_fop_engine
!
!
   subroutine run_fop_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(fop_engine) :: engine
      class(ccs)         :: wf
!
      real(dp) :: energy_threshold, residual_threshold
!
      logical, dimension(:), allocatable :: skip_states
!
!     Cholesky decomposition
!
      call engine%do_cholesky(wf)
!
      call wf%mo_preparations()
!
!     Ground state solution
!
      call engine%do_ground_state(wf)
!
!     Prepare for excited state calculation
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      if (wf%integrals%get_eri_t1_mem()) call wf%integrals%update_g_pqrs_t1_in_memory()
!
      if(wf%integrals%get_eri_t1_mem()) &
         call output%printf('Note: All T1-integrals are stored in memory',fs='(/t3, a)',pl='normal')
!
!     Determine multipliers
!
      call engine%do_multipliers(wf)
!
!     Excited state solutions
!
      call engine%do_excited_state(wf, 'right')
      call engine%do_excited_state(wf, 'left')
!
      call mem%alloc(skip_states, wf%n_singlet_states)
!
      call engine%get_thresholds(energy_threshold, residual_threshold)
!
      call wf%biorthonormalize_L_and_R(energy_threshold, residual_threshold, skip_states)
!
!     EOM properties if requested
!
      if (engine%eom) then
!
         call wf%prepare_for_density()
         call wf%initialize_gs_density()
         call wf%construct_gs_density()
!
         call wf%initialize_transition_densities()
         call engine%do_eom(wf, skip_states)
         call wf%destruct_transition_densities()
!
      end if
!
      call mem%dealloc(skip_states, wf%n_singlet_states)
!
   end subroutine run_fop_engine
!
!
   subroutine read_settings_fop_engine(engine)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(fop_engine) :: engine
!
      call engine%read_gs_settings()
      call engine%read_es_settings()
      call engine%read_fop_settings()
!
   end subroutine read_settings_fop_engine
!
!
   subroutine read_fop_settings_fop_engine(engine)
!!
!!    Read FOP settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(fop_engine) :: engine
!
      if (input%requested_keyword_in_section('dipole length','cc fop')) engine%dipole_length = .true.
      if (input%requested_keyword_in_section('eom','cc fop')) engine%eom   = .true.
      if (input%requested_keyword_in_section('lr','cc fop')) engine%lr     = .true.
!
!     Sanity checks
!
      if (engine%eom .and. engine%lr) call output%error_msg('can not run lr and eom in same calculation.')
      if (engine%lr) call output%error_msg('lr not yet available.')
      if (.not. engine%eom .and. .not. engine%lr) call output%error_msg('specify either eom og lr for fop.')
!
   end subroutine read_fop_settings_fop_engine


   subroutine do_eom_fop_engine(engine, wf, skip_states)
!!
!!    Do EOM
!!    Written by Written by Josefine H. Andersen, Sarai D. Folkestad 
!!    and Alexander C. Paul, June 2019
!!
      implicit none
!
      class(fop_engine) :: engine
      class(ccs)        :: wf
!
      logical, dimension(wf%n_singlet_states), intent(in) :: skip_states
!
      real(dp), dimension(:,:,:), allocatable :: operator
!
      real(dp), dimension(3) :: transition_strength, transition_moment_left, transition_moment_right
!
      integer :: state, component
!
      real(dp) :: trace_r_tdm, trace_l_tdm
      integer :: p
!
      type(sequential_file) :: left_file, right_file
      character(len=15)     :: name_left, name_right
!
      type(timings) :: EOM_timer
!
      EOM_timer = timings('Total time for EOM FOP')
!
      call EOM_timer%turn_on()
!
      call output%long_string_print('EOM first order properties calculation', fs='(/t3,a)', colons=.true.)
      call output%long_string_print(engine%author, fs='(t3,a)', colons=.true.)
!
      if (engine%dipole_length) then
!
         call mem%alloc(operator, wf%n_mo, wf%n_mo, 3)
!
!        Constructs dipole operator in t1-transformed basis.
         call wf%construct_mu(operator)
!
!        Loop over excited states, construct transition density
!        and calculate transition strength
!
         call output%printf('- Summary of EOM first order properties calculation:', fs='(/t3,a)')
!
         do state = 1, wf%n_singlet_states
!
            if(skip_states(state)) then
               call output%printf('Warning: Skipped state (i0) because it is &
               & parallel to the previous state', ints=[state], fs='(/t3,a)')
               cycle
            end if
!
            call wf%construct_right_transition_density(state)
!
            call wf%construct_left_transition_density(state)
!
!           Initialize and write the files to store the transition densities
!
            write(name_left, '(a, i3.3)') 'left_tdm_', state
            write(name_right, '(a, i3.3)') 'right_tdm_', state
!
            left_file  = sequential_file(trim(name_left))
            right_file = sequential_file(trim(name_right))
!
            call left_file%open_('write')
            call right_file%open_('write')
!
            call left_file%write_(wf%left_transition_density, wf%n_mo**2)
            call right_file%write_(wf%right_transition_density, wf%n_mo**2)
!
            call left_file%close_()
            call right_file%close_()
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
            call engine%print_summary_eom(transition_strength, transition_moment_left, &
                                          transition_moment_right, state,              &
                                          wf%right_excitation_energies(state))
!
!           Print traces of the density matrices for Debugging
!
            trace_r_tdm = zero
            trace_l_tdm = zero
!
            do p = 1, wf%n_mo
               trace_r_tdm = trace_r_tdm + wf%right_transition_density(p,p)
               trace_l_tdm = trace_l_tdm + wf%left_transition_density(p,p)
            end do
            call output%printf('Trace left transition density:  (f15.10)',    &
                                fs='(t6,a)', reals=[trace_l_tdm], pl='v')
            call output%printf('Trace right transition density: (f15.10)',    &
                                fs='(t6,a/)', reals=[trace_r_tdm], pl='v')
!
         enddo
!
      endif
!
      call EOM_timer%turn_off()
!
   end subroutine do_eom_fop_engine
!
!
   subroutine get_thresholds_fop_engine(energy_threshold, residual_threshold)
!!
!!    Get thresholds from input
!!    written by Alexander C. Paul, Oct 2019
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
   end subroutine get_thresholds_fop_engine
!
!
   subroutine print_summary_eom_fop_engine(transition_strength, transition_moment_left, &
                                           transition_moment_right, state, excitation_energy)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      real(dp), dimension(3), intent(in) :: transition_strength
      real(dp), dimension(3), intent(in) :: transition_moment_left 
      real(dp), dimension(3), intent(in) :: transition_moment_right
      real(dp), intent(in) :: excitation_energy
      integer, intent(in)  :: state
!
      character(len=1), dimension(3) :: components = ['X', 'Y', 'Z']
!
      integer :: component
!
      real(dp) :: sum_strength
!
      call output%printf('State (i0):', pl='m', fs='(/t6,a)', ints=[state])
      call output%print_separator('m', 10, '-', fs='(t6,a)')
!
      call output%printf('Excitation energy [E_h]:      (f19.12)', &
                         pl='m', reals=[excitation_energy], fs='(t6,a)')
      call output%printf('Excitation energy [eV]:       (f19.12)', &
                         pl='m', reals=[excitation_energy*Hartree_to_eV], fs='(t6,a)')
!
      call output%printf('Hartree-to-eV (CODATA 2014):  (f19.8)', &
                         pl='m', reals=[Hartree_to_eV], fs='(t6,a)')
!
      call output%printf('Transition moments [a.u.]         Transition strength [a.u.]', &
                         pl='m', fs='(/t6,14X,a)', ll=74)
      call output%print_separator('m', 74, '-', fs='(t6,a)')
      call output%printf('Comp. q     < k |q| 0 >       < 0 |q| k >        < 0 |q| k > < k |q| 0 >  ', &
                         pl='m', fs='(t6,a)', ll=74)
      call output%print_separator('m', 74, '-', fs='(t6,a)')
!
      sum_strength = zero
!
      do component = 1, 3
!
         call output%printf('(a0)      (f17.10) (f17.10)       (f17.10)',  &
                           reals=[transition_moment_left(component),       &
                           transition_moment_right(component),             &
                           transition_strength(component)],                &
                           chars=[components(component)], pl='m', fs='(t6,a)')
!
         sum_strength = sum_strength + transition_strength(component)
!
      enddo
!
      call output%print_separator('m', 74, '-', fs='(t6,a)')
!
      call output%printf('Oscillator strength [a.u.]: (f19.12)', &
                        reals=[(two/three)*excitation_energy*sum_strength], &
                        pl='m', fs='(t6,a)')
!
      flush(output%unit)
!
   end subroutine print_summary_eom_fop_engine
!
!
   subroutine set_printables_fop_engine(engine)
!!
!!    Set printables
!!    Written by sarai D. Folkestad, May 2019
!!
      implicit none
!
      class(fop_engine) :: engine
!
      character(len=5) :: fop_type
!
      engine%name_       = 'First order coupled cluster properties engine'
      engine%author      = 'J. H. Andersen, S. D. Folkestad, E. F. Kjønstad, A. Paul 2019'
!
      engine%tag = 'first order properties'
!
      if (engine%eom) then
!
         fop_type = 'EOM'
!
      else
!
         fop_type = 'LR'
!
      endif
!
      engine%tasks = [character(len=150) ::                                                                 &
      'Cholesky decomposition of the ERI-matrix',                                                           &
      'Calculation of the ground state amplitudes and energy ('//trim(engine%gs_algorithm)//'-algorithm)',  &
      'Calculation of the multipliers ('//trim(engine%multipliers_algorithm)//'-algorithm)',                &
      'Calculation of the ' //trim(engine%es_type) //' excitation vectors and&
      & energies ('//trim(engine%es_algorithm)//'-algorithm)',                                              &
      'Calculation of the first order property ('//trim(fop_type)//')']
!
      engine%description = 'Calculates dipole transition moments and oscillator strengths between &
                           &the ground state and the excited states.'
!
   end subroutine set_printables_fop_engine
!
!
end module fop_engine_class
