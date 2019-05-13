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
   use es_engine_class
!
   type, extends(es_engine) :: fop_engine
!
      character(len=100) :: tag           = 'First order coupled cluster properties'
      character(len=100) :: author        = 'J. H. Andersen, S. D. Folkestad, E. F. Kjønstad, 2019'
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
      procedure :: prepare             => prepare_fop_engine
      procedure :: run                 => run_fop_engine
!
      procedure :: read_settings       => read_settings_fop_engine
      procedure :: read_fop_settings   => read_fop_settings_fop_engine
!
      procedure :: do_eom              => do_eom_fop_engine
!
      procedure, nopass :: print_summary_eom   => print_summary_eom_fop_engine
!
   end type fop_engine
!
contains
!
!
   subroutine prepare_fop_engine(engine)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(fop_engine) :: engine
!
      engine%name_       = 'First order properties engine'
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
   end subroutine prepare_fop_engine
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
!     Cholesky decomposition
!
      call engine%do_cholesky(wf, wf%orbital_coefficients)
!
!     Ground state solution
!
      call engine%do_ground_state(wf)
!
!     Prepare for excited state calculation
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      call wf%integrals%can_we_keep_g_pqrs_t1()
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
!     EOM properties if requested
!
      if (engine%eom) call engine%do_eom(wf)
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
!
!
   subroutine do_eom_fop_engine(engine, wf)
!!
!!    Do EOM
!!    Written by Josefine H. Andersen and Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(fop_engine) :: engine
      class(ccs)         :: wf
!
      real(dp), dimension(:,:,:), allocatable :: operator
!
      character(len=1), dimension(3) :: components = ['X', 'Y', 'Z']
!
      real(dp), dimension(:,:), allocatable :: transition_strength, transition_moment_left, transition_moment_right
!
      real(dp), dimension(:), allocatable :: etaX, csiX, excitation_energies
!
      integer :: component, n_states, state
!
      call wf%prepare_for_eom_fop()
!
      call mem%alloc(etaX, wf%n_es_amplitudes)
      call mem%alloc(csiX, wf%n_es_amplitudes)
!
!     Print banner
!
      engine%description = 'Calculates dipole transition moments and oscillator strengths between &
                              &the ground state and the excited states.'
!
      call long_string_print(engine%tag,'(//t3,a)',.true.)
      call long_string_print(engine%author,'(t3,a/)',.true.)
      call long_string_print(engine%description,'(t3,a)',.false.,'(t3,a)','(t3,a)')
!
      call input%get_required_keyword_in_section('singlet states', 'solver cc es', n_states)
      call mem%alloc(excitation_energies, n_states)
!
      call wf%read_excitation_energies(n_states, excitation_energies)
!
      call mem%alloc(transition_strength, 3, n_states)
      call mem%alloc(transition_moment_left, 3, n_states)
      call mem%alloc(transition_moment_right, 3, n_states)
!
      transition_strength     = zero
      transition_moment_right = zero
      transition_moment_left  = zero
!
      if (engine%dipole_length) then
!
         call mem%alloc(operator, wf%n_mo, wf%n_mo, 3)
!
         call wf%construct_mu(operator)  ! Constructs dipole operator in t1-transformed basis.
!
         do component = 1, size(components)
!
            call wf%construct_csiX(operator(:,:,component), csiX)      
!
            call wf%construct_eom_etaX(operator(:,:,component), csiX, etaX)
!
!           Loop over excited states and calculate transition strength
!  
            do state = 1, n_states
!
               call wf%calculate_transition_strength(transition_strength(component, state), etaX, &
                   csiX, state, transition_moment_left(component, state), transition_moment_right(component, state))
!    
            enddo
!
         enddo
!
         call engine%print_summary_eom(transition_strength, transition_moment_left, &
                                    transition_moment_right, n_states, excitation_energies)
!
      endif
!
      call mem%dealloc(transition_strength, 3, n_states)
      call mem%dealloc(transition_moment_left, 3, n_states)
      call mem%dealloc(transition_moment_right, 3, n_states)      
!
   end subroutine do_eom_fop_engine
!
!
   subroutine print_summary_eom_fop_engine(transition_strength, transition_moment_left, &
                                    transition_moment_right, n_states, excitation_energies)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      integer, intent(in)  :: n_states
!
      real(dp), dimension(3, n_states), intent(in) :: transition_strength, transition_moment_left, transition_moment_right
      real(dp), dimension(n_states), intent(in) :: excitation_energies
!
      character(len=1), dimension(3) :: components = ['X', 'Y', 'Z']
!
      integer :: component, state
!
      real(dp) :: sum_strength
!              
      write(output%unit, '(/t3,a)') '- Summary of EOM first order properties calculation:'
!
      do state = 1, n_states
!
         write(output%unit, '(/t6, a6 ,i3, a1)') 'State ', state, ':'
         write(output%unit, '(t6, a)') '----------'
         write(output%unit, '(t6, a30, f19.12)') 'Excitation energy [E_h]:      ', excitation_energies(state)
         write(output%unit, '(t6, a30, f19.12)') 'Excitation energy [eV]:       ', excitation_energies(state)*Hartree_to_eV
!
         write(output%unit, '(t6, a30, f19.8)')  'Hartree-to-eV (CODATA 2014):  ', Hartree_to_eV
!
         write(output%unit, '(/t6,a)')  '                 Transition moments               Transition strength   '    
         write(output%unit, '(t6,a)')   '------------------------------------------------------------------------'
         write(output%unit, '(t6,a)')   'Comp. q     < k |q| 0 >       < 0 |q| k >       < k |q| 0 > < 0 |q| k > '
         write(output%unit, '(t6,a)')   '------------------------------------------------------------------------'
!
         sum_strength = zero
!
         do component = 1, 3
!
            write(output%unit, '(t6,a1,6x,f17.10,1x,f17.10,6x,f17.10)') components(component),                       &
                                                                        transition_moment_left(component, state),    &
                                                                        transition_moment_right(component, state),   &
                                                                        transition_strength(component, state)
!
            sum_strength = sum_strength + transition_strength(component, state)
!
         enddo
!   
         write(output%unit, '(t6,a)')   '-------------------------------------------------------------------------'
!
         write(output%unit, '(t6, a21, f19.12)') 'Oscillator strength: ', (two/three)*excitation_energies(state)*sum_strength
!
!
      enddo
!
   end subroutine print_summary_eom_fop_engine
!
!
end module fop_engine_class
