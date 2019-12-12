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
module es_engine_class
!!
!!    Coupled cluster ground state engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
   use kinds
   use global_in,       only: input
   use global_out,      only: output
   use timings_class,   only: timings
   use task_list_class, only: task_list
!
   use gs_engine_class, only: gs_engine
   use ccs_class,       only: ccs
!
   type, extends(gs_engine) :: es_engine
!
      character(len=200) :: es_algorithm
      character(len=200) :: es_type
      character(len=200) :: es_transformation
!
      logical :: es_restart
!
   contains
!
      procedure :: run                       => run_es_engine
!
      procedure :: read_settings             => read_settings_es_engine
      procedure :: read_es_settings          => read_es_settings_es_engine
!
      procedure :: do_excited_state          => do_excited_state_es_engine
!
      procedure :: set_printables            => set_printables_es_engine
!
      procedure :: restart_handling          => restart_handling_es_engine
!
   end type es_engine
!
!
   interface es_engine
!
      procedure :: new_es_engine 
!
   end interface es_engine
!
!
contains
!
!
   function new_es_engine() result(engine)
!!
!!    New ES engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(es_engine) :: engine
!
!     Set standards and then read if nonstandard
!
      engine%es_algorithm           = 'davidson'
      engine%gs_algorithm           = 'diis'
      engine%es_type                = 'valence'
      engine%es_transformation      = 'right'
      engine%multipliers_algorithm  = 'davidson'
!
      engine%gs_restart            = .false.
      engine%multipliers_restart   = .false.
      engine%es_restart            = .false.
!
      call engine%read_settings()
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
   end function new_es_engine
!
!
   subroutine read_settings_es_engine(engine)
!!
!!    Read es settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(es_engine) :: engine
!
      call engine%read_gs_settings()
      call engine%read_es_settings()
!
   end subroutine read_settings_es_engine
!
!
   subroutine read_es_settings_es_engine(engine)
!!
!!    Read es settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(es_engine) :: engine
!
      call input%get_keyword_in_section('algorithm', 'solver cc es', engine%es_algorithm)
!
      if (input%requested_keyword_in_section('core excitation', 'solver cc es') .and. .not. &
          input%requested_keyword_in_section('ionization', 'solver cc es')) engine%es_type = 'core'
!
      if (input%requested_keyword_in_section('ionization', 'solver cc es') .and. .not. &
         input%requested_keyword_in_section('core excitation', 'solver cc es')) engine%es_type = 'ionize'
!
      if (input%requested_keyword_in_section('ionization', 'solver cc es') .and.    &
         input%requested_keyword_in_section('core excitation', 'solver cc es'))     &
            call output%error_msg('XPS still not implemented.')
!
      if (input%requested_keyword_in_section('left eigenvectors', 'solver cc es')) &
                                                engine%es_transformation = 'left'
!
      if (input%requested_keyword_in_section('right eigenvectors', 'solver cc es')) &
                                                engine%es_transformation = 'right'
!
      engine%es_restart = input%requested_keyword_in_section('restart', 'solver cc es')
!
   end subroutine read_es_settings_es_engine
!
!
   subroutine run_es_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
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
!     Excited state solutions
!
      call engine%do_excited_state(wf, engine%es_transformation)
!
   end subroutine run_es_engine
!
!
   subroutine do_excited_state_es_engine(engine, wf, transformation)
!!
!!    Do excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Solves the excited state (valence or cvs) using
!!    either a DIIS or Davidson solver
!!
!!    or calculates the excitation energies and oscillator strengths 
!!    using the asymmetric Lanczos algorithm. 
!!
!!    Modified by Torsha Moitra, S. Coriani and Sarai D. Folkestad, Sep-Nov 2019
!!
!!       Added the asymmetric Lanczos solver 
!!
      use abstract_cc_es_class, only: abstract_cc_es
      use davidson_cc_es_class, only: davidson_cc_es
      use diis_cc_es_class, only: diis_cc_es
      use asymmetric_lanczos_cc_es_class, only: asymmetric_lanczos_cc_es
!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
!
      character(len=*), intent(in) :: transformation
!
      class(asymmetric_lanczos_cc_es), allocatable :: cc_es_solver_asymmetric_lanczos

      class(abstract_cc_es), allocatable :: cc_es_solver
!
      call engine%tasks%print_('es solver', &
            append_string='Calculating ' // trim(transformation) //' vectors', append_fs='(t6,a)')
!
!     Prepare for excited state
!
      if (engine%es_algorithm == 'asymmetric lanczos') then
!
         call engine%do_multipliers(wf)       
         cc_es_solver_asymmetric_lanczos = asymmetric_lanczos_cc_es(wf)
         call cc_es_solver_asymmetric_lanczos%run(wf)
         call cc_es_solver_asymmetric_lanczos%cleanup(wf)
!
      else
!
         if (engine%es_algorithm == 'diis') then
!
            cc_es_solver = diis_cc_es(transformation, wf, engine%es_restart)

! 
         elseif (engine%es_algorithm == 'davidson') then
!
            if (trim(wf%name_) == 'low memory cc2' .or. trim(wf%name_) == 'cc3') then
!
                call output%error_msg('Davidson not implemented for CC3 and lowmem CC2')
!
            end if
!
            cc_es_solver = davidson_cc_es(transformation, wf, engine%es_restart)
!
         else
               call output%error_msg('Could not start excited state solver. It may be that the &
                                    &algorithm is not implemented for the method specified.')
         endif
!
         call cc_es_solver%run(wf)
         call cc_es_solver%cleanup(wf)
!
      endif
!
   end subroutine do_excited_state_es_engine
!
!
   subroutine set_printables_es_engine(engine)
!!
!!    Set printables
!!    Written by sarai D. Folkestad, May 2019
!!
!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(es_engine) :: engine
!
      engine%name_  = 'Excited state coupled cluster engine'
!
      engine%tag = 'excited state'
!
!     Prepare the list of tasks
!
      engine%tasks = task_list()
!
      call engine%tasks%add(label='mo preparations',                             &
                            description='Preparation of MO basis and integrals')
!
      call engine%tasks%add(label='gs solver',                                &
                           description='Calculation of the ground state ('//  &
                           trim(engine%gs_algorithm)//' algorithm)')
!
      if (engine%es_algorithm == 'asymmetric lanczos') then
!
         call engine%tasks%add(label='multipliers solver',                    &
                           description='Calculation of the multipliers ('     &
                           //trim(engine%multipliers_algorithm)&
                           //' algorithm)')
!
         call engine%tasks%add(label='es solver',                             &
                           description='Calculation of the excited state ('// &
                           trim(engine%es_algorithm)//' algorithm)')
!
      else
!
         call engine%tasks%add(label='es solver',                              &
                           description='Calculation of the excited state ('//  &
                           trim((engine%es_algorithm))//' algorithm)')
!
      endif
!
      engine%description  = 'Calculates the coupled cluster excitation vectors and excitation energies'
!
   end subroutine set_printables_es_engine
!
!
   subroutine restart_handling_es_engine(engine, wf)
!!
!!    Restart handling
!!    Written by Sarai D. Folkestad, Nov 2019
!!
!!    Writes the restart information 
!!    if restart is not requested.
!!
!!    If restart is requested performs safety 
!!    checks for restart
!!
      implicit none
!
      class(es_engine), intent(in) :: engine
      class(ccs), intent(in) :: wf
!
      if (.not. engine%restart) then
!
         call wf%write_cc_restart()
!
      else
!
         if (engine%gs_restart .or. engine%multipliers_restart) &
                                 call wf%is_restart_safe('ground state')
!
         if (engine%es_restart) call wf%is_restart_safe('excited state')
!
      endif
!
   end subroutine restart_handling_es_engine
!
!
end module es_engine_class
