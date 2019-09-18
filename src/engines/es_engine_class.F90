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
   use gs_engine_class
!
   type, extends(gs_engine) :: es_engine
!
      character(len=200) :: es_algorithm
      character(len=200) :: es_type
      character(len=200) :: es_transformation
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
      engine%name_  = 'Excited state engine'
      engine%author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
!     Set standards and then read if nonstandard
!
      engine%es_algorithm        = 'davidson'
      engine%gs_algorithm        = 'diis'
      engine%es_type             = 'valence'
      engine%es_transformation   = 'right'
!
      call engine%read_settings()
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
      if (input%requested_keyword_in_section('left eigenvectors', 'solver cc es')) engine%es_transformation = 'left'
!
      if (input%requested_keyword_in_section('right eigenvectors', 'solver cc es')) engine%es_transformation = 'right'
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
!     Cholesky decomposition
!
      call engine%do_cholesky(wf, wf%orbital_coefficients)
!
!     Ground state solution
!
      call engine%do_ground_state(wf)
!
      call wf%integrals%write_t1_cholesky(wf%t1)
!      
      if (wf%need_g_abcd()) call wf%integrals%can_we_keep_g_pqrs_t1()
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
      use davidson_cc_es_class, only: davidson_cc_es
      use diis_cc_gs_class, only: diis_cc_gs
      use diis_cc_es_class, only: diis_cc_es
      use diis_A_inv_cc_es_class, only: diis_A_inv_cc_es
!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
!
      character(len=*), intent(in) :: transformation
!
      class(diis_cc_es), allocatable :: cc_es_solver_diis
!
      type(diis_A_inv_cc_es), allocatable :: cc_es_solver_A_inv
!
      class(davidson_cc_es), allocatable :: cc_es_solver_davidson
!
!     Prepare for excited state
!
      if (engine%es_algorithm == 'diis' .or. trim(wf%name_) == 'low memory cc2' .or. trim(wf%name_) == 'cc3') then
!
        if (trim(engine%es_algorithm) == 'davidson' .and. (trim(wf%name_) == 'low memory cc2' .or. trim(wf%name_) == 'cc3')) then
!            
            call output%warning_msg("for " // trim(wf%name_) // "excited states, the DIIS algorithm will be used, " // &
                                    "even though 'davidson' is default or was specified.")
!
        endif
!
        cc_es_solver_diis = diis_cc_es(transformation, wf)
        call cc_es_solver_diis%run(wf)
        call cc_es_solver_diis%cleanup(wf)
!
      elseif (engine%es_algorithm == 'diis a inverse') then
!
         cc_es_solver_A_inv = diis_A_inv_cc_es(transformation, wf)
         call cc_es_solver_A_inv%run(wf)
         call cc_es_solver_A_inv%cleanup(wf)
!
      elseif (engine%es_algorithm == 'davidson') then
!
         cc_es_solver_davidson = davidson_cc_es(transformation, wf)
         call cc_es_solver_davidson%run(wf)
         call cc_es_solver_davidson%cleanup(wf)
!
      else
!
         call output%error_msg('Could not start excited state solver. It may be that the &
                                 &algorithm is not implemented for the method specified.')
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
      implicit none
!
      class(es_engine) :: engine
!
      engine%tag = 'excited state'
!
      engine%tasks = [character(len=150) ::                                                              &
      'Cholesky decomposition of the ERI-matrix',                                                        &
      'Calculation of the ground state amplitudes ('//trim(engine%gs_algorithm)//'-algorithm)',          &
      'Calculation of the ground state energy',                                                          &
      'Calculation of the ' // trim(engine%es_transformation) //' hand side '// trim(engine%es_type) //  &
      ' excitation vectors ('//trim(engine%es_algorithm)//'-algorithm)',                                 &
      'Calculation of the excitation energies ('//trim(engine%es_algorithm)//'-algorithm)']
!
      engine%description  = 'Calculates the coupled cluster excitation vectors and excitation energies'
!
   end subroutine set_printables_es_engine
!
!
end module es_engine_class
