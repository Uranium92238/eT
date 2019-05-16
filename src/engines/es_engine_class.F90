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
      procedure :: prepare                   => prepare_es_engine
      procedure :: run                       => run_es_engine
!
      procedure :: read_settings             => read_settings_es_engine
      procedure :: read_es_settings          => read_es_settings_es_engine
!
      procedure :: do_excited_state          => do_excited_state_es_engine
!
      procedure, private :: set_printables   => set_printables_es_engine
!
   end type es_engine
!
contains
!
   subroutine prepare_es_engine(engine)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(es_engine) :: engine
!
      engine%name_  = 'Excited state engine'
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
   end subroutine prepare_es_engine
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
      if (input%requested_keyword_in_section('core excitation', 'solver cc es')) engine%es_type = 'core'
      if (input%requested_keyword_in_section('left eigenvectors', 'solver cc es')) engine%es_transformation = 'left'
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
      call wf%integrals%can_we_keep_g_pqrs_t1()
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
      use davidson_cc_es_class
      use davidson_cvs_cc_es_class
      use diis_cc_gs_class
      use diis_cc_es_class
!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
!
      character(len=*), intent(in) :: transformation
!
      type(diis_cc_es), allocatable                  :: cc_es_solver_diis
!
      type(davidson_cc_es), allocatable, target      ::  cc_valence_es
      type(davidson_cvs_cc_es), allocatable, target  ::  cc_core_es
!
      class(davidson_cc_es), pointer :: cc_es_solver
!
!     Prepare for excited state
!
      if (engine%es_algorithm == 'diis' .or. wf%name_ == 'low memory cc2' .or. wf%name_ == 'cc3') then
!
         allocate(cc_es_solver_diis)
!
         call cc_es_solver_diis%prepare(transformation)
         call cc_es_solver_diis%run(wf)
         call cc_es_solver_diis%cleanup()
!
         deallocate(cc_es_solver_diis)
!
      elseif (engine%es_algorithm == 'davidson') then
!
         if (engine%es_type == 'core') then
!
            allocate(cc_core_es)
            cc_es_solver => cc_core_es
!
            call cc_es_solver%prepare(transformation)
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_core_es)
!
         elseif(engine%es_type == 'valence ionized') then
!
            call output%error_msg('valence ionized not implemented yet')
!
         elseif(engine%es_type == 'core ionized') then
!
            call output%error_msg('core ionized not implemented yet')
!
         else ! es_type = valence
!
            allocate(cc_valence_es)
            cc_es_solver => cc_valence_es
!
            call cc_es_solver%prepare(transformation)
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_valence_es)
!
         endif
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
      engine%tag    = 'excited state'
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
