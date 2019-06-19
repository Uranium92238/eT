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
module hf_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use hf_class 
!
   use scf_diis_hf_class
   use scf_hf_class
!
   type hf_engine 
!
      character(len=200) :: algorithm 
      logical :: restart 
!
   contains 
!
      procedure :: ignite                    => ignite_hf_engine
!
      procedure, private :: run              => run_hf_engine
!
      procedure, private :: read_settings    => read_settings_hf_engine
!
   end type hf_engine 
!
!
   interface hf_engine
!
      procedure :: new_hf_engine 
!
   end interface hf_engine
!
!
contains
!
!
   function new_hf_engine() result(engine)
!!
!!    New HF engine 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      type(hf_engine) :: engine
!
      engine%algorithm = 'scf-diis'
      engine%restart = .false.
!
      call engine%read_settings()
!
   end function new_hf_engine
!
!
   subroutine run_hf_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine)  :: engine 
      class(hf)         :: wf 
!
      type(scf_hf), allocatable :: scf
      type(scf_diis_hf), allocatable :: scf_diis
!
      if (trim(engine%algorithm) == 'scf-diis') then
!
         scf_diis = scf_diis_hf(wf, engine%restart)
         call scf_diis%run(wf)
         call scf_diis%cleanup(wf)
!
      elseif (trim(engine%algorithm) == 'scf') then 
!
         if (engine%restart) call output%error_msg('SCF does not support restart.')
!
         scf = scf_hf(wf)
         call scf%run(wf)
         call scf%cleanup(wf)
!
      else
!
         call output%error_msg('did not recognize hf algorithm: '// engine%algorithm)
!
      endif
!
   end subroutine run_hf_engine
!
!
   subroutine read_settings_hf_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none
!
      class(hf_engine) :: engine 
!
      call input%get_keyword_in_section('algorithm', 'solver hf', engine%algorithm)
      if (input%requested_keyword_in_section('restart', 'solver hf')) engine%restart = .true.
!
   end subroutine read_settings_hf_engine
!
!
   subroutine ignite_hf_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none 
!
      class(hf_engine) :: engine 
      class(hf)        :: wf 
!
      call engine%run(wf)
!
   end subroutine ignite_hf_engine
!
!
end module hf_engine_class
