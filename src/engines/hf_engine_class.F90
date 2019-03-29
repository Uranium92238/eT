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
!
   contains 
!
      procedure         :: prepare        => prepare_hf_engine
      procedure         :: run            => run_hf_engine
      procedure, nopass :: cleanup        => cleanup_hf_engine
!
      procedure         :: read_settings => read_settings_hf_engine
!
   end type hf_engine 
!
contains
!
!
   subroutine prepare_hf_engine(engine)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine) :: engine 
!
      engine%algorithm = 'scf-diis'
!
      call engine%read_settings()
!
   end subroutine prepare_hf_engine
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
      type(scf_diis_hf), allocatable :: scf_diis
      type(scf_hf), allocatable      :: scf
!
      if (trim(engine%algorithm) == 'scf-diis') then
!
         allocate(scf_diis)
!
         call scf_diis%prepare(wf)
         call scf_diis%run(wf)
         call scf_diis%cleanup(wf)
!
         deallocate(scf_diis)
!
      elseif (trim(engine%algorithm) == 'scf') then 
!
         allocate(scf)
!
         call scf%prepare(wf)
         call scf%run(wf)
         call scf%cleanup(wf)
!
         deallocate(scf)
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
   subroutine cleanup_hf_engine()
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
   end subroutine cleanup_hf_engine
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
      if (input%requested_section('hf')) then 
!
         call input%read_keyword_in_section('algorithm', 'hf', engine%algorithm)
!
      endif
!
   end subroutine read_settings_hf_engine
!
!
end module hf_engine_class
