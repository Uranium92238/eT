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
module hf_geoopt_engine_class
!!
!!    Hartree-Fock geometry optimization engine class module 
!!    Written by Eirik F. Kjønstad, 2019 
!!
   use hf_class 
!
   use bfgs_geoopt_hf_class
!
   type hf_geoopt_engine
!
      character(len=200) :: algorithm 
      logical :: restart 
!
   contains 
!
      procedure :: ignite                    => ignite_hf_geoopt_engine
!
      procedure, private :: run              => run_hf_geoopt_engine
!
      procedure, private :: read_settings    => read_settings_hf_geoopt_engine
!
   end type hf_geoopt_engine 
!
!
   interface hf_geoopt_engine
!
      procedure :: new_hf_geoopt_engine 
!
   end interface hf_geoopt_engine
!
!
contains
!
!
   function new_hf_geoopt_engine() result(engine)
!!
!!    New HF engine 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      type(hf_geoopt_engine) :: engine
!
      engine%algorithm = 'bfgs'
      engine%restart = .false.
!
      call engine%read_settings()
!
   end function new_hf_geoopt_engine
!
!
   subroutine run_hf_geoopt_engine(engine, wf)
!!
!!    Run 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      class(hf_geoopt_engine) :: engine 
      class(hf)               :: wf 
!
      type(bfgs_geoopt_hf) :: bfgs_geoopt
!
      if (trim(engine%algorithm) == 'bfgs') then
!
         bfgs_geoopt = bfgs_geoopt_hf(engine%restart)
         call bfgs_geoopt%run(wf)
         call bfgs_geoopt%cleanup()
!
      else
!
         call output%error_msg('did not recognize hf geoopt algorithm: '// engine%algorithm)
!
      endif
!
   end subroutine run_hf_geoopt_engine
!
!
   subroutine read_settings_hf_geoopt_engine(engine)
!!
!!    Read settings 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none
!
      class(hf_geoopt_engine) :: engine 
!
      call input%get_keyword_in_section('algorithm', 'solver hf geoopt', engine%algorithm)
      if (input%requested_keyword_in_section('restart', 'solver hf geoopt')) engine%restart = .true.
!
   end subroutine read_settings_hf_geoopt_engine
!
!
   subroutine ignite_hf_geoopt_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      class(hf_geoopt_engine) :: engine 
      class(hf)        :: wf 
!
      call engine%run(wf)
!
   end subroutine ignite_hf_geoopt_engine
!
!
end module hf_geoopt_engine_class
