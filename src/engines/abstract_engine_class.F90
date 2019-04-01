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
module abstract_engine_class
!
!!
!!    Abstract engine class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use ccs_class
!
   implicit none
!
   type, abstract :: abstract_engine
!
      character(len=200) :: name_
!
   contains
!
      procedure :: ignite => ignite_abstract_engine
!
      procedure(essential_engine), deferred      :: prepare 
      procedure(essential_engine), deferred      :: cleanup   
      procedure(essential_engine_w_wf), deferred :: run        
!
   end type abstract_engine
!
!
   abstract interface
!
      subroutine essential_engine(engine)
!
         import :: abstract_engine
!
         implicit none 
!
         class(abstract_engine) :: engine
!
      end subroutine essential_engine
!
!
      subroutine essential_engine_w_wf(engine, wf)
!
         import :: abstract_engine, ccs
!
         implicit none 
!
         class(abstract_engine) :: engine
!
         class(ccs) :: wf
!
      end subroutine essential_engine_w_wf
!
!
   end interface
!
contains
!
!
   subroutine ignite_abstract_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
!!    Prepare, run, cleanup
!!
      implicit none
!
      call engine%prepare()
      call engine%run(wf)
      call engine%cleanup()
!
   end subroutine ignite_abstract_engine
!
!
end module abstract_engine_class
!
