!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module cc_engine_factory_class
!
!!
!!    CC engine factory class
!!    Written by Eirik F. Kjønstad, May 2022
!!
!
   use cc_engine_class, only: cc_engine
!
   implicit none 
!
   type :: cc_engine_factory
!
   contains
!
      procedure, public, nopass :: create & 
                                => create_cc_engine_factory
!
   end type cc_engine_factory
!
!
contains 
!
!
   subroutine create_cc_engine_factory(engine)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use cc_engine_class,                   only: cc_engine
      use cc_gs_engine_class,                only: cc_gs_engine
      use cc_es_engine_class,                only: cc_es_engine
      use cc_gs_mean_value_engine_class,     only: cc_gs_mean_value_engine
      use td_cc_engine_class,                only: td_cc_engine
!
      use cc_response_engine_factory_class,  only: cc_response_engine_factory
!
      use global_out,                        only: output 
      use global_in,                         only: input
!
      implicit none 
!
      class(cc_engine), allocatable, intent(out) :: engine  
!
      type(cc_response_engine_factory), allocatable :: response_engine_factory
!
      if (input%is_keyword_present('response', 'do')) then
!
         allocate(cc_response_engine_factory::response_engine_factory)
!
         call response_engine_factory%create(engine)
!
      elseif (input%is_keyword_present('excited state', 'do')) then
!
         allocate(cc_es_engine::engine)
!
      elseif (input%is_keyword_present('mean value', 'do')) then
!
         allocate(cc_gs_mean_value_engine::engine)
!
      elseif (input%is_keyword_present('ground state', 'do')) then
!
         allocate(cc_gs_engine::engine)
!
      elseif (input%is_keyword_present('time dependent state', 'do')) then
!
         allocate(td_cc_engine::engine)
!
      else
!
         call output%error_msg('could not recognize coupled cluster task')
!
      endif
!
   end subroutine create_cc_engine_factory
!
!
end module cc_engine_factory_class
