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
module reference_engine_factory_class
!
!!
!!    Reference engine factory class
!!    Written by Eirik F. Kjønstad, May 2022
!!
!
   implicit none 
!
   type :: reference_engine_factory
!
   contains
!
      procedure, public, nopass :: create & 
                                => create_reference_engine_factory
!
   end type reference_engine_factory
!
!
contains 
!
!
   subroutine create_reference_engine_factory(engine)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use reference_engine_class, only: reference_engine
      use hf_geoopt_engine_class, only: hf_geoopt_engine
      use tdhf_engine_class,      only: tdhf_engine
!
      use global_in,              only: input
!
      implicit none 
!
      class(reference_engine), allocatable, intent(out) :: engine  
!
      if (input%is_keyword_present('ground state geoopt', 'do')) then
!
         engine = hf_geoopt_engine()
!
      elseif (input%is_keyword_present('time dependent hf', 'do')) then
!
         engine = tdhf_engine()
!
      else
!
         engine = reference_engine()
!
      endif
!
   end subroutine create_reference_engine_factory
!
!
end module reference_engine_factory_class
