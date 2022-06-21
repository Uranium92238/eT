!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module abstract_G_tool_factory_class
!
!!
!! Abstract G tool factory class
!! Written by Sarai D. Folkestad, 2021
!!
!
   use kinds
!
   use abstract_G_adder_class,    only: abstract_G_adder
   use abstract_G_screener_class, only: abstract_G_screener
!
   implicit none
!
   type, abstract :: abstract_G_tool_factory
!
   contains
!
      procedure(create_abstract), deferred :: create
!
   end type abstract_G_tool_factory
!
   abstract interface
!
      subroutine create_abstract(this, screener, adder)
!
         import abstract_G_tool_factory
         import abstract_G_adder
         import abstract_G_screener
!
         implicit none
!
         class(abstract_G_tool_factory)           , intent(in)    :: this
         class(abstract_G_screener),   allocatable, intent(inout) :: screener
         class(abstract_G_adder),      allocatable, intent(inout) :: adder
!
      end subroutine create_abstract
!
   end interface
!
!
end module abstract_G_tool_factory_class
