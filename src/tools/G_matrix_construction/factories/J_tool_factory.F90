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
module J_tool_factory_class
!
!!
!! J tool factory class
!! Written by Sarai D. Folkestad, 2021
!!
!
   use kinds
!
   use abstract_G_tool_factory_class, only: abstract_G_tool_factory
!
   use abstract_G_adder_class,    only: abstract_G_adder
   use abstract_G_screener_class, only: abstract_G_screener
!
   use J_adder_class,    only: J_adder
   use J_screener_class, only: J_screener
!
   implicit none
!
   type, extends(abstract_G_tool_factory) :: J_tool_factory
!
      real(dp) :: J_threshold
!
   contains
!
      procedure ::  create => create_J_tool_factory
!
   end type J_tool_factory
!
   interface J_tool_factory
      procedure :: new_J_tool_factory
   end interface J_tool_factory
!
contains
!
!
   pure function new_J_tool_factory(J_threshold) result(this)
!!
!!    New J tool Factory
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      type(J_tool_factory) :: this
!
      real(dp), intent(in) :: J_threshold
!
      this%J_threshold = J_threshold
!
   end function new_J_tool_factory
!
   subroutine create_J_tool_factory(this, screener, adder)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(J_tool_factory),                     intent(in)    :: this
      class(abstract_G_screener),   allocatable, intent(inout) :: screener
      class(abstract_G_adder),      allocatable, intent(inout) :: adder
!
      screener = J_screener(this%J_threshold)
      adder = J_adder()
!
   end subroutine create_J_tool_factory
!
!
end module J_tool_factory_class
