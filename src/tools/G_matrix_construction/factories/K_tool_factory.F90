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
module K_tool_factory_class
!
!!
!! K tool factory class
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
   use K_adder_class,    only: K_adder
   use K_screener_class, only: K_screener
!
   implicit none
!
   type, extends(abstract_G_tool_factory) :: K_tool_factory
!
      real(dp) :: K_threshold
!
   contains
!
      procedure ::  create => create_K_tool_factory
!
   end type K_tool_factory
!
   interface K_tool_factory
      procedure :: new_K_tool_factory
   end interface K_tool_factory
!
contains
!
!
   pure function new_K_tool_factory(K_threshold) result(this)
!!
!!    New K tool factory
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      type(K_tool_factory) :: this
!
      real(dp), intent(in) :: K_threshold
!
      this%K_threshold = K_threshold
!
   end function new_K_tool_factory
!
   subroutine create_K_tool_factory(this, screener, adder)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(K_tool_factory),                     intent(in)    :: this
      class(abstract_G_screener),   allocatable, intent(inout) :: screener
      class(abstract_G_adder),      allocatable, intent(inout) :: adder
!
      screener = K_screener(this%K_threshold)
      adder = K_adder()
!
   end subroutine create_K_tool_factory
!
!
end module K_tool_factory_class
