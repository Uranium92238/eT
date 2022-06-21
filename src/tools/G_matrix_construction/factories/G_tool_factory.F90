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
module G_tool_factory_class
!
!!
!! G tool factory class
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
   use G_adder_class,    only: G_adder
   use G_screener_class, only: G_screener
!
   implicit none
!
   type, extends(abstract_G_tool_factory) :: G_tool_factory
!
      real(dp) :: J_threshold, K_threshold
!
   contains
!
      procedure ::  create => create_G_tool_factory
!
   end type G_tool_factory
!
   interface G_tool_factory
      procedure :: new_G_tool_factory
   end interface G_tool_factory
!
contains
!
!
   pure function new_G_tool_factory(J_threshold, K_threshold) result(this)
!!
!!    New G-tool factory
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      type(G_tool_factory) :: this
!
      real(dp), intent(in) :: J_threshold, K_threshold
!
      this%J_threshold = J_threshold
      this%K_threshold = K_threshold
!
   end function new_G_tool_factory
!
   subroutine create_G_tool_factory(this, screener, adder)
!!
!!    Create G tool
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(G_tool_factory),                     intent(in)    :: this
      class(abstract_G_screener),   allocatable, intent(inout) :: screener
      class(abstract_G_adder),      allocatable, intent(inout) :: adder
!
      screener = G_screener(this%J_threshold, this%K_threshold)
      adder = G_adder()
!
   end subroutine create_G_tool_factory
!
!
end module G_tool_factory_class
