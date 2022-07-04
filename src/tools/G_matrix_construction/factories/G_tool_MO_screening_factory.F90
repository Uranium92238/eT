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
module G_tool_MO_screening_factory_class
!
!!
!! G tool MO screening factory class
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
   use G_MO_screener_class, only: G_MO_screener
!
   implicit none
!
   type, extends(abstract_G_tool_factory) :: G_tool_MO_screening_factory
!
      real(dp) :: J_threshold, K_threshold
!
      integer :: n_mo, n_ao
!
      real(dp), dimension(:,:), pointer      :: C     
!
   contains
!
      procedure ::  create => create_G_tool_MO_screening_factory
!
   end type G_tool_MO_screening_factory
!
   interface G_tool_MO_screening_factory
      procedure :: new_G_tool_MO_screening_factory
   end interface G_tool_MO_screening_factory
!
contains
!
!
   function new_G_tool_MO_screening_factory(J_threshold, K_threshold, n_ao, n_mo, C) result(this)
!!
!!    New G tool MO screening factory
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      type(G_tool_MO_screening_factory) :: this
!
!
      real(dp),                        intent(in)         :: J_threshold, K_threshold
      integer,                         intent(in)         :: n_ao, n_mo
      real(dp), dimension(n_ao, n_mo), intent(in), target :: C
!
      this%J_threshold = J_threshold
      this%K_threshold = K_threshold
!
      this%n_mo = n_mo
      this%n_ao = n_ao
!
      this%C(1:this%n_ao, 1:this%n_mo) => C(:,:)
!
   end function new_G_tool_MO_screening_factory
!
   subroutine create_G_tool_MO_screening_factory(this, screener, adder)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2021
!!
!!
      implicit none
!
      class(G_tool_MO_screening_factory),        intent(in)    :: this
      class(abstract_G_screener),   allocatable, intent(inout) :: screener
      class(abstract_G_adder),      allocatable, intent(inout) :: adder
!
      screener = G_MO_screener(this%J_threshold, this%K_threshold, this%n_ao, this%n_mo, this%C)
      adder = G_adder()
!
   end subroutine create_G_tool_MO_screening_factory
!
!
end module G_tool_MO_screening_factory_class
