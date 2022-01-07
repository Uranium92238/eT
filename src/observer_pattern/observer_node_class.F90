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
module observer_node_class
!
!!
!! Observer node class
!! Written by Sarai D. Folkestad, Oct 2021
!!
!! Part of the general implementation of the
!! Observer pattern
!!
!! A concrete observer_node must implement an 'update'
!!
!
   use parameters
!
   use global_out, only: output
!
   use observer_class, only: observer
!
   implicit none
!
   type :: observer_node
!
      class(observer_node), pointer :: next => null()
      class(observer_node), pointer :: previous => null()
!
      character(len=200), private :: tag
!
      class(observer), pointer :: observer_
!
   contains
!
      procedure :: get_tag => get_tag_observer_node
      procedure :: update => update_observer_node
!
   end type observer_node   
!
   interface observer_node
!
      procedure :: new_observer_node
!
   end interface observer_node
!
contains
!
   function new_observer_node(observer_, tag) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observer), intent(inout), target :: observer_
!
      character(len=*), intent(in) :: tag
!
      type(observer_node) :: this
!
      this%observer_ => observer_
      this%tag = trim(tag)
!
   end function new_observer_node
!
   function get_tag_observer_node(this) result(tag)
!!
!!    Get tag
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observer_node), intent(in) :: this
!
      character(len=200) :: tag
!
      tag = this%tag
!
   end function get_tag_observer_node
!
!
   subroutine update_observer_node(this)
!!
!!    Update
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observer_node) :: this
!
      call this%observer_%update()
!
   end subroutine update_observer_node
!
end module observer_node_class
