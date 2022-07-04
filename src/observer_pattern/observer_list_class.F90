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
module observer_list_class
!
!!
!! Observer list class
!! Written by Sarai D. Folkestad, Oct 2021
!!
!! Part of the general implementation of the
!! Observer design pattern
!!
!
   use parameters
!
   use observer_node_class, only: observer_node
!
   implicit none
!
   type :: observer_list
!
      integer :: n = 0
!
      class(observer_node), pointer :: head => null()
      class(observer_node), pointer :: tail => null()
!
   contains
!
      procedure :: update  => update_observer_list
      procedure :: add     => add_observer_list
      procedure :: remove  => remove_observer_list
!
      procedure, private :: tag_taken
      procedure, private :: get_node
!
      procedure :: remove_all_observers
!
      final :: cleanup
!
   end type observer_list
!
!
contains
!
!
   subroutine update_observer_list(this)
!!
!!    Update
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observer_list) :: this
!
      integer :: i
      class(observer_node), pointer :: node
!
      do i = 1, this%n
!
         call this%get_node(i, node)
         call node%update()
!
      enddo
!
   end subroutine update_observer_list
!
!
   subroutine add_observer_list(this, tag, o)
!!
!!    Add observer
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      use global_out, only : output
      use observer_class, only: observer
!
      implicit none
!
      class(observer_list) :: this
!
      class(observer), target, intent(inout) :: o
      character(len=*) :: tag
!
      type(observer_node), pointer :: node
!
      if (this%tag_taken(tag)) &
         call output%error_msg('tried to add observer with same tag as another observer!')
!
      allocate(node)
      node = observer_node(o, tag)
!
      if (associated(this%head)) then
!
         this%tail%next => node
         this%tail%next%previous => this%tail
!
         this%tail => this%tail%next
         this%tail%next => null()
!
      else
!
         this%head => node
         this%tail => this%head
!
         this%head%previous => null()
         this%head%next => null()
!
         this%tail%previous => null()
         this%tail%next => null()
!
      endif
!
      this%n = this%n + 1
!
   end subroutine add_observer_list
!
!
   subroutine remove_observer_list(this, tag)
!!
!!    Remove observer
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      use global_out, only : output
!
      implicit none
!
      class(observer_list) :: this
!
      integer :: i
      character(len=*) :: tag
!
      class(observer_node), pointer :: node, node_previous, node_next
!
      if (.not. this%tag_taken(tag)) &
         call output%error_msg('tried to remove non-existing observer!')
!
      do i = 1, this%n
!
         call this%get_node(i, node)

         if (trim(node%get_tag()) == trim(tag)) then
!
            if (i == this%n) then ! remove at the end
!
!
               this%tail => null()
               this%head => null()
!
            else ! remove inside list
!
               node_next => node%next
               node_previous => node%previous
!
               node_previous%next => node_next
               node_next%previous => node_previous
!
            endif
!
            node%observer_ => null()
            deallocate(node)
!
            this%n = this%n - 1
!
         endif
      enddo
!
   end subroutine remove_observer_list
!
!
   subroutine remove_all_observers(this)
!!
!!    Remove all observers
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observer_list), intent(inout) :: this
!
      class(observer_node), pointer :: node
!
      do while (this%n > 0)
!
         call this%get_node(this%n, node)
!
         if (this%n > 1) then
!
            this%tail =>  this%tail%previous
            this%tail%next => null()
!
         else ! remove last observer
!
            this%tail => null()
            this%head => null()
!
         endif
!
         node%observer_ => null()
         deallocate(node)
!
         this%n = this%n - 1
!
      enddo
!
   end subroutine remove_all_observers


   subroutine cleanup(this)
!!
!!    Remove all observers
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      type(observer_list) :: this
!
      call this%remove_all_observers()
!
   end subroutine cleanup
!
!
   function tag_taken(this, tag) result(taken)
!!
!!    Tag taken
!!    Written by Sarai D. Folkestad, Oct 2021
!!
!!
      implicit none
!
      class(observer_list) :: this
      character(len=*) :: tag
!
      character(len=200) :: current_tag
!
      logical :: taken
!
      integer :: i
      class(observer_node), pointer :: node
!
      taken = .false.
      do i = 1, this%n
!
         call this%get_node(i, node)
         current_tag = node%get_tag()
!
         if (trim(current_tag) == trim(tag)) taken = .true.
!
      enddo
!
   end function tag_taken
!
!
   subroutine get_node(this, i, node)
!!
!!    Get node
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observer_list), intent(in) :: this
!
      integer, intent(in) :: i
!
      class(observer_node), pointer, intent(out) :: node
!
      integer :: j
!
      node => this%head
!
      do j = 1, i - 1
         node => node%next
      enddo
!
   end subroutine get_node
!
!
end module observer_list_class
