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
module cholesky_block_list_class
!!
!! Cholesky block list class
!! Written by Sarai D. Folkestad, 2021
!!
   use kinds
   use memory_manager_class, only : mem
   use cholesky_block_node_class, only : cholesky_block_node
   use global_out, only : output
   use range_class, only: range_
!
   implicit none
!
   type :: cholesky_block_list
!
      integer :: n_nodes = 0
      integer :: n_J
!
      type(cholesky_block_node), pointer :: head => null()
      type(cholesky_block_node), pointer :: tail => null()
!
   contains
!
      procedure, public :: finalize   => finalize_cholesky_block_list
      procedure, public :: remove     => remove_cholesky_block_list
      procedure, public :: push_back  => push_back_cholesky_block_list
      procedure, public :: get_array  => get_array_cholesky_block_list
!
      procedure, private :: pop_back
      procedure, private :: pop_front
      procedure, private :: get_node
!
      procedure, public :: node_in_list
!
   end type cholesky_block_list
!
   interface cholesky_block_list
!
      procedure :: new_cholesky_block_list
!
   end interface cholesky_block_list
!
!
contains
!
!
   function new_cholesky_block_list(n_J) result(list)
!!
!!    New Cholesky block list
!!    Written by Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      type(cholesky_block_list) :: list
      integer, intent(in) :: n_J
!
      list%n_J = n_J
!
      list%n_nodes = 0
      list%head => null()
      list%tail => null()
!
   end function new_cholesky_block_list
!
!
   subroutine finalize_cholesky_block_list(list)
!!
!!    Finalize
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Delete all nodes, by deleting the last node until the
!!    list is empty. Then set head and tail to null()
!!
      implicit none
!
      class(cholesky_block_list) :: list
!
      do while (list%n_nodes .gt. 0)
!
         call list%pop_back()
!
      enddo
!
   end subroutine finalize_cholesky_block_list
!
!
   subroutine push_back_cholesky_block_list(list, range_p, range_q)
!!
!!    Push back
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Insert a node at the end of the list
!!
      implicit none
!
      class(cholesky_block_list), intent(inout) :: list
!
      type(range_), intent(in) :: range_p, range_q
!
      if (list%node_in_list(range_p, range_q)) &
         call output%error_msg('cannot add node corresponding to existing node!')
!
      if (associated(list%head)) then 
!
         allocate(list%tail%next)
         list%tail%next = cholesky_block_node(list%n_J, range_p, range_q)
         call list%tail%next%initialize()
!
         list%tail%next%previous => list%tail
!
         list%tail => list%tail%next
!
      else 
!
         allocate(list%head)
!
         list%head = cholesky_block_node(list%n_J, range_p, range_q)
         call list%head%initialize()
!
         list%tail => list%head
!
      endif
!
      list%n_nodes = list%n_nodes + 1
!
   end subroutine push_back_cholesky_block_list
!
!
   subroutine pop_back(list)
!!
!!    Pop back
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Remove at the end of the list
!!
      implicit none
!
      class(cholesky_block_list), intent(inout) :: list
!
      if (associated(list%head)) then
!
         call list%tail%destruct()
!
         if (associated(list%tail%previous)) then
!
            list%tail =>  list%tail%previous
            list%tail%next => null()
!
         else ! List is now emty
!
            list%head => null()
            list%tail => null()
!
         endif
!
      else 
!
         call output%error_msg('Attempted to remove node from empty list.')
!
      endif
!
      list%n_nodes = list%n_nodes - 1
!
   end subroutine pop_back
!
!
   subroutine pop_front(list)
!!
!!    Pop front
!!    Written by Sarai D. Folkestad, Jan 2022
!!
!!    Remove at the start of the list
!!
      implicit none
!
      class(cholesky_block_list), intent(inout) :: list
!
      if (associated(list%head)) then
!
         call list%head%destruct()
!
         if (associated(list%head%next)) then
!
            list%head =>  list%head%next
            list%head%previous => null()
!
         else ! List is now emty
!
            list%head => null()
            list%tail => null()
!
         endif
!
      else
!
         call output%error_msg('Attempted to remove node from empty list.')
!
      endif
!
      list%n_nodes = list%n_nodes - 1
!
   end subroutine pop_front
!
!
   subroutine get_array_cholesky_block_list(list, array_pointer, range_p, range_q, first, last)
!!
!!    Get array
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      class(cholesky_block_list), intent(in) :: list
!
      real(dp), dimension(:,:,:), pointer, intent(out) :: array_pointer
!
      type(range_), intent(in) :: range_p, range_q
!
      integer, intent(in), optional :: first, last
!
      integer :: first_local, last_local
!
      type(cholesky_block_node), pointer :: node_pointer
!
      call list%get_node(node_pointer, range_p, range_q)
!
      first_local = 1
      if (present(first)) first_local = first
!
      last_local = node_pointer%range_q%length
      if (present(last)) last_local = last
!
      if (last_local > node_pointer%range_q%length .or. last_local < 1) &
           call output%error_msg('in get_array in the Cholesky block list. last is out of range')
      if (first_local > node_pointer%range_q%length .or. first_local < 1) &
           call output%error_msg('in get_array in the Cholesky block list. first is out of range')
      if (first_local > last_local) &
           call output%error_msg('in get_array in the Cholesky block list. last is before first')
!
      array_pointer => node_pointer%array(:,:, first_local:last_local)
!
   end subroutine get_array_cholesky_block_list
!
!
   subroutine get_node(list, node_pointer, range_p, range_q)
!!
!!    Get node
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Set a pointer to the n'th node in the list
!!
      implicit none
!
      class(cholesky_block_list), intent(in) :: list
      type(range_), intent(in) :: range_p, range_q
!
      type(cholesky_block_node), pointer, intent(out) :: node_pointer
!
      integer :: element
!
      node_pointer => list%head
!
      do element = 1, list%n_nodes
!
         if(node_pointer%range_p == range_p .and. &
            node_pointer%range_q == range_q) return
!
         node_pointer => node_pointer%next
!
      enddo
!
      call output%error_msg('in cholesky_block_list. Tried to get non-existing node')
!
   end subroutine get_node
!
!
   subroutine remove_cholesky_block_list(list, range_p, range_q)
!!
!!    Remove
!!    Written by Sarai D. Folkestad, Jan 2022
!!
!!    Remove at position n
!!
      implicit none
!
      class(cholesky_block_list), intent(inout) :: list
!
      type(range_), intent(in) :: range_p, range_q
!
      type(cholesky_block_node), pointer :: node_pointer, node_prev_pointer, node_next_pointer
!
      if (associated(list%head)) then ! List has nodes
!
         if (list%tail%range_p == range_p .and. &
             list%tail%range_q == range_q) then
!
            call list%pop_back()
!
         else if (list%head%range_p == range_p .and. &
                  list%head%range_q == range_q) then
!
            call list%pop_front()
!
         else

            call list%get_node(node_pointer, range_p, range_q)
!
            node_next_pointer => node_pointer%next
            node_prev_pointer => node_pointer%previous
!
            call node_pointer%destruct()
            node_pointer => null()
!
            node_prev_pointer%next => node_next_pointer
            node_next_pointer%previous => node_prev_pointer
!
            list%n_nodes = list%n_nodes - 1
!
         endif
      else
!
         call output%error_msg('in cholesky_block_list. Tried to remove node from empty list')
!
      endif
!
   end subroutine remove_cholesky_block_list
!
!
   function node_in_list(list, range_p, range_q) result(found)
!!
!!    Node in list
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      class(cholesky_block_list), intent(inout) :: list
      type(range_), intent(in) :: range_p, range_q
!
      logical :: found
      integer :: n
!
      type(cholesky_block_node), pointer :: node_pointer
!
      found = .false.
!
      node_pointer => list%head
!
      do n = 1, list%n_nodes
!
         if(node_pointer%range_p == range_p .and. &
            node_pointer%range_q == range_q) then
            found = .true.
            return
         endif
         node_pointer => node_pointer%next
!
      enddo
!
   end function node_in_list
!
!
end module cholesky_block_list_class
