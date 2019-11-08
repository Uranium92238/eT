!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
module array_list_class
!
!!
!!    Array list class module
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, July 2018
!!
!!    A linked list of arrays. The nodes are defined in array_node class.
!!
!!    The class contains all functionality needed to use the linked list
!!
!!       - Different methods of insertion of a new array node
!!       - Removal of an array node
!!       - Get and set routines to access information of a specific node
!!
!
   use kinds
   use memory_manager_class, only : mem
   use array_node_class, only : array_node
   use global_out, only : output
!
   implicit none
!
   type :: array_list
!
      integer :: n_nodes = 0
!
      type(array_node), pointer :: head => null()
      type(array_node), pointer :: tail => null()
!
   contains
!
      procedure :: initialize    => initialize_array_list
      procedure :: finalize      => finalize_array_list
!
      procedure :: push_back     => push_back_array_list
      procedure :: pop_back      => pop_back_array_list
      procedure :: push_front    => push_front_array_list
      procedure :: insert        => insert_array_list
      procedure :: remove        => remove_array_list
!
      procedure :: get_array     => get_array_array_list
      procedure :: set_array     => set_array_array_list
!
      procedure :: get_node      => get_node_array_list
!
      procedure :: get_n_rows_element     => get_n_rows_element_array_list
      procedure :: get_n_columns_element  => get_n_columns_element_array_list
!
      procedure :: keep_columns           => keep_columns_array_list
!
   end type array_list
!
!
   interface array_list
!
      procedure :: new_array_list
!
   end interface array_list
!
!
contains
!
!
   function new_array_list() result(list)
!!
!!    New array list
!!    Written by Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      type(array_list) :: list
!
      list%n_nodes = 0
      list%head => null()
      list%tail => null()
!
   end function new_array_list
!
!
   subroutine initialize_array_list(list)
!!
!!    Initialize 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Initialize the empty array list:
!! 
!!       Set both head and tail poiters to null()
!!       Set the number of nodes to 0
!!
      implicit none
!
      class(array_list) :: list
!
      list%n_nodes = 0
      list%head => null()
      list%tail => null()
!
   end subroutine initialize_array_list
!
!
   subroutine finalize_array_list(list)
!!
!!    Finalize
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Delete all nodes, by deleting the last node until the
!!    list is empty. Then set head and tail to null()
!!
      implicit none
!
      class(array_list) :: list
!
      do while (list%n_nodes .gt. 0)
!
         call list%pop_back()
!
      enddo
!
      list%head => null()
      list%tail => null()
!
   end subroutine finalize_array_list
!
!
   subroutine keep_columns_array_list(list, n, first, last)
!!
!!    Keep columns
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Keep only a specific range (first, last) of columns of 
!!    array node n.
!!
!!    Redefine the node with these columns only.
!!
      implicit none
!
      class(array_list), intent(inout) :: list
!
      integer, intent(in) :: n, first, last
!
      integer :: element, n_rows
!
      real(dp), dimension(:,:), allocatable :: temp_array
!
      type(array_node), pointer :: node_pointer
!
!     Point at the correct node (n)
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
!     Construct a temporary array using only the requested columns
!
      call mem%alloc(temp_array, node_pointer%n_rows, last - first + 1)
!
      call dcopy((last - first + 1)*node_pointer%n_rows, &
                  node_pointer%array(1, first), 1, temp_array, 1)
!
!     Redefine the node
!
      call node_pointer%destruct()
!
      n_rows = node_pointer%n_rows
!
      call node_pointer%initialize(n_rows, last - first + 1)
!
      call dcopy(node_pointer%n_rows*node_pointer%n_columns, temp_array, 1, node_pointer%array,1)
!      
      call mem%dealloc(temp_array, node_pointer%n_rows, last - first + 1)
!
   end subroutine keep_columns_array_list
!
!
   subroutine push_back_array_list(list, n_rows, n_columns)
!!
!!    Push back
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Insert a node at the end of the list
!!
      implicit none
!
      class(array_list), intent(inout) :: list
      integer, intent(in) :: n_rows, n_columns
!
      if (associated(list%head)) then ! List has nodes already
!
!        Allocate the next node of the last node (tail)
!
         allocate(list%tail%next)
!
!        Set the node specifics
!
         call list%tail%next%initialize(n_rows, n_columns)
!
!        The previous node for the new node is tail
!
         list%tail%next%previous => list%tail
!
!        Reset tail to the new node
!
         list%tail => list%tail%next
!
!        The next node of the last node does not exist
!
         list%tail%next => null()
!
      else ! List is empty
!
!        Allocate head
!
         allocate(list%head)
!
!        Set the node specifics
!
         call list%head%initialize(n_rows, n_columns) 
!
!        This is the first node so head = tail
!
         list%tail => list%head
!
!        Only one node, so previous and next do not exist
!
         list%head%next => null()
         list%head%previous => null()
!
      endif
!
      list%n_nodes = list%n_nodes + 1
!
   end subroutine push_back_array_list
!
!
   subroutine push_front_array_list(list, n_rows, n_columns)
!!
!!    Push front
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Insert a node at the beginning of the list
!!
      implicit none
!
      class(array_list), intent(inout) :: list
      integer, intent(in) :: n_rows, n_columns
!
      if (associated(list%head)) then ! List has nodes already
!
!        Allocate new node, it is before head
!
         allocate(list%head%previous)
!
!        Set node specifics
!
         call list%head%previous%initialize(n_rows, n_columns)
!
!        The new node is before (old) head
!
         list%head%previous%next => list%head
!
!        The new node is the new head 
!
         list%head => list%head%previous
         list%head%previous => null()
!
      else ! List is empty
!
!        Allocate head
!
         allocate(list%head)
!
!        Set node specifics
!
         call list%head%initialize(n_rows, n_columns)
!
!        This is the first node so head = tail
!
         list%tail => list%head
!
!        Only one node, so previous and next do not exist
!
         list%head%next => null()
         list%head%previous => null()
!
      endif
!
      list%n_nodes = list%n_nodes + 1
!
   end subroutine push_front_array_list
!
!
   subroutine pop_back_array_list(list)
!!
!!    Pop back
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Remove at the end of the list
!!
      implicit none
!
      class(array_list), intent(inout) :: list
!
      if (associated(list%head)) then ! List has nodes already
!
!        Deallocate the array
!
         call list%tail%destruct()
!
         if (associated(list%tail%previous)) then
!
!           Reset the tail
!
            list%tail =>  list%tail%previous
!
!           Deallocate the node
!
            deallocate(list%tail%next)
!
         else
!
            deallocate(list%tail)
!
         endif

!
      else ! List is empty
!
         call output%error_msg('Attempted to remove node from empty list.')
!
      endif
!
      list%n_nodes = list%n_nodes - 1
!
   end subroutine pop_back_array_list
!
!
   subroutine get_array_array_list(list, array_pointer, n)
!!
!!    Get array
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Set a pointer to the n'th array in the list
!!
      implicit none
!
      class(array_list), intent(inout) :: list
!
      real(dp), dimension(:,:), pointer, intent(out) :: array_pointer
!
      integer, intent(in) :: n
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         call output%error_msg('Attempted to get an element not in the list.')
!
      endif
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
      array_pointer => node_pointer%array
!
   end subroutine get_array_array_list
!
!
   subroutine get_node_array_list(list, node_pointer, n)
!!
!!    Get node
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Set a pointer to the n'th node in the list
!!
      implicit none
!
      class(array_list), intent(in) :: list
!
      integer, intent(in) :: n
!
      type(array_node), pointer, intent(out) :: node_pointer
!
      integer :: element
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         call output%error_msg('Attempted to get a node not in the list.')
!
      endif
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
   end subroutine get_node_array_list
!
!
   subroutine set_array_array_list(list, n, array, n_rows, n_columns)
!!
!!    Set array
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Sets the array (and dimensions) of n'th node in list
!!
      implicit none
!
      class(array_list), intent(inout) :: list
!
      integer, intent(in) :: n
      integer, intent(in) :: n_rows
      integer, intent(in) :: n_columns
!
      real(dp), dimension(n_rows, n_columns) :: array
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         call output%error_msg('Attempted to set an element not in the list.')
!
      endif
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
      if (n_rows .ne. node_pointer%n_rows .or. n_columns .ne. node_pointer%n_columns) then
!
         call output%error_msg('Attempted to set an element with wrong dimensions.')
!
      endif
!
      node_pointer%array = array
!
   end subroutine set_array_array_list
!
!
   integer function get_n_rows_element_array_list(list, n)
!!
!!    Get number of rows
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Get the number of rows of the n'th node in the list
!!
      implicit none
!
      class(array_list), intent(in) :: list
!
      integer, intent(in) :: n
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         call output%error_msg('Attempted to get an element not in the list.')
!
      endif
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
      get_n_rows_element_array_list = node_pointer%n_rows
!
   end function get_n_rows_element_array_list
!
!
   integer function get_n_columns_element_array_list(list, n)
!!
!!    Get number of columns
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Get the number of columns of the n'th node in the list
!!
      implicit none
!
      class(array_list), intent(in) :: list
!
      integer, intent(in) :: n
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         call output%error_msg('Attempted to get an element not in the list.')
!
      endif
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
      get_n_columns_element_array_list = node_pointer%n_columns
!
   end function get_n_columns_element_array_list
!
!
   subroutine insert_array_list(list, n_rows, n_columns, n)
!!
!!    Insert
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Insert at position n
!!
      implicit none
!
      class(array_list), intent(inout) :: list
      integer, intent(in) :: n_rows, n_columns, n
!
      type(array_node), pointer :: nth_pointer, nth_prev_pointer
!
      integer :: element
!
       if (n .gt. (list%n_nodes + 1) .or. n .lt. 1) then
!
         call output%error_msg('Attempted to get an element not in the list.')
!
      endif
!
      if (associated(list%head)) then ! List has nodes already
!
         if (n == (list%n_nodes + 1) ) then ! insert at the end
!
            call list%push_back(n_rows, n_columns)
            return
!
         elseif (n == 1) then ! insert at the beginning
!
            call list%push_front(n_rows, n_columns)
            return
! 
         else
!
!           point to the correct place in the list
!
            nth_pointer => list%head
!
            do element = 1, n - 1
!
               nth_pointer => nth_pointer%next
!
            enddo
!
!              Point to the node previous to the insertion point
!
               nth_prev_pointer => nth_pointer%previous
               nth_prev_pointer%next => null()
!
!              Allocate a node following that previous node
!
               allocate(nth_prev_pointer%next)
!
!              Set node specifics
!
               call nth_prev_pointer%next%initialize(n_rows, n_columns)
!
!              The old n'th node is after the new node
!              and the node previous to the old n'th node 
!              is previous to the new n'th node
!
               nth_pointer%previous => nth_prev_pointer%next
               nth_prev_pointer%next%previous => nth_prev_pointer
               nth_prev_pointer%next%next => nth_pointer

         endif
!
      else ! List is empty
!
         if (n == 1) then
!
!           Allocate head
!
            allocate(list%head)
!
!           Set node specifics
!
            call list%head%initialize(n_rows, n_columns)
!
!           This is the first node so head = tail
!
            list%tail => list%head
!
!           Only one node, so previous and next do not exist
!
            list%head%next => null()
            list%head%previous => null()
!
         else
!
            call output%error_msg('Attempted to insert element into empty list '// &
                                  'at position (i0)', ints=[n])
!
         endif
!
      endif
!
      list%n_nodes = list%n_nodes + 1
!
   end subroutine insert_array_list
!
!
   subroutine remove_array_list(list, n)
!!
!!    Remove
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Remove at position n
!!
      implicit none
!
      class(array_list), intent(inout) :: list
      integer, intent(in) :: n
!
      type(array_node), pointer :: nth_pointer, nth_prev_pointer, nth_next_pointer
!
      integer :: element
!
       if (n .gt. (list%n_nodes + 1) .or. n .lt. 1) then
!
         call output%error_msg('Attempted to remove an element not in the list.')
!
      endif
!
      if (associated(list%head)) then ! List has nodes already
!
         if (n == list%n_nodes) then ! Remove at the end
!
            call list%pop_back()
            return
!
         else
!
            nth_pointer => list%head
!
            do element = 1, n - 1
!
               nth_pointer => nth_pointer%next
!
            enddo
!
!           Reset pointers
!
            nth_next_pointer => nth_pointer%next 
            nth_prev_pointer => nth_pointer%previous
!
            nth_prev_pointer%next => nth_next_pointer
            nth_next_pointer%previous => nth_prev_pointer
!
!           Deallocate array
!           
            call nth_pointer%destruct()

         endif
!
      endif
!
      list%n_nodes = list%n_nodes - 1
!
   end subroutine remove_array_list
!
!
end module array_list_class
