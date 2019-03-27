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
!
   use kinds
   use memory_manager_class
   use array_node_class
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
      procedure :: initialize => initialize_array_list
      procedure :: finalize   => finalize_array_list
!
      procedure :: push_back  => push_back_array_list
      procedure :: pop_back   => pop_back_array_list
      procedure :: push_front  => push_front_array_list
      procedure :: insert     => insert_array_list
      procedure :: remove     => remove_array_list
!
      procedure :: get_element => get_element_array_list
      procedure :: set_element => set_element_array_list
!
      procedure :: get_node => get_node_array_list
!
      procedure :: get_n_rows_element => get_n_rows_element_array_list
      procedure :: get_n_columns_element => get_n_columns_element_array_list
!
      procedure :: keep_columns => keep_columns_array_list
!
   end type array_list
!
contains
!
!
   subroutine initialize_array_list(list)
!!
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
!!    Keep only a specific range of columns of n
!!
      implicit none
!
      class(array_list) :: list
!
      integer :: n, first, last
!
      integer :: element
!
      real(dp), dimension(:,:), allocatable :: temp_array
!
      type(array_node), pointer :: node_pointer
!
      node_pointer => list%head
!
      do element = 1, n - 1
!
         node_pointer => node_pointer%next
!
      enddo
!
      call mem%alloc(temp_array, node_pointer%n_rows, last - first + 1)
!
      temp_array(:, :) = node_pointer%array(:, first:last)
!
      call mem%dealloc(node_pointer%array, node_pointer%n_rows, node_pointer%n_columns)
      call mem%alloc(node_pointer%array, node_pointer%n_rows, last - first + 1)
      node_pointer%array = temp_array
!
      node_pointer%n_columns = last - first + 1
      call mem%dealloc(temp_array, node_pointer%n_rows, last - first + 1)
!
   end subroutine keep_columns_array_list
!
!
   subroutine push_back_array_list(list, n_rows, n_columns)
!!
!!    Insert at the end of the list
!!
      implicit none
!
      class(array_list) :: list
      integer, intent(in) :: n_rows, n_columns
!
      if (associated(list%head)) then ! List has nodes already
!
         allocate(list%tail%next)
         call mem%alloc(list%tail%next%array, n_rows, n_columns)
         list%tail%next%n_rows = n_rows
         list%tail%next%n_columns = n_columns
!
         list%tail%next%previous => list%tail
         list%tail => list%tail%next

         list%tail%next => null()
!
      else ! List is empty
!
         allocate(list%head)
         call mem%alloc(list%head%array, n_rows, n_columns)
         list%head%n_rows = n_rows
         list%head%n_columns = n_columns
!
         list%tail => list%head

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
!!    Insert at the beginning of the list
!!
      implicit none
!
      class(array_list) :: list
      integer, intent(in) :: n_rows, n_columns
!
      if (associated(list%head)) then ! List has nodes already
!
         allocate(list%head%previous)
         call mem%alloc(list%head%previous%array, n_rows, n_columns)
         list%head%previous%n_rows = n_rows
         list%head%previous%n_columns = n_columns
!
         list%head%previous%next => list%head
         list%head => list%head%previous

         list%head%previous => null()
!
      else ! List is empty
!
         allocate(list%head)
         call mem%alloc(list%head%array, n_rows, n_columns)
         list%head%n_rows = n_rows
         list%head%n_columns = n_columns
!
         list%tail => list%head

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
!!    Remove at the end of the list
!!
      implicit none
!
      class(array_list) :: list
!
      if (associated(list%head)) then ! List has nodes already
!
         call mem%dealloc(list%tail%array, list%tail%n_rows, list%tail%n_columns)
!
         if (associated(list%tail%previous)) then
            list%tail =>  list%tail%previous
            deallocate(list%tail%next)
         else
            deallocate(list%tail)
         endif

!
      else ! List is empty
!
         write(output%unit, '(a)') 'Error: attempted to remove node from empty list.'
         stop
!
      endif
!
      list%n_nodes = list%n_nodes - 1
!
   end subroutine pop_back_array_list
!
!
   subroutine get_element_array_list(list, array_pointer, n)
!!
!!    Set a pointer to the nth array in the list
!!
      implicit none
!
      class(array_list) :: list
!
      real(dp), dimension(:,:), pointer :: array_pointer
!
      integer, intent(in) :: n
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to get an element not in the list.'
         stop
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
   end subroutine get_element_array_list
!
!
   subroutine get_node_array_list(list, node_pointer, n)
!!
!!    Set a pointer to the nth node in the list
!!
      implicit none
!
      class(array_list) :: list
!
      integer :: n
      integer :: element
!
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to get a node not in the list.'
         stop
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
   subroutine set_element_array_list(list, n, array, n_rows, n_columns)
!!
!!    Sets value of nth element in list
!!
      implicit none
!
      class(array_list) :: list
!
      integer :: n
!
      integer :: n_rows
      integer :: n_columns
!
      real(dp), dimension(n_rows, n_columns) :: array
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to set an element not in the list.'
         flush(output%unit)
         stop
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
         write(output%unit, '(t3,a)') 'Error: attempted to set an element with wrong dimensions.'
         flush(output%unit)
         stop
!
      endif
!
      node_pointer%array = array
!
   end subroutine set_element_array_list
!
!
   integer function get_n_rows_element_array_list(list, n)
!!
!!    Set a pointer to the nth array in the list
!!
      implicit none
!
      class(array_list) :: list
!
      integer, intent(in) :: n
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to get an element not in the list.'
         stop
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
!!    Set a pointer to the nth array in the list
!!
      implicit none
!
      class(array_list) :: list
!
      integer, intent(in) :: n
!
      integer :: element
      type(array_node), pointer :: node_pointer
!
      if (n .gt. list%n_nodes .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to get an element not in the list.'
         stop
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
!!    Insert at position n
!!
      implicit none
!
      class(array_list) :: list
      integer, intent(in) :: n_rows, n_columns, n
!
      type(array_node), pointer :: nth_pointer, nth_prev_pointer
!
      integer :: element
!
       if (n .gt. (list%n_nodes + 1) .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to get an element not in the list.'
         stop
!
      endif
!
      if (associated(list%head)) then ! List has nodes already
!
         if (n == (list%n_nodes + 1) ) then
!
            call list%push_back(n_rows, n_columns)
            return
!
         elseif (n == 1) then
!
            call list%push_front(n_rows, n_columns)
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
               nth_prev_pointer => nth_pointer%previous
               nth_prev_pointer%next => null()
!
               allocate(nth_prev_pointer%next)
               call mem%alloc(nth_prev_pointer%next%array, n_rows, n_columns)
               nth_prev_pointer%next%n_rows = n_rows
               nth_prev_pointer%next%n_columns = n_columns
               
               nth_pointer%previous => nth_prev_pointer%next
               nth_prev_pointer%next%previous => nth_prev_pointer
               nth_prev_pointer%next%next => nth_pointer

         endif
!
      else ! List is empty
!
         if (n == 1) then
!
            allocate(list%head)
            call mem%alloc(list%head%array, n_rows, n_columns)
!
            list%tail => list%head
            list%head%n_rows = n_rows
            list%head%n_columns = n_columns

            list%head%next => null()
            list%head%previous => null()
!
         else
!
            write(output%unit, '(a47, i3, a15)') 'Error: attempted to insert element at position ', n, ' of empty list!'
            stop
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
!!    Remove at position n
!!
      implicit none
!
      class(array_list) :: list
      integer, intent(in) :: n
!
      type(array_node), pointer :: nth_pointer, nth_prev_pointer, nth_next_pointer
!
      integer :: element
!
       if (n .gt. (list%n_nodes + 1) .or. n .lt. 1) then
!
         write(output%unit, '(t3,a)') 'Error: attempted to remove an element not in the list.'
         stop
!
      endif
!
      if (associated(list%head)) then ! List has nodes already
!
         if (n == list%n_nodes) then
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
            nth_next_pointer => nth_pointer%next 
            nth_prev_pointer => nth_pointer%previous
!
            nth_prev_pointer%next => nth_next_pointer
            nth_next_pointer%previous => nth_prev_pointer
               
            call mem%dealloc(nth_pointer%array, nth_pointer%n_rows, nth_pointer%n_columns)

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
