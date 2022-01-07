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
module array_3D_list_class
!!
!! Array 3D list class
!! Written by Sarai D. Folkestad, 2021
!!
   use kinds
   use memory_manager_class, only : mem
   use array_3D_node_class, only : array_3D_node
   use global_out, only : output
!
   implicit none
!
   type :: array_3D_list
!
      integer :: n_nodes = 0
!
      type(array_3D_node), pointer :: head => null()
      type(array_3D_node), pointer :: tail => null()
!
   contains
!
      procedure, public :: finalize => finalize_array_3D_list
      procedure, public :: push_back => push_back_array_3D_list
      procedure, public :: pop_back => pop_back_array_3D_list
      procedure, public :: get_array_n => get_array_n_array_3D_list
      procedure, public :: set_array_n => set_subarray_n_array_3D_list
!
      procedure, private :: get_node
!
   end type array_3D_list
!
!
   interface array_3D_list
!
      procedure :: new_array_3D_list
!
   end interface array_3D_list
!
!
contains
!
!
   function new_array_3D_list() result(list)
!!
!!    New array list
!!    Written by Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      type(array_3D_list) :: list
!
      list%n_nodes = 0
      list%head => null()
      list%tail => null()
!
   end function new_array_3D_list
!
!
   subroutine finalize_array_3D_list(list)
!!
!!    Finalize
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Delete all nodes, by deleting the last node until the
!!    list is empty. Then set head and tail to null()
!!
      implicit none
!
      class(array_3D_list) :: list
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
   end subroutine finalize_array_3D_list
!
!
   subroutine push_back_array_3D_list(list, dim_1, dim_2, dim_3)
!!
!!    Push back
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Insert a node at the end of the list
!!
      implicit none
!
      class(array_3D_list), intent(inout) :: list
!
      integer, intent(in) :: dim_1, dim_2, dim_3
!
      if (associated(list%head)) then 
!
         allocate(list%tail%next)
!
         list%tail%next = array_3D_node(dim_1, dim_2, dim_3)
         call list%tail%next%initialize()
!
         list%tail%next%previous => list%tail
!
         list%tail => list%tail%next
!
         list%tail%next => null()
!
      else 
!
         allocate(list%head)
!
         list%head = array_3D_node(dim_1, dim_2, dim_3)
         call list%head%initialize()
!
         list%tail => list%head
!
         list%head%next => null()
         list%head%previous => null()
!
      endif
!
      list%n_nodes = list%n_nodes + 1
!
   end subroutine push_back_array_3D_list
!
!
   subroutine pop_back_array_3D_list(list)
!!
!!    Pop back
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Remove at the end of the list
!!
      implicit none
!
      class(array_3D_list), intent(inout) :: list
!
      if (associated(list%head)) then
!
         call list%tail%destruct()
!
         if (associated(list%tail%previous)) then
!
            list%tail =>  list%tail%previous
            deallocate(list%tail%next)
!
         else
!
            deallocate(list%tail)
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
   end subroutine pop_back_array_3D_list
!
!
   subroutine get_array_n_array_3D_list(list, array_pointer, n)
!!
!!    Get array
!!
      implicit none
!
      class(array_3D_list), intent(in) :: list
!
      real(dp), dimension(:,:,:), pointer, intent(out) :: array_pointer
!
      integer, intent(in) :: n
!
      type(array_3D_node), pointer :: node_pointer
!
      call list%get_node(node_pointer, n)
!
      array_pointer => node_pointer%array
!
   end subroutine get_array_n_array_3D_list
!
!
   subroutine get_node(list, node_pointer, n)
!!
!!    Get node
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Set a pointer to the n'th node in the list
!!
      implicit none
!
      class(array_3D_list), intent(in) :: list
!
      integer, intent(in) :: n
!
      type(array_3D_node), pointer, intent(out) :: node_pointer
!
      integer :: element
!
      if (n .gt. list%n_nodes .or. n .lt. 1) &
         call output%error_msg('Attempted to get a node not in the list.')
!
      node_pointer => list%head
!
      do element = 1, n - 1
         node_pointer => node_pointer%next
      enddo
!
   end subroutine get_node
!
!
   subroutine set_subarray_n_array_3D_list(list, n, array, first_1, last_1, &
                                             first_2, last_2, first_3, last_3)
!!
!!    Set array
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jul 2018
!!
!!    Sets the array (and dimensions) of n'th node in list
!!
      implicit none
!
      class(array_3D_list), intent(inout) :: list
!
      integer, intent(in) :: n
      integer, intent(in) :: first_1, last_1, first_2, last_2, first_3, last_3
!
      real(dp), dimension(first_1:last_1, first_2:last_2, first_3:last_3) :: array
!
      type(array_3D_node), pointer :: node_pointer
!
      integer :: i, j, k
!
      call list%get_node(node_pointer, n)
!
!$omp parallel do private(i, j, k)
      do i = first_3, last_3
         do j = first_2, last_2
            do k = first_1, last_1
!
               node_pointer%array(k, j, i) = array(k - first_1 + 1, &
                                                   j - first_2 + 1, &
                                                   i - first_3 + 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine set_subarray_n_array_3D_list
!
!
end module array_3D_list_class
