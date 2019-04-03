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
module cholesky_array_list_class
!
!!
!!    Cholesky array list class module
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, July 2018
!!
!
   use kinds
   use memory_manager_class
   use array_list_class
   use array_utilities
!
   implicit none
!
   type, extends(array_list) :: cholesky_array_list
!
!     Nothing here
!
   contains
!
      procedure :: reduce => reduce_cholesky_array_list
!
   end type cholesky_array_list
!
contains
!
!
   subroutine reduce_cholesky_array_list(cholesky_array, block_firsts, block_significant, n_blocks, dim_reduced)
!!
!!    Reduce Cholesky array by deleting rows
!!
!!    block_firsts: integer array of offsets to first element in each block
!!    block_significant: logical array whose elements are true if a block is significant and is to be kept
!!    n_blocks: total number of blocks
!!    dim_reduced: number of rows after reduction
!!
      implicit none
!
      class(cholesky_array_list) :: cholesky_array
!
      integer :: n_blocks
!
      integer, dimension(n_blocks + 1) :: block_firsts
      logical, dimension(n_blocks)     :: block_significant
!
      integer :: dim_reduced
!
      integer :: list_element
!
      type(array_node), pointer :: node_ptr, node_ptr_new
!
      do list_element = 1, cholesky_array%n_nodes
!
!        Get pointer to node associated with list element
!
         call cholesky_array%get_node(node_ptr, list_element)
!
!        Insert new reduced node at position 'list_element',
!        then get the pointer to this new element 
!
         call cholesky_array%insert(dim_reduced, node_ptr%n_columns, list_element)
         call cholesky_array%get_node(node_ptr_new, list_element)
!
!        Determine reduced array by cutting away insignificant blocks
!
         call reduce_array(node_ptr%array,     &
                           node_ptr_new%array, &
                           block_firsts,       &
                           block_significant,  &
                           n_blocks,           &
                           node_ptr%n_rows,    &
                           dim_reduced,        &
                           node_ptr%n_columns)
!
!        Remove old unreduced node from list 
!
         call cholesky_array%remove(list_element + 1)
!
      enddo
!
   end subroutine reduce_cholesky_array_list
!
!
end module cholesky_array_list_class
