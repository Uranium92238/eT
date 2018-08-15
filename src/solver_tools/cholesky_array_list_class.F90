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
      integer(i15), dimension(n_blocks + 1, 1) :: block_firsts
      logical, dimension(n_blocks, 1)          :: block_significant
!
      integer(i15) :: n_blocks
      integer(i15) :: dim_reduced
!
      integer(i15) :: list_element
!
      type(array_node), pointer :: node_ptr
!
      real(dp), dimension(:,:), allocatable :: temp_reduced_array
!
      integer(i15) :: n_rows, n_columns
!
      do list_element = 1, cholesky_array%n_nodes
!
!        Get pointer to node associated with list element
!
         call cholesky_array%get_node(node_ptr, list_element)
!
!        Determine reduced array by cutting away insignificant blocks
!
         call mem%alloc(temp_reduced_array, dim_reduced, node_ptr%n_columns)
!
         call reduce_array(node_ptr%array,     &
                           temp_reduced_array, &
                           block_firsts,       &
                           block_significant,  &
                           n_blocks,           &
                           node_ptr%n_rows,    &
                           dim_reduced,        &
                           node_ptr%n_columns)
!
!        Reallocate and set reduced node
!
         call mem%dealloc(node_ptr%array, node_ptr%n_rows, node_ptr%n_columns)
!
         node_ptr%n_rows = dim_reduced
!
         call mem%alloc(node_ptr%array, node_ptr%n_rows, node_ptr%n_columns)
!
         node_ptr%array = temp_reduced_array
!
         call mem%dealloc(temp_reduced_array, dim_reduced, node_ptr%n_columns)
!
      enddo
!
   end subroutine reduce_cholesky_array_list
!
!
end module cholesky_array_list_class
