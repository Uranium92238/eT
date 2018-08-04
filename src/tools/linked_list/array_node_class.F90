module array_node_class
!!
!!    Array node class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
   use kinds
   use memory_manager_class
!
   implicit none
!
   type :: array_node
!
      type(array_node), pointer :: next => null()
      type(array_node), pointer :: previous => null()
!
      real(dp), dimension(:,:), allocatable :: array
!
      integer(i15) :: n_rows
      integer(i15) :: n_columns
!
   contains
!
      procedure :: initialize => initialize_array_node
      procedure :: finalize   => finalize_array_node
!
   end type array_node
!
!
contains
!
!
   subroutine initialize_array_node(node, n_rows, n_columns)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
      implicit none
!
      class(array_node) :: node
!
      integer(i15), intent(in) :: n_rows
      integer(i15), intent(in) :: n_columns
!
      node%n_rows    = n_rows
      node%n_columns = n_columns
!
      call mem%alloc(node%array, node%n_rows, node%n_columns)
!
   end subroutine initialize_array_node
!
!
   subroutine finalize_array_node(node)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
      implicit none
!
      class(array_node) :: node
!
      call mem%dealloc(node%array, node%n_rows, node%n_columns)
!
   end subroutine finalize_array_node
!
!
end module array_node_class
