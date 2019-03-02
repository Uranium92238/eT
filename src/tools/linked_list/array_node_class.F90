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
      integer :: n_rows
      integer :: n_columns
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
      integer, intent(in) :: n_rows
      integer, intent(in) :: n_columns
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
