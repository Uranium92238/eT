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
module array_3D_node_class
!!
!! Array 3D node class
!! Written by Sarai D. Folkestad, 2021
!!
!! 
   use kinds
   use memory_manager_class, only : mem
!
   implicit none
!
   type :: array_3D_node
!
      type(array_3D_node), pointer :: next => null()
      type(array_3D_node), pointer :: previous => null()
!
      real(dp), dimension(:,:,:), allocatable :: array
!
      integer :: dim_1
      integer :: dim_2
      integer :: dim_3
!
   contains
!
      procedure :: initialize => initialize_array_3D_node
      procedure :: destruct   => destruct_array_3D_node
!
   end type array_3D_node
!
!
   interface array_3D_node
      procedure :: new_array_3D_node
   end interface
!
!
contains
!
!
   function new_array_3D_node(dim_1, dim_2, dim_3) result(node)
!!
!!    New array node
!!    Written by Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      type(array_3D_node) :: node
!
      integer, intent(in) :: dim_1, dim_2, dim_3
!
      node%dim_1 = dim_1
      node%dim_2 = dim_2
      node%dim_3 = dim_3
!
      node%next => null()
      node%previous => null()
!
   end function new_array_3D_node
!
!
   subroutine initialize_array_3D_node(node)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!! 
!!    Set the number of rows and collumns of the array.
!!    Allocate the array
!!
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(array_3D_node) :: node
!
      call mem%alloc(node%array, node%dim_1, node%dim_2, node%dim_3)
      call zero_array(node%array, node%dim_1*node%dim_2*node%dim_3)
!
   end subroutine initialize_array_3D_node
!
!
   subroutine destruct_array_3D_node(node)
!!
!!    Destruct
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
!!    Deallocates the array
!!
      implicit none
!
      class(array_3D_node) :: node
!
      call mem%dealloc(node%array, node%dim_1, node%dim_2, node%dim_3)
!
   end subroutine destruct_array_3D_node
!
!
end module array_3D_node_class
