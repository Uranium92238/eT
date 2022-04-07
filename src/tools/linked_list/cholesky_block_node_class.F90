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
module cholesky_block_node_class
!!
!! Cholesky block node class
!! Written by Sarai D. Folkestad, Jan 2022
!!
   use parameters
   use memory_manager_class, only: mem
   use range_class, only: range_

!
   implicit none
!
   type :: cholesky_block_node
!
!
      type(cholesky_block_node), pointer :: next => null()
      type(cholesky_block_node), pointer :: previous => null()
!
      real(dp), dimension(:,:,:), allocatable :: array
!
      type(range_) :: range_p, range_q
!
      integer :: n_J
!
   contains
!
      procedure :: initialize => initialize_cholesky_block_node
      procedure :: destruct   => destruct_cholesky_block_node
!
   end type cholesky_block_node
!
   interface cholesky_block_node
!
      procedure :: new_cholesky_block_node
!
   end interface
!
contains
!
!
   function new_cholesky_block_node(n_J, range_p, range_q) result(node)
!!
!!    New Cholesky block node
!!    Written by Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      type(cholesky_block_node) :: node
!
      type(range_), intent(in) :: range_p, range_q
      integer, intent(in) :: n_J
!
      node%n_J = n_J
      node%range_p = range_p
      node%range_q = range_q
!
      node%next => null()
      node%previous => null()
!
   end function new_cholesky_block_node
!
!
   subroutine initialize_cholesky_block_node(node)
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
      class(cholesky_block_node) :: node
!
      call mem%alloc(node%array, node%n_J, node%range_p%length, node%range_q%length)
      call zero_array(node%array, node%n_J*node%range_p%length*node%range_q%length)
!
   end subroutine initialize_cholesky_block_node
!
!
   subroutine destruct_cholesky_block_node(node)
!!
!!    Destruct
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
!!    Deallocates the array
!!
      implicit none
!
      class(cholesky_block_node) :: node
!
      call mem%dealloc(node%array, node%n_J, node%range_p%length, node%range_q%length)
!
   end subroutine destruct_cholesky_block_node
!
!
end module cholesky_block_node_class
