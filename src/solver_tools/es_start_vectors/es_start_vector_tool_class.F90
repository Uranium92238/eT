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
module es_start_vector_tool_class
!
!!
!!    Excited state start vector tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    Sets start vectors with ones at particular excitations given by "indices".
!!    These depend on the calculation (valence/CVS/IP) and are set accordingly in 
!!    the constructors of the descendants. 
!!
!!    A solver can request the nth start vector R as follows:
!!
!!       call start_vector_tool%get_vector(R, n).
!!
!!    Initialization depends on the constructor; see descendants. 
!!
!
   use parameters
   use array_utilities, only: zero_array
!
   implicit none
!
   type, abstract :: es_start_vector_tool
!
      integer, dimension(:), allocatable :: indices
!
      integer :: n_vectors 
      integer :: vector_length
!
   contains
!
      procedure :: get_vector => get_vector_es_start_vector_tool
!
   end type es_start_vector_tool
!
!
contains 
!
!
   subroutine get_vector_es_start_vector_tool(tool, R, n)
!!
!!    Get vector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    Sets the nth vector R equal to the starting guess according to the 
!!    nth index in indices (see constructors for the initialization of tool%indices)
!!
      implicit none 
!
      class(es_start_vector_tool), intent(in) :: tool
!
      real(dp), dimension(tool%vector_length), intent(inout) :: R 
!
      integer, intent(in) :: n 
!
      call zero_array(R, tool%vector_length)
      R(tool%indices(n)) = one
!
   end subroutine get_vector_es_start_vector_tool
!
!
end module es_start_vector_tool_class
