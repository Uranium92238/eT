!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
!!       call start_vector_tool%get(R, n).
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
      procedure :: get => get_es_start_vector_tool
!
   end type es_start_vector_tool
!
!
contains 
!
!
   subroutine get_es_start_vector_tool(tool, wf, n, vector, energy, side, restart)
!!
!!    Get
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    Modified by Alexander C. Paul, Sep 2020 - now also handles restart
!!
!!    Returns the nth start vector
!!    If restart is requested, it checks for the vector on file
!!    otherwise the start vector are determined from tool%indices.
!!
      use array_utilities, only: get_l2_norm
      use ccs_class, only: ccs
!
      implicit none 
!
      class(es_start_vector_tool), intent(in) :: tool
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n
!
      real(dp), dimension(tool%vector_length), intent(inout) :: vector
      real(dp), intent(inout) :: energy
!
      logical, intent(in) :: restart
!
      character(len=*), intent(in) :: side
!
      logical  :: vector_found
      real(dp) :: norm_
!
      if (.not. restart) then
!
         energy = zero
         call zero_array(vector, tool%vector_length)
         vector(tool%indices(n)) = one
!
      else ! restarting
!
         call wf%check_and_get_restart_vector(vector, energy, n, side, vector_found)
!
         if (vector_found) then
!
            norm_ = get_l2_norm(vector, tool%vector_length)
            call dscal(tool%vector_length, one/norm_, vector, 1)
!               
         else
!
            energy = zero
            call zero_array(vector, tool%vector_length)
            vector(tool%indices(n)) = one
!
         end if
      end if
!
   end subroutine get_es_start_vector_tool
!
!
end module es_start_vector_tool_class
