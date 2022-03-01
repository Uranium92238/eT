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
!!    A solver can request the nth start start_vector R as follows:
!!
!!       call start_vector_this%get(R, n).
!!
!!    Initialization depends on the constructor; see descendants. 
!!
!
   use parameters
   use array_utilities, only: zero_array
   use start_vector_tool_class, only: start_vector_tool
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(start_vector_tool), abstract :: es_start_vector_tool
!
      integer, dimension(:), allocatable :: indices
      integer :: n_vectors
      class(ccs), pointer :: wf
      logical :: restart
      character(len=200) :: side
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
   subroutine get_es_start_vector_tool(this, start_vector, I, energy)
!!
!!    Get
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    Modified by Alexander C. Paul, Sep 2020 - now also handles restart
!!
!!    Returns the nth start start_vector
!!    If restart is requested, it checks for the start_vector on file
!!    otherwise the start start_vector are determined from this%indices.
!!
      use array_utilities, only: get_l2_norm
!
      implicit none 
!
      class(es_start_vector_tool), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: start_vector
      integer, intent(in) :: I
!
      real(dp), intent(out) :: energy
!
      logical  :: vector_found
      real(dp) :: norm_
!
      if (.not. this%restart) then
!
         energy = zero
         call zero_array(start_vector, this%n_parameters)
         start_vector(this%indices(I)) = one
!
      else ! restarting
!
         call this%wf%check_and_get_restart_vector(start_vector, energy, &
                                                   I, trim(this%side), vector_found)
!
         if (vector_found) then
!
            norm_ = get_l2_norm(start_vector, this%n_parameters)
            call dscal(this%n_parameters, one/norm_, start_vector, 1)
!               
         else
!
            energy = zero
            call zero_array(start_vector, this%n_parameters)
            start_vector(this%indices(I)) = one
!
         end if
      end if
!
   end subroutine get_es_start_vector_tool
!
!
end module es_start_vector_tool_class
