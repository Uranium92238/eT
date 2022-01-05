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
module tdhf_start_vector_tool_class
!
!!
!!    TDHF start vector tool class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!
   use kinds
   use start_vector_tool_class, only: start_vector_tool
   use hf_class,                only: hf
!
   implicit none
!
   type, extends(start_vector_tool) :: tdhf_start_vector_tool
!
      logical, private :: restart
      class(hf), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_tdhf_start_vector_tool
!
   end type tdhf_start_vector_tool
!
   interface  tdhf_start_vector_tool
!
      procedure :: new_tdhf_start_vector_tool
!
   end interface  tdhf_start_vector_tool
!
contains
!
!
   function new_tdhf_start_vector_tool(wf, restart) result(this)
!!
!!    New TDHF start vector tool
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(in), target :: wf
      logical, intent(in) :: restart
      type(tdhf_start_vector_tool) :: this
!
      this%restart = restart
      this%wf => wf
      this%n_parameters = wf%n_o*wf%n_v
!
   end function new_tdhf_start_vector_tool
!
!
   subroutine get_tdhf_start_vector_tool(this, start_vector, I)
!!
!!    Get
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Gets the I'th start vector for some purpose
!!
      implicit none
!
      class(tdhf_start_vector_tool), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: start_vector
      integer, intent(in) :: I
!
      call this%wf%get_tdhf_start_vector(I, start_vector, this%restart)
!
   end subroutine
!
!
end module tdhf_start_vector_tool_class
