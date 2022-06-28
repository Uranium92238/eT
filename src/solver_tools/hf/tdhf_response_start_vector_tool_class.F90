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
module tdhf_response_start_vector_tool_class
!
!!
!!    TDHF response start vector tool class module
!!    Written by Sarai D. Folkestad, Mar 2022
!!
!
   use parameters
   use linear_equation_start_vector_tool_class, only: linear_equation_start_vector_tool
   use hf_class,                only: hf
!
   implicit none
!
   type, extends(linear_equation_start_vector_tool) :: tdhf_response_start_vector_tool
!
      class(hf), pointer, private :: wf
      integer, private :: component
!
   contains
!
      procedure, public :: get => get_tdhf_response_start_vector_tool
!
   end type tdhf_response_start_vector_tool
!
   interface  tdhf_response_start_vector_tool
!
      procedure :: new_tdhf_response_start_vector_tool
!
   end interface  tdhf_response_start_vector_tool
!
contains
!
!
   function new_tdhf_response_start_vector_tool(wf, component) result(this)
!!
!!    New TDHF response start vector tool
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(in), target :: wf
      integer, intent(in) :: component
      type(tdhf_response_start_vector_tool) :: this
!
      this%wf => wf
      this%n_parameters = 2*wf%n_o*wf%n_v
      this%component = component
!
   end function new_tdhf_response_start_vector_tool
!
!
   subroutine get_tdhf_response_start_vector_tool(this, start_vector, I)
!!
!!    Get
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Gets the I'th start vector for some purpose
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(tdhf_response_start_vector_tool), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: start_vector
      integer,    intent(in) :: I
!
      call do_nothing(I)
!
      call this%wf%get_dipole_gradient(start_vector, this%component)
      call dscal(this%n_parameters, -one, start_vector, 1)
!
   end subroutine
!
!
end module tdhf_response_start_vector_tool_class
