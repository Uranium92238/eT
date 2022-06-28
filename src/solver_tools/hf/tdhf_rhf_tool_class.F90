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
module tdhf_rhs_tool_class
!
!!
!!    TDHF right-hand side tool class module
!!    Written by Regina Matveeva, Sept 2021
!!
!
   use parameters
   use rhs_linear_equation_tool_class,     only : rhs_linear_equation_tool
   use ccs_class,                          only : ccs
!
   implicit none
!
   type, extends(rhs_linear_equation_tool) :: tdhf_rhs_tool
!
      class(hf), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_tdhf_rhs_tool
!
   end type tdhf_rhs_tool
!
!
   interface tdhf_rhs_tool
!
      procedure :: new_tdhf_rhs_tool
!
   end interface tdhf_rhs_tool
!
!
contains
!
   function new_tdhf_rhs_tool(wf) result(this)
!!
!!    New
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      type(tdhf_rhs_tool)  :: this
!
      this%wf => wf
      this%n_parameters = 2*wf%n_o*wf%n_v
      this%n_rhs        = 1
!
   end function new_tdhf_rhs_tool
!
!
   subroutine get_tdhf_rhs_tool(this, rhs)
!!
!!    Get
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(tdhf_rhs_tool),  intent(inout)    :: this
      real(dp), dimension(this%n_parameters, this%n_rhs), intent(inout) :: rhs
!
      call this%wf%get_dipole_gradient(rhs)
!
   end subroutine get_tdhf_rhs_tool
!
!
end module tdhf_rhs_tool_class
