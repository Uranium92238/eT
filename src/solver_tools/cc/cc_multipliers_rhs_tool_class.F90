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
module cc_multipliers_rhs_tool_class
!
!!
!!    Coupled Cluster multipliers right-hand side (of a linear equation) class module
!!    Written by Regina Matveeva, Sept 2021
!!
!
   use parameters
   use rhs_linear_equation_tool_class,     only : rhs_linear_equation_tool
   use ccs_class,                          only : ccs
!
   implicit none
!
   type, extends(rhs_linear_equation_tool) :: cc_multipliers_rhs_tool
!
      class(ccs), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_cc_multipliers_rhs_tool
!
   end type cc_multipliers_rhs_tool
!
!
   interface cc_multipliers_rhs_tool
!
      procedure :: new_cc_multipliers_rhs_tool
!
   end interface cc_multipliers_rhs_tool
!
!
contains
!
   function new_cc_multipliers_rhs_tool(wf) result(this)
!!
!!    New
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      type(cc_multipliers_rhs_tool)  :: this
!
      this%wf => wf
      this%n_parameters = wf%n_gs_amplitudes
      this%n_rhs        = 1
!
   end function new_cc_multipliers_rhs_tool
!
!
   subroutine get_cc_multipliers_rhs_tool(this, rhs)
!!
!!    Get
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(cc_multipliers_rhs_tool),  intent(inout)    :: this
      real(dp), dimension(this%n_parameters, this%n_rhs), intent(inout) :: rhs
!
      call this%wf%construct_eta(rhs)
      call dscal(this%wf%n_gs_amplitudes, -one, rhs, 1)
!
   end subroutine get_cc_multipliers_rhs_tool
!
!
end module cc_multipliers_rhs_tool_class
