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
module cc_multipliers_linear_storage_tool_class
!
!!
!!    Coupled Cluster multipliers linear storage tool class module
!!    Written by Regina Matveeva, Sept 2021
!!
!
   use global_out, only: output
   use parameters
   use linear_storage_tool_class, only: linear_storage_tool
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(linear_storage_tool) :: cc_multipliers_linear_storage_tool
!
      class(ccs), pointer, private :: wf
!
   contains
!
      procedure, public :: store => store_cc_multipliers_linear_storage_tool
!
   end type  cc_multipliers_linear_storage_tool
!
!
   interface  cc_multipliers_linear_storage_tool
!
      procedure :: new_cc_multipliers_linear_storage_tool
!
   end interface  cc_multipliers_linear_storage_tool
!
!
contains
!
   function new_cc_multipliers_linear_storage_tool(wf) result(this)
!!
!!    New Coupled Cluster multipliers storage tool
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ccs), intent(in), target  :: wf
      type(cc_multipliers_linear_storage_tool)   :: this
!
      this%wf => wf
      this%n_parameters = wf%n_gs_amplitudes
!
   end function new_cc_multipliers_linear_storage_tool
!
!
   subroutine store_cc_multipliers_linear_storage_tool(this, solution_vector, n)
!!
!!    Store
!!    Written by Regina Matveeva, Sept 2021
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(cc_multipliers_linear_storage_tool),   intent(inout)     :: this
      integer,                                  intent(in) :: n
      real(dp), dimension(this%n_parameters),   intent(in) :: solution_vector
!
      call do_nothing(n)
!
      call this%wf%set_multipliers(solution_vector)
      call this%wf%save_multipliers()
!
   end subroutine store_cc_multipliers_linear_storage_tool
!
!
end module cc_multipliers_linear_storage_tool_class
