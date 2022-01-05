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
module cc_multipliers_start_vector_tool_class
!
!!
!!    Coupled Cluster multipliers start vector tool class module
!!    Written by Regina Matveeva, Sept 2021
!!
!!
!
   use kinds
   use start_vector_tool_class,             only: start_vector_tool
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(start_vector_tool) :: cc_multipliers_start_vector_tool
!
   logical, private :: restart
   class(ccs), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_cc_multipliers_start_vector_tool
!
   end type cc_multipliers_start_vector_tool
!
!
   interface  cc_multipliers_start_vector_tool
!
      procedure :: new_cc_multipliers_start_vector_tool
!
   end interface  cc_multipliers_start_vector_tool
!
!
contains
!
!
   function new_cc_multipliers_start_vector_tool(wf, restart) result(this)
!!
!!    New Coupled Cluster multipliers start vector tool
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      logical, intent(in) :: restart
      type(cc_multipliers_start_vector_tool) :: this
!
      this%restart = restart
      this%wf => wf
      this%n_parameters = wf%n_gs_amplitudes
!
   end function new_cc_multipliers_start_vector_tool
!
!
   subroutine get_cc_multipliers_start_vector_tool(this, start_vector, I)
!!
!!    Get
!!    Written by Regina Matveeva, Sept 2021
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(cc_multipliers_start_vector_tool), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: start_vector
      integer, intent(in) :: I
!
      call do_nothing(I)
      call this%wf%get_initial_cc_multipliers(start_vector, this%restart)
!
   end subroutine
!
!
end module cc_multipliers_start_vector_tool_class
