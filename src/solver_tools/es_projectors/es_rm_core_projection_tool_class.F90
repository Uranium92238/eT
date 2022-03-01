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
module es_rm_core_projection_tool_class
!
!!
!!    Excited state remove core projection tool class module
!!    Written by Sarai D. Folkestad, 2021
!!
!
   use abstract_projection_tool_class
   use ccs_class, only: ccs 
   use eigen_davidson_tool_class, only: eigen_davidson_tool
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(abstract_projection_tool) :: es_rm_core_projection_tool
!
      class(ccs), pointer, private :: wf
!
   contains
!
      procedure :: get_projector => get_projector_es_rm_core_projection_tool
!
   end type es_rm_core_projection_tool
!
!
   interface es_rm_core_projection_tool
!
      procedure :: new_es_rm_core_projection_tool
!
   end interface es_rm_core_projection_tool
!
!
contains 
!
!
   function new_es_rm_core_projection_tool(wf) result(tool)
!!
!!    New es remove core projection tool 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      type(es_rm_core_projection_tool) :: tool 
!
      class(ccs), target :: wf
!
      tool%n_parameters = wf%n_es_amplitudes
!
      tool%active = .true.
      tool%wf => wf
!
   end function new_es_rm_core_projection_tool
!
!
   subroutine get_projector_es_rm_core_projection_tool(tool, projector)
!!
!!    Get projector
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(es_rm_core_projection_tool), intent(in) :: tool
      real(dp), dimension(tool%n_parameters), intent(out) :: projector
!
      call tool%wf%get_rm_core_projector(projector, tool%wf%n_core_MOs, tool%wf%core_MOs)
!
   end subroutine get_projector_es_rm_core_projection_tool
!
!
end module es_rm_core_projection_tool_class
