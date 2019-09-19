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
module es_ip_projection_tool_class
!
!!
!!    Excited state IP projection tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!
   use es_projection_tool_class
   use ccs_class, only: ccs 
   use eigen_davidson_tool_class, only: eigen_davidson_tool
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(es_projection_tool) :: es_ip_projection_tool
!
   contains
!
      final :: destructor 
!
   end type es_ip_projection_tool
!
!
   interface es_ip_projection_tool
!
      procedure :: new_es_ip_projection_tool
!
   end interface es_ip_projection_tool
!
!
contains 
!
!
   function new_es_ip_projection_tool(wf) result(tool)
!!
!!    New ES IP projection tool 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
      implicit none 
!
      type(es_ip_projection_tool) :: tool 
!
      class(ccs) :: wf 
!
      tool%vector_length = wf%n_es_amplitudes
      call mem%alloc(tool%projector, tool%vector_length)
!
      tool%active = .true.
      call wf%get_ip_projector(tool%projector)
!
   end function new_es_ip_projection_tool
!
!
   subroutine destructor(tool)
!!
!!    Destructor 
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none 
!
      type(es_ip_projection_tool) :: tool 
!
      if (allocated(tool%projector)) call mem%dealloc(tool%projector, tool%vector_length)
!
   end subroutine destructor
!
!
end module es_ip_projection_tool_class
