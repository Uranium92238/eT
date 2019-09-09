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
module es_ip_start_vector_tool_class
!
!!
!!    Excited state IP start vector tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!
   use ccs_class, only: ccs 
   use es_start_vector_tool_class
!
   use memory_manager_class, only: mem 
   use global_out, only: output 
!
   implicit none
!
   type, extends(es_start_vector_tool) :: es_ip_start_vector_tool
!
   contains
!
      final :: destructor
!
   end type es_ip_start_vector_tool
!
!
   interface es_ip_start_vector_tool
!
      procedure :: new_es_ip_start_vector_tool
!
   end interface es_ip_start_vector_tool
!
!
contains 
!
!
   function new_es_ip_start_vector_tool(wf) result(tool)
!!
!!    New ES IP start vector tool 
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      type(es_ip_start_vector_tool) :: tool 
!
      if (.not. wf%bath_orbital) call output%error_msg('Calculation of IPs requires bath orbitals.')
!
      tool%n_vectors = wf%n_excited_states
      tool%vector_length = wf%n_es_amplitudes
!
      call mem%alloc(tool%indices, tool%n_vectors)
      call wf%set_ip_start_indices(tool%indices)
!
   end function new_es_ip_start_vector_tool
!
!
   subroutine destructor(tool)
!!
!!    Destructor 
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none 
!
      type(es_ip_start_vector_tool) :: tool 
!
      call mem%dealloc(tool%indices, tool%n_vectors)
!
   end subroutine destructor
!
!
end module es_ip_start_vector_tool_class
