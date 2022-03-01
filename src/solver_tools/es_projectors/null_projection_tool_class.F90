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
module null_projection_tool_class
!
!!
!!    Null projection tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!
   use parameters
   use ccs_class, only: ccs
   use abstract_projection_tool_class, only: abstract_projection_tool
!
   implicit none
!
   type, extends(abstract_projection_tool) :: null_projection_tool
!
   contains
!
      procedure :: get_projector => get_projector_null_projection_tool
!
   end type null_projection_tool
!
!
   interface null_projection_tool
!
      procedure :: new_null_projection_tool
!
   end interface null_projection_tool
!
!
contains 
!
!
   function new_null_projection_tool(n_parameters) result(tool)
!!
!!    New
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none 
!
      type(null_projection_tool) :: tool
      integer, intent(in) :: n_parameters
!
      tool%active = .false. ! No projection needed
      tool%n_parameters = n_parameters
!
   end function new_null_projection_tool
!
!
!
   subroutine get_projector_null_projection_tool(tool, projector)
!!
!!    Get projector
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(null_projection_tool), intent(in) :: tool
      real(dp), dimension(tool%n_parameters), intent(out) :: projector
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, tool%n_parameters
!
         projector(i) = one
!
      enddo
!$omp end parallel do
!
   end subroutine get_projector_null_projection_tool
!
!
end module null_projection_tool_class
