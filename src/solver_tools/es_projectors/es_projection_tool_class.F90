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
module es_projection_tool_class
!
!!
!!    Excited state projector tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    Vectors can be projected as follows:
!!
!!       call projection_tool%project(R)     (*)
!!
!!    Initialization depends on the calculation type; see descendants.
!!
!!    Note that the call (*) is usually not done for Davidson
!!    algorithms, where the actual projection in done inside the Davidson tool.
!!    Initialization is still handled by the projection tool.
!!
!
   use kinds 
   use global_out, only: output
!
   implicit none
!
   type, abstract :: es_projection_tool
!
      logical :: active ! Do projection if true, don't if false 
!
      real(dp), dimension(:), allocatable :: projector 
!
      integer :: vector_length
!
   contains
!
      procedure :: do_ => do_es_projection_tool 
!
   end type es_projection_tool
!
!
contains 
!
!
   subroutine do_es_projection_tool(tool, R)
!!
!!    Do
!!    Written by Sarai D. Folkestad, 2018-2019
!!
!!    Does the projection operation:
!!
!!       R(i) = R(i)*projector(i).
!!
!!    The projector is set by the constructor.
!!
      implicit none
!
      class(es_projection_tool), intent(in) :: tool
!
      real(dp), dimension(tool%vector_length), intent(inout) :: R 
!
      integer :: i
!
      if (.not. tool%active) call output%error_msg('asked to project, but projection tool is set to inactive') 
      if (.not. allocated(tool%projector)) call output%error_msg('cannot project without projector allocated')
!
!$omp parallel do private(i)
      do i = 1, tool%vector_length
!
         R(i) = R(i)*tool%projector(i)
!
      enddo  
!$omp end parallel do
!
   end subroutine do_es_projection_tool
!
!
end module es_projection_tool_class
