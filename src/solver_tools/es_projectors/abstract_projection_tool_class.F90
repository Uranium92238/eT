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
module abstract_projection_tool_class
!
!!
!!    Abstract projector tool class module
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
   use memory_manager_class, only: mem
!
   implicit none
!
   type, abstract :: abstract_projection_tool
!
      logical :: active ! Do projection if true, don't if false
!
      integer :: n_parameters
      real(dp), dimension(:), allocatable :: projector

!
   contains
!
      procedure :: project => project_abstract_projection_tool
      procedure :: initialize
      procedure :: cleanup
      procedure (get_projector_abstract), deferred :: get_projector
!
   end type abstract_projection_tool
!
!
   abstract interface
!
      subroutine get_projector_abstract(tool, projector)
!
         use parameters
         import abstract_projection_tool
!
         implicit none
!
         class(abstract_projection_tool), intent(in) :: tool
         real(dp), dimension(tool%n_parameters), intent(out) :: projector
!
      end subroutine get_projector_abstract
!
   end interface
!
!
contains 
!
!
   subroutine initialize(tool)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, 2022
!!
      implicit none
!
      class(abstract_projection_tool), intent(inout) :: tool
!
      if (.not. tool%active) return
!
      call mem%alloc(tool%projector, tool%n_parameters)
      call tool%get_projector(tool%projector)
!
   end subroutine initialize
!
!
   subroutine cleanup(tool)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, 2022
!!
      implicit none
!
      class(abstract_projection_tool), intent(inout) :: tool
!
      if (.not. tool%active) return
      if (allocated(tool%projector)) call mem%dealloc(tool%projector, tool%n_parameters)
!
   end subroutine cleanup
!
!
   subroutine project_abstract_projection_tool(tool, R)
!!
!!    Project
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
      class(abstract_projection_tool), intent(in) :: tool
!
      real(dp), dimension(tool%n_parameters), intent(inout) :: R
!
      integer :: i

      if (.not. tool%active) return
!
!$omp parallel do private(i)
      do i = 1, tool%n_parameters
!
         R(i) = R(i)*tool%projector(i)
!
      enddo  
!$omp end parallel do
!
   end subroutine project_abstract_projection_tool
!
!
end module abstract_projection_tool_class
