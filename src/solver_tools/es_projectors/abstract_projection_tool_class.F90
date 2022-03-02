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
      integer :: n_parameters
!
   contains
!
      procedure (project_abstract), deferred :: project
!
   end type abstract_projection_tool
!
!
   abstract interface
!
      subroutine project_abstract(tool, vector)
!
         use parameters
         import abstract_projection_tool
!
         implicit none
!
         class(abstract_projection_tool),        intent(in)    :: tool
         real(dp), dimension(tool%n_parameters), intent(inout) :: vector
!
      end subroutine project_abstract
!
   end interface
!
!
contains 
!
!
end module abstract_projection_tool_class
