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
module start_vector_tool_class
!
!!
!!    Start vector tool class module
!!    Written by Sarai D. Folkestad, May 2021
!!
   implicit none
!
   type, abstract :: start_vector_tool
!
      integer :: n_parameters
!
   contains
!
      procedure(get_start_vector_tool), deferred, public :: get
!
   end type start_vector_tool
!
!
   abstract interface
!
      subroutine get_start_vector_tool(this, start_vector, I, energy)
!!
!!       Get
!!       Written by Sarai D. Folkestad, May 2021
!!
!!       Gets the I'th start vector for some purpose
!!
         use kinds
!
         import start_vector_tool
!
         implicit none
!
         class(start_vector_tool), intent(in) :: this
         real(dp), dimension(this%n_parameters), intent(out) :: start_vector
         integer, intent(in) :: I
         real(dp), intent(out) :: energy
!
      end subroutine
!
   end interface
!
end module start_vector_tool_class
