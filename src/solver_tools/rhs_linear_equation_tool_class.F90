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
module rhs_linear_equation_tool_class
!
!!
!!    Abstract right-hand side linear equation class module
!!    Written by Regina Matveeva, Sept 2021
!!
!!    Gets the right-hand side of a linear equation
!!
!
   implicit none
!
   type, abstract :: rhs_linear_equation_tool
!
      integer :: n_parameters
!
      integer :: n_rhs
!
   contains
!
      procedure (get_rhs_linear_equation_tool), deferred, public :: get
!
   end type  rhs_linear_equation_tool
!
   abstract interface
!
   subroutine get_rhs_linear_equation_tool(this, rhs)
!!
!!    Get
!!    Written by Regina Matveeva, Sept 2021
!!
      use kinds
!
      import rhs_linear_equation_tool
!
      implicit none
!
      class(rhs_linear_equation_tool),                    intent(inout) :: this
      real(dp), dimension(this%n_parameters, this%n_rhs), intent(inout)  :: rhs
!
   end subroutine get_rhs_linear_equation_tool
!
   end interface
!
end module rhs_linear_equation_tool_class
