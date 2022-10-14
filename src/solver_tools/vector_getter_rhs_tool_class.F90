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
module vector_getter_rhs_tool_class
!
!!
!!    Vector getter RHS tool
!!    Written by Eirik Kjønstad, Oct 2022
!!
!
   use parameters
   use rhs_linear_equation_tool_class, only: rhs_linear_equation_tool
!
   implicit none
!
   type, extends(rhs_linear_equation_tool) :: vector_getter_rhs_tool
!
      real(dp), dimension(:,:), pointer :: rhs 
!
   contains
!
      procedure, public :: get => get_vector_getter_rhs_tool
!
   end type vector_getter_rhs_tool
!
!
   interface vector_getter_rhs_tool
!
      procedure :: new_vector_getter_rhs_tool
!
   end interface vector_getter_rhs_tool
!
!
contains
!
   function new_vector_getter_rhs_tool(n_parameters, n_rhs, rhs) result(this)
!!
!!    New
!!    Written by Eirik Kjønstad, Oct 2022
!!
      implicit none
!
      type(vector_getter_rhs_tool) :: this
!
      integer, intent(in) :: n_parameters, n_rhs
!
      real(dp), dimension(n_parameters * n_rhs), target :: rhs 
!
      this%n_parameters = n_parameters
      this%n_rhs        = n_rhs 
!
      this%rhs(1 : n_parameters, 1 : n_rhs) => rhs(1 : n_parameters * n_rhs)
!
   end function new_vector_getter_rhs_tool
!
!
   subroutine get_vector_getter_rhs_tool(this, rhs)
!!
!!    Get
!!    Written by Eirik Kjønstad, Oct 2022
!!
      implicit none
!
      class(vector_getter_rhs_tool),  intent(inout) :: this
!
      real(dp), dimension(this%n_parameters, this%n_rhs), intent(inout) :: rhs
!
      call dcopy(this%n_parameters * this%n_rhs, this%rhs, 1, rhs, 1)
!
   end subroutine get_vector_getter_rhs_tool
!
!
end module vector_getter_rhs_tool_class
