!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module null_linear_storage_tool_class
!
!!
!!    Null linear storage tool class module
!!    Written by Regina Matveeva and Sarai D. Folkestad, Oct 2021
!!
!!    Does nothing
!!
!
   use parameters
   use linear_equation_storage_tool_class, only: linear_equation_storage_tool
!
   implicit none
!
      type, extends(linear_equation_storage_tool) :: null_linear_storage_tool
!
   contains
!
      procedure, public :: store => store_null
!
   end type  null_linear_storage_tool
!
   interface null_linear_storage_tool
!
      procedure :: new_null_linear_storage_tool
!
   end interface
!
contains
!
!
   pure function new_null_linear_storage_tool() result(this)
!!
!!    New
!!    Written by Regina Matveeva and Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      type(null_linear_storage_tool) :: this
!
      this%n_parameters = 0
!
   end function new_null_linear_storage_tool
!
!
   subroutine store_null(this, solution_vector, n)
!!
!!    Store
!!    Written by Regina Matveeva and Sarai D. Folkestad, Oct 2021
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(null_linear_storage_tool),          intent(inout) :: this
      integer,                                  intent(in) :: n
      real(dp), dimension(this%n_parameters),   intent(in) :: solution_vector
!
      call do_nothing(this)
      call do_nothing(solution_vector)
      call do_nothing(n)
!
   end subroutine store_null
!
!
end module null_linear_storage_tool_class
