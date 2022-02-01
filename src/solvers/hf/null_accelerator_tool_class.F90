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
module null_accelerator_tool_class
!
!!
!! Null accelerator tool,
!! Written by Sarai D. Folkestad, 2020
!!
!! Performs the null acceleration (no acceleration)
!!
!
   use parameters
!
   use accelerator_tool_class, only: accelerator_tool
!
   implicit none
!
   type, extends(accelerator_tool) :: null_accelerator_tool
!
   contains
!
      procedure :: do_        => do_null_accelerator_tool
      procedure :: initialize => initialize_null_accelerator_tool
      procedure :: finalize   => finalize_null_accelerator_tool
!
   end type null_accelerator_tool
!
   interface null_accelerator_tool
!
      procedure :: new_null_accelerator_tool
!
   end interface null_accelerator_tool
!
contains
!
   function new_null_accelerator_tool(x_dimension, e_dimension) result(this)
!!
!!    New null accelerator tool
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      type(null_accelerator_tool)   :: this
      integer, intent(in)           :: x_dimension
      integer, intent(in)           :: e_dimension
!
      this%x_dimension = x_dimension
      this%e_dimension = e_dimension
!
   end function new_null_accelerator_tool
!
!
   subroutine do_null_accelerator_tool(this, x, e)
!!
!!    Do
!!    Written by Sarai D. Folkestad, 2020
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(null_accelerator_tool), intent(inout)           :: this
      real(dp), dimension(this%x_dimension), intent(inout)  :: x
      real(dp), dimension(this%e_dimension), intent(in)     :: e
!
      call do_nothing(this)
      call do_nothing(x)
      call do_nothing(e)
!
   end subroutine do_null_accelerator_tool
!
!
   subroutine initialize_null_accelerator_tool(this)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, 2020
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(null_accelerator_tool), intent(inout) :: this
!
      call do_nothing(this)
!
   end subroutine initialize_null_accelerator_tool
!
!
   subroutine finalize_null_accelerator_tool(this)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad, 2020
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(null_accelerator_tool), intent(inout) :: this
!
      call do_nothing(this)
!
   end subroutine finalize_null_accelerator_tool
!
!
end module null_accelerator_tool_class
