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
module accelerator_tool_class
!
!!
!!    Accelerator tool
!!    Written by Sarai D. Folkestad, 2020-2021
!!
!!    Abstract class for handling of
!!    convergence acceleration to solve equations
!!
!
   use kinds
   use parameters
!
   implicit none
!
   type, abstract :: accelerator_tool
!
      integer :: x_dimension
      integer :: e_dimension
!
   contains
!
      procedure (do_acceleration),                  deferred :: do_
      procedure (initialize_finalize_acceleration), deferred :: initialize
      procedure (initialize_finalize_acceleration), deferred :: finalize
!
   end type accelerator_tool
!
   abstract interface
!
      subroutine do_acceleration(this, x, e)
!
         use kinds
         import accelerator_tool
!
         implicit none
!
         class(accelerator_tool),               intent(inout)  :: this
         real(dp), dimension(this%x_dimension), intent(inout)  :: x
         real(dp), dimension(this%e_dimension), intent(in)     :: e
!
      end subroutine
!
      subroutine initialize_finalize_acceleration(this)
!
         import accelerator_tool
!
         implicit none
!
         class(accelerator_tool), intent(inout) :: this
!
      end subroutine
!
  end interface
!
end module accelerator_tool_class
