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
module diis_accelerator_tool_class
!
!!
!! DIIS accelerator tool,
!! Written by Sarai D. Folkestad, 2020
!!
!! Handles DIIS (or CROP) acceleration by use of the 
!! DIIS tool.
!!
!
   use kinds
   use parameters
!
   use global_out, only: output
!
   use accelerator_tool_class, only: accelerator_tool
   use diis_tool_class,        only: diis_tool
!
   implicit none
!
   type, extends(accelerator_tool) :: diis_accelerator_tool
!
      type(diis_tool), private :: diis
      logical,         private :: records_in_memory
!
   contains
!
      procedure :: do_        => do_diis_accelerator_tool
      procedure :: initialize => initialize_diis_accelerator_tool
      procedure :: finalize   => finalize_diis_accelerator_tool
!
      procedure, private :: print_accelerator_info => print_accelerator_info_diis_accelerator_tool
!
   end type diis_accelerator_tool
!
   interface diis_accelerator_tool
!
     procedure :: new_diis_accelerator_tool
!
   end interface diis_accelerator_tool
!
contains
!
   function new_diis_accelerator_tool(x_dimension,       &
                                      e_dimension,       &
                                      records_in_memory, &
                                      diis) result(this)
!!
!!    New diis accelerator tool
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      type(diis_accelerator_tool)        :: this
      integer,                intent(in) :: x_dimension
      integer,                intent(in) :: e_dimension
      logical,                intent(in) :: records_in_memory
      type(diis_tool),        intent(in) :: diis
!
      this%x_dimension = x_dimension
      this%e_dimension = e_dimension
!
      this%records_in_memory = records_in_memory
!
      this%diis = diis
!
      call this%print_accelerator_info()
!
   end function new_diis_accelerator_tool
!
!
   subroutine do_diis_accelerator_tool(this, x, e)
!!
!!    Do
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(diis_accelerator_tool),          intent(inout)  :: this
      real(dp), dimension(this%x_dimension), intent(inout)  :: x
      real(dp), dimension(this%e_dimension), intent(in)     :: e
!
      call this%diis%update(e, x)
!
   end subroutine do_diis_accelerator_tool
!
!
   subroutine initialize_diis_accelerator_tool(this)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Initializes records (on disk or memory) 
!!    of the DIIS history.
!!
      implicit none
!
      class(diis_accelerator_tool), intent(inout) :: this
!
      call this%diis%initialize_storers(this%records_in_memory)
!
   end subroutine initialize_diis_accelerator_tool
!
!
   subroutine finalize_diis_accelerator_tool(this)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(diis_accelerator_tool), intent(inout) :: this
!
      call this%diis%finalize_storers()
!
   end subroutine finalize_diis_accelerator_tool
!
!
   subroutine print_accelerator_info_diis_accelerator_tool(this)
!!
!!    Print
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(diis_accelerator_tool), intent(inout)  :: this
!
      call this%diis%print_settings()
!
   end subroutine print_accelerator_info_diis_accelerator_tool
!
end module diis_accelerator_tool_class
