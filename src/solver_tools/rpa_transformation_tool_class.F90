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
module rpa_transformation_tool_class
!
!!
!!    RPA transformation class
!!    Written by Sarai D. Folkestad, 2021
!!
!
   use kinds
!
   use hf_class,                    only: hf
   use transformation_tool_class,   only: transformation_tool
!
   implicit none
!
   type, extends(transformation_tool) :: rpa_transformation_tool
!
      class(hf), pointer, private :: wf
!
   contains
!
      procedure, public :: transform => transform_rpa
      procedure, public :: initialize  => initialize_rpa
!
   end type  rpa_transformation_tool
!
   interface  rpa_transformation_tool
!
      procedure :: new_rpa_transformation_tool
!
   end interface  rpa_transformation_tool
!
contains
!
   function new_rpa_transformation_tool(wf) result(this)
!!
!!    New RPA transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(in), target :: wf
      type(rpa_transformation_tool) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_o*wf%n_v
!
   end function new_rpa_transformation_tool
!
!
   subroutine transform_rpa(this, trial, transform, frequency)
!!
!!    Transform
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Modified by Regina Matveeva, Sep 2021
!!    Added frequency (necessary due to a modification of the transformation_tool)
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(rpa_transformation_tool), intent(in)   :: this
      real(dp), dimension(this%n_parameters) :: trial, transform
      real(dp), intent(in) :: frequency
!
      call do_nothing(frequency)
      call this%wf%rpa_transformation(trial, transform)
!
   end subroutine transform_rpa
!
!
   subroutine initialize_rpa(this)
!!
!!    Initialize
!!    Written by Regina Matveeva, Sept 2021
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(rpa_transformation_tool), intent(in)   :: this
!
      call do_nothing(this)
!
   end subroutine initialize_rpa
!
!
end module rpa_transformation_tool_class
