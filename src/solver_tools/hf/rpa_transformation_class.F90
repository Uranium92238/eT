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
module rpa_transformation_class
!
!!
!!    RPA transformation class
!!    Written by Sarai D. Folkestad, 2021
!!
!
   use kinds
!
   use hf_class,                    only: hf
   use transformation_class,   only: transformation
!
   implicit none
!
   type, extends(transformation) :: rpa_transformation
!
      class(hf), pointer, private :: wf
!
   contains
!
      procedure, public :: transform => transform_rpa
      procedure, public :: initialize  => initialize_rpa
!
   end type  rpa_transformation
!
   interface  rpa_transformation
!
      procedure :: new_rpa_transformation
!
   end interface  rpa_transformation
!
contains
!
   function new_rpa_transformation(wf) result(this)
!!
!!    New RPA transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(in), target :: wf
      type(rpa_transformation) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_o*wf%n_v
!
   end function new_rpa_transformation
!
!
   subroutine transform_rpa(this, trial, transform)
!!
!!    Transform
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Modified by Regina Matveeva, Sep 2021
!!    Added frequency (necessary due to a modification of the transformation)
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(rpa_transformation), intent(in)   :: this
      real(dp), dimension(this%n_parameters) :: trial, transform
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
      class(rpa_transformation), intent(in)   :: this
!
      call do_nothing(this)
!
   end subroutine initialize_rpa
!
!
end module rpa_transformation_class
