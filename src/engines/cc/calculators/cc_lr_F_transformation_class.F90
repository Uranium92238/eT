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
module cc_lr_F_transformation_class
!
!!
!! CC LR F transformation class
!! Written by Eirik F. Kjønstad, Jan 2022
!!
!
   use parameters
   use ccs_class, only: ccs
   use cc_F_transformation_class, only: cc_F_transformation
!
   implicit none
!
   type, extends(cc_F_transformation) :: cc_lr_F_transformation
!
   contains
!
      procedure, public :: transform &
                        => transform_cc_lr_F_transformation
!
   end type cc_lr_F_transformation
!
!
   interface cc_lr_F_transformation
!
      procedure :: new_cc_lr_F_transformation
!
   end interface cc_lr_F_transformation
!
!
contains
!
!
   function new_cc_lr_F_transformation(wf) result(this)
!!
!!    New CC LR F transformation
!!    Written by Eirik F. Kjønstad, Jan 2022
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
!
      type(cc_lr_F_transformation) :: this
!
      this%wf => wf
!
   end function new_cc_lr_F_transformation
!
!
   subroutine transform_cc_lr_F_transformation(this, t, Ft)
!!
!!    Transform
!!    Written by Eirik F. Kjønstad, Jan 2022
!!
      implicit none
!
      class(cc_lr_F_transformation) :: this
!
      real(dp), dimension(this%wf%n_es_amplitudes), intent(in)  :: t
      real(dp), dimension(this%wf%n_es_amplitudes), intent(out) :: Ft
!
      call this%wf%F_transformation(t, Ft)
!
   end subroutine transform_cc_lr_F_transformation
!
!
end module cc_lr_F_transformation_class
