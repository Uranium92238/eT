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
module cc_jacobian_transformation_class
!
!!
!!    Coupled Cluster Jacobian transformation class
!!    Written by Regina Matveeva, Sept 2021
!!
!!    Based on abstract_jacobian_transformer_class by Eirik F. Kjønstad, 2021
!
   use global_out, only: output
   use kinds
!
   use ccs_class,                   only: ccs
   use transformation_class,   only: transformation
!
   implicit none
!
   type, extends(transformation) :: cc_jacobian_transformation
!
      character(len=:), allocatable :: side ! 'right' or 'left'
      class(ccs), pointer, private  :: wf
!
   contains
!
      procedure, public :: initialize  => initialize_cc_jacobian
      procedure, public :: transform   => transform_cc_jacobian
!
   end type  cc_jacobian_transformation
!
   interface  cc_jacobian_transformation
!
      procedure :: new_cc_jacobian_transformation
!
   end interface  cc_jacobian_transformation
!
   contains
!
   function new_cc_jacobian_transformation(wf, side, n_parameters) result(this)
!!
!!    New Coupled Cluster Jacobian transformation
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      character(len=*), intent(in)   :: side
      integer, intent(in) :: n_parameters
      type(cc_jacobian_transformation) :: this
!
      this%wf => wf
      this%n_parameters = n_parameters
!
      this%side = side
!
   end function new_cc_jacobian_transformation
!
!
   subroutine transform_cc_jacobian(this, trial, transform)
!!
!!    Transform Coupled Cluster Jacobian
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(cc_jacobian_transformation), intent(in)   :: this
      real(dp), dimension(this%n_parameters) :: trial, transform
!
      if (this%side == 'left') then
         call this%wf%jacobian_transpose_transformation(trial, transform)
      elseif (this%side == 'right') then
         call this%wf%jacobian_transformation(trial, transform)
      endif
!
   end subroutine transform_cc_jacobian
!
!
   subroutine initialize_cc_jacobian(this)
!!
!!    Prepare for Coupled Cluster Jacobian
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(cc_jacobian_transformation), intent(in) :: this
!
      call output%printf('v', '- Prepare for multiplier equation')
!
      call this%wf%prepare_for_Jacobians(this%side)
!
end subroutine initialize_cc_jacobian
!
!
end module cc_jacobian_transformation_class
