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
module cc_triplet_jacobian_transformation_class
!
!!
!!    CC triplet Jacobian transformation class
!!    Written by Regina Matveeva, Sept 2021
!!
!!    Based on abstract_jacobian_transformer_class by Eirik F. Kjønstad, 2021
!
   use global_out, only: output
   use kinds
!
   use ccs_class,             only: ccs
   use transformation_class,  only: transformation
!
   implicit none
!
   type, extends(transformation) :: cc_triplet_jacobian_transformation
!
      character(len=:), allocatable :: side ! 'right' or 'left'
      class(ccs), pointer, private  :: wf
!
   contains
!
      procedure, public :: initialize  => initialize_cc_triplet_jacobian
      procedure, public :: transform   => transform_cc_triplet_jacobian
!
   end type  cc_triplet_jacobian_transformation
!
   interface  cc_triplet_jacobian_transformation
!
      procedure :: new_cc_triplet_jacobian_transformation
!
   end interface  cc_triplet_jacobian_transformation
!
   contains
!
   function new_cc_triplet_jacobian_transformation(wf, side) result(this)
!!
!!    New
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      character(len=*), intent(in)   :: side
      type(cc_triplet_jacobian_transformation) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_triplet_amplitudes
!
      this%side = side
!
   end function new_cc_triplet_jacobian_transformation
!
!
   subroutine transform_cc_triplet_jacobian(this, trial, transform)
!!
!!    Transform
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(cc_triplet_jacobian_transformation), intent(in)   :: this
      real(dp), dimension(this%n_parameters) :: trial, transform
!
      call this%wf%construct_triplet_jacobian_transform(this%side, trial, transform)
!
   end subroutine transform_cc_triplet_jacobian
!
!
   subroutine initialize_cc_triplet_jacobian(this)
!!
!!    Initialize
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(cc_triplet_jacobian_transformation), intent(in) :: this
!
      call this%wf%prepare_for_triplet_jacobians(this%side)
!
end subroutine initialize_cc_triplet_jacobian
!
!
end module cc_triplet_jacobian_transformation_class
