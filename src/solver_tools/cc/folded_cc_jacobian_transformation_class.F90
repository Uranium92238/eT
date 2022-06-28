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
module folded_cc_jacobian_transformation_class
!
!!
!!    Perturbative CC Jacobian transformation class
!!    Written by Regina Matveeva, Sept 2021
!!
!
   use global_out, only: output
   use kinds
!
   use ccs_class,                   only: ccs
   use transformation_class,   only: transformation
!
   implicit none
!
   type, extends(transformation) :: folded_cc_jacobian_transformation
!
      real(dp), private :: frequency
!
      character(len=:), allocatable, private :: side ! 'right' or 'left'
      class(ccs), pointer, private  :: wf
!
   contains
!
      procedure, public :: initialize  => initialize_cc_jacobian
      procedure, public :: transform   => transform_cc_jacobian
!
   end type  folded_cc_jacobian_transformation
!
   interface  folded_cc_jacobian_transformation
!
      procedure :: new_folded_cc_jacobian_transformation
!
   end interface  folded_cc_jacobian_transformation
!
contains
!
   function new_folded_cc_jacobian_transformation(wf, side, &
                                                       n_parameters, frequency) result(this)
!!
!!    New perturbative CC Jacobian transformation
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      character(len=*), intent(in)   :: side
      integer, intent(in) :: n_parameters
      real(dp), intent(in) :: frequency

      type(folded_cc_jacobian_transformation) :: this
!
      this%wf => wf
      this%n_parameters = n_parameters
!
      this%side = side
      this%frequency = frequency
!
   end function new_folded_cc_jacobian_transformation
!
!
   subroutine transform_cc_jacobian(this, trial, transform)
!!
!!    Transform CC Jacobian
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(folded_cc_jacobian_transformation), intent(in)   :: this
      real(dp), dimension(this%n_parameters) :: trial, transform
!
      call this%wf%construct_Jacobian_transform(this%side, trial, transform, this%frequency)
!
   end subroutine transform_cc_jacobian
!
!
   subroutine initialize_cc_jacobian(this)
!!
!!    initialize CC Jacobian
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
      class(folded_cc_jacobian_transformation), intent(in) :: this
!
      call output%printf('v', '- Prepare for multiplier equation')
!
      call this%wf%prepare_for_Jacobians(this%side)
!
end subroutine initialize_cc_jacobian
!
!
end module folded_cc_jacobian_transformation_class
