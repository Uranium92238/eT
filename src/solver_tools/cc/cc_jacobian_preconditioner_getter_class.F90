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
module cc_jacobian_preconditioner_getter_class
!
!!
!!    Coupled Cluster jacobian preconditioner getter class module
!!    Written by Regina Matveeva, Sept 2021
!!
!
   use parameters
   use preconditioner_getter_class, only: preconditioner_getter
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(preconditioner_getter) :: cc_jacobian_preconditioner_getter
!
      class(ccs), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_cc_jacobian_preconditioner_getter
!
   end type cc_jacobian_preconditioner_getter
!
!
   interface  cc_jacobian_preconditioner_getter
!
      procedure :: new_cc_jacobian_preconditioner_getter
!
   end interface  cc_jacobian_preconditioner_getter
!
!
contains
!
   function new_cc_jacobian_preconditioner_getter(wf, n_parameters) result(this)
!!
!!    New Coupled Cluster jacobian preconditioner getter
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ccs), intent(in), target :: wf
      integer, intent(in) :: n_parameters
      type(cc_jacobian_preconditioner_getter) :: this
!
      this%wf => wf
      this%n_parameters = n_parameters
!
   end function new_cc_jacobian_preconditioner_getter
!
!
   subroutine get_cc_jacobian_preconditioner_getter(this, preconditioner)
!!
!!    Get
!!    Written by Regina Matveeva, Sep 2021
!!
      implicit none
!
      class(cc_jacobian_preconditioner_getter), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: preconditioner
!
      call this%wf%get_orbital_differences(preconditioner, this%n_parameters)
!
   end subroutine get_cc_jacobian_preconditioner_getter
!
!
end module cc_jacobian_preconditioner_getter_class
