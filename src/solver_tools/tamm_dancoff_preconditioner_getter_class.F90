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
module tamm_dancoff_preconditioner_getter_class
!
!!
!!    Tamm-dancoff preconditioner getter class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!
   use parameters
   use preconditioner_getter_class, only: preconditioner_getter
   use hf_class, only: hf
!
   implicit none
!
   type, extends(preconditioner_getter) :: tamm_dancoff_preconditioner_getter
!
      class(hf), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_tamm_dancoff_preconditioner_getter
!
   end type tamm_dancoff_preconditioner_getter
!
!
   interface  tamm_dancoff_preconditioner_getter
!
      procedure :: new_tamm_dancoff_preconditioner_getter
!
   end interface  tamm_dancoff_preconditioner_getter
!
contains
!
   function new_tamm_dancoff_preconditioner_getter(wf) result(this)
!!
!!    New Tamm Dancoff preconditioner getter
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(in), target :: wf
      type(tamm_dancoff_preconditioner_getter) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_o*wf%n_v
!
   end function new_tamm_dancoff_preconditioner_getter
!
!
   subroutine get_tamm_dancoff_preconditioner_getter(this, preconditioner)
!!
!!    Get
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(tamm_dancoff_preconditioner_getter), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: preconditioner
!
      call this%wf%get_tamm_dancoff_preconditioner(preconditioner)
!
   end subroutine get_tamm_dancoff_preconditioner_getter
!
end module tamm_dancoff_preconditioner_getter_class
