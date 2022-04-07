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
module fci_preconditioner_getter_class
!
!!
!! FCI preconditioner getter class module
!! Written by Enrico Ronca, 2022
!!
!
   use parameters
   use preconditioner_getter_class, only: preconditioner_getter
   use fci_class, only: fci
!
   implicit none
!
   type, extends(preconditioner_getter) :: fci_preconditioner_getter
!
      class(fci), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_fci_preconditioner_getter
!
   end type fci_preconditioner_getter
!
   interface  fci_preconditioner_getter
!
      procedure :: new_fci_preconditioner_getter
!
   end interface  fci_preconditioner_getter
!
contains
!
   function new_fci_preconditioner_getter(wf) result(this)
!!
!!    New FCI preconditioner tool
!!    Written by Enrico Ronca, 2022
!!
      implicit none
!
      class(fci), intent(in), target :: wf
      type(fci_preconditioner_getter) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_determinants
!
   end function new_fci_preconditioner_getter
!
!
   subroutine get_fci_preconditioner_getter(this, preconditioner)
!!
!!    Get
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(fci_preconditioner_getter), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: preconditioner
!
      call this%wf%construct_h_diagonal(preconditioner)
!
   end subroutine get_fci_preconditioner_getter
!
end module fci_preconditioner_getter_class
