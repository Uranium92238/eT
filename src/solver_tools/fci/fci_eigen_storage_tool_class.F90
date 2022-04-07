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
module fci_eigen_storage_tool_class
!
!!
!! FCI eigen storage tool class module
!! Written by Enrico Ronca, 2022
!!
!
   use parameters
   use eigen_storage_tool_class, only: eigen_storage_tool
   use fci_class, only: fci
!
   implicit none
!
   type, extends(eigen_storage_tool) :: fci_eigen_storage_tool
!
      class(fci), pointer, private :: wf
!
   contains
!
      procedure, public :: store => store_fci_eigen_storage_tool
!
   end type  fci_eigen_storage_tool
!
   interface  fci_eigen_storage_tool
!
      procedure :: new_fci_eigen_storage_tool
!
   end interface  fci_eigen_storage_tool
!
contains
!
   function new_fci_eigen_storage_tool(wf) result(this)
!!
!!    New FCI transformation
!!    Written by Enrico Ronca, 2022
!!
      implicit none
!
      class(fci), intent(in), target :: wf
      type(fci_eigen_storage_tool) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_determinants
!
   end function new_fci_eigen_storage_tool
!
   subroutine store_fci_eigen_storage_tool(this, eigenvalue, eigenvector, I)
!!
!!    Store
!!    Written by Enrico Ronca, 2022
!!
      implicit none
!
      class(fci_eigen_storage_tool),            intent(inout) :: this
      integer,                                  intent(in) :: I
      real(dp), dimension(this%n_parameters),   intent(in) :: eigenvector
      real(dp),                                 intent(in) :: eigenvalue
!
      this%wf%energies(I) = eigenvalue
      call dcopy(this%n_parameters, eigenvector, 1, this%wf%ci_coefficients(:,:,I), 1)
!
      call this%wf%save_fci_state(eigenvector, eigenvalue, I)
!
   end subroutine store_fci_eigen_storage_tool
!
!
end module fci_eigen_storage_tool_class
