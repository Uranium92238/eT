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
module cc_eigen_storage_tool_class
!
!!
!!    CC eigen storage tool class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!
   use parameters
   use eigen_storage_tool_class, only: eigen_storage_tool
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(eigen_storage_tool) :: cc_eigen_storage_tool
!
      class(ccs), pointer, private :: wf
      character(len=200) :: side
!
   contains
!
      procedure, public :: store => store_cc_eigen_storage_tool
!
   end type  cc_eigen_storage_tool
!
   interface  cc_eigen_storage_tool
!
      procedure :: new_cc_eigen_storage_tool
!
   end interface  cc_eigen_storage_tool
!
contains
!
   function new_cc_eigen_storage_tool(wf, side) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(ccs), intent(inout), target :: wf
      character(len=*), intent(in) :: side
      type(cc_eigen_storage_tool) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_es_amplitudes
      this%side = side
!
   end function new_cc_eigen_storage_tool
!
   subroutine store_cc_eigen_storage_tool(this, eigenvalue, eigenvector, I)
!!
!!    Store
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(cc_eigen_storage_tool),             intent(inout) :: this
      integer,                                  intent(in) :: I
      real(dp), dimension(this%n_parameters),   intent(in) :: eigenvector
      real(dp),                                 intent(in) :: eigenvalue
!
      call this%wf%save_excited_state(eigenvector, I, I, this%side, [eigenvalue])
!
      if (trim(this%side) == 'right') this%wf%right_excitation_energies(I) = eigenvalue
      if (trim(this%side) == 'left') this%wf%left_excitation_energies(I) = eigenvalue
!
   end subroutine store_cc_eigen_storage_tool
!
!
end module cc_eigen_storage_tool_class
