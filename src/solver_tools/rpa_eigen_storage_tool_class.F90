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
module rpa_eigen_storage_tool_class
!
!!
!!    RPA eigen storage tool class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!
   use parameters
   use eigen_storage_tool_class, only: eigen_storage_tool
   use hf_class, only: hf
!
   implicit none
!
   type, extends(eigen_storage_tool) :: rpa_eigen_storage_tool
!
      class(hf), pointer, private :: wf
!
   contains
!
      procedure, public :: store => store_rpa_eigen_storage_tool
!
   end type  rpa_eigen_storage_tool
!
   interface  rpa_eigen_storage_tool
!
      procedure :: new_rpa_eigen_storage_tool
!
   end interface  rpa_eigen_storage_tool
!
contains
!
   function new_rpa_eigen_storage_tool(wf) result(this)
!!
!!    New RPA transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(in), target :: wf
      type(rpa_eigen_storage_tool) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_o*wf%n_v
!
   end function new_rpa_eigen_storage_tool
!
   subroutine store_rpa_eigen_storage_tool(this, eigenvalue, eigenvector, I)
!!
!!    Store
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(rpa_eigen_storage_tool),            intent(inout) :: this
      integer,                                  intent(in) :: I
      real(dp), dimension(this%n_parameters),   intent(in) :: eigenvector
      real(dp),                                 intent(in) :: eigenvalue
!
      call this%wf%save_tdhf_vector(eigenvector, eigenvalue, I)
!
      this%wf%tdhf_excitation_energies(I) = sqrt(eigenvalue) ! Note that RPA excitation energies 
                                                             ! are the square roots of the 
                                                             ! eigenvalues. 
                                                             ! (see tdhf_transformation_hf.F90)
!
   end subroutine store_rpa_eigen_storage_tool
!
!
end module rpa_eigen_storage_tool_class
