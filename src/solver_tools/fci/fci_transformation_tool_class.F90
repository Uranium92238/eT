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
module fci_transformation_tool_class
!
!!
!! FCI transformation class
!! Written by Enrico Ronca, 2022
!!
!
   use kinds
!
   use fci_class,                   only: fci
   use transformation_tool_class,   only: transformation_tool
!
   implicit none
!
   type, extends(transformation_tool) :: fci_transformation_tool
!
      class(fci), pointer, private :: wf
!
   contains
!
      procedure, public :: transform => transform_fci
      procedure, public :: initialize  => initialize_fci
!
   end type  fci_transformation_tool
!
   interface  fci_transformation_tool
!
      procedure :: new_fci_transformation_tool
!
   end interface  fci_transformation_tool
!
contains
!
   function new_fci_transformation_tool(wf) result(this)
!!
!!    New FCI transformation
!!    Written by Enrico Ronca, 2022
!!
      implicit none
!
      class(fci), intent(in), target :: wf
      type(fci_transformation_tool) :: this
!
      this%wf => wf
      this%n_parameters = wf%n_determinants
!
   end function new_fci_transformation_tool
!
!
   subroutine transform_fci(this, trial, transform)
!!
!!    Transform
!!    Written by Enrico Ronca, 2022
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(fci_transformation_tool), intent(in)   :: this
      real(dp), dimension(this%n_parameters) :: trial, transform
      call this%wf%hamiltonian_transformation_no_spin_symmetry(trial, transform)
!
   end subroutine transform_fci
!
!
   subroutine initialize_fci(this)
!!
!!    Initialize
!!    Written by Enrico Ronca, 2022
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(fci_transformation_tool), intent(in)   :: this
!
      call do_nothing(this)
!
   end subroutine initialize_fci
!
!
end module fci_transformation_tool_class
