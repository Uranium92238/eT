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
module fci_start_vector_tool_class
!
!!
!! FCI start vector tool class module
!! Written by Enrico Ronca, 2022
!!
!
   use kinds
   use start_vector_tool_class, only: start_vector_tool
   use fci_class,               only: fci
!
   implicit none
!
   type, extends(start_vector_tool) :: fci_start_vector_tool
!
      logical, private :: restart
      class(fci), pointer, private :: wf
!
   contains
!
      procedure, public :: get => get_fci_start_vector_tool
!
   end type fci_start_vector_tool
!
   interface  fci_start_vector_tool
!
      procedure :: new_fci_start_vector_tool
!
   end interface  fci_start_vector_tool
!
contains
!
!
   function new_fci_start_vector_tool(wf, restart) result(this)
!!
!!    New FCI start vector tool
!!    Written by Enrico Ronca, 2022
!!
      implicit none
!
      class(fci), intent(in), target :: wf
      logical, intent(in) :: restart
      type(fci_start_vector_tool) :: this
!
      this%restart = restart
      this%wf => wf
      this%n_parameters = wf%n_determinants
!
   end function new_fci_start_vector_tool
!
!
   subroutine get_fci_start_vector_tool(this, start_vector, I, energy)
!!
!!    Get
!!    Written by Enrico Ronca, 2022
!!
!!    Gets the I'th start vector for some purpose
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(fci_start_vector_tool), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: start_vector
      integer, intent(in) :: I
      real(dp), intent(out) :: energy
!
      call do_nothing(energy)
      call this%wf%get_fci_start_vector(I, start_vector, this%restart)
!
   end subroutine
!
!
end module fci_start_vector_tool_class
