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
module cc_es_eigen_davidson_print_tool_class
!
!!
!!    CC ES eigen davidson print tool class module
!!    Written by Regina Matveeva, Dec 2021
!!
!
   use ccs_class, only: ccs
   use eigen_davidson_print_tool_class, only: eigen_davidson_print_tool
!
   implicit none
!
   type, extends(eigen_davidson_print_tool) :: cc_es_eigen_davidson_print_tool
!
      character(len=200) :: side
      class(ccs), pointer, private :: wf
!
   contains
!
      procedure, public :: print_summary &
                                => print_summary_cc_es_eigen_davidson_print_tool
!
   end type cc_es_eigen_davidson_print_tool
!
   interface cc_es_eigen_davidson_print_tool
!
      procedure :: new_cc_es_eigen_davidson_print_tool
!
   end interface cc_es_eigen_davidson_print_tool
!
contains
!
!
   function new_cc_es_eigen_davidson_print_tool(wf, side) result(this)
!!
!!    New cc es eigen davidson print tool
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(ccs), target, intent(in) :: wf
      character(len=*), intent(in) :: side
!
      type(cc_es_eigen_davidson_print_tool) :: this
!
      this%side = side
      this%wf => wf
!
   end function new_cc_es_eigen_davidson_print_tool
!
!
   subroutine print_summary_cc_es_eigen_davidson_print_tool(this, n_solutions, converged, iteration)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(cc_es_eigen_davidson_print_tool), intent(in) :: this
      integer, intent(in) :: n_solutions, iteration
      logical, dimension(n_solutions), intent(in) :: converged
!
      call this%eigen_davidson_print_tool%print_summary(n_solutions, converged, iteration)
!
      if ( all(converged)) call this%wf%print_es_summary(trim(this%side))
!
   end subroutine print_summary_cc_es_eigen_davidson_print_tool
!
!
end module cc_es_eigen_davidson_print_tool_class
