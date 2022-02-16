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
module cc_mean_values_task_class
!
!!
!! CC mean value task class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class,     only: ccs
   use cc_task_class, only: cc_task
   use global_in,     only: input
!
   implicit none
!
   type, extends(cc_task) :: cc_mean_values_task
!
      logical, private :: dipole
      logical, private :: quadrupole
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_mean_values_task
!
   end type cc_mean_values_task
!
!
   interface cc_mean_values_task
!
      procedure :: new_cc_mean_values_task
!
   end interface cc_mean_values_task
!
!
contains
!
!
   function new_cc_mean_values_task() result(this)
!!
!!    New CC mean values task
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(cc_mean_values_task) :: this
!
      this%name_ = 'Determining CC mean values'
!
      this%dipole     = input%is_keyword_present('dipole', 'cc mean value')
      this%quadrupole = input%is_keyword_present('quadrupole', 'cc mean value')
!
   end function new_cc_mean_values_task
!
!
   subroutine execute_cc_mean_values_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_mean_values_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      call this%print_header()
      call this%start_timer()
!
      call wf%initialize_gs_density()
      call wf%prepare_for_properties()
      call wf%construct_gs_density()
!
      if (this%dipole) call wf%calculate_and_print_dipole()
      if (this%quadrupole) call wf%calculate_and_print_quadrupole()
!
      call this%end_timer()
!
   end subroutine execute_cc_mean_values_task
!
!
end module cc_mean_values_task_class
