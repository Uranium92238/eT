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
module biorthonormalization_task_class
!
!!
!! Biorthonormalization task class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use parameters
   use ccs_class,     only: ccs
   use cc_task_class, only: cc_task
   use global_in,     only: input
!
   implicit none
!
   type, extends(cc_task) :: biorthonormalization_task
!
      real(dp), private :: residual_threshold, energy_threshold
      character(len=200), private :: es_algorithm
!
   contains
!
      procedure, public :: execute &
                        => execute_biorthonormalization_task
!
   end type biorthonormalization_task
!
!
   interface biorthonormalization_task
!
      procedure :: new_biorthonormalization_task
!
   end interface biorthonormalization_task
!
!
contains
!
!
   function new_biorthonormalization_task() result(this)
!!
!!    New biorthonormalization task
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      type(biorthonormalization_task) :: this
!
      this%name_ = 'Biorthonormalization of excited CC states'
!
      call input%get_keyword('algorithm', 'solver cc es', this%es_algorithm)
!
      this%energy_threshold = 1.0d-3
      call input%get_keyword('energy threshold', 'solver cc es', this%energy_threshold)
!
      this%residual_threshold = this%energy_threshold
      call input%get_keyword('residual threshold', 'solver cc es', this%residual_threshold)
!
   end function new_biorthonormalization_task
!
!
   subroutine execute_biorthonormalization_task(this, wf)
!!
!!    Execute
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      class(biorthonormalization_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      call this%print_header()
      call this%start_timer()
!
      if (this%es_algorithm .eq. 'diis') then
         call wf%remove_parallel_states(this%residual_threshold, 'both')
      end if
!
      call wf%biorthonormalize_L_and_R(this%energy_threshold, this%residual_threshold)
!
      call this%end_timer()
!
   end subroutine execute_biorthonormalization_task
!
!
end module biorthonormalization_task_class
