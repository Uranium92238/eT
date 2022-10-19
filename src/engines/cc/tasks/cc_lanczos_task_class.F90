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
module cc_lanczos_task_class
!
!!
!! CC Lanczos task class
!! Written by Alexander C. Paul
!!
!
   use ccs_class,             only: ccs
   use abstract_solver_class, only: abstract_solver
   use cc_task_class,         only: cc_task
!
   implicit none
!
   type, extends(cc_task) :: cc_lanczos_task
!
      class(abstract_solver), allocatable, private :: solver
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_lanczos_task
!
   end type cc_lanczos_task
!
!
   interface cc_lanczos_task
!
      procedure :: new_cc_lanczos_task
!
   end interface cc_lanczos_task
!
!
contains
!
!
   function new_cc_lanczos_task() result(this)
!!
!!    New
!!    Written by Alexander C. Paul, Jun 2022
!!
      implicit none
!
      type(cc_lanczos_task) :: this
!
      this%name_ = 'Determining CC excited states using Lanczos solver'
!
   end function new_cc_lanczos_task
!
!
   subroutine execute_cc_lanczos_task(this, wf)
!!
!!    Execute
!!    Written by Alexander C. Paul, Jun 2022
!!
      use asymmetric_lanczos_cc_es_class, only: asymmetric_lanczos_cc_es
!
      implicit none
!
      class(cc_lanczos_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      call this%print_header()
      call this%start_timer()
!
      this%solver = asymmetric_lanczos_cc_es(wf)
!
      call wf%construct_fock()
!
      call this%solver%run()
!
      call this%solver%cleanup()
!
      call this%end_timer()
!
   end subroutine execute_cc_lanczos_task
!
!
end module cc_lanczos_task_class
