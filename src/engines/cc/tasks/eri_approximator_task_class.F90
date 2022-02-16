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
module eri_approximator_task_class
!
!!
!! ERI approximator task class
!! Written by Alex C. Paul, 2022
!!
!
   use ccs_class,     only: ccs
   use cc_task_class, only: cc_task
   use eri_approximator_factory_class, only: eri_approximator_factory
!
   implicit none
!
   type, extends(cc_task) :: eri_approximator_task

      type(eri_approximator_factory), allocatable :: factory
      class(cc_task), allocatable :: eri_approximator
!
   contains
!
      procedure, public :: execute &
                        => execute_eri_approximator_task
!
      procedure, public :: print_header &
                        => print_header_eri_approximator_task
!
   end type eri_approximator_task
!
!
contains
!
!
   subroutine execute_eri_approximator_task(this, wf)
!!
!!    Execute
!!    Written by Alex C. Paul, 2022
!!
      implicit none
!
      class(eri_approximator_task), intent(inout) :: this
      class(ccs), target,           intent(inout) :: wf
!
      call this%print_header()
!
      this%factory = eri_approximator_factory()
      this%eri_approximator = this%factory%create()
      call this%eri_approximator%execute(wf)
!
   end subroutine execute_eri_approximator_task
!
!
   subroutine print_header_eri_approximator_task(this)
!!
!!    Print header
!!    Written by Eirik F. Kjønstad, 2022
!!
      use warning_suppressor
!
      implicit none
!
      class(eri_approximator_task) :: this
!
      call do_nothing(this) ! wrapper object for a task that will print banner
!
   end subroutine print_header_eri_approximator_task
!
!
end module eri_approximator_task_class
