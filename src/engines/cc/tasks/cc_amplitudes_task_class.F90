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
module cc_amplitudes_task_class
!
!!
!! CC amplitudes task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use ccs_class,                          only: ccs
   use abstract_solver_class,              only: abstract_solver
   use cc_task_class,                      only: cc_task
   use cc_amplitudes_solver_factory_class, only: cc_amplitudes_solver_factory
!
   implicit none
!
   type, extends(cc_task) :: cc_amplitudes_task
!
      class(abstract_solver),             allocatable, private :: solver
      type(cc_amplitudes_solver_factory),              private :: solver_factory
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_amplitudes_task
!
   end type cc_amplitudes_task
!
!
   interface cc_amplitudes_task
!
      procedure :: new_cc_amplitudes_task
!
   end interface cc_amplitudes_task
!
!
contains
!
!
   function new_cc_amplitudes_task() result(this)
!!
!!    New CC amplitudes task
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(cc_amplitudes_task) :: this
!
      this%name_ = 'Determining CC cluster amplitudes'
!
   end function new_cc_amplitudes_task
!
!
   subroutine execute_cc_amplitudes_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_amplitudes_task), intent(inout) :: this
      class(ccs), target,        intent(inout) :: wf
!
      call this%print_header()
      call this%start_timer()
!
      call this%solver_factory%create(wf, this%solver)
      call this%solver%run()
      call this%solver%cleanup()
!
      call this%end_timer()
!
   end subroutine execute_cc_amplitudes_task
!
!
end module cc_amplitudes_task_class
