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
module cc_es_amplitudes_task_class
!
!!
!! CC excited state amplitudes task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use ccs_class,                               only: ccs
   use abstract_solver_class,                   only: abstract_solver
   use cc_task_class,                           only: cc_task
   use cc_es_amplitudes_solver_factory_class,   only: cc_es_amplitudes_solver_factory
!
   implicit none
!
   type, extends(cc_task) :: cc_es_amplitudes_task
!
      class(abstract_solver),                allocatable, private :: solver
      type(cc_es_amplitudes_solver_factory), allocatable, private :: solver_factory
!
      character(len=200), private :: transformation
      logical, private :: restart
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_es_amplitudes_task
!
   end type cc_es_amplitudes_task
!
!
   interface cc_es_amplitudes_task
!
      procedure :: new_cc_es_amplitudes_task
!
   end interface cc_es_amplitudes_task
!
!
contains
!
!
   function new_cc_es_amplitudes_task(transformation, restart) result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      character(len=*), intent(in) :: transformation
!
      logical, intent(in) :: restart
!
      type(cc_es_amplitudes_task) :: this
!
      this%name_ = 'Determining CC excited state amplitudes'
!
      this%transformation = transformation
      this%restart = restart
!
   end function new_cc_es_amplitudes_task
!
!
   subroutine execute_cc_es_amplitudes_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_es_amplitudes_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      call this%print_header()
      call this%start_timer()
!
      this%solver_factory = cc_es_amplitudes_solver_factory(wf%name_,            &
                                                            this%transformation, &
                                                            this%restart)
!
      call this%solver_factory%create(wf, this%solver)
!
      call wf%construct_fock(task='es')
!
      call this%solver%run()
!
      call this%solver%cleanup()
!
      call this%end_timer()
!
   end subroutine execute_cc_es_amplitudes_task
!
!
end module cc_es_amplitudes_task_class
