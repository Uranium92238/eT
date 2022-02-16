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
module cc_propagation_task_class
!
!!
!! CC propagation task class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class,                   only: ccs
   use cc_propagation_class,        only: cc_propagation
   use electric_field_class,        only: electric_field
   use cc_task_class,               only: cc_task
   use cc_propagator_factory_class, only: cc_propagator_factory
!
   implicit none
!
   type, extends(cc_task) :: cc_propagation_task
!
      class(cc_propagation), allocatable  :: propagator
      type(electric_field), allocatable   :: field
!
      type(cc_propagator_factory), allocatable :: propagator_factory
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_propagation_task
!
   end type cc_propagation_task
!
!
   interface cc_propagation_task
!
      procedure :: new_cc_propagation_task
!
   end interface cc_propagation_task
!
!
contains
!
!
   function new_cc_propagation_task() result(this)
!!
!!    New
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      type(cc_propagation_task) :: this
!
      this%name_ = 'Propagating CC amplitudes'
!
   end function new_cc_propagation_task
!
!
   subroutine execute_cc_propagation_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_propagation_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      call this%print_header()
      call this%start_timer()
!
      this%propagator_factory = cc_propagator_factory()
      call this%propagator_factory%create(wf, this%propagator)
!
      this%field = electric_field()
!
      call this%propagator%initializations()
      call this%propagator%run(wf, this%field)
      call this%propagator%cleanup(wf)
!
      call this%field%cleanup()
!
      call this%end_timer()
!
   end subroutine execute_cc_propagation_task
!
!
end module cc_propagation_task_class
