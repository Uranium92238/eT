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
module tdhf_task_class
!
!!
!! TDHF task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, May 2022
!!
!
   use hf_class,      only: hf
   use hf_task_class, only: hf_task
!
   implicit none
!
   type, extends(hf_task) :: tdhf_task
!
   contains
!
      procedure, public :: execute &
                        => execute_tdhf_task
!
   end type tdhf_task
!
!
   interface tdhf_task
!
      procedure :: new_tdhf_task
!
   end interface tdhf_task
!
!
contains
!
!
   function new_tdhf_task() result(this)
!!
!!    New
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      type(tdhf_task) :: this
!
      this%name_ = 'Determining TDHF excitation energies'
!
   end function new_tdhf_task
!
!
   subroutine execute_tdhf_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad, May 2021
!!
      use eigen_davidson_solver_class, only: eigen_davidson_solver
      use tdhf_solver_factory_class,   only: tdhf_solver_factory
!
      implicit none
!
      class(tdhf_task), intent(inout) :: this
      class(hf), target, intent(inout) :: wf
!
      class(eigen_davidson_solver), allocatable :: solver
      type(tdhf_solver_factory) :: solver_factory
!
      call this%print_header()
      call this%start_timer()
!
      call solver_factory%create(wf, solver)
!
      call wf%initialize_tdhf_quantities(solver%get_n_solutions())
!
      call solver%run()
      call wf%tdhf_summary()
!
      call this%end_timer()
!
   end subroutine execute_tdhf_task
!
!
end module tdhf_task_class
   