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
module fci_eigenproblem_task_class
!
!!
!! FCI eigenproblem task class
!! Written by Alexander C. Paul, May 2022
!!
!
   use fci_class,      only: fci
   use fci_task_class, only: fci_task
!
   implicit none
!
   type, extends(fci_task) :: fci_eigenproblem_task
!
      logical :: restart
!
   contains
!
      procedure, public :: execute &
                        => execute_fci_eigenproblem_task
!
      procedure, private :: read_settings
!
   end type fci_eigenproblem_task
!
!
   interface fci_eigenproblem_task
!
      procedure :: new_fci_eigenproblem_task
!
   end interface fci_eigenproblem_task
!
!
contains
!
!
   function new_fci_eigenproblem_task() result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(fci_eigenproblem_task) :: this
!
      this%name_ = 'Determining FCI eigenstates'
!
   end function new_fci_eigenproblem_task
!
!
   subroutine execute_fci_eigenproblem_task(this, wf)
!!
!!    Execute
!!    Written by Enrico Ronca and Alexander C. Paul, 2021-2022
!!
      use eigen_davidson_solver_class, only: eigen_davidson_solver
      use fci_solver_factory_class,   only: fci_solver_factory
!
      implicit none
!
      class(fci_eigenproblem_task), intent(inout) :: this
      class(fci), target, intent(inout) :: wf
!
      class(eigen_davidson_solver), allocatable :: solver
      type(fci_solver_factory) :: factory
!
      call this%print_header()
      call this%start_timer()
!
      call this%read_settings()
!
      call factory%create(wf, solver)
      call solver%run()
      call wf%print_fci_summary()
!
      call this%end_timer()
!
   end subroutine execute_fci_eigenproblem_task
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Alexander C. Paul, 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(fci_eigenproblem_task), intent(inout) :: this
!
      this%restart = (input%is_keyword_present('restart', 'solver fci') .or. &
                      input%is_keyword_present('restart', 'do'))
!
   end subroutine read_settings
!
!
end module fci_eigenproblem_task_class
   