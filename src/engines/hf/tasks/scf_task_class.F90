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
module scf_task_class
!
!!
!! SCF task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use hf_class,      only: hf
   use hf_task_class, only: hf_task
!
   implicit none
!
   type, extends(hf_task) :: scf_task
!
      logical :: restart, skip_scf, write_mo_info, write_molden_file
!
   contains
!
      procedure, public :: execute &
                        => execute_scf_task
!
      procedure, private :: read_settings
!
   end type scf_task
!
!
   interface scf_task
!
      procedure :: new_scf_task
!
   end interface scf_task
!
!
contains
!
!
   function new_scf_task() result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(scf_task) :: this
!
      this%name_ = 'Determining reference state'
!
   end function new_scf_task
!
!
   subroutine execute_scf_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use scf_solver_factory_class, only: scf_solver_factory
      use scf_solver_class, only: scf_solver
!
      implicit none
!
      class(scf_task), intent(inout) :: this
      class(hf), target, intent(inout) :: wf
!
      type(scf_solver_factory), allocatable :: factory
      class(scf_solver), allocatable :: solver
!
      call this%print_header()
      call this%start_timer()
!
      call this%read_settings(wf)
!
      call wf%prepare_for_scf(this%restart, this%skip_scf)
!
      factory = scf_solver_factory()
      call factory%create(wf, solver, this%skip_scf)
!
      call solver%run(wf)
!
      call wf%finalize_gs()
!
      if (.not. this%skip_scf) call wf%flip_final_orbitals()
      call wf%print_summary()
!
      if (this%write_mo_info) call wf%write_orbital_info()
      if (this%write_molden_file) call wf%write_molden_file()
!
      call this%end_timer()
!
   end subroutine execute_scf_task
!
!
   subroutine read_settings(this, wf)
!!
!!    Read settings
!!    Written by Alexander C. Paul, 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(scf_task), intent(inout) :: this
      class(hf), intent(in) :: wf
!
      this%restart = .false.
!
      if (input%is_keyword_present('restart', 'solver scf') .or. &
          input%is_keyword_present('restart', 'do')) then
!
         this%restart = wf%is_restart_possible()
!
      end if
!
      this%skip_scf = input%is_keyword_present('skip', 'solver scf')
!
      this%write_mo_info = input%is_keyword_present('print orbitals', 'solver scf')
      if (input%is_keyword_present('print orbitals', 'system')) then
         this%write_mo_info = .true.
      end if
!
      this%write_molden_file = input%is_keyword_present('write molden', 'solver scf')
!
   end subroutine read_settings
!
!
end module scf_task_class
