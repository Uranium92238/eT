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
module hf_mean_value_task_class
!
!!
!! HF Mean value task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use hf_class,      only: hf
   use hf_task_class, only: hf_task
!
   implicit none
!
   type, extends(hf_task) :: hf_mean_value_task
!
      logical :: dipole, quadrupole, task_requested
!
   contains
!
      procedure, public :: execute &
                        => execute_hf_mean_value_task
!
      procedure, private :: read_settings
!
   end type hf_mean_value_task
!
!
   interface hf_mean_value_task
!
      procedure :: new_hf_mean_value_task
!
   end interface hf_mean_value_task
!
!
contains
!
!
   function new_hf_mean_value_task() result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(hf_mean_value_task) :: this
!
      this%name_ = 'Determining HF properties'
!
   end function new_hf_mean_value_task
!
!
   subroutine execute_hf_mean_value_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use scf_solver_factory_class, only: scf_solver_factory
      use scf_solver_class, only: scf_solver
!
      implicit none
!
      class(hf_mean_value_task), intent(inout) :: this
      class(hf), target, intent(inout) :: wf
!
      call this%read_settings()
      if (.not. this%task_requested) return
!
      call this%print_header()
      call this%start_timer()
!
      if(this%dipole) call wf%calculate_and_print_dipole()
!
      if (this%quadrupole) call wf%calculate_and_print_quadrupole()
!
      call this%end_timer()
!
   end subroutine execute_hf_mean_value_task
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
      class(hf_mean_value_task), intent(inout) :: this
!
      this%task_requested = input%is_section_present('hf mean value')
      this%dipole = input%is_keyword_present('dipole','hf mean value')
      this%quadrupole = input%is_keyword_present('quadrupole','hf mean value')
!
   end subroutine read_settings
!
!
end module hf_mean_value_task_class
   