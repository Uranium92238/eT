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
module hf_gs_engine_class
!
!!
!! HF ground state engine class
!! Written by Sarai D. Folkestad, Eirik F. Kjønstad, and Alexander C. Paul, 2018-2022
!!
!
   use hf_engine_class, only: hf_engine
   use hf_class,        only: hf
!
   use scf_task_class, only: scf_task
   use sad_generation_task_class, only: sad_generation_task
   use hf_mean_value_task_class, only: hf_mean_value_task
   use hf_visualization_task_class, only: hf_visualization_task
!
   implicit none
!
   type, extends(hf_engine) :: hf_gs_engine
!
      type(sad_generation_task), allocatable, private :: sad_generation
      type(scf_task), allocatable, private :: scf
      type(hf_mean_value_task), allocatable, private :: mean_value
      type(hf_visualization_task), allocatable, private :: visualization
!
   contains
!
      procedure, public :: ignite => ignite_hf_gs_engine
!
      procedure, private, nopass :: do_sad
!
   end type hf_gs_engine
!
contains
!
!
   subroutine ignite_hf_gs_engine(this, wf)
!!
!!    Ignite
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(hf_gs_engine), intent(inout) :: this
      class(hf),           intent(inout) :: wf
!
      if (this%do_sad(wf)) then
!
         this%sad_generation = sad_generation_task()
         call this%sad_generation%execute(wf)
!
      end if
!
      this%scf = scf_task()
      call this%scf%execute(wf)
!
      this%mean_value = hf_mean_value_task()
      call this%mean_value%execute(wf)
!
!     If we don't prepare here, we cannot plot active HF densities
      call wf%prepare_for_post_HF_method()
!
      this%visualization = hf_visualization_task()
      call this%visualization%execute(wf)
!
   end subroutine ignite_hf_gs_engine
!
!
   function do_sad(wf) result(do_)
!!
!!    Do SAD
!!    Written by Alexander C. Paul, May 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(hf), intent(inout) :: wf
      logical :: do_
      character(len=200) :: start_guess
!
      do_ = .true.
!
      call input%get_keyword('ao density guess', 'solver scf', start_guess)
!
      if (start_guess == "core") then
!
         do_ = .false.
         return
!
      end if
!
      if (input%is_keyword_present('restart', 'solver scf') .or. &
          input%is_keyword_present('restart', 'do') .or.         &
          input%is_keyword_present('skip', 'solver scf')) then
!
         do_ = .not. wf%is_restart_possible()
         return
!
      end if
!
   end function do_sad
!
!
   end module hf_gs_engine_class

