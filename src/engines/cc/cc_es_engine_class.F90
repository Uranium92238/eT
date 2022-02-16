!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module cc_es_engine_class
!
!!
!! CC excited state engine class
!! Written by Eirik F. Kjønstad, Alexander Paul, and Sarai D. Folkestad, 2021-2022
!!
!
   use global_in,                   only: input
   use ccs_class,                   only: ccs
   use cc_engine_class,             only: cc_engine
   use cc_amplitudes_task_class,    only: cc_amplitudes_task
   use cc_es_amplitudes_task_class, only: cc_es_amplitudes_task
   use cc_visualization_task_class, only: cc_visualization_task
   use eri_approximator_task_class, only: eri_approximator_task
   use cc_wavefunctions_class,      only: cc_wavefunctions
!
   implicit none
!
   type, extends(cc_engine) :: cc_es_engine
!
      type(eri_approximator_task),              private :: eri_approximator
      type(cc_amplitudes_task), allocatable,    private :: ground_state_amplitudes
      type(cc_es_amplitudes_task), allocatable, private :: excited_state_amplitudes
      type(cc_visualization_task), allocatable, private :: visualization
!
   contains
!
      procedure, public :: ignite => ignite_cc_es_engine
      procedure, public :: set_allowed_wfs => set_allowed_wfs_cc_es_engine
!
      procedure, private, nopass :: do_right
      procedure, private, nopass :: do_left
!
   end type cc_es_engine
!
contains
!
!
   subroutine ignite_cc_es_engine(this, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2021-2022
!!
      implicit none
!
      class(cc_es_engine), intent(inout) :: this
      class(ccs),          intent(inout) :: wf
!
      logical :: restart
!
      call this%set_allowed_wfs()
      call this%check_wavefunctions(wf)
!
      call this%eri_approximator%execute(wf)
!
      this%ground_state_amplitudes = cc_amplitudes_task()
      call this%ground_state_amplitudes%execute(wf)
!
      restart = input%is_keyword_present('restart', 'solver cc es') .or. &
                input%is_keyword_present('restart', 'do')
!
      if (this%do_right()) then
!
         this%excited_state_amplitudes = cc_es_amplitudes_task(transformation='right', &
                                                               restart=restart)
         call this%excited_state_amplitudes%execute(wf)
!
      endif
!
      if (this%do_left()) then
!
         if (this%do_right()) restart = .true.
         this%excited_state_amplitudes = cc_es_amplitudes_task(transformation='left', &
                                                               restart=restart)
         call this%excited_state_amplitudes%execute(wf)
!
      endif
!
      this%visualization = cc_visualization_task()
      call this%visualization%execute(wf)
!
   end subroutine ignite_cc_es_engine
!
!
   function do_right() result(do_)
!!
!!    Do right
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      logical :: do_
!
      do_ = .false.
!
      if (input%is_keyword_present('right eigenvectors', 'solver cc es')) then
         do_ = .true.
      else
         if (.not. input%is_keyword_present('left eigenvectors', 'solver cc es')) do_ = .true.
      endif
!
   end function do_right
!
!
   function do_left() result(do_)
!!
!!    Do left
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      logical :: do_
!
      do_ = input%is_keyword_present('left eigenvectors', 'solver cc es')
!
   end function do_left
!
!
   subroutine set_allowed_wfs_cc_es_engine(this)
!!
!!    Set allowed wavefunctions
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(cc_es_engine), intent(inout) :: this
!
      this%allowed_cc_wfs = cc_wavefunctions()
!
      call this%allowed_cc_wfs%set('ccs',            allowed=.true.)
      call this%allowed_cc_wfs%set('cc2',            allowed=.true.)
      call this%allowed_cc_wfs%set('low memory cc2', allowed=.true.)
      call this%allowed_cc_wfs%set('ccsd',           allowed=.true.)
      call this%allowed_cc_wfs%set('cc3',            allowed=.true.)
      call this%allowed_cc_wfs%set('mlcc2',          allowed=.true.)
      call this%allowed_cc_wfs%set('mlccsd',         allowed=.true.)
!
   end subroutine set_allowed_wfs_cc_es_engine
!
end module cc_es_engine_class
