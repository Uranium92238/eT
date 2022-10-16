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
      logical, dimension(2) :: do_side
      logical, dimension(2) :: do_spin
!
   contains
!
      procedure, public :: ignite => ignite_cc_es_engine
      procedure, public :: set_allowed_wfs => set_allowed_wfs_cc_es_engine
!
      procedure, private :: determine_side
      procedure, private :: determine_spin
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
      integer :: side, spin
!
      character(len=10), dimension(2), parameter :: side_string = [character(len=10) :: 'right', 'left']
      character(len=10), dimension(2), parameter :: spin_string = [character(len=10) :: 'singlet', 'triplet']
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
      call this%determine_side()
      call this%determine_spin()
!
      do side = 1, 2
!
         if (.not. this%do_side(side)) cycle
!
         if (side == 2 .and. this%do_side(1)) restart = .true.
!
         do spin = 1, 2
!
            if (.not. this%do_spin(spin)) cycle
!
            this%excited_state_amplitudes = cc_es_amplitudes_task(trim(side_string(side)), &
                                                                  trim(spin_string(spin)), &
                                                                  restart)
            call this%excited_state_amplitudes%execute(wf)
!
         enddo
      enddo
!
      this%visualization = cc_visualization_task()
      call this%visualization%execute(wf)
!
   end subroutine ignite_cc_es_engine
!
!
   subroutine determine_side(this)
!!
!!    Determine side
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Determines the transformation
!!
!!    this%do_side(k)
!!
!!    k = 1 -> right
!!    k = 2 -> left
!!
      implicit none
!
      class(cc_es_engine), intent(inout) :: this
!
!     Default is only right
      this%do_side(1) = .true.
      this%do_side(2) = .false.
!
      if (input%is_keyword_present('left eigenvectors', 'solver cc es')) then
         this%do_side(1) = .false.
         this%do_side(2) = .true.
      endif
!
      if (input%is_keyword_present('right eigenvectors', 'solver cc es')) then
         this%do_side(1) = .true.
      endif
!
   end subroutine determine_side
!
!
   subroutine determine_spin(this)
!!
!!    Determine spin
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Determines the spin symmetry to run
!!
!!    this%do_spin(k)
!!
!!    k = 1 -> singlet
!!    k = 2 -> triplet
!!
      implicit none
!
      class(cc_es_engine), intent(inout) :: this
!
!     Default is only singlet - This is needed as long as Lanczos is treated in the es engine
      this%do_spin(1) = .true.
      this%do_spin(2) = .false.
!
      if (input%is_keyword_present('triplet states', 'solver cc es')) then
         this%do_spin(2) = .true.
         this%do_spin(1) = .false.
      endif
!
      if (input%is_keyword_present('singlet states', 'solver cc es')) then
         this%do_spin(1) = .true.
      endif
!
   end subroutine determine_spin
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
