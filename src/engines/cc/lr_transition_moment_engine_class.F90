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
module lr_transition_moment_engine_class
!
!!
!! Linear response transition moment engine class
!! Written by Alexander C. Paul, Dec 2021
!!
!
   use ccs_class,                            only: ccs
   use cc_engine_class,                      only: cc_engine
   use cc_amplitudes_task_class,             only: cc_amplitudes_task
   use cc_multipliers_task_class,            only: cc_multipliers_task
   use cc_es_amplitudes_task_class,          only: cc_es_amplitudes_task
   use biorthonormalization_task_class,      only: biorthonormalization_task
   use cc_lr_transition_moments_task_class,  only: cc_lr_transition_moments_task
   use eri_approximator_task_class,          only: eri_approximator_task
   use cc_wavefunctions_class,               only: cc_wavefunctions
!
   implicit none
!
   type, extends(cc_engine) :: lr_transition_moment_engine
!
      type(eri_approximator_task),                        private :: eri_approximator
      type(cc_amplitudes_task),              allocatable, private :: ground_state_amplitudes
      type(cc_multipliers_task),             allocatable, private :: ground_state_multipliers
      type(cc_es_amplitudes_task),           allocatable, private :: excited_state_amplitudes
      type(biorthonormalization_task),       allocatable, private :: biorthonormalization
      type(cc_lr_transition_moments_task),   allocatable, private :: lr_transition_moments
!
   contains
!
      procedure, public :: ignite => ignite_lr_transition_moment_engine
      procedure, public :: set_allowed_wfs => set_allowed_wfs_lr_transition_moment_engine
!
   end type lr_transition_moment_engine
!
!
contains
!
!
   subroutine ignite_lr_transition_moment_engine(this, wf)
!!
!!    Ignite
!!    Written by Alexander C. Paul, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(lr_transition_moment_engine), intent(inout) :: this
!
      class(ccs), intent(inout) :: wf
!
      logical :: restart
!
      call this%set_allowed_wfs()
      call this%check_wavefunctions(wf)
!
      restart = input%is_keyword_present('restart', 'solver cc es') .or. &
                input%is_keyword_present('restart', 'do')
!
      call this%eri_approximator%execute(wf)
!
      this%ground_state_amplitudes = cc_amplitudes_task()
      call this%ground_state_amplitudes%execute(wf)
!
      this%ground_state_multipliers = cc_multipliers_task()
      call this%ground_state_multipliers%execute(wf)
!
      call wf%save_tbar_intermediates()
!
      this%excited_state_amplitudes = cc_es_amplitudes_task(transformation='right', &
                                                            restart=restart)
!
      call this%excited_state_amplitudes%execute(wf)
!
      this%excited_state_amplitudes = cc_es_amplitudes_task(transformation='left', &
                                                            restart=.true.)
!
      call this%excited_state_amplitudes%execute(wf)
!
      this%biorthonormalization = biorthonormalization_task()
      call this%biorthonormalization%execute(wf)
!
      this%lr_transition_moments = cc_lr_transition_moments_task()
      call this%lr_transition_moments%execute(wf)
!
   end subroutine ignite_lr_transition_moment_engine
!
!
   subroutine set_allowed_wfs_lr_transition_moment_engine(this)
!!
!!    Set allowed wavefunctions
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(lr_transition_moment_engine), intent(inout) :: this
!
      this%allowed_cc_wfs = cc_wavefunctions()
!
      call this%allowed_cc_wfs%set('ccs',  allowed=.true.)
      call this%allowed_cc_wfs%set('cc2',  allowed=.true.)
      call this%allowed_cc_wfs%set('ccsd', allowed=.true.)
!
   end subroutine set_allowed_wfs_lr_transition_moment_engine
!
!
end module lr_transition_moment_engine_class
