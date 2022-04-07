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
module cc_gs_mean_value_engine_class
!
!!
!! CC ground state engine class
!! Written by Tor S. Haugland, Eirik F. Kjønstad,
!! Alexander Paul, and Sarai D. Folkestad, 2019-2022
!!
!
   use ccs_class,                   only: ccs
   use cc_engine_class,             only: cc_engine
   use cc_amplitudes_task_class,    only: cc_amplitudes_task
   use cc_multipliers_task_class,   only: cc_multipliers_task
   use cc_mean_values_task_class,   only: cc_mean_values_task
   use cc_visualization_task_class, only: cc_visualization_task
   use eri_approximator_task_class, only: eri_approximator_task
   use cc_wavefunctions_class,      only: cc_wavefunctions
!
   implicit none
!
   type, extends(cc_engine) :: cc_gs_mean_value_engine
!
      type(cc_amplitudes_task), allocatable,    private :: ground_state_amplitudes
      type(cc_multipliers_task), allocatable,   private :: ground_state_multipliers
      type(cc_mean_values_task),   allocatable, private :: mean_values
      type(cc_visualization_task), allocatable, private :: visualization
      type(eri_approximator_task),              private :: eri_approximator
!
   contains
!
      procedure, public :: ignite => ignite_cc_gs_mean_value_engine
      procedure, public :: set_allowed_wfs => set_allowed_wfs_cc_gs_mean_value_engine

!
   end type cc_gs_mean_value_engine
!
contains
!
!
   subroutine ignite_cc_gs_mean_value_engine(this, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_gs_mean_value_engine), intent(inout) :: this
      class(ccs),                     intent(inout) :: wf
!
      call this%set_allowed_wfs()
      call this%check_wavefunctions(wf)
!
      call this%eri_approximator%execute(wf)
!
      this%ground_state_amplitudes = cc_amplitudes_task()
      call this%ground_state_amplitudes%execute(wf)
!
      this%ground_state_multipliers = cc_multipliers_task()
      call this%ground_state_multipliers%execute(wf)
!
      this%mean_values = cc_mean_values_task()
      call this%mean_values%execute(wf)
!
      this%visualization = cc_visualization_task()
      call this%visualization%execute(wf)
!
   end subroutine ignite_cc_gs_mean_value_engine
!
!
   subroutine set_allowed_wfs_cc_gs_mean_value_engine(this)
!!
!!    Set allowed wavefunctions
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(cc_gs_mean_value_engine), intent(inout) :: this
!
      this%allowed_cc_wfs = cc_wavefunctions()
!
      call this%allowed_cc_wfs%set('ccs',  allowed=.true.)
      call this%allowed_cc_wfs%set('cc2',  allowed=.true.)
      call this%allowed_cc_wfs%set('ccsd', allowed=.true.)
      call this%allowed_cc_wfs%set('cc3',  allowed=.true.)
!
   end subroutine set_allowed_wfs_cc_gs_mean_value_engine
!
end module cc_gs_mean_value_engine_class
