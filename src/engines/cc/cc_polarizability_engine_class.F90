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
module cc_polarizability_engine_class
!
!!
!! Coupled cluster polarizability engine class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2021
!!
!
   use ccs_class,                      only: ccs
   use cc_engine_class,                only: cc_engine
   use cc_amplitudes_task_class,       only: cc_amplitudes_task
   use cc_multipliers_task_class,      only: cc_multipliers_task
   use cc_polarizability_task_class,   only: cc_polarizability_task
   use eri_approximator_task_class,    only: eri_approximator_task
   use cc_wavefunctions_class,         only: cc_wavefunctions
!
   implicit none
!
   type, extends(cc_engine) :: cc_polarizability_engine
!
      type(eri_approximator_task),               private :: eri_approximator
      type(cc_amplitudes_task), allocatable,     private :: t_task
      type(cc_multipliers_task), allocatable,    private :: tbar_task
      type(cc_polarizability_task), allocatable, private :: alpha_task
!
   contains
!
      procedure, public :: ignite => ignite_cc_polarizability_engine
      procedure, public :: set_allowed_wfs => set_allowed_wfs_cc_polarizability_engine
!
      final :: destructor
!
   end type cc_polarizability_engine
!
contains
!
!
   subroutine ignite_cc_polarizability_engine(this, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_polarizability_engine), intent(inout) :: this
!
      class(ccs), intent(inout) :: wf
!
      call this%set_allowed_wfs()
      call this%check_wavefunctions(wf)
!
      call this%eri_approximator%execute(wf)
!
      this%t_task = cc_amplitudes_task()
      call this%t_task%execute(wf)
!
      this%tbar_task = cc_multipliers_task()
      call this%tbar_task%execute(wf)
!
      this%alpha_task = cc_polarizability_task()
      call this%alpha_task%execute(wf)
!
   end subroutine ignite_cc_polarizability_engine
!
!
   subroutine destructor(this)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(cc_polarizability_engine) :: this
!
      if (allocated(this%alpha_task)) deallocate(this%alpha_task)
!
   end subroutine destructor
!
!
   subroutine set_allowed_wfs_cc_polarizability_engine(this)
!!
!!    Set allowed wavefunctions
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(cc_polarizability_engine), intent(inout) :: this
!
      this%allowed_cc_wfs = cc_wavefunctions()
!
      call this%allowed_cc_wfs%set('ccs',  allowed=.true.)
      call this%allowed_cc_wfs%set('cc2',  allowed=.true.)
      call this%allowed_cc_wfs%set('ccsd',  allowed=.true.)
!
   end subroutine set_allowed_wfs_cc_polarizability_engine
!
!
end module cc_polarizability_engine_class
