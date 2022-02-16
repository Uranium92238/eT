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
module cc_multipliers_solver_factory_class_
!
!!
!! CC multiplier solver factory class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class, only: ccs
   use global_in, only: input
!
   implicit none
!
   type :: cc_multipliers_solver_factory_
!
      character(len=200), private :: algorithm
!
      logical, private :: restart
!
   contains
!
      procedure, public :: create
!
   end type cc_multipliers_solver_factory_
!
!
   interface cc_multipliers_solver_factory_
!
      procedure :: new_cc_multipliers_solver_factory
!
   end interface cc_multipliers_solver_factory_
!
!
contains
!
!
   function new_cc_multipliers_solver_factory(method) result(this)
!!
!!    New CC multiplier solver factory
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(cc_multipliers_solver_factory_) :: this
!
      character(len=*), intent(in) :: method
!
      if (trim(method) .eq. 'ccs'            .or. &
          trim(method) .eq. 'cc2'            .or. &
          trim(method) .eq. 'cc3'            .or. &
          trim(method) .eq. 'low memory cc2' .or. &
          trim(method) .eq. 'mlcc2') then
!
         this%algorithm = 'diis'
!
      else
!
         this%algorithm = 'davidson'
!
      end if
!
      call input%get_keyword('algorithm', 'solver cc multipliers', this%algorithm)
!
      this%restart = input%is_keyword_present('restart', 'solver cc multipliers')
!
      if (input%is_keyword_present('restart', 'do')) then
!
         this%restart = .true.
!
      end if
!
   end function new_cc_multipliers_solver_factory
!
!
   subroutine create(this, wf, solver)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, 2021
!!
      use global_out, only: output
      use abstract_solver_class,         only: abstract_solver
      use diis_cc_multipliers_class,      only: diis_cc_multipliers
      use cc_multipliers_solver_factory_class, only: cc_multipliers_solver_factory
      use amplitude_updater_class, only: amplitude_updater
      use quasi_newton_updater_class, only: quasi_newton_updater
!
      implicit none
!
      class(cc_multipliers_solver_factory_), intent(in) :: this
!
      class(ccs), intent(in) :: wf
!
      class(abstract_solver), allocatable :: solver
!
      type(cc_multipliers_solver_factory) :: davidson_factory
!
      class(amplitude_updater), allocatable :: t_updater
!
      if (this%algorithm == 'diis') then
!
         t_updater = quasi_newton_updater(n_amplitudes     = wf%n_gs_amplitudes, &
                                          scale_amplitudes = .false.,   &
                                          scale_residual   = .false.)
!
         solver = diis_cc_multipliers(wf, this%restart, t_updater)
!
      else if (this%algorithm == 'davidson') then
!
         call davidson_factory%create(wf, solver)
!
      else
!
         call output%error_msg('Did not recognize multipliers algorithm: (a0)', &
                               chars=[this%algorithm])
!
      end if
!
   end subroutine create
!
!
end module cc_multipliers_solver_factory_class_
