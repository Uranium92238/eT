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
module cc_es_amplitudes_solver_factory_class
!
!!
!! CC excited state amplitudes solver factory class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class,    only: ccs
   use global_in,    only: input
   use global_out,   only: output
!
   implicit none
!
!
   type :: cc_es_amplitudes_solver_factory
!
      character(len=200), private :: transformation
!
      logical, private :: restart
!
      character(len=200), private :: algorithm
!
   contains
!
      procedure, public :: create
!
   end type cc_es_amplitudes_solver_factory
!
!
   interface cc_es_amplitudes_solver_factory
!
      procedure :: new_cc_es_amplitudes_solver_factory
!
   end interface cc_es_amplitudes_solver_factory
!
!
contains
!
!
   function new_cc_es_amplitudes_solver_factory(method, transformation, restart) result(this)
!!
!!    New CC excited state amplitudes solver factory
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(cc_es_amplitudes_solver_factory) :: this
!
      character(len=*), intent(in) :: transformation
!
      character(len=*), intent(in) :: method
!
      logical, intent(in) :: restart
!
      this%transformation = trim(transformation)
      this%restart = restart
!
      if (method .eq. 'cc3' .or. &
          method .eq. 'low memory cc2') then
!
         this%algorithm = 'non-linear davidson'
!
      else
!
         this%algorithm = 'davidson'
!
      end if
!
      call input%get_keyword('algorithm', 'solver cc es', this%algorithm)
!
   end function new_cc_es_amplitudes_solver_factory
!
!
   subroutine create(this, wf, solver)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, 2021
!!
      use abstract_solver_class, only: abstract_solver
!
      use asymmetric_lanczos_cc_es_class, only: asymmetric_lanczos_cc_es
      use diis_cc_es_class,               only: diis_cc_es
      use nonlinear_davidson_cc_es_class, only: nonlinear_davidson_cc_es
!
      use cc_multipliers_task_class,           only: cc_multipliers_task
      use davidson_cc_es_solver_factory_class, only: davidson_cc_es_solver_factory
!
      implicit none
!
      class(cc_es_amplitudes_solver_factory), intent(in) :: this
!
      class(ccs), intent(inout) :: wf
!
      class(abstract_solver), allocatable, intent(out) :: solver
!
      class(cc_multipliers_task), allocatable :: multipliers
      type(davidson_cc_es_solver_factory), allocatable :: davidson_factory
!
      if (this%algorithm == 'asymmetric lanczos') then
!
         multipliers = cc_multipliers_task()
         call multipliers%execute(wf)
!
         solver = asymmetric_lanczos_cc_es(wf)
!
      else
!
         if (this%algorithm == 'diis') then
!
            solver = diis_cc_es(this%transformation, wf, this%restart)
!
         elseif (this%algorithm == 'davidson') then
!
            if (trim(wf%name_) == 'low memory cc2' .or. &
                trim(wf%name_) == 'cc3') then
!
               call output%error_msg('Davidson not implemented for CC3 and lowmem CC2.')
!
            end if
!
            davidson_factory = davidson_cc_es_solver_factory(this%transformation, this%restart)
            call davidson_factory%create(wf, solver)
!
         elseif (this%algorithm == 'non-linear davidson') then
!
            solver = nonlinear_davidson_cc_es(this%transformation, wf, this%restart)
!
         else
!
            call output%error_msg('Could not create excited state solver. It may be that the &
                                    &algorithm is not implemented for the method specified.')
         endif
!
         call wf%initialize_excited_state_files()
         call wf%initialize_excitation_energies()
!
      endif
!
   end subroutine create
!
!
end module cc_es_amplitudes_solver_factory_class
