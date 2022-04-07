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
module cc_amplitudes_solver_factory_class
!
!!
!! CC amplitudes solver factory class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use kinds
   use ccs_class,                            only: ccs
   use global_in,                            only: input
   use abstract_jacobian_transformer_class,  only: abstract_jacobian_transformer
   use global_in,                            only: input
   use global_out,                           only: output
   use amplitude_updater_class,              only: amplitude_updater
!
   implicit none
!
!
   type :: newton_raphson_settings
!
      real(dp), public           :: relative_threshold
      integer, public            :: max_iterations
      character(len=200), public :: storage
      logical, public            :: records_in_memory
!
   end type newton_raphson_settings
!
!
   type :: cc_amplitudes_solver_factory
!
      character(len=200), private :: algorithm
      character(len=200), private :: multimodel_newton
      logical, private :: restart
      character(len=200), private :: storage
!
      class(amplitude_updater), allocatable, private :: t_updater
      class(abstract_jacobian_transformer), allocatable, private :: transformer
!
      type(newton_raphson_settings), private :: microiterations
!
   contains
!
      procedure, public :: create
!
      procedure, private :: read_settings
!
      procedure, private :: create_amplitude_updater
!
      procedure, private :: create_quasi_newton_updater
      procedure, private :: create_newton_raphson_updater
!
      procedure, private :: create_jacobian_transformer
!
      procedure, private :: determine_microiteration_settings
!
   end type cc_amplitudes_solver_factory
!
!
contains
!
!
   subroutine create(this, wf, solver)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, 2021
!!
      use abstract_solver_class,    only: abstract_solver
      use diis_cc_gs_class,   only: diis_cc_gs
!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
!
      class(ccs), intent(in) :: wf
!
      class(abstract_solver), allocatable :: solver
!
      call this%read_settings(wf)
!
      call this%create_amplitude_updater(wf)
!
      solver = diis_cc_gs(wf, this%restart, this%t_updater, this%storage)
!
   end subroutine create
!
!
   subroutine create_amplitude_updater(this, wf)
!!
!!    Create amplitude updater
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
!
      class(ccs), intent(in) :: wf
!
      if (trim(this%algorithm) == 'diis') then
!
         call this%create_quasi_newton_updater(wf)
!
      elseif (trim(this%algorithm) == 'newton-raphson') then
!
         call this%create_newton_raphson_updater(wf)
!
      else
!
         call output%error_msg('Could not recognize algorithm (a0) in &
                              &cc_amplitudes_solver_factory.F90.', &
                              chars=[trim(this%algorithm)])
!
      endif
!
   end subroutine create_amplitude_updater
!
!
   subroutine create_quasi_newton_updater(this, wf)
!!
!!    Create Quasi-Newton updater
!!    Written by Eirik F. Kjønstad, 2021
!!
      use quasi_newton_updater_class, only: quasi_newton_updater
!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
!
      class(ccs), intent(in) :: wf
!
      this%t_updater = quasi_newton_updater(n_amplitudes     = wf%n_gs_amplitudes, &
                                            scale_amplitudes = .true.,             &
                                            scale_residual   = .true.)
!
   end subroutine create_quasi_newton_updater
!
!
   subroutine create_newton_raphson_updater(this, wf)
!!
!!    Create Newton-Raphson updater
!!    Written by Eirik F. Kjønstad, 2021
!!
      use newton_raphson_updater_class, only: newton_raphson_updater
!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
!
      class(ccs), intent(in) :: wf
!
      if (trim(wf%name_) == 'cc2') &
         call output%error_msg('Newton-Raphson not implemented for CC2')
!
      call this%determine_microiteration_settings()
!
      call this%create_jacobian_transformer()
!
      this%t_updater = &
         newton_raphson_updater(n_amplitudes       = wf%n_gs_amplitudes,                        &
                                scale_amplitudes   = .true.,                                    &
                                scale_residual     = .false.,                                   &
                                relative_threshold = this%microiterations%relative_threshold,   &
                                records_in_memory  = this%microiterations%records_in_memory,    &
                                max_iterations     = this%microiterations%max_iterations,       &
                                transformer        = this%transformer)
!
   end subroutine create_newton_raphson_updater
!
!
   subroutine determine_microiteration_settings(this)
!!
!!    Determine microiteration settings
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
!
      this%microiterations%relative_threshold = 1.0d-2
      this%microiterations%max_iterations     = 100
      this%microiterations%storage            = 'disk'
!
      call input%get_keyword('rel micro threshold',      &
                             'solver cc gs',             &
                             this%microiterations%relative_threshold)
!
      call input%get_keyword('max micro iterations',     &
                             'solver cc gs',             &
                             this%microiterations%max_iterations)
!
      call input%get_keyword('micro iteration storage',  &
                             'solver cc gs',             &
                             this%microiterations%storage)
!
      this%microiterations%records_in_memory = .false.
!
      if (trim(this%storage) == 'memory') then ! If macroiterations use memory storage, then use
                                               ! this also for the microiterations
!
         this%microiterations%records_in_memory = .true.
!
      end if
 !
      if (trim(this%microiterations%storage) == 'memory') then
!
         this%microiterations%records_in_memory = .true.
!
      end if
!
   end subroutine determine_microiteration_settings
!
!
   subroutine create_jacobian_transformer(this)
!!
!!    Create Jacobian transformer
!!    Written by Eirik F. Kjønstad, 2021
!!
      use citation_printer_class,                 only: eT_citations
      use jacobian_transformer_class,             only: jacobian_transformer
      use approximate_jacobian_transformer_class, only: approximate_jacobian_transformer
!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
!
      if (trim(this%multimodel_newton) == 'on') then
!
         this%transformer = approximate_jacobian_transformer('right')
!
         call eT_citations%add('Multimodel Newton algorithm')
!
      else if (trim(this%multimodel_newton) == 'off') then
!
         this%transformer = jacobian_transformer('right')
!
      else
!
         call output%error_msg('did not recognize transformer type in cc_amplitudes_solver_factory')
!
      endif
!
   end subroutine create_jacobian_transformer
!
!
   subroutine read_settings(this, wf)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_amplitudes_solver_factory), intent(inout) :: this
      class(ccs),                          intent(in) :: wf
!
      this%algorithm = 'diis'
      if (trim(wf%name_) == 'cc3') this%algorithm = 'newton-raphson'
!
      this%multimodel_newton = 'off'
      if (trim(wf%name_) == 'cc3') this%multimodel_newton = 'on'

      call input%get_keyword('multimodel newton', 'solver cc gs', this%multimodel_newton)
      if (trim(this%multimodel_newton) == 'on') this%algorithm = 'newton-raphson'
!
      call input%get_keyword('algorithm', 'solver cc gs', this%algorithm)
!
      this%storage   = 'disk'
      this%restart = input%is_keyword_present('restart', 'solver cc gs')
!
      if (input%is_keyword_present('restart', 'do')) then
!
         this%restart = .true.
!
      end if
!
      call input%get_keyword('storage', 'solver cc gs', this%storage)
!
   end subroutine read_settings
!
!
end module cc_amplitudes_solver_factory_class
