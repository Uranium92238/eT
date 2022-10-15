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
   use kinds
   use ccs_class,  only: ccs
   use global_in,  only: input
   use global_out, only: output
!
   use abstract_solver_class, only: abstract_solver
   use convergence_tool_class, only: convergence_tool
   use start_vector_tool_class, only: start_vector_tool
   use abstract_projection_tool_class, only: abstract_projection_tool
!
   implicit none
!
!
   type :: cc_es_amplitudes_solver_factory
!
      class(convergence_tool), allocatable :: convergence_checker
      class(start_vector_tool), allocatable :: start_vectors
      class(abstract_projection_tool), allocatable :: projector
!
      character(len=200), private :: algorithm, transformation, es_type
      integer, private :: n_states, max_iterations
      logical, private :: energy_convergence, restart, records_in_memory
      real(dp), private :: energy_threshold, residual_threshold
!
   contains
!
      procedure, public :: create
!
      procedure, private :: read_settings
      procedure, private :: determine_es_type
      procedure, private :: create_common_es_solver_components
      procedure, private :: create_davidson_es_solver
!
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
      use diis_cc_es_class,               only: diis_cc_es
      use nonlinear_davidson_cc_es_class, only: nonlinear_davidson_cc_es
!
      implicit none
!
      class(cc_es_amplitudes_solver_factory), intent(inout) :: this
!
      class(ccs), intent(inout) :: wf
!
      class(abstract_solver), allocatable, intent(out) :: solver
!
      call this%read_settings()
      call wf%prepare_for_excited_states(this%n_states, this%es_type)
      call this%create_common_es_solver_components(wf)
!
      if (this%algorithm == 'diis') then
!
         solver = diis_cc_es(this%transformation, wf, &
                             this%start_vectors, &
                             this%projector, &
                             this%convergence_checker, &
                             this%n_states, &
                             this%max_iterations, &
                             this%records_in_memory, &
                             this%es_type)
!
      elseif (this%algorithm == 'davidson') then
!
         call this%create_davidson_es_solver(wf, solver)
!
      elseif (this%algorithm == 'non-linear davidson') then
!
         solver = nonlinear_davidson_cc_es(this%transformation, wf,  &
                                           this%start_vectors,       &
                                           this%projector,           &
                                           this%convergence_checker, &
                                           this%n_states,            &
                                           this%max_iterations,      &
                                           this%records_in_memory)
!
      else
!
         call output%error_msg('Could not create excited state solver. It may be that the &
                                 &algorithm is not implemented for the method specified.')
      endif
!
   end subroutine create
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad, May 2021
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_es_amplitudes_solver_factory), intent(inout) :: this
!
      this%n_states           = 0
      this%max_iterations     = 100
      this%energy_threshold   = 1.0d-3
      this%residual_threshold = 1.0d-3
      this%records_in_memory  = .false.
!
      call input%get_keyword('singlet states', 'solver cc es', this%n_states)
      if (this%n_states == 0) call output%error_msg('can not solve for 0 excited states')
!
      call input%get_keyword('max iterations', 'solver cc es', this%max_iterations)

      call input%get_keyword('residual threshold', 'solver cc es', this%residual_threshold)

      this%energy_convergence = input%is_keyword_present('energy threshold', 'solver cc es')
!
      if (this%energy_convergence) &
         call input%get_keyword('energy threshold', 'solver cc es', this%energy_threshold)
!
      call input%place_records_in_memory('solver cc es', this%records_in_memory)
!
      call this%determine_es_type()
!
   end subroutine read_settings
!
!
   subroutine determine_es_type(this)
!!
!!    Determine es type
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(cc_es_amplitudes_solver_factory), intent(inout)  :: this
!
      this%es_type = 'valence'
!
      if (input%is_keyword_present('core excitation', 'solver cc es') .and. .not. &
          input%is_keyword_present('ionization', 'solver cc es')) then
!
         this%es_type = 'core'
!
      end if
!
      if (input%is_keyword_present('ionization', 'solver cc es') .and. .not. &
          input%is_keyword_present('core excitation', 'solver cc es')) then
!
         this%es_type = 'ionize'
!
      end if
!
      if (input%is_keyword_present('remove core', 'solver cc es')) then
!
         this%es_type = 'remove core'
!
      end if
!
   end subroutine determine_es_type
!
!
   subroutine create_common_es_solver_components(this, wf)
!!
!!    Create common es solver components
!!    Written by Alexander C. Paul, May 2022
!!
      use cc_es_start_vector_factory_class, only: cc_es_start_vector_factory
      use cc_es_projector_factory_class, only: cc_es_projector_factory
!
      implicit none
!
      class(cc_es_amplitudes_solver_factory), intent(inout) :: this
      class(ccs), intent(in) :: wf
!
      class(cc_es_start_vector_factory), allocatable :: start_vector_factory
      class(cc_es_projector_factory), allocatable :: projector_factory
!
      start_vector_factory = cc_es_start_vector_factory(this%es_type,        &
                                                        this%transformation, &
                                                        this%restart)
      call start_vector_factory%create(wf, this%start_vectors)
!
      projector_factory = cc_es_projector_factory(this%es_type)
      call projector_factory%create(wf, this%projector)
!
      this%convergence_checker = convergence_tool(this%energy_threshold,   &
                                                  this%residual_threshold, &
                                                  this%energy_convergence)
!
   end subroutine create_common_es_solver_components
!
!
   subroutine create_davidson_es_solver(this, wf, solver)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Initializes the davidson solver with the correct instances of the following types:
!!
!!       - Solver tool (eigen_davidson)
!!       - Convergence tool
!!       - Transformation tool,
!!       - Storage tool
!!       - Start vector tool
!!       - Preconditioner tool
!!
!!    for the determination of CC excitation energies
!!
      use eigen_davidson_solver_class, only: eigen_davidson_solver
!
      use eigen_davidson_tool_class, only: eigen_davidson_tool
      use cc_eigen_storage_tool_class, only: cc_eigen_storage_tool
      use cc_jacobian_preconditioner_getter_class, only: cc_jacobian_preconditioner_getter
      use cc_jacobian_transformation_class, only: cc_jacobian_transformation
!
      implicit none
!
      class(cc_es_amplitudes_solver_factory), intent(inout) :: this
      class(ccs),                             intent(inout) :: wf
      class(abstract_solver), allocatable,    intent(out)   :: solver
!
      class(eigen_davidson_tool),               allocatable :: davidson
      class(cc_eigen_storage_tool),             allocatable :: storer
      class(cc_jacobian_transformation),        allocatable :: transformer
      class(cc_jacobian_preconditioner_getter), allocatable :: preconditioner
!
      real(dp) :: lindep_threshold
      integer :: max_dim_red
!
      if (trim(wf%name_) == 'low memory cc2' .or. trim(wf%name_) == 'cc3') then
!
         call output%error_msg('Davidson not implemented for (a0).', &
                               chars=[trim(wf%name_)])
!
      end if
!
      max_dim_red = max(100, 10 * this%n_states)
      call input%get_keyword('max reduced dimension', 'solver cc es', max_dim_red)
!
      preconditioner = cc_jacobian_preconditioner_getter(wf, wf%n_es_amplitudes)
      storer = cc_eigen_storage_tool(wf, trim(this%transformation))
!
      transformer = cc_jacobian_transformation(wf, trim(this%transformation), wf%n_es_amplitudes)
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
!
      davidson = eigen_davidson_tool(name_             = 'cc_es_davidson',    &
                                     n_parameters      = wf%n_es_amplitudes,  &
                                     n_solutions       = this%n_states,       &
                                     lindep_threshold  = lindep_threshold,    &
                                     max_dim_red       = max_dim_red,         &
                                     records_in_memory = this%records_in_memory)
!
      solver = eigen_davidson_solver(transformer         = transformer,              &
                                     davidson            = davidson,                 &
                                     convergence_checker = this%convergence_checker, &
                                     storer              = storer,                   &
                                     start_vector        = this%start_vectors,       &
                                     preconditioner      = preconditioner,           &
                                     projector           = this%projector,           &
                                     n_solutions         = this%n_states,            &
                                     max_iterations      = this%max_iterations)
!
   end subroutine create_davidson_es_solver
!
!
end module cc_es_amplitudes_solver_factory_class
