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
module cc_triplet_es_amplitudes_solver_factory_class
!
!!
!! CC triplet excited state amplitudes solver factory class
!! Written by Sarai D. Folkestad, 2022
!!
!
   use kinds
   use ccs_class,  only: ccs
   use global_in,  only: input
   use global_out, only: output
!
   use abstract_solver_class, only: abstract_solver
!
   implicit none
!
!
   type :: cc_triplet_es_amplitudes_solver_factory
!
      character(len=200), private :: transformation
      integer,            private :: n_states, max_iterations
      logical,            private :: energy_convergence, restart, records_in_memory
      real(dp),           private :: energy_threshold, residual_threshold
!
   contains
!
      procedure, public :: create
!
      procedure, private :: read_settings
      procedure, private :: create_davidson_es_solver
!
   end type cc_triplet_es_amplitudes_solver_factory
!
!
   interface cc_triplet_es_amplitudes_solver_factory
!
      procedure :: new_cc_triplet_es_amplitudes_solver_factory
!
   end interface cc_triplet_es_amplitudes_solver_factory
!
!
contains
!
!
   function new_cc_triplet_es_amplitudes_solver_factory(method, transformation, restart) result(this)
!!
!!    New triplet CC excited state amplitudes solver factory
!!    Written by Sarai D. Folkestad, 2022
!!
      implicit none
!
      type(cc_triplet_es_amplitudes_solver_factory) :: this
!
      character(len=*), intent(in) :: transformation
      character(len=*), intent(in) :: method
!
      logical, intent(in) :: restart
      character(len=200)  :: algorithm
!
      if (method .ne. 'ccs' .and. method .ne. 'cc2') &
          call output%error_msg('Triplets not implemented for '// method)
!
      this%transformation = trim(transformation)
      this%restart = restart
!
      algorithm = 'davidson'
      call input%get_keyword('algorithm', 'solver cc es', algorithm)
!
      if (trim(algorithm) .ne. 'davidson') &
         call output%error_msg('Only the Davidson algorithm can be used for triplet excited states')
!
   end function new_cc_triplet_es_amplitudes_solver_factory
!
!
   subroutine create(this, wf, solver)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2022
!!
      implicit none
!
      class(cc_triplet_es_amplitudes_solver_factory), intent(inout) :: this
!
      class(ccs), intent(inout) :: wf
!
      class(abstract_solver), allocatable, intent(out) :: solver
!
      call this%read_settings()
      call wf%prepare_for_triplet_excited_states(this%n_states)
!
      call this%create_davidson_es_solver(wf, solver)
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
      class(cc_triplet_es_amplitudes_solver_factory), intent(inout) :: this
!
      this%n_states           = 0
      this%max_iterations     = 100
      this%energy_threshold   = 1.0d-3
      this%residual_threshold = 1.0d-3
      this%records_in_memory  = .false.
!
      call input%get_keyword('triplet states', 'solver cc es', this%n_states)
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
!
   end subroutine read_settings
!
!
   subroutine create_davidson_es_solver(this, wf, solver)
!!
!!    Create davidson es solver
!!    Written by Sarai D. Folkestad, 2020
!!
      use eigen_davidson_solver_class, only: eigen_davidson_solver
      use eigen_davidson_tool_class,   only: eigen_davidson_tool
!
      use transformation_class,                             only: transformation
      use convergence_tool_class,                           only: convergence_tool
      use start_vector_tool_class,                          only: start_vector_tool
      use eigen_storage_tool_class,                         only: eigen_storage_tool
      use null_projection_tool_class,                       only: null_projection_tool
      use preconditioner_getter_class,                      only: preconditioner_getter
      use abstract_projection_tool_class,                   only: abstract_projection_tool
      use cc_triplet_es_storage_tool_class,                 only: cc_triplet_es_storage_tool
      use cc_triplet_jacobian_transformation_class,         only: cc_triplet_jacobian_transformation
      use triplet_es_valence_start_vector_tool_class,       only: triplet_es_valence_start_vector_tool
      use cc_triplet_jacobian_preconditioner_getter_class,  only: cc_triplet_jacobian_preconditioner_getter
!
      implicit none
!
!
      class(cc_triplet_es_amplitudes_solver_factory), intent(inout) :: this
      class(ccs),                                     intent(inout) :: wf
      class(abstract_solver), allocatable,            intent(out)   :: solver
!
      class(eigen_davidson_tool),         allocatable :: davidson
      class(eigen_storage_tool),          allocatable :: storer
      class(transformation),              allocatable :: transformer
      class(preconditioner_getter),       allocatable :: preconditioner
      class(convergence_tool),            allocatable :: convergence_checker
      class(start_vector_tool),           allocatable :: start_vectors
      class(abstract_projection_tool),    allocatable :: projector
!
      real(dp) :: lindep_threshold
      integer  :: max_dim_red
!
      max_dim_red = max(100, 10 * this%n_states)
      call input%get_keyword('max reduced dimension', 'solver cc es', max_dim_red)
!
      start_vectors = triplet_es_valence_start_vector_tool(wf, this%transformation, this%restart)
      projector = null_projection_tool(wf%n_triplet_amplitudes)
!
      preconditioner = cc_triplet_jacobian_preconditioner_getter(wf)
      storer = cc_triplet_es_storage_tool(wf, trim(this%transformation))
!
      transformer = cc_triplet_jacobian_transformation(wf, trim(this%transformation))
!
      convergence_checker = convergence_tool(this%energy_threshold,   &
                                             this%residual_threshold, &
                                             this%energy_convergence)
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
!

      davidson = eigen_davidson_tool(name_             = 'cc_es_davidson',    &
                                     n_parameters      = wf%n_triplet_amplitudes,  &
                                     n_solutions       = this%n_states,       &
                                     lindep_threshold  = lindep_threshold,    &
                                     max_dim_red       = max_dim_red,         &
                                     records_in_memory = this%records_in_memory)
!
      solver = eigen_davidson_solver(transformer         = transformer,              &
                                     davidson            = davidson,                 &
                                     convergence_checker = convergence_checker, &
                                     storer              = storer,                   &
                                     start_vector        = start_vectors,       &
                                     preconditioner      = preconditioner,           &
                                     projector           = projector,           &
                                     n_solutions         = this%n_states,            &
                                     max_iterations      = this%max_iterations)
!
   end subroutine create_davidson_es_solver
!
!
end module cc_triplet_es_amplitudes_solver_factory_class
