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
module fci_solver_factory_class
!
!!
!! FCI solver factory class
!! Written by Enrico Ronca, 2022
!!
!
   use kinds
   use global_out,                      only: output
   use transformation_class,            only: transformation
   use eigen_storage_tool_class,        only: eigen_storage_tool
   use fci_class,                       only: fci
   use eigen_davidson_tool_class,       only: eigen_davidson_tool
   use convergence_tool_class,          only: convergence_tool
   use eigen_davidson_solver_class,     only: eigen_davidson_solver
   use start_vector_tool_class,         only: start_vector_tool
   use fci_start_vector_tool_class,     only: fci_start_vector_tool
   use preconditioner_getter_class,     only: preconditioner_getter
!
   implicit none
!
   type:: fci_solver_factory
!
      integer  :: n_states, max_dim_red, max_iterations
      real(dp) :: energy_threshold, residual_threshold
      logical  :: records_in_memory, energy_convergence, tamm_dancoff, restart
!
   contains
!
      procedure, public :: create => create_fci_solver_factory
!
      procedure, private :: read_settings
!
   end type fci_solver_factory
!
contains
!
!
   subroutine create_fci_solver_factory(this, wf, solver)
!!
!!    Create
!!    Written by Enrico Ronca, 2022
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
!!    for the determination of FCI energies and eigenfunctions.
!!
      use null_projection_tool_class,      only: null_projection_tool
      use fci_transformation_class,        only: fci_transformation
      use fci_eigen_storage_tool_class,    only: fci_eigen_storage_tool
      use fci_preconditioner_getter_class, only: fci_preconditioner_getter
!
      implicit none
!
      class(fci_solver_factory),                 intent(inout)  :: this
      class(fci),                                intent(in)     :: wf
      class(eigen_davidson_solver), allocatable, intent(out)    :: solver
!
      class(eigen_davidson_tool),   allocatable :: davidson
      class(convergence_tool),      allocatable :: convergence_checker
      class(transformation),        allocatable :: transformer
      class(eigen_storage_tool),    allocatable :: storer
      class(start_vector_tool),     allocatable :: start_vector
      class(preconditioner_getter), allocatable :: preconditioner
      class(null_projection_tool),  allocatable :: projector
!
      real(dp) :: lindep_threshold
!
      call this%read_settings()
!
      convergence_checker = convergence_tool(energy_threshold   = this%energy_threshold,     &
                                             residual_threshold = this%residual_threshold,   &
                                             energy_convergence = this%energy_convergence)
!
      start_vector = fci_start_vector_tool(wf, restart=this%restart)
      projector = null_projection_tool(wf%n_determinants)
!
      transformer    = fci_transformation(wf)
      storer         = fci_eigen_storage_tool(wf)
      preconditioner = fci_preconditioner_getter(wf)
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
!
      davidson = eigen_davidson_tool(name_             = 'fci_davidson',      &
                                     n_parameters      = wf%n_determinants,   &
                                     n_solutions       = this%n_states,       &
                                     lindep_threshold  = lindep_threshold,    &
                                     max_dim_red       = this%max_dim_red,    &
                                     records_in_memory = this%records_in_memory)
!
      solver = eigen_davidson_solver(transformer           = transformer,         &
                                     davidson              = davidson,            &
                                     convergence_checker   = convergence_checker, &
                                     storer                = storer,              &
                                     start_vector          = start_vector,        &
                                     preconditioner        = preconditioner,      &
                                     projector             = projector,           &
                                     n_solutions           = this%n_states,       &
                                     max_iterations        = this%max_iterations)
!
   end subroutine create_fci_solver_factory
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Enrico Ronca, 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(fci_solver_factory),   intent(inout)  :: this
      character(len=200) :: storage
!
!     Defaults
!
      this%n_states             = 1
      this%max_iterations       = 100
      this%max_dim_red          = 100
      this%energy_threshold     = 1.0d-6
      this%residual_threshold   = 1.0d-6
      this%energy_convergence   = .false.
      this%records_in_memory    = .false.
      this%restart              = .false.
!
!
      call input%get_keyword('residual threshold', 'solver fci', this%residual_threshold)
      call input%get_keyword('max iterations', 'solver fci', this%max_iterations)
      call input%get_keyword('max reduced dimension', 'solver fci', this%max_dim_red)
!
      call input%get_keyword('states', 'solver fci', this%n_states)
!
      this%restart = input%is_keyword_present('restart', 'solver fci') &
                    .or. input%is_keyword_present('restart', 'do')
!
      storage = 'disk'
      call input%get_keyword('storage', 'solver fci', storage)
      if (trim(storage) == 'memory') this%records_in_memory = .true.
!
      if (input%is_keyword_present('energy threshold', 'solver fci')) then
!
         call input%get_keyword('energy threshold', 'solver fci', this%energy_threshold)
         this%energy_convergence = .true.
!
      endif
!
   end subroutine read_settings
!
!
end module fci_solver_factory_class
