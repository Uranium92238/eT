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
module tdhf_solver_factory_class
!
!!
!!    TDHF solver factory class
!!    Written by Sarai D. Folkestad, May 2021
!!
!

   use kinds
   use global_out,                   only: output
   use transformation_tool_class,    only: transformation_tool
   use eigen_storage_tool_class,     only: eigen_storage_tool
   use hf_class,                     only: hf
   use eigen_davidson_tool_class,    only: eigen_davidson_tool
   use convergence_tool_class,       only: convergence_tool
   use general_eigen_davidson_class, only: general_eigen_davidson
   use start_vector_tool_class,      only: start_vector_tool
   use tdhf_start_vector_tool_class, only: tdhf_start_vector_tool
   use preconditioner_getter_class,  only: preconditioner_getter
!
   implicit none
!
   type:: tdhf_solver_factory
!
      integer :: n_states, max_dim_red, max_iterations
      real(dp) :: energy_threshold, residual_threshold
      logical :: records_in_memory, energy_convergence, tamm_dancoff, restart
!
   contains
!
      procedure, public :: create => create_tdhf_solver_factory
!
      procedure, private :: read_settings
!
      procedure, private, nopass :: create_rpa
      procedure, private, nopass :: create_tamm_dancoff
!
   end type tdhf_solver_factory
!
contains
!
!
    subroutine create_tdhf_solver_factory(this, wf, solver)
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
!!    for the determination of TDHF (RPA or Tamm-Dancoff)
!!    excitation energies.
!!
      implicit none
!
      class(tdhf_solver_factory),                 intent(inout)  :: this
      class(hf),                                  intent(in)  :: wf
      class(general_eigen_davidson), allocatable, intent(out) :: solver
!
      class(eigen_davidson_tool),  allocatable :: davidson
      class(convergence_tool),     allocatable :: convergence_checker
      class(transformation_tool),  allocatable :: transformer
      class(eigen_storage_tool),   allocatable :: storer
      class(start_vector_tool),    allocatable :: start_vector
      class(preconditioner_getter), allocatable :: preconditioner
!
      real(dp) :: lindep_threshold
!
      call this%read_settings()
!
      convergence_checker = convergence_tool(energy_threshold   = this%energy_threshold,     &
                                             residual_threshold = this%residual_threshold,   &
                                             energy_convergence = this%energy_convergence)
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
      davidson = eigen_davidson_tool(name_            = 'tdhf_davidson',      &
                                     n_parameters      = wf%n_v*wf%n_o,       &
                                     n_solutions       = this%n_states,       &
                                     lindep_threshold  = lindep_threshold,    &
                                     max_dim_red       = this%max_dim_red,    &
                                     records_in_memory = this%records_in_memory)
!
      start_vector = tdhf_start_vector_tool(wf, restart=this%restart)
!
      if (this%tamm_dancoff) then
!
         call output%printf('m','Tamm-Dancoff approximations (CIS/CCS) enabled!', ffs='(/t3, a)')
!
         call this%create_tamm_dancoff(wf, transformer, storer, preconditioner)
!
      else
!
         call this%create_rpa(wf, transformer, storer, preconditioner)
!
      endif
!
      solver = general_eigen_davidson(transformer           = transformer,         &
                                      davidson              = davidson,            &
                                      convergence_checker   = convergence_checker, &
                                      storer                = storer,              &
                                      start_vector          = start_vector,        &
                                      preconditioner        = preconditioner,      &
                                      n_solutions           = this%n_states,       &
                                      max_iterations        = this%max_iterations)
!
    end subroutine create_tdhf_solver_factory
!
!
   subroutine create_rpa(wf, transformer, storer, preconditioner)
!!
!!    Create RPA
!!    Written by Sarai D. Folkestad, May 2021
!!
!
      use rpa_transformation_tool_class,           only: rpa_transformation_tool
      use rpa_eigen_storage_tool_class,            only: rpa_eigen_storage_tool
      use rpa_preconditioner_getter_class,           only: rpa_preconditioner_getter
!
      implicit none
!
      class(transformation_tool),        allocatable, intent(inout)  :: transformer
      class(eigen_storage_tool),         allocatable, intent(inout)  :: storer
      class(hf),                                      intent(in)     :: wf
      class(preconditioner_getter),      allocatable, intent(inout)  :: preconditioner
!
      transformer    = rpa_transformation_tool(wf)
      storer         = rpa_eigen_storage_tool(wf)
      preconditioner = rpa_preconditioner_getter(wf)
!
   end subroutine create_rpa
!
!
   subroutine create_tamm_dancoff(wf, transformer, storer, preconditioner)
!!
!!    Create Tamm-Dancoff
!!    Written by Sarai D. Folkestad, May 2021
!!
!
      use tamm_dancoff_transformation_tool_class,  only: tamm_dancoff_transformation_tool
      use tamm_dancoff_eigen_storage_tool_class,   only: tamm_dancoff_eigen_storage_tool
      use tamm_dancoff_preconditioner_getter_class,only: tamm_dancoff_preconditioner_getter
!
      implicit none
!
      class(transformation_tool),        allocatable, intent(inout)  :: transformer
      class(eigen_storage_tool),         allocatable, intent(inout)  :: storer
      class(hf),                                      intent(in)     :: wf
      class(preconditioner_getter),      allocatable, intent(inout)  :: preconditioner
!
      transformer    = tamm_dancoff_transformation_tool(wf)
      storer         = tamm_dancoff_eigen_storage_tool(wf)
      preconditioner = tamm_dancoff_preconditioner_getter(wf)
!
   end subroutine create_tamm_dancoff
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad, May 2021
!!
!
      use global_in, only: input
!
      implicit none
!
      class(tdhf_solver_factory),   intent(inout)  :: this
!
      character(len=200) :: storage
!
!     Defaults
!
      this%n_states             = 0
      this%max_dim_red          = 50
      this%max_iterations       = 100
      this%energy_threshold     = 1.0d-3
      this%residual_threshold   = 1.0d-3
      this%energy_convergence   = .false.
      this%records_in_memory    = .true.
      this%tamm_dancoff         = .false.
      this%restart              = .false.
!
!
      call input%get_keyword('states', 'solver tdhf', this%n_states)
      if (this%n_states == 0) call output%error_msg('can not solve for 0 TDHF eigenstates')
!
      call input%get_keyword('max reduced dimension', 'solver tdhf', this%max_dim_red)
      call input%get_keyword('max iterations', 'solver tdhf', this%max_iterations)

      call input%get_keyword('residual threshold', 'solver tdhf', this%residual_threshold)

      if (input%is_keyword_present('energy threshold', 'solver tdhf')) then
!
         call input%get_keyword('energy threshold', 'solver tdhf', this%energy_threshold)
         this%energy_convergence = .true.
!
      endif
!
      call input%get_keyword('storage', 'solver tdhf', storage)
      if (trim(storage) == 'file') this%records_in_memory = .false.
!
      this%tamm_dancoff = input%is_keyword_present('tamm-dancoff', 'solver tdhf')
!
      this%restart = input%is_keyword_present('restart', 'solver tdhf') &
                    .or. input%is_keyword_present('restart', 'do')
!
   end subroutine read_settings
!
end module tdhf_solver_factory_class
