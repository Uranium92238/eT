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
module davidson_cc_es_solver_factory_class
!
!!
!!    Davidson CC es solver factory class
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!
   use kinds
   use global_out,                      only: output
   use global_in,                       only: input
   use abstract_solver_class,           only: abstract_solver
   use transformation_tool_class,       only: transformation_tool
   use eigen_storage_tool_class,        only: eigen_storage_tool
   use ccs_class,                       only: ccs
   use eigen_davidson_tool_class,       only: eigen_davidson_tool
   use convergence_tool_class,          only: convergence_tool
   use eigen_davidson_solver_class,     only: eigen_davidson_solver
   use start_vector_tool_class,         only: start_vector_tool
   use preconditioner_getter_class,     only: preconditioner_getter
   use cc_es_eigen_davidson_print_tool_class, only: cc_es_eigen_davidson_print_tool
   use cc_jacobian_transformation_tool_class, only: cc_jacobian_transformation_tool
   use abstract_projection_tool_class, only : abstract_projection_tool
!
   implicit none
!
   type:: davidson_cc_es_solver_factory
!
      integer :: n_states, max_dim_red, max_iterations
      real(dp) :: energy_threshold, residual_threshold
      logical :: records_in_memory, energy_convergence, tamm_dancoff, restart
      character(len=200) :: transformation, es_type
!
   contains
!
      procedure, public :: create => create_davidson_cc_es_solver_factory
!
      procedure, private :: read_settings
      procedure, private :: get_start_vector_tool
      procedure, private :: get_projection_tool
!
   end type davidson_cc_es_solver_factory
!
   interface davidson_cc_es_solver_factory
!
      procedure :: new_davidson_cc_es_solver_factory
!
   end interface
!
contains
!
!
   function new_davidson_cc_es_solver_factory(transformation, restart) result(this)
!!
!!    New
!!
      implicit none
!
      character(len=*), intent(in) :: transformation
!
      type(davidson_cc_es_solver_factory) :: this
      logical :: restart
!
      this%transformation = transformation
      this%restart = restart
!
   end function new_davidson_cc_es_solver_factory
!
!
   subroutine create_davidson_cc_es_solver_factory(this, wf, solver)
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
      use cc_es_eigen_davidson_print_tool_class, only: cc_es_eigen_davidson_print_tool
      use cc_eigen_storage_tool_class, only: cc_eigen_storage_tool
      use cc_jacobian_preconditioner_getter_class, only: cc_jacobian_preconditioner_getter
!
      implicit none
!
      class(davidson_cc_es_solver_factory),                      intent(inout)  :: this
      class(ccs),                                  intent(inout)     :: wf
      class(abstract_solver), allocatable,        intent(out)    :: solver
!
      class(eigen_davidson_tool),       allocatable :: davidson
      class(convergence_tool),          allocatable :: convergence_checker
      class(transformation_tool),       allocatable :: transformer
      class(eigen_storage_tool),        allocatable :: storer
      class(start_vector_tool),         allocatable :: start_vector
      class(preconditioner_getter),     allocatable :: preconditioner
      type(cc_es_eigen_davidson_print_tool), allocatable :: printer
      class(abstract_projection_tool),        allocatable :: projector
!
      real(dp) :: lindep_threshold
!
      call this%read_settings()
!
      wf%n_singlet_states = this%n_states
!
      if (trim(this%es_type) == 'core') then
!
         call wf%read_cvs_settings()

      elseif (trim(this%es_type) == 'remove core') then
!
         call wf%read_rm_core_settings()
!
      endif
!
      convergence_checker = convergence_tool(energy_threshold   = this%energy_threshold,     &
                                             residual_threshold = this%residual_threshold,   &
                                             energy_convergence = this%energy_convergence)
!
!
      call this%get_start_vector_tool(wf, start_vector)
      call this%get_projection_tool(wf, projector)
!
      preconditioner = cc_jacobian_preconditioner_getter(wf, wf%n_es_amplitudes)
      storer = cc_eigen_storage_tool(wf, trim(this%transformation))
!
      transformer = cc_jacobian_transformation_tool(wf, trim(this%transformation), wf%n_es_amplitudes)
!
      printer = cc_es_eigen_davidson_print_tool(wf, trim(this%transformation))
      call printer%print_banner()
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
      davidson = eigen_davidson_tool(name_            = 'cc_es_davidson',     &
                                     n_parameters      = wf%n_es_amplitudes,  &
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
                                     printer               = printer,             &
                                     projector             = projector,           &
                                     n_solutions           = this%n_states,       &
                                     max_iterations        = this%max_iterations)
!
   end subroutine create_davidson_cc_es_solver_factory
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
      class(davidson_cc_es_solver_factory),   intent(inout)  :: this
!
      character(len=200) :: storage
!
      this%n_states             = 0
      this%max_iterations       = 100
      this%energy_threshold     = 1.0d-3
      this%residual_threshold   = 1.0d-3
      this%energy_convergence   = .false.
      this%records_in_memory    = .true.
!
!
      call input%get_keyword('singlet states', 'solver cc es', this%n_states)
      if (this%n_states == 0) call output%error_msg('can not solve for 0 excited states')
!
      this%max_dim_red = max(100, 10*this%n_states)
!
      call input%get_keyword('max reduced dimension', 'solver cc es', this%max_dim_red)
      call input%get_keyword('max iterations', 'solver cc es', this%max_iterations)

      call input%get_keyword('residual threshold', 'solver cc es', this%residual_threshold)

      if (input%is_keyword_present('energy threshold', 'solver cc es')) then
!
         call input%get_keyword('energy threshold', 'solver cc es', this%energy_threshold)
         this%energy_convergence = .true.
!
      endif
!
      call input%get_keyword('storage', 'solver cc es', storage)
      if (trim(storage) == 'file') this%records_in_memory = .false.
!
      this%es_type = 'valence'
      if (input%is_keyword_present('core excitation', 'solver cc es') .and. .not. &
          input%is_keyword_present('ionization', 'solver cc es')) this%es_type = 'core'
!
      if (input%is_keyword_present('ionization', 'solver cc es') .and. .not. &
          input%is_keyword_present('core excitation', 'solver cc es')) this%es_type = 'ionize'
!
      if (input%is_keyword_present('remove core', 'solver cc es')) this%es_type = 'remove core'
!
   end subroutine read_settings
!
!
   subroutine get_start_vector_tool(this, wf, start_vector)
!!
!!    Initialize start vector tool
!!    Written by Eirik F. Kjønstad, Sep 2019
!!
      use global_in, only: input
!
      use start_vector_tool_class,              only: start_vector_tool
      use es_manual_start_vector_tool_class,    only: es_manual_start_vector_tool
      use es_valence_start_vector_tool_class,   only: es_valence_start_vector_tool
      use es_cvs_start_vector_tool_class,       only: es_cvs_start_vector_tool
      use es_ip_start_vector_tool_class,        only: es_ip_start_vector_tool
!
      implicit none
!
      class(davidson_cc_es_solver_factory), intent(in)  :: this
      class(ccs) :: wf
      class(start_vector_tool), allocatable, intent(out) :: start_vector
!
      if (input%is_keyword_present('state guesses', 'solver cc es')) then
!
         start_vector = es_manual_start_vector_tool(wf, &
                                                    this%transformation, &
                                                    this%restart)
!
      else
!
         if (trim(this%es_type) == 'valence') then
!
            start_vector = es_valence_start_vector_tool(wf, &
                                                        this%transformation, &
                                                        this%restart)
!
         elseif (trim(this%es_type) == 'core') then
!
            start_vector = es_cvs_start_vector_tool(wf, &
                                                    this%transformation, &
                                                    this%restart)
!
         elseif (trim(this%es_type) == 'ionize') then
!
            start_vector = es_ip_start_vector_tool(wf, &
                                                   this%transformation, &
                                                   this%restart)
!
         elseif (trim(this%es_type) == 'remove core') then
!
            start_vector = es_valence_start_vector_tool(wf, &
                                                        this%transformation, &
                                                        this%restart)
!
         else
!
            call output%error_msg('could not recognize excited state type in abstract_cc_es')
!
         endif
!
      endif
!
   end subroutine get_start_vector_tool
!
!
   subroutine get_projection_tool(this, wf, projector)
!!
!!    Get projection tool
!!    Written by Eirik F. Kjønstad, Sep 2019
!!
!
      use null_projection_tool_class, only: null_projection_tool
      use es_cvs_projection_tool_class, only: es_cvs_projection_tool
      use es_ip_projection_tool_class, only: es_ip_projection_tool
      use es_rm_core_projection_tool_class, only: es_rm_core_projection_tool
!
      implicit none
!
      class(davidson_cc_es_solver_factory) :: this
      class(ccs) :: wf
      class(abstract_projection_tool), allocatable, intent(out) :: projector
!
      if (trim(this%es_type) == 'valence') then
!
         projector = null_projection_tool(wf%n_es_amplitudes)
!
      elseif (trim(this%es_type) == 'core') then
!
         projector = es_cvs_projection_tool(wf)
!
      elseif (trim(this%es_type) == 'ionize') then
!
         projector = es_ip_projection_tool(wf)
!
      elseif (trim(this%es_type) == 'remove core') then
!
         projector = es_rm_core_projection_tool(wf)
!
      else
!
         call output%error_msg('could not recognize excited state type in abstract_cc_es')
!
      endif
!
   end subroutine get_projection_tool
!
!
end module davidson_cc_es_solver_factory_class
