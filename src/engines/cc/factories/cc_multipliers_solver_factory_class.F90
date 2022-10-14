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
module cc_multipliers_solver_factory_class
!
!!
!!    Coupled Cluster multipliers solver factory class
!!    Written by Regina Matveeva, Sept 2021
!!
!!    Creates all the tools needed for the solution of the multiplier or left
!!    Coupled Cluster ground state equation
!!
!!                tbar^T A = - eta^T
!!
!!    for t-bar. Here, A is the Coupled Cluster Jacobian
!!
!!    A_mu,nu = < mu |[H-bar, tau_nu] | HF >,    H-bar = e-T H eT,
!!
!!    and
!!
!!    eta_mu = < HF | [H-bar, tau_nu] | HF >.
!!
!!    The multipliers are tbar and give the left CC ground state as
!!
!!    < Lambda| = < HF| e-T + sum_mu tbar_mu < mu| e-T
!!
!!    The solutions are determined using the Davidson
!!    reduced space algorithm, where the eigenvalue problem
!!    is solved in a subspace generated from the residuals* obtained
!!    in the preceding iterations. See E. R. Davidson, J. Comput. Phys.
!!    17, 87 (1975) for more details.
!!
!!    * Preconditioned using orbital differences approximation of A
!!     (See davidson_cc_es solver for more details.)
!

   use parameters
   use abstract_solver_class, only: abstract_solver
   use global_out,            only: output
   use ccs_class,             only: ccs
!
   implicit none
!
   type :: cc_multipliers_solver_factory
!
      integer            :: max_iterations
      real(dp)           :: residual_threshold
      logical            :: records_in_memory, restart
      character(len=200) :: algorithm
!
   contains
!
      procedure, public :: create => create_cc_multipliers_solver_factory
!
      procedure, private :: create_davidson
      procedure, private :: create_diis
      procedure, private :: create_nr
!
      procedure, private, nopass :: read_davidson_settings
      procedure, private, nopass :: read_nr_settings
      procedure, private :: read_algorithm
      procedure, private :: read_general_settings
!
   end type cc_multipliers_solver_factory
!
!
   interface cc_multipliers_solver_factory
!
      procedure :: new_cc_multipliers_solver_factory
!
   end interface cc_multipliers_solver_factory
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
      type(cc_multipliers_solver_factory) :: this
!
      character(len=*), intent(in) :: method
!
      call this%read_algorithm(method)
      call this%read_general_settings()
!
   end function new_cc_multipliers_solver_factory
!
!
   subroutine create_cc_multipliers_solver_factory(this, wf, solver)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Initializes the Davidson solver with the correct instances of the following types:
!!
!!       - Solver tool (linear_davidson)
!!       - Transformation tool
!!       - RHS of a linear equation tool
!!       - Storage tool
!!       - Start vector tool
!!       - Preconditioner tool
!!
!!    for the determination of Coupled Cluster multipliers
!!
      implicit none
!
      class(cc_multipliers_solver_factory), intent(inout) :: this
      class(ccs),                           intent(in)    :: wf
      class(abstract_solver), allocatable,  intent(out)   :: solver
!
      if (trim(this%algorithm) == 'davidson') call this%create_davidson(wf, solver)
      if (trim(this%algorithm) == 'diis') call this%create_diis(wf, solver)
      if (trim(this%algorithm) == 'newton-raphson') call this%create_nr(wf, solver)
!
   end subroutine create_cc_multipliers_solver_factory
!
!
   subroutine create_davidson(this, wf, solver)
!!
!!    Create Davidson
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Initializes the Davidson solver with the correct instances of the following types:
!!
!!       - Solver tool (linear_davidson)
!!       - Transformation tool
!!       - RHS of a linear equation tool
!!       - Storage tool
!!       - Start vector tool
!!       - Preconditioner tool
!!
!!    for the determination of Coupled Cluster multipliers
!!

      use transformation_class,                             only: transformation
      use cc_jacobian_transformation_class,                 only: cc_jacobian_transformation
      use cc_multipliers_linear_equation_storage_tool_class,only: cc_multipliers_linear_equation_storage_tool
      use linear_davidson_tool_class,                       only: linear_davidson_tool
      use linear_davidson_solver_class,                     only: linear_davidson_solver
      use cc_multipliers_rhs_tool_class,                    only: cc_multipliers_rhs_tool
      use cc_multipliers_start_vector_tool_class,           only: cc_multipliers_start_vector_tool
      use cc_jacobian_preconditioner_getter_class,          only: cc_jacobian_preconditioner_getter
      use linear_davidson_single_solution_print_tool_class, only: linear_davidson_single_solution_print_tool
      use folded_cc_jacobian_transformation_class, only: folded_cc_jacobian_transformation
!
      implicit none
!
      class(cc_multipliers_solver_factory), intent(inout) :: this
      class(ccs),                           intent(in)    :: wf
      class(abstract_solver), allocatable,  intent(out)   :: solver
!
      class(linear_davidson_tool),                       allocatable :: davidson
      class(transformation),                             allocatable :: transformer
      class(cc_multipliers_linear_equation_storage_tool),allocatable :: storer
      class(cc_multipliers_rhs_tool),                    allocatable :: rhs_getter
      class(cc_multipliers_start_vector_tool),           allocatable :: start_vector
      class(cc_jacobian_preconditioner_getter),          allocatable :: preconditioner
      class(linear_davidson_single_solution_print_tool), allocatable :: printer
!
      real(dp)               :: lindep_threshold
      real(dp), dimension(1) :: frequencies
      integer                :: max_reduced_dimension
!
      call this%read_davidson_settings(max_reduced_dimension)
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
!
      printer = linear_davidson_single_solution_print_tool()
!
      davidson = linear_davidson_tool(name_             = 'multipliers_davidson',&
                                      n_parameters      = wf%n_gs_amplitudes,    &
                                      lindep_threshold  = lindep_threshold,      &
                                      max_dim_red       = max_reduced_dimension, &
                                      n_equations       = 1,                     &
                                      records_in_memory = this%records_in_memory)
!
      start_vector   = cc_multipliers_start_vector_tool(wf, restart=this%restart)
!
      if (wf%name_ == 'lowmem cc2' .or. wf%name_ == 'cc3' ) then
!
         transformer = folded_cc_jacobian_transformation(wf,          &
                                                         side='left', &
                                                         n_parameters=wf%n_gs_amplitudes,&
                                                         frequency=zero)
!
      else
!
         transformer = cc_jacobian_transformation(wf,          &
                                                  side='left', &
                                                  n_parameters=wf%n_gs_amplitudes)
      endif
!
      storer         = cc_multipliers_linear_equation_storage_tool(wf)
      preconditioner = cc_jacobian_preconditioner_getter(wf, wf%n_gs_amplitudes)
      rhs_getter     = cc_multipliers_rhs_tool(wf)
!
      frequencies = zero
      solver = linear_davidson_solver(transformer        = transformer,             &
                                      davidson           = davidson,                &
                                      storer             = storer,                  &
                                      preconditioner     = preconditioner,          &
                                      rhs_getter         = rhs_getter,              &
                                      printer            = printer,                 &
                                      n_rhs              = 1,                       &
                                      n_solutions        = 1,                       &
                                      max_iterations     = this%max_iterations,     &
                                      residual_threshold = this%residual_threshold, &
                                      frequencies        = frequencies,             &
                                      start_vector       = start_vector)
!
   end subroutine create_davidson
!
!
   subroutine create_diis(this, wf, solver)
!!
!!    Create diis
!!    Written by Sarai D. Folkestad, 2022
!
      use quasi_newton_updater_class, only: quasi_newton_updater
      use diis_cc_multipliers_class,    only: diis_cc_multipliers
!
      implicit none
!
      class(cc_multipliers_solver_factory), intent(inout) :: this
      class(ccs),                           intent(in)    :: wf
      class(abstract_solver), allocatable,  intent(out)   :: solver
!
      class(quasi_newton_updater), allocatable :: t_updater
!
      t_updater = quasi_newton_updater(n_amplitudes=wf%n_gs_amplitudes, &
                                       scale_amplitudes=.false.,        &
                                       scale_residual=.false.)
!
      solver = diis_cc_multipliers(wf, this%restart, t_updater)
!
   end subroutine create_diis
!
!
   subroutine create_nr(this, wf, solver)
!!
!!    Create nr
!!    Written by Sarai D. Folkestad, 2022
!!
      use citation_printer_class, only: eT_citations
!
      use newton_raphson_updater_class, only: newton_raphson_updater
      use diis_cc_multipliers_class,    only: diis_cc_multipliers
!
      use abstract_jacobian_transformer_class,     only: abstract_jacobian_transformer
      use jacobian_transformer_class,              only: jacobian_transformer
      use approximate_jacobian_transformer_class,  only: approximate_jacobian_transformer
!
      implicit none
!
      class(cc_multipliers_solver_factory), intent(inout) :: this
      class(ccs),                           intent(in)    :: wf
      class(abstract_solver), allocatable,  intent(out)   :: solver
!
      class(newton_raphson_updater), allocatable :: t_updater
      class(abstract_jacobian_transformer), allocatable :: transformer
!
      character(len=200) :: multimodel_newton
      real(dp) :: relative_threshold
      logical :: micro_records_in_memory
      integer :: max_micro_iterations
!
      call this%read_nr_settings(wf,                        &
                                 multimodel_newton,         &
                                 relative_threshold,        &
                                 micro_records_in_memory,   &
                                 max_micro_iterations)
!
      if (trim(multimodel_newton) == 'on') then
!
         transformer = approximate_jacobian_transformer('left')
         call eT_citations%add('Multimodel Newton algorithm')
!
      else
!
         transformer = jacobian_transformer('left')
!
      end if
!
      t_updater = newton_raphson_updater(n_amplitudes          = wf%n_gs_amplitudes,      &
                                            scale_amplitudes   = .false.,                 &
                                            scale_residual     = .false.,                 &
                                            relative_threshold = relative_threshold,      &
                                            records_in_memory  = micro_records_in_memory, &
                                            max_iterations     = max_micro_iterations,    &
                                            transformer        = transformer)
!
      solver = diis_cc_multipliers(wf, this%restart, t_updater)
!
   end subroutine create_nr
!
!
   subroutine read_davidson_settings(max_reduced_dimension)
!!
!!    Read davidson settings
!!    Written by Regina Matveeva, Sept 2021
!!
      use global_in, only: input
!
      implicit none
!
      integer, intent(out) :: max_reduced_dimension
!
      max_reduced_dimension          = 50
      call input%get_keyword('max reduced dimension', &
                             'solver cc multipliers', &
                             max_reduced_dimension)
!
   end subroutine read_davidson_settings
!
!
   subroutine read_algorithm(this, method)
!!
!!    Read algorithm
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_multipliers_solver_factory), intent(inout) :: this
      character(len=*), intent(in) :: method
!
      if (method == 'cc3') then
!
         this%algorithm = 'newton-raphson'
!
      else if (method .eq. 'ccs'            .or. &
               method .eq. 'cc2'            .or. &
               method .eq. 'low memory cc2' .or. &
               method .eq. 'mlcc2') then
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
   end subroutine read_algorithm
!
!
   subroutine read_general_settings(this)
!!
!!    Read general settings
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_multipliers_solver_factory),   intent(inout)  :: this
!
      character(len=200) :: storage
!
!     Defaults
!
      this%max_iterations       = 100
      this%residual_threshold   = 1.0d-5
      this%records_in_memory    = .true.
      this%restart              = .false.
!
      call input%get_keyword('max iterations', 'solver cc multipliers', this%max_iterations)
!
      call input%get_keyword('threshold', 'solver cc multipliers', this%residual_threshold)
!
      storage = 'disk'
      call input%get_keyword('storage', 'solver cc multipliers', storage)
!
      if (trim(storage) == 'disk') this%records_in_memory = .false.
!
      this%restart = input%is_keyword_present('restart', 'solver cc multipliers') .or. &
                     input%is_keyword_present('restart', 'do')
!
   end subroutine read_general_settings
!
!
   subroutine read_nr_settings(wf,                        &
                               multimodel_newton,         &
                               relative_threshold,        &
                               micro_records_in_memory,   &
                               max_micro_iterations)
!!
!!    Read NR settings
!!    Written by Sarai D. Folkestad, Jam 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      character(len=200), intent(out) :: multimodel_newton
      real(dp), intent(out) :: relative_threshold
      logical, intent(out) :: micro_records_in_memory
      integer, intent(out) :: max_micro_iterations
!
      character(len=200) :: storage
!
      relative_threshold   = 1.0d-2
      max_micro_iterations = 100
      storage              = 'disk'
!
      multimodel_newton = 'off'
      if (trim(wf%name_) == 'cc3') multimodel_newton = 'on'
!
      call input%get_keyword('rel micro threshold',     'solver cc multipliers', relative_threshold)
      call input%get_keyword('max micro iterations',    'solver cc multipliers', max_micro_iterations)
      call input%get_keyword('micro iteration storage', 'solver cc multipliers', storage)
      call input%get_keyword('multimodel newton', 'solver cc multipliers', multimodel_newton)
!
      micro_records_in_memory = .false.
      if (trim(storage) == 'memory') micro_records_in_memory = .true.
!
   end subroutine read_nr_settings
!
end module cc_multipliers_solver_factory_class
