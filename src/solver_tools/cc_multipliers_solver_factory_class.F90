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
!!    A_mu,nu = < mu | [H-bar, tau_nu] | HF >,    H-bar = e-T H eT,
!!
!!    and
!!
!!    eta_mu = < HF | [H-bar, tau_nu] | HF >.
!!
!!    The multipliers are tbar and give the left CC ground state as
!!
!!    < Lambda | = < HF | e-T + sum_mu tbar_mu < mu | e-T
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

   use kinds
   use global_out,                                     only: output
   use transformation_tool_class,                      only: transformation_tool
   use cc_jacobian_transformation_tool_class,          only: cc_jacobian_transformation_tool
   use linear_storage_tool_class,                      only: linear_storage_tool
   use cc_multipliers_linear_storage_tool_class,       only: cc_multipliers_linear_storage_tool
   use ccs_class,                                      only: ccs
   use linear_davidson_tool_class,                     only: linear_davidson_tool
   use general_linear_davidson_class,                  only: general_linear_davidson
   use rhs_linear_equation_tool_class,                 only: rhs_linear_equation_tool
   use cc_multipliers_rhs_tool_class,                  only: cc_multipliers_rhs_tool
   use start_vector_tool_class,                        only: start_vector_tool
   use cc_multipliers_start_vector_tool_class,         only: cc_multipliers_start_vector_tool
   use preconditioner_getter_class,                    only: preconditioner_getter
   use cc_multipliers_preconditioner_getter_class,     only: cc_multipliers_preconditioner_getter
!
   implicit none
!
   type:: cc_multipliers_solver_factory
!
      integer            :: n_frequencies, max_dim_red, max_iterations
      real(dp)           :: residual_threshold
      logical            :: records_in_memory, restart
      character(len=200) :: side
!
   contains
!
      procedure, public :: create => create_cc_multipliers_solver_factory
!
      procedure, private :: read_settings
!
   end type cc_multipliers_solver_factory
!
contains
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
      class(cc_multipliers_solver_factory),        intent(inout)  :: this
      class(ccs),                                  intent(in)  :: wf
      class(general_linear_davidson), allocatable, intent(out) :: solver
!
      class(linear_davidson_tool),        allocatable  :: davidson
      class(transformation_tool),         allocatable  :: transformer
      class(linear_storage_tool),         allocatable  :: storer
      class(rhs_linear_equation_tool),    allocatable  :: rhs_getter
      class(start_vector_tool),           allocatable  :: start_vector
      class(preconditioner_getter),       allocatable  :: preconditioner
!
      real(dp) :: lindep_threshold
!
      call this%read_settings()
!
      lindep_threshold = min(1.0d-11, 0.1d0 * this%residual_threshold)
      davidson = linear_davidson_tool(name_            = 'multipliers_davidson',      &
                                     n_parameters      = wf%n_gs_amplitudes,          &
                                     lindep_threshold  = lindep_threshold,    &
                                     max_dim_red       = this%max_dim_red,    &
                                     n_equations       = this%n_frequencies,  &
                                     records_in_memory = this%records_in_memory)
!
      start_vector   = cc_multipliers_start_vector_tool(wf, restart=this%restart)
      transformer    = cc_jacobian_transformation_tool(wf, side=this%side)
      storer         = cc_multipliers_linear_storage_tool(wf)
      preconditioner = cc_multipliers_preconditioner_getter(wf)
      rhs_getter     = cc_multipliers_rhs_tool(wf)
!
      solver = general_linear_davidson(transformer           = transformer,         &
                                       davidson              = davidson,            &
                                       storer                = storer,              &
                                       start_vector          = start_vector,        &
                                       preconditioner        = preconditioner,      &
                                       rhs_getter            = rhs_getter,          &
                                       n_rhs                 = 1,                   &
                                       n_solutions           = this%n_frequencies,  &
                                       max_iterations        = this%max_iterations, &
                                       residual_threshold    = this%residual_threshold)
!
    end subroutine create_cc_multipliers_solver_factory
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Regina Matveeva, Sept 2021
!!
!!    Based on read_setting in tdhf_solver_factory by Sarai D. Folkestad
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
      this%n_frequencies        = 1
      this%max_dim_red          = 50
      this%max_iterations       = 100
      this%residual_threshold   = 1.0d-5
      this%records_in_memory    = .true.
      this%restart              = .false.
      this%side                 = 'left'
!
!
      call input%get_keyword('max reduced dimension', 'solver cc multipliers', this%max_dim_red)
      call input%get_keyword('max iterations', 'solver cc multipliers', this%max_iterations)

      call input%get_keyword('threshold', 'solver cc multipliers', this%residual_threshold)
!
      call input%get_keyword('storage', 'solver cc multipliers', storage)
      if (trim(storage) == 'disk') this%records_in_memory = .false.
!
      this%restart = input%is_keyword_present('restart', 'solver cc multipliers') &
                    .or. input%is_keyword_present('restart', 'do')
!
   end subroutine read_settings
!
end module cc_multipliers_solver_factory_class
