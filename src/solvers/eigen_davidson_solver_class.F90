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
module eigen_davidson_solver_class
!
!!
!! General Davidson solver class
!! Written by Sarai D. Folkestad, May 2021
!!
!
   use parameters
   use memory_manager_class, only: mem
   use eigen_davidson_tool_class, only: eigen_davidson_tool
   use timings_class, only: timings
!
   use abstract_solver_class,           only: abstract_solver
   use convergence_tool_class,          only: convergence_tool
   use transformation_class,            only: transformation
   use eigen_storage_tool_class,        only: eigen_storage_tool
   use start_vector_tool_class,         only: start_vector_tool
   use preconditioner_getter_class,     only: preconditioner_getter
   use eigen_davidson_print_tool_class, only: eigen_davidson_print_tool
   use abstract_projection_tool_class,  only: abstract_projection_tool
!
   implicit none
!
   type, extends(abstract_solver) :: eigen_davidson_solver
!
      integer, private :: n_solutions
      integer, private :: max_iterations
!
      class(eigen_davidson_tool), allocatable, private  :: davidson
!
      class(convergence_tool),          allocatable, private :: convergence_checker
      class(transformation),            allocatable, private :: transformer
      class(eigen_storage_tool),        allocatable, private :: storer
      class(start_vector_tool),         allocatable, private :: start_vector
      class(preconditioner_getter),     allocatable, private :: preconditioner
      class(eigen_davidson_print_tool), allocatable, private :: printer
      class(abstract_projection_tool),  allocatable, private :: projector
!
   contains
!
      procedure, public :: run => run_eigen_davidson_solver
      procedure, public :: get_n_solutions => get_n_solutions_eigen_davidson_solver
!
      procedure, private :: test_convergence_and_add_trials
!
end type eigen_davidson_solver
!
   interface eigen_davidson_solver
!
      procedure :: new_eigen_davidson_solver
!
end interface eigen_davidson_solver
!
contains
!
   function new_eigen_davidson_solver(transformer,         &
                                      davidson,            &
                                      convergence_checker, &
                                      storer,              &
                                      start_vector,        &
                                      preconditioner,      &
                                      projector,           &
                                      n_solutions,         &
                                      max_iterations) result(this)
!!
!!    New eigen Davidson
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      type(eigen_davidson_solver) :: this
!
      class(transformation),            intent(in) :: transformer
      type(eigen_davidson_tool),        intent(in) :: davidson
      type(convergence_tool),           intent(in) :: convergence_checker
      class(eigen_storage_tool),        intent(in) :: storer
      class(start_vector_tool),         intent(in) :: start_vector
      class(preconditioner_getter),     intent(in) :: preconditioner
      class(abstract_projection_tool),  intent(in) :: projector
      integer,                          intent(in) :: max_iterations, n_solutions
!
      this%davidson            = davidson
      this%transformer         = transformer
      this%convergence_checker = convergence_checker
      this%storer              = storer
      this%start_vector        = start_vector
      this%preconditioner      = preconditioner
      this%projector           = projector
      this%n_solutions         = n_solutions
      this%max_iterations      = max_iterations
!
      this%printer = eigen_davidson_print_tool()
!
      call this%printer%print_banner()
      call this%davidson%print_settings
!
      this%total_timer = timings("Eigen davidson solver total", pl='m')
      this%iteration_timer = timings("Eigen davidson solver iteration", pl='m')
!
   end function new_eigen_davidson_solver
!
!
   subroutine run_eigen_davidson_solver(this)
!!
!!    Run
!!    Written by Sarai D. Folkestad, May 2021
!!
      use array_initialization, only: zero_array, zero_array_complex
!
      implicit none
!
      class(eigen_davidson_solver), intent(inout) :: this

      logical, dimension(:), allocatable :: converged
      logical, dimension(:), allocatable :: residual_lt_lindep
!
      integer :: iteration, trial
!
      real(dp), dimension(:), allocatable :: residual_norm
      real(dp), dimension(:), allocatable :: c, residual, omega
      complex(dp), dimension(:), allocatable :: new_omega
!
      call this%total_timer%turn_on()
!
      call mem%alloc(omega, this%n_solutions)
      call mem%alloc(new_omega, this%n_solutions)
!
      call this%davidson%initialize()
      call this%transformer%initialize()
!
      call mem%alloc(c, this%davidson%n_parameters)
      call this%preconditioner%get(c)
      call this%davidson%set_preconditioner(c)
!
      do trial = 1,this%n_solutions
!
         call this%start_vector%get(c, trial, omega(trial))
         call this%davidson%set_trial(c, trial)
!
      enddo
!
      call this%printer%print_settings(this%n_solutions, this%max_iterations)
!
      call zero_array(omega, this%n_solutions)
      call zero_array_complex(new_omega, this%n_solutions)
!
      call mem%alloc(converged, this%n_solutions, set_to=.false.)
      call mem%alloc(residual_lt_lindep, this%n_solutions, set_to=.false.)
!
      call mem%alloc(residual_norm, this%n_solutions)
      call mem%alloc(residual, this%davidson%n_parameters)
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. this%max_iterations))
!
         call this%iteration_timer%turn_on()
!
         iteration = iteration + 1
!
         call this%davidson%update_reduced_space()
!
         do trial = this%davidson%first_new_trial(), this%davidson%last_new_trial()
!
            call this%davidson%get_trial(c, trial)
            call this%transformer%transform(c, residual)
!
            call this%projector%project(residual)
            call this%davidson%set_transform(residual, trial)
!
         enddo
!
         call this%davidson%solve_reduced_problem()
!
         new_omega = this%davidson%get_omega(this%n_solutions)
!
         call this%test_convergence_and_add_trials(omega,               &
                                                   new_omega,           &
                                                   residual_norm,       &
                                                   converged,           &
                                                   residual_lt_lindep,  &
                                                   iteration)
!
         call this%printer%print_iteration(this%n_solutions, omega, new_omega, &
                                           residual_norm, iteration, this%davidson%dim_red)
!
         call dcopy(this%n_solutions, real(new_omega), 1, omega, 1)
!
         call this%iteration_timer%turn_off()
         call this%iteration_timer%reset()
!
      enddo
!
      call mem%dealloc(c, this%davidson%n_parameters)
      call mem%dealloc(residual, this%davidson%n_parameters)
      call mem%dealloc(residual_norm, this%n_solutions)
!
      call this%printer%print_summary(this%n_solutions, converged, iteration)
!
      call mem%dealloc(converged, this%n_solutions)
      call mem%dealloc(residual_lt_lindep, this%n_solutions)
!
      call mem%dealloc(omega, this%n_solutions)
      call mem%dealloc(new_omega, this%n_solutions)
!
      call this%davidson%cleanup()
!
      call this%total_timer%turn_off()
!
   end subroutine run_eigen_davidson_solver
!
!
   subroutine test_convergence_and_add_trials(this,               &
                                              omega,              &
                                              new_omega,          &
                                              residual_norm,      &
                                              converged,          &
                                              residual_lt_lindep, &
                                              iteration)
!!
!!    Test convergence and add trial
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    - Get residual and solution
!!    - Test convergence
!!    - Add new trial if needed and residual norm is greater
!!      than linear dependence threshold
!!
      use global_out, only: output
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(eigen_davidson_solver), intent(inout) :: this
      real(dp), dimension(this%n_solutions), intent(in)  :: omega
      complex(dp), dimension(this%n_solutions), intent(in)  :: new_omega
      real(dp), dimension(this%n_solutions), intent(out) :: residual_norm
      logical, dimension(this%n_solutions), intent(out)  :: converged, residual_lt_lindep
      integer, intent(in) :: iteration
!
      real(dp), dimension(:), allocatable :: residual
      real(dp), dimension(:), allocatable :: solution
!
      integer :: n
!
      call mem%alloc(residual, this%davidson%n_parameters)
      call mem%alloc(solution, this%davidson%n_parameters)
!
      do n = 1, this%n_solutions
!
         call this%davidson%construct_solution(solution, n)
         call this%davidson%construct_residual(residual, solution, n)
!
         call this%storer%store(real(new_omega(n)), solution, n)
!
         residual_norm(n) = get_l2_norm(residual, this%davidson%n_parameters)
         converged(n)     = this%convergence_checker%has_converged(residual_norm(n), &
                            real(new_omega(n))-omega(n), iteration)
!
         residual_lt_lindep(n) = (residual_norm(n) .le. this%davidson%lindep_threshold)
!
         if (.not. converged(n)) then
!
            if (residual_lt_lindep(n)) then
!
               call output%warning_msg('Residual norm for root (i0) smaller than linear ' // &
                                       'dependence threshold, but energy and residual  '  // &
                                       'thresholds have not yet been met. No new trial ' // &
                                       'added for this root.', ints=[n])
!
            else
!
!
               call this%davidson%add_new_trial(residual, n)
!
            endif
!
         endif
!
      enddo
!
!     Special case when residuals are below lindep_threshold,
!     but energy has not converged
!
      if (all(residual_lt_lindep) .and. .not. all(converged)) then
!
         converged = .true.
         call output%printf('m', 'Note: all residuals are converged and &
                            &have norms that are smaller than linear &
                            &dependence threshold. ' // 'Energy &
                            &convergence therefore not tested in this &
                            &calculation.')
!
      endif
!
      call mem%dealloc(residual, this%davidson%n_parameters)
      call mem%dealloc(solution, this%davidson%n_parameters)
!
   end subroutine test_convergence_and_add_trials
!
!
   pure function get_n_solutions_eigen_davidson_solver(this) result(n_solutions)
!!
!!    Get n solutions
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(eigen_davidson_solver), intent(in) :: this
      integer :: n_solutions
!
      n_solutions = this%n_solutions
!
   end function get_n_solutions_eigen_davidson_solver
!
!
   end module eigen_davidson_solver_class
