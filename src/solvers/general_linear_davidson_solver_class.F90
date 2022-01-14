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
module general_linear_davidson_solver_class
!
!!
!!    General linear Davidson solver class
!!    Written by Regina Matveeva and Ida-Marie Høyvik, Sept 2021
!!
!!    Based on eigen_davidson_class by Sarai D. Folkestad
!!
!!    This solver uses the Davidson algorithm to solve the sets of linear equations
!!
!!       (A - omega_k I) X_k = G_k,    or     (A^T - omega_k I) X_k = G_k,
!!
!!    where omega_k is a set of real numbers referred to as "frequencies",
!!    and the G_k can either be one vector, meaning that we solve
!!
!!       (A - omega_k I) X_k = G,      or     (A^T - omega_k I) X_k = G,
!!
!!    or there can be different right-hand-sides, G_k, for different frequencies, omega_k.
!!
!!    The solver builds a shared reduced space for the set of linear equations.
!!
!!    This solver is a generalization of the ground state multipliers Davidson
!!    solver (authors Eirik F. Kjønstad, Sarai D. Folkestad) as well as the
!!    Davidson coupled cluster linear equations solver (authors Eirik F. Kjønstad,
!!    Sarai D. Folkestad, Josefine H. Andersen)
!!
!
   use kinds
   use hf_class
   use linear_davidson_tool_class
!
   use transformation_tool_class,              only: transformation_tool
   use linear_storage_tool_class,              only: linear_storage_tool
   use start_vector_tool_class,                only: start_vector_tool
   use preconditioner_getter_class,            only: preconditioner_getter
   use rhs_linear_equation_tool_class,         only: rhs_linear_equation_tool
   use linear_davidson_print_tool_class,       only: linear_davidson_print_tool
!
   implicit none
!
   type :: general_linear_davidson_solver
!
      integer  :: n_solutions
      integer  :: n_rhs
      integer, private  :: max_iterations
      real(dp), private :: residual_threshold
!
      class(linear_davidson_tool), allocatable, private :: davidson
!
      class(transformation_tool),         allocatable          :: transformer
      class(linear_storage_tool),         allocatable, private :: storer
      class(start_vector_tool),           allocatable, private :: start_vector
      class(preconditioner_getter),       allocatable, private :: preconditioner
      class(rhs_linear_equation_tool),    allocatable, private :: rhs_getter
      class(linear_davidson_print_tool),  allocatable, private :: printer
!
   contains
!
      procedure, public :: run => run_general_linear_davidson_solver
      procedure, public :: get_n_solutions => get_n_solutions_general_linear_davidson_solver
!
      procedure, private :: test_convergence_and_add_trials
!
end type general_linear_davidson_solver
!
   interface general_linear_davidson_solver
!
      procedure :: new_general_linear_davidson_solver
!
end interface general_linear_davidson_solver
!
contains
!
   function new_general_linear_davidson_solver(transformer,         &
                                               davidson,            &
                                               storer,              &
                                               start_vector,        &
                                               preconditioner,      &
                                               rhs_getter,          &
                                               printer,             &
                                               n_rhs,               &
                                               n_solutions,         &
                                               max_iterations,      &
                                               residual_threshold) result(this)
!!
!!    New general linear Davidson
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      type(general_linear_davidson_solver) :: this
!
      class(transformation_tool),        intent(in) :: transformer
      type(linear_davidson_tool),        intent(in) :: davidson
      class(linear_storage_tool),        intent(in) :: storer
      class(start_vector_tool),          intent(in) :: start_vector
      class(preconditioner_getter),      intent(in) :: preconditioner
      class(rhs_linear_equation_tool),   intent(in) :: rhs_getter
      class(linear_davidson_print_tool), intent(in) :: printer
      integer,                           intent(in) :: max_iterations, n_solutions, n_rhs
      real(dp),                          intent(in) :: residual_threshold
!

      this%transformer         = transformer
      this%davidson            = davidson
      this%storer              = storer
      this%start_vector        = start_vector
      this%preconditioner      = preconditioner
      this%rhs_getter          = rhs_getter
      this%printer             = printer
!
      this%n_solutions         = n_solutions
      this%n_rhs               = n_rhs
      this%max_iterations      = max_iterations
      this%residual_threshold  = residual_threshold
!
   end function new_general_linear_davidson_solver
!
!
   subroutine run_general_linear_davidson_solver(this, frequencies)
!!
!!    Run
!!    Written by Regina Matveeva and Ida-Marie Høyvik, Sept 2021
!!
      implicit none
!
      class(general_linear_davidson_solver) :: this
!
      real(dp), dimension(this%n_solutions), intent(in) :: frequencies
!
      logical, dimension(:), allocatable  :: converged
!
      integer  :: n, iteration, trial
!
      real(dp), dimension(:),    allocatable     :: residual_norm
      real(dp), dimension(:),    allocatable     :: c, residual
      real(dp), dimension(:),    allocatable     :: solution
      real(dp), dimension(:,:),  allocatable     :: rhs
      real(dp) :: dummy = zero
!
      call mem%alloc(rhs, this%davidson%n_parameters, this%n_rhs)
!
      call this%rhs_getter%get(rhs)
!
      call mem%alloc(converged, this%n_rhs)
      converged            = .false.
!
      do n = 1, this%n_rhs
!
         if (get_l2_norm(rhs(:, n), this%davidson%n_parameters*this%n_rhs) < this%residual_threshold) then
!
            converged(n) = .true.
!
         end if
!
      end do
!
      if (all(converged)) then
!
         call output%printf('m', 'Right hand side is zero to within threshold (e6.1).', &
                            reals=[this%residual_threshold], fs='(/t3,a)')
!
         call mem%alloc(solution, this%davidson%n_parameters)
!
         call zero_array(solution, this%davidson%n_parameters)
!
         call this%storer%store(solution, 0)
!
         call mem%dealloc(solution, this%davidson%n_parameters)
         call mem%dealloc(converged, this%n_rhs)
         call mem%dealloc(rhs, this%davidson%n_parameters, this%n_rhs)
!
         return
!
      else if (any(converged)) then
!
         call output%error_msg('Some of right hand side vectors are zero &
                              &to within threshold (e6.1).', &
                              reals=[this%residual_threshold], fs='(/t3,a)')
!
      end if
!
      call mem%dealloc(converged, this%n_rhs)
!
      call this%davidson%initialize()
!
      call this%davidson%set_rhs(rhs)
      call this%davidson%set_frequencies(frequencies)
!
      call mem%alloc(c, this%davidson%n_parameters)
      call this%preconditioner%get(c)
      call this%davidson%set_preconditioner(c)
!
      do trial = 1, this%n_solutions
         call this%start_vector%get(c, trial)
         call this%davidson%set_trial(c, trial)
      end do
!
      call this%printer%print_settings(this%residual_threshold, this%max_iterations)
!
      call this%transformer%initialize()
!
!     :: Iterative loop ::
!
      call mem%alloc(converged, this%n_solutions)
      call mem%alloc(residual_norm, this%n_solutions)
      call mem%alloc(residual, this%davidson%n_parameters)
!
      converged = .false.
!
      iteration = 0
!
      call this%printer%print_iteration_header()
!
      do while (.not. all(converged) .and. (iteration .le. this%max_iterations))
!
         iteration = iteration + 1
!
!        Reduced space preparations
!
         call this%davidson%update_reduced_space()
!
!        Transform new trial vectors
!
         do trial = this%davidson%first_new_trial(), this%davidson%last_new_trial()
!
            call this%davidson%get_trial(c, trial)
            call this%transformer%transform(c, residual, dummy)
            call this%davidson%set_transform(residual, trial)
!
         enddo
!
!        Solve linear equation(s) in reduced space
!
         call this%davidson%solve_reduced_problem()
!
!        Loop over roots and check residuals,
!        then generate new trial vectors for roots not yet converged
!
         call this%test_convergence_and_add_trials(residual_norm, converged)
!
         call this%printer%print_iteration(this%n_solutions,     &
                                           residual_norm,        &
                                           iteration,            &
                                           this%davidson%dim_red)
!
         call mem%alloc(solution, this%davidson%n_parameters)
!
         do n = 1, this%n_solutions

            call this%davidson%construct_solution(solution, n)
!
            call this%storer%store(solution, n)

         end do
!
         call mem%dealloc(solution, this%davidson%n_parameters)
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
!
      call this%davidson%cleanup()
      call mem%dealloc(rhs, this%davidson%n_parameters, this%n_rhs)
!
   end subroutine run_general_linear_davidson_solver
!
!
   subroutine test_convergence_and_add_trials(this,                     &
                                              residual_norm,            &
                                              converged)
!!
!!    Test convergence and add trial
!!    Written by Regina Matveeva, Sept 2021
!!
!!    Based on test_convergence_and_add_trials in general_eigen_davison_solver by
!!    Sarai D. Folkestad
!!
!!    - Get residual
!!    - Test convergence
!!    - Add new trial if needed and residual norm is greater
!!      than linear dependence threshold
!!
      implicit none
!
      class(general_linear_davidson_solver), intent(inout) :: this
      real(dp), dimension(this%n_solutions), intent(out)   :: residual_norm
      logical, dimension(this%n_solutions), intent(out)    :: converged
!
      real(dp), dimension(:), allocatable :: residual
!
      integer :: n
!
      call mem%alloc(residual, this%davidson%n_parameters)
!
      converged = .true.
!
      do n = 1, this%n_solutions
!
         call this%davidson%construct_residual(residual, n)
!
         residual_norm(n) = get_l2_norm(residual, this%davidson%n_parameters)
!
         if (residual_norm(n) > this%residual_threshold) then
!
            converged(n) = .false.
!
            call this%davidson%add_new_trial(residual, n)
!
         endif
!
      enddo
!
      call mem%dealloc(residual, this%davidson%n_parameters)
!
!
   end subroutine test_convergence_and_add_trials
!
!
   pure function get_n_solutions_general_linear_davidson_solver(this) result(n_solutions)
!!
!!    Get n solutions
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(general_linear_davidson_solver), intent(in) :: this
      integer :: n_solutions
!
      n_solutions = this%n_solutions
!
   end function get_n_solutions_general_linear_davidson_solver
!
!
end module general_linear_davidson_solver_class
