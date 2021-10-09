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
module general_eigen_davidson_class
!
!!
!! General Davidson solver class
!! Written by Sarai D. Folkestad, May 2021
!!
!
   use kinds
   use hf_class
   use eigen_davidson_tool_class
!
   use convergence_tool_class,      only: convergence_tool
   use transformation_tool_class,   only: transformation_tool
   use eigen_storage_tool_class,    only: eigen_storage_tool
   use start_vector_tool_class,     only: start_vector_tool
   use preconditioner_getter_class, only: preconditioner_getter
!
   implicit none
!
   type :: general_eigen_davidson
!
      integer, private :: n_solutions
      integer, private :: max_iterations
!
      class(eigen_davidson_tool), allocatable, private :: davidson
!
      class(convergence_tool),     allocatable, private :: convergence_checker
      class(transformation_tool),  allocatable, private :: transformer
      class(eigen_storage_tool),   allocatable, private :: storer
      class(start_vector_tool),    allocatable, private :: start_vector
      class(preconditioner_getter),allocatable, private :: preconditioner
!
   contains
!
      procedure, public :: run => run_general_eigen_davidson
      procedure, public :: get_n_solutions => get_n_solutions_general_eigen_davidson
!
      procedure, private :: print_iteration
      procedure, private :: print_summary
      procedure, private :: print_settings
      procedure, private :: test_convergence_and_add_trials
!
   end type general_eigen_davidson
!
   interface general_eigen_davidson
!
      procedure :: new_general_eigen_davidson
!
   end interface general_eigen_davidson
!
contains
!
   function new_general_eigen_davidson(transformer,         &
                                       davidson,            &
                                       convergence_checker, &
                                       storer,              &
                                       start_vector,        &
                                       preconditioner,      &
                                       n_solutions,         &
                                       max_iterations) result(this)
!!
!!    New general Davidson
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      type(general_eigen_davidson) :: this
!
      class(transformation_tool),   intent(in) :: transformer
      type(eigen_davidson_tool),    intent(in) :: davidson
      type(convergence_tool),       intent(in) :: convergence_checker
      class(eigen_storage_tool),    intent(in) :: storer
      class(start_vector_tool),     intent(in) :: start_vector
      class(preconditioner_getter), intent(in) :: preconditioner
      integer,                      intent(in) :: max_iterations, n_solutions
!

      this%davidson            = davidson
      this%transformer         = transformer
      this%convergence_checker = convergence_checker
      this%storer              = storer
      this%start_vector        = start_vector
      this%preconditioner      = preconditioner
      this%n_solutions         = n_solutions
      this%max_iterations      = max_iterations
!
   end function new_general_eigen_davidson
!
!
   subroutine run_general_eigen_davidson(this)
!!
!!    Run
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(general_eigen_davidson) :: this

      logical, dimension(:), allocatable :: converged
      logical, dimension(:), allocatable :: residual_lt_lindep
!
      integer :: iteration, trial
!
      real(dp), dimension(:), allocatable :: residual_norm
      real(dp), dimension(:), allocatable :: c, residual
      real(dp), dimension(:), allocatable :: omega
!
      call this%davidson%initialize()
!
      call mem%alloc(c, this%davidson%n_parameters)
      call this%preconditioner%get(c)
      call this%davidson%set_preconditioner(c)
!
      do trial = 1,this%n_solutions
!
         call this%start_vector%get(c, trial)
         call this%davidson%set_trial(c, trial)
!
      enddo
!
      call this%print_settings()
!
      call mem%alloc(omega, this%n_solutions)
      omega = zero
!
      call mem%alloc(converged, this%n_solutions)
      call mem%alloc(residual_lt_lindep, this%n_solutions)
!
      converged            = .false.
      residual_lt_lindep   = .false.
!
      call mem%alloc(residual_norm, this%n_solutions)
      call mem%alloc(residual, this%davidson%n_parameters)
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. this%max_iterations))
!
         iteration = iteration + 1
!
         call this%davidson%update_reduced_space()
!
         do trial = this%davidson%first_new_trial(), this%davidson%last_new_trial()
!
            call this%davidson%get_trial(c, trial)
            call this%transformer%transform(c, residual)
            call this%davidson%set_transform(residual, trial)
!
         enddo
!
         call this%davidson%solve_reduced_problem()
!
         call this%test_convergence_and_add_trials(omega,               &
                                                   residual_norm,       &
                                                   converged,           &
                                                   residual_lt_lindep,  &
                                                   iteration)
!
         call this%print_iteration(omega, residual_norm, iteration)
         omega = this%davidson%get_omega_re(this%n_solutions)
!
      enddo
!
      call mem%dealloc(c, this%davidson%n_parameters)
      call mem%dealloc(residual, this%davidson%n_parameters)
      call mem%dealloc(residual_norm, this%n_solutions)
!
      call this%print_summary(converged, iteration)
!
      call mem%dealloc(converged, this%n_solutions)
      call mem%dealloc(residual_lt_lindep, this%n_solutions)
      call mem%dealloc(omega, this%n_solutions)
!
      call this%davidson%cleanup()
!
   end subroutine run_general_eigen_davidson
!
!
   subroutine print_iteration(this, omega, residual_norm, iteration)
!!
!!    Print iteration
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none

      class(general_eigen_davidson), intent(in) :: this
      real(dp), dimension(this%n_solutions), intent(in) :: residual_norm, omega
      integer, intent(in) :: iteration
!
      integer :: n
      real(dp), dimension(:), allocatable :: new_omega_re, new_omega_im
!
      call mem%alloc(new_omega_re, this%n_solutions)
      new_omega_re = this%davidson%get_omega_re(this%n_solutions)
!
      call mem%alloc(new_omega_im, this%n_solutions)
      new_omega_im = this%davidson%get_omega_im(this%n_solutions)
!
      call output%printf('n', 'Iteration:               (i4)', &
                         ints=[iteration], fs='(/t3,a)')
      call output%printf('n', 'Reduced space dimension: (i4)', ints=[this%davidson%dim_red])
!
      call output%printf('n', ' Root  omega (Re)        omega &
                            &(Im)         |residual|   Delta omega (Re)', fs='(/t3,a)', ll=120)
!
      call output%print_separator('n', 72,'-')
!
      do n = 1, this%n_solutions
!
        call output%printf('n', '(i4) (f16.12)  (f16.12)    (e11.4)  (e11.4) ',   &
                           ints=[n],                                              &
                           reals=[new_omega_re(n),                                &
                                  new_omega_im(n),                                &
                                  residual_norm(n),                               &
                                  dabs(this%davidson%omega_re(n) - omega(n))],    &
                                  ll=120)
!
      enddo
!
      call output%print_separator('n', 72,'-')
      call mem%dealloc(new_omega_re, this%n_solutions)
      call mem%dealloc(new_omega_im, this%n_solutions)
!
   end subroutine print_iteration
!
!
   subroutine print_summary(this, converged, iteration)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none

      class(general_eigen_davidson), intent(in) :: this
      logical, dimension(this%n_solutions), intent(in) :: converged
      integer, intent(in) :: iteration
!
      if (.not. all(converged)) then
!
         call output%error_msg("Did not converge in the max number of &
                               &iterations in Davidson solver.")
!
      else
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
      endif
!
   end subroutine print_summary
!
!
   subroutine print_settings(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none

      class(general_eigen_davidson), intent(in) :: this
!
      call output%printf('m', '- Davidson solver settings', fs='(/t3,a)')
!
      call output%printf('m', 'Number of singlet states:     (i11)', &
                         ints=[this%n_solutions], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations:     (i11)', &
                         ints=[this%max_iterations], fs='(t6,a)')
!
   end subroutine print_settings
!
!
   subroutine test_convergence_and_add_trials(this,               &
                                              omega,              &
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
      implicit none
!
      class(general_eigen_davidson), intent(inout) :: this
      real(dp), dimension(this%n_solutions), intent(in)  :: omega
      real(dp), dimension(this%n_solutions), intent(out) :: residual_norm
      logical, dimension(this%n_solutions), intent(out)  :: converged, residual_lt_lindep
      integer, intent(in) :: iteration
!
      real(dp), dimension(:), allocatable :: residual
      real(dp), dimension(:), allocatable :: solution
!
      real(dp), dimension(:), allocatable :: new_omega_re
!
      integer :: n
!
      call mem%alloc(residual, this%davidson%n_parameters)
      call mem%alloc(solution, this%davidson%n_parameters)
!
      call mem%alloc(new_omega_re, this%n_solutions)
      new_omega_re = this%davidson%get_omega_re(this%n_solutions)
!
      do n = 1, this%n_solutions
!
         call this%davidson%construct_solution(solution, n)
         call this%davidson%construct_residual(residual, solution, n)
!
         call this%storer%store(new_omega_re(n), solution, n)
!
         residual_norm(n) = get_l2_norm(residual, this%davidson%n_parameters)
         converged(n)     = this%convergence_checker%has_converged(residual_norm(n), &
                            new_omega_re(n)-omega(n), iteration)
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
               call this%davidson%add_new_trial(residual, n)
!
            endif
!
         endif
!
      enddo
!
      call mem%dealloc(residual, this%davidson%n_parameters)
      call mem%dealloc(solution, this%davidson%n_parameters)
      call mem%dealloc(new_omega_re, this%n_solutions)
!
   end subroutine test_convergence_and_add_trials
!
!
   pure function get_n_solutions_general_eigen_davidson(this) result(n_solutions)
!!
!!    Get n solutions
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(general_eigen_davidson), intent(in) :: this
      integer :: n_solutions
!
      n_solutions = this%n_solutions
!
   end function get_n_solutions_general_eigen_davidson
!
!
end module general_eigen_davidson_class
