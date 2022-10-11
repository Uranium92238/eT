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
module davidson_cc_linear_equations_class
!
!!
!!    Davidson coupled cluster linear equations solver class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Josefine H. Andersen, 2018-2019
!!
!!    Uses the Davidson algorithm to solve the sets of linear equations
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
!!    The solver builds a shared reduced space for the set of linear equations,
!!    and stores the solution(s) on file(s) provided on input.
!!
!!    This solver was set up by Eirik F. Kjønstad, Nov 2019, and is a simple generalization
!!    of the ground state multipliers Davidson solver (authors Eirik F. Kjønstad, Sarai D. Folkestad)
!!    as well as the solvers for response made by Josefine H. Andersen, spring 2019.
!!
!
   use parameters
   use ccs_class, only: ccs
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use linear_davidson_tool_class, only: linear_davidson_tool
!
   implicit none
!
   type :: davidson_cc_linear_equations
!
      character(len=100) :: name_   = 'Davidson coupled cluster linear equations solver'
!
      character(len=500) :: description = 'A Davidson solver that solves linear equations &
                                          &involving the Jacobian transformation. When solving &
                                          &multiple linear equations, it builds a shared reduced &
                                          &space in which to express each of the solutions to the &
                                          &equations.'
!
      character(len=500) :: eq_description ! provided by the programmer and printed to the user
!
      integer :: max_iterations, max_dim_red
!
      real(dp) :: residual_threshold
!
      character(len=200) :: storage
      character(len=200) :: section
!
      logical :: restart, records_in_memory
!
      type(timings) :: timer
!
      type(linear_davidson_tool), allocatable :: davidson
!
   contains
!
      procedure :: cleanup                         => cleanup_davidson_cc_linear_equations
!
      procedure :: print_banner                    => print_banner_davidson_cc_linear_equations
      procedure :: print_settings                  => print_settings_davidson_cc_linear_equations
!
      procedure :: read_settings                   => read_settings_davidson_cc_linear_equations
!
      procedure :: run                             => run_davidson_cc_linear_equations
!
      procedure :: set_precondition_vector => set_precondition_vector_davidson_cc_linear_equations
!
   end type davidson_cc_linear_equations
!
!
   interface davidson_cc_linear_equations
!
      procedure :: new_davidson_cc_linear_equations
!
   end interface davidson_cc_linear_equations
!
!
contains
!
!
   function new_davidson_cc_linear_equations(wf, section, eq_description, n_frequencies, n_rhs) result(solver)
!!
!!    New Davidson CC linear equations
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    wf:               The coupled cluster wavefunction to use for the Jacobian transformation
!!
!!    section:          Thresholds will be read from input in "solver section", i.e. if
!!                      section == "cc response", then thresholds are read from the section
!!                      in the input named "solver cc response".
!!
!!    eq_description:   Equation description. String describing the particular problem the
!!                      solver is meant to solve. This is printed to output so that the user
!!                      can tell which kinds of linear equations are solved.
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      type(davidson_cc_linear_equations) :: solver
!
      class(ccs) :: wf
!
      character(len=*), intent(in) :: section
      character(len=*), intent(in) :: eq_description !
      integer,          intent(in) :: n_rhs
      integer,          intent(in) :: n_frequencies
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' linear equations solver')
      call solver%timer%turn_on()
!
!     Set default settings
!
      solver%max_iterations      = 100
      solver%residual_threshold  = 1.0d-3
      solver%restart             = .false.
      solver%max_dim_red         = 100
      solver%records_in_memory   = .false.
      solver%storage             = 'disk'
      solver%section             = trim(section)
      solver%eq_description      = trim(eq_description)
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Read non-default settings and print settings to output
!
      call solver%read_settings()
      call solver%print_settings()
!
!     Determine whether to store records in memory or on file
!
      if (trim(solver%storage) == 'memory') then
!
         solver%records_in_memory = .true.
!
      elseif (trim(solver%storage) == 'disk') then
!
         solver%records_in_memory = .false.
!
      else
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                 trim(solver%storage))
!
      endif
!
!     :: Initialize solver tool and set preconditioner and start vectors ::
!
      solver%davidson = linear_davidson_tool('davidson_t_response',                 &
                                       wf%n_es_amplitudes,                          &
                                       min(1.0d-11, solver%residual_threshold),     &
                                       solver%max_dim_red, n_rhs,                   &
                                       n_frequencies, solver%records_in_memory)
!
   end function new_davidson_cc_linear_equations
!
!
   subroutine print_settings_davidson_cc_linear_equations(solver)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(davidson_cc_linear_equations) :: solver
!
      call output%printf('m', '- Davidson CC linear equations solver settings:', fs='(/t3,a)')
!
      call output%printf('m', 'Residual threshold:       (e9.2)', &
                         reals=[solver%residual_threshold], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations: (i9)', &
                         ints=[solver%max_iterations], fs='(t6,a)')
!
   end subroutine print_settings_davidson_cc_linear_equations
!
!
   subroutine run_davidson_cc_linear_equations(solver, wf, G,  &
                                 frequencies, solution_files,  &
                                 transformation)
!!
!!    Run
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Josefine H. Andersen, 2018-2019
!!
!!    Solves the equations
!!
!!       (A - omega_k I) X_k = G_k,    or     (A^T - omega_k I) X_k = G_k,
!!
!!    where A is the coupled cluster Jacobian matrix and omega_k, G_k is specified
!!    by the programmer. The omega_k are the frequencies. G is referred to as the
!!    right-hand-side(s).
!!
!!    List of arguments:
!!
!!       G:                The right-hand-side vector or vectors. If there are multiple
!!                         right-hand-sides, then G_k is placed in column k of the array G.
!!
!!       n_rhs:            Number of columns of G.
!!
!!       frequencies:      The vector containing the frequencies omega_k.
!!
!!       n_frequencies:    k = 1, 2, ..., n_frequencies
!!
!!       solution_files:   Stream file array to store the solutions X_k
!!
!!       transformation:   Whether to use A or A^T ("right" or "left", respectively).
!!
      use stream_file_class, only: stream_file
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(davidson_cc_linear_equations) :: solver
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes, solver%davidson%n_rhs), intent(in) :: G
!
      real(dp), dimension(solver%davidson%n_solutions), intent(in) :: frequencies
!
      character(len=*), intent(in) :: transformation
!
      type(stream_file), dimension(solver%davidson%n_solutions) :: solution_files
!
      logical :: converged_residual
!
      real(dp), dimension(:), allocatable :: c, X, solution
!
      integer :: iteration, root, trial
!
      real(dp) :: residual_norm
!
      call solver%davidson%initialize()
      call solver%davidson%set_rhs(G)
      call solver%davidson%set_frequencies(frequencies)
!
      call solver%set_precondition_vector(wf)
!
      call solver%davidson%set_trials_to_preconditioner_guess()
!
!     :: Iterative loop ::
!
      iteration = 0
      converged_residual = .false.
!
      call output%printf('n', '- Entering iterative loop to solve equations.', fs='(/t3,a)')
!
      do while (.not. converged_residual .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1
!
!        Reduced space preparations
!
         if (solver%davidson%red_dim_exceeds_max()) call solver%davidson%set_trials_to_solutions()
!
         call solver%davidson%update_reduced_dim()
!
         call solver%davidson%orthonormalize_trial_vecs()
!
         call output%printf('n', 'Iteration:               (i0)', &
                            ints=[iteration], fs='(/t3,a)')
         call output%printf('n', 'Reduced space dimension: (i0)', ints=[solver%davidson%dim_red])
!
!        Transform new trial vectors
!
         call mem%alloc(X, wf%n_es_amplitudes)
         call mem%alloc(c, wf%n_es_amplitudes)
!
         do trial = solver%davidson%first_new_trial(), solver%davidson%last_new_trial()
!
            call solver%davidson%get_trial(c, trial)
!
            call wf%construct_Jacobian_transform(trim(transformation), c, X)
!
            call solver%davidson%set_transform(X, trial)
!
         enddo
!
         call mem%dealloc(c, wf%n_es_amplitudes)
!
!        Solve linear equation(s) in reduced space
!
         call solver%davidson%solve_reduced_problem()
!
!        Loop over roots and check residuals,
!        then generate new trial vectors for roots not yet converged
!
         call output%printf('n', 'Frequency      Residual norm', fs='(/t3,a)')
         call output%print_separator('n', 30, '-')
!
         converged_residual = .true.
!
         do root = 1, solver%davidson%n_solutions
!
            call solver%davidson%construct_residual(X, root)
!
            residual_norm = get_l2_norm(X, wf%n_es_amplitudes)
!
            if (residual_norm > solver%residual_threshold) then
!
               converged_residual = .false.
               call solver%davidson%add_new_trial(X, root)
!
            endif
!
            call output%printf('n', '(e9.2)     (e9.2)', &
                               reals=[frequencies(root), residual_norm])
!
         enddo
!
         call mem%dealloc(X, wf%n_es_amplitudes)
!
         call output%print_separator('n', 30, '-')
!
      enddo ! End of iterative loop
!
!     :: Summary of calculation ::
!
      if (converged_residual) then
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
         call mem%alloc(solution, wf%n_es_amplitudes)
!
         do root = 1, solver%davidson%n_solutions
!
            call solver%davidson%construct_solution(solution, root)
!
            call solution_files(root)%open_('rewind')
            call solution_files(root)%write_(solution, wf%n_es_amplitudes)
            call solution_files(root)%close_()
!
         enddo
!
         call mem%dealloc(solution, wf%n_es_amplitudes)
!
      else
!
         call output%error_msg('was not able to converge the equations in the given ' //&
                                 'number of iterations in run_davidson_cc_linear_equations.')
!
      endif
!
      call solver%davidson%cleanup()
!
   end subroutine run_davidson_cc_linear_equations
!
!
   subroutine cleanup_davidson_cc_linear_equations(solver, wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(ccs) :: wf
      class(davidson_cc_linear_equations) :: solver
!
      call solver%timer%turn_off()
!
      call output%printf('m', '- Finished solving the ' //  &
                         trim(convert_to_uppercase(wf%name_)) // ' linear equations', &
                         fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec):  (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6,a)')
!
      call output%printf('m', 'Total cpu time (sec):   (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_davidson_cc_linear_equations
!
!
   subroutine print_banner_davidson_cc_linear_equations(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_linear_equations) :: solver
!
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('m', solver%description, ffs='(t3,a)', fs='(t3,a)')
      call output%printf('m', solver%eq_description, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_davidson_cc_linear_equations
!
!
   subroutine set_precondition_vector_davidson_cc_linear_equations(solver, wf)
!!
!!    Set precondition vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets precondition vector to orbital differences
!!
      implicit none
!
      class(davidson_cc_linear_equations) :: solver
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_es_amplitudes)
      call wf%get_orbital_differences(preconditioner, wf%n_es_amplitudes)
      call solver%davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_es_amplitudes)
!
   end subroutine set_precondition_vector_davidson_cc_linear_equations
!
!
   subroutine read_settings_davidson_cc_linear_equations(solver)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      use global_in, only: input
!
      implicit none
!
      class(davidson_cc_linear_equations) :: solver
!
      call input%get_keyword('threshold',                                  &
                                        'solver ' // trim(solver%section), &
                                        solver%residual_threshold)
!
      call input%get_keyword('max iterations',                             &
                                        'solver ' // trim(solver%section), &
                                        solver%max_iterations)
!
      call input%get_keyword('storage',                                    &
                                        'solver ' // trim(solver%section), &
                                         solver%storage)
!
   end subroutine read_settings_davidson_cc_linear_equations
!
!
end module davidson_cc_linear_equations_class
