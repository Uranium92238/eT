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
module davidson_cc_es_class
!
!!
!! Davidson coupled cluster excited state solver class module
!! Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2019
!!
!! Solves the CC excited state eigenvalue equation
!!
!!    A R = omega R   or  L^T A = omega L^T
!!
!! for a set of right or left states (R, L) and the associated
!! excitation energies omega. Here, A is the coupled cluster
!! Jacobian matrix:
!!
!!       A_mu,nu = < mu | [H-bar, tau_nu] | HF >,    H-bar = e-T H eT.
!!
!! The solutions are determined using the Davidson
!! reduced space algorithm, where the eigenvalue problem
!! is solved in a subspace generated from the residuals* obtained
!! in the preceding iterations. See E. R. Davidson, J. Comput. Phys.
!! 17, 87 (1975) for more details.
!!
!! * More precisely, using preconditioned residuals based on the
!!   orbital differences approximation of A,
!!
!!    A_mu,nu = < mu | [F, tau_nu] | HF > = epsilon_mu delta_mu,nu
!!
!!   Here,
!!
!!               epsilon_ai = orbital_energy(a) - orbital_energy(i)
!!
!!             epsilon_aibj = orbital_energy(a) + orbital_energy(b)
!!                          - orbital_energy(i) - orbital_energy(j),
!!
!!   and so on.
!!
!
   use kinds
   use ccs_class
   use eigen_davidson_tool_class
   use abstract_cc_es_class, only: abstract_cc_es
!
   use es_projection_tool_class
   use es_valence_projection_tool_class, only: es_valence_projection_tool
   use es_cvs_projection_tool_class, only: es_cvs_projection_tool
   use es_ip_projection_tool_class, only: es_ip_projection_tool
!
   use convergence_tool_class, only: convergence_tool
!
   implicit none
!
   type, extends(abstract_cc_es) :: davidson_cc_es
!
      integer :: max_dim_red
      type(eigen_davidson_tool), allocatable :: davidson
!
   contains
!
      procedure, non_overridable :: run => run_davidson_cc_es
!
      procedure :: set_precondition_vector  &
                => set_precondition_vector_davidson_cc_es
!
      procedure :: read_settings &
                => read_settings_davidson_cc_es
!
   end type davidson_cc_es
!
!
   interface davidson_cc_es
!
      procedure :: new_davidson_cc_es
!
   end interface davidson_cc_es
!
!
contains
!
!
   function new_davidson_cc_es(transformation, wf, restart) result(solver)
!!
!!    New Davidson CC ES
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Modified by Andreas S. Skeidsvoll, Oct 2020
!!    Added max_dim_red and n_singlet_states consistency checks.
!!
      implicit none
!
      type(davidson_cc_es) :: solver
!
      class(ccs),       intent(inout)  :: wf
      character(len=*), intent(in)     :: transformation
      logical,          intent(in)     :: restart
!
      real(dp) :: lindep_threshold
      integer  :: max_dim_red
      logical  :: records_in_memory
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' excited state (' &
                     // trim(transformation) //')')
      call solver%timer%turn_on()
!
      solver%name_ = 'Davidson coupled cluster excited state solver'
      solver%tag   = 'Davidson'
!
      solver%description1 = 'A Davidson solver that calculates the lowest eigenvalues and &
               &the right or left eigenvectors of the Jacobian matrix, A. The eigenvalue &
               &problem is solved in a reduced space, the dimension of which is &
               &expanded until the convergence criteria are met.'
!
      solver%description2 = 'A complete description of the algorithm can be found in &
                            &E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%n_singlet_states     = 0
      solver%max_iterations       = 100
      solver%restart              = restart
      solver%transformation       = trim(transformation)
      solver%es_type              = 'valence'
!
      records_in_memory           = .false.
      max_dim_red                 = max(100, 10*solver%n_singlet_states)
!
!     Initialize convergence checker with default threshols
!
      solver%convergence_checker = convergence_tool(energy_threshold   = 1.0d-3,   &
                                                    residual_threshold = 1.0d-3,   &
                                                    energy_convergence = .false.)
!
      call solver%read_settings(max_dim_red, records_in_memory)
!
      call solver%print_es_settings()
!
      wf%n_singlet_states = solver%n_singlet_states
!
      lindep_threshold = min(1.0d-11, 0.1d0*solver%convergence_checker%residual_threshold)
      solver%davidson = eigen_davidson_tool('cc_es_davidson',        &
                                             wf%n_es_amplitudes,     &
                                             solver%n_singlet_states,&
                                             lindep_threshold,       &
                                             max_dim_red,            &
                                             records_in_memory)
!
   end function new_davidson_cc_es
!
!
   subroutine run_davidson_cc_es(solver, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_es) :: solver
!
      class(ccs) :: wf
!
      logical, dimension(:), allocatable :: converged
      logical, dimension(:), allocatable :: residual_lt_lindep
!
      integer :: iteration, trial, n, state
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:), allocatable :: c
      real(dp), dimension(:,:), allocatable :: r
!
      real(dp), dimension(:), allocatable :: residual
      real(dp), dimension(:), allocatable :: solution
!
      type(timings) :: construct_new_trial, iteration_time
      character(len=30) :: timer_name
!
      construct_new_trial = timings("Davidson: construct new trial", pl="v")
!
!     :: Preparations ::
!
      call solver%initialize_start_vector_tool(wf)
      call solver%initialize_projection_tool(wf)
!
      call solver%prepare_wf_for_excited_state(wf)
!
      call solver%davidson%initialize()
!
      call solver%set_precondition_vector(wf)
!
!     :: Set start vectors and handle restart ::
!
      call solver%initialize_energies()
!
      call mem%alloc(c, wf%n_es_amplitudes)
!
      do trial = 1, solver%n_singlet_states
!
         call solver%set_initial_guesses(wf, c, trial, trial)
         call solver%davidson%set_trial(c, trial)
!
      enddo
!
      call mem%dealloc(c, wf%n_es_amplitudes)
!
!     :: Iterative loop ::
!
      call mem%alloc(converged, solver%n_singlet_states)
      call mem%alloc(residual_lt_lindep, solver%n_singlet_states)
!
      converged          = .false.
      residual_lt_lindep = .false.
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1
!
         write(timer_name, '(a,i0)') "Davidson: CC ES iteration ", iteration
         iteration_time = timings(trim(timer_name), pl="n")
         call iteration_time%turn_on()
!
!        Reduced space preparations
!
         if (solver%davidson%red_dim_exceeds_max()) call solver%davidson%set_trials_to_solutions()
!
         call solver%davidson%update_reduced_dim()
!
         call solver%davidson%orthonormalize_trial_vecs()
!
!        Print iteration information
!
         call output%printf('n', 'Iteration:               (i4)', &
                            ints=[iteration], fs='(/t3,a)')
         call output%printf('n', 'Reduced space dimension: (i4)', ints=[solver%davidson%dim_red])
!
!        Transform new trial vectors and write them to file
!
         call mem%alloc(c, wf%n_es_amplitudes)
         call mem%alloc(residual, wf%n_es_amplitudes)
!
         do trial = solver%davidson%first_new_trial(), solver%davidson%last_new_trial()
!
            call solver%davidson%get_trial(c, trial)
!
            call wf%construct_jacobian_transform(solver%transformation, c, residual)
!
            if (solver%projector%active) call solver%projector%do_(residual)
!
            call solver%davidson%set_transform(residual, trial)
!
         enddo ! Done transforming new trials
!
         call mem%dealloc(c, wf%n_es_amplitudes)
!
!        Solve problem in the reduced space
!
         call solver%davidson%solve_reduced_problem()
!
!        Construct the full space solutions, test for convergence,
!        and use the residuals to construct the next trial vectors
!
         call mem%alloc(solution, wf%n_es_amplitudes)
!
         call output%printf('n', ' Root  Eigenvalue (Re)   Eigenvalue &
                            &(Im)    |residual|   Delta E (Re)', fs='(/t3,a)', ll=120)
         call output%print_separator('n', 72,'-')
!
         call construct_new_trial%turn_on()
!
         do n = 1, solver%n_singlet_states
!
            call solver%davidson%construct_solution(solution, n)
!
            call wf%save_excited_state(solution, n, n, solver%transformation, solver%energies(n))
!
            call solver%davidson%construct_residual(residual, solution, n)
!
            if (solver%projector%active) call solver%projector%do_(residual) ! CVS projection,
                                                                             ! for instance
!
            residual_norm = get_l2_norm(residual, wf%n_es_amplitudes)
!
            converged(n) = solver%convergence_checker%has_converged(residual_norm, &
                                 solver%davidson%omega_re(n)-solver%energies(n), iteration)
!
            residual_lt_lindep(n) = (residual_norm <= solver%davidson%lindep_threshold)
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
                  call solver%davidson%add_new_trial(residual, n)
!
               endif
!
            endif
!
            call output%printf('n', '(i4) (f16.12)  (f16.12)    (e11.4)  (e11.4) ',             &
                               ints=[n], ll=120,                                                &
                               reals=[solver%davidson%omega_re(n),                              &
                                      solver%davidson%omega_im(n),                              &
                                      residual_norm,                                            &
                                      abs(solver%davidson%omega_re(n) - solver%energies(n))])
!
         enddo ! Done constructing new trials from residuals
!
         call construct_new_trial%turn_off()
         call construct_new_trial%reset()
!
         call output%print_separator('n', 72,'-')
!
         call mem%dealloc(residual, wf%n_es_amplitudes)
         call mem%dealloc(solution, wf%n_es_amplitudes)
!
!        Update energies and save them
!
         solver%energies = solver%davidson%omega_re
!
!        Special case when residuals converge in first iteration, e.g. on restart
!        => Exit without testing energy convergence
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
         call iteration_time%turn_off()
         call iteration_time%reset()
!
      enddo ! End of iterative loop
!
!     :: Calculation summary ::
!
      if (.not. all(converged)) then
!
         call output%error_msg("Did not converge in the max number of &
                               &iterations in run_davidson_cc_es.")
!
      else
!
         call wf%set_excitation_energies(solver%energies, solver%transformation)
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
         call mem%alloc(r, wf%n_es_amplitudes, solver%n_singlet_states)
!
         do state = 1, solver%n_singlet_states
!
            call solver%davidson%construct_solution(r(:,state), state)
!
         enddo
!
         call solver%print_summary(wf, r)
!
         call mem%dealloc(r, wf%n_es_amplitudes, solver%n_singlet_states)
!
      endif
!
      call solver%davidson%cleanup()
!
      call mem%dealloc(converged, solver%n_singlet_states)
      call mem%dealloc(residual_lt_lindep, solver%n_singlet_states)
!
   end subroutine run_davidson_cc_es
!
!
   subroutine set_precondition_vector_davidson_cc_es(solver, wf)
!!
!!    Set precondition vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets precondition vector to orbital differences.
!!
      implicit none
!
      class(davidson_cc_es) :: solver
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_es_amplitudes)
      call wf%get_orbital_differences(preconditioner, wf%n_es_amplitudes)
      call solver%davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_es_amplitudes)
!
   end subroutine set_precondition_vector_davidson_cc_es
!
!
   subroutine read_settings_davidson_cc_es(solver, max_dim_red, records_in_memory)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      implicit none
!
      class(davidson_cc_es) :: solver
      integer, intent(inout) :: max_dim_red
      logical, intent(inout) :: records_in_memory
!
      call solver%read_es_settings(records_in_memory)
      call input%get_keyword('max reduced dimension', 'solver cc es', max_dim_red)
!
   end subroutine read_settings_davidson_cc_es
!
!
end module davidson_cc_es_class
