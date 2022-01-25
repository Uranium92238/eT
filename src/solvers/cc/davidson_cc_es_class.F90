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
   use ccs_class, only: ccs
   use memory_manager_class, only: mem
   use timings_class, only: timings
   use global_out, only: output
   use abstract_cc_es_class, only: abstract_cc_es
   use eigen_davidson_tool_class, only: eigen_davidson_tool
!
   use convergence_tool_class, only: convergence_tool
!
   implicit none
!
   type, extends(abstract_cc_es) :: davidson_cc_es
!
      integer :: max_dim_red
!
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
      procedure :: read_davidson_settings &
                => read_davidson_settings_davidson_cc_es
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
   function new_davidson_cc_es(transformation, wf, restart) result(this)
!!
!!    New Davidson CC ES
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Modified by Andreas S. Skeidsvoll, Oct 2020
!!    Added max_dim_red and n_singlet_states consistency checks.
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      type(davidson_cc_es) :: this
!
      class(ccs),       intent(inout), target  :: wf
      character(len=*), intent(in)     :: transformation
      logical,          intent(in)     :: restart
!
      real(dp) :: lindep_threshold
      logical  :: records_in_memory
!
      this%wf => wf 
!
      this%timer = timings(trim(convert_to_uppercase(this%wf%name_)) // ' excited state (' &
                     // trim(transformation) //')')
      call this%timer%turn_on()
!
      this%name_ = 'Davidson coupled cluster excited state solver'
      this%tag   = 'Davidson'
!
      this%description1 = 'A Davidson solver that calculates the lowest eigenvalues and &
               &the right or left eigenvectors of the Jacobian matrix, A. The eigenvalue &
               &problem is solved in a reduced space, the dimension of which is &
               &expanded until the convergence criteria are met.'
!
      this%description2 = 'A complete description of the algorithm can be found in &
                            &E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call this%print_banner()
!
!     Set defaults
!
      this%n_singlet_states     = 0
      this%max_dim_red          = 0
      this%max_iterations       = 100
      this%restart              = restart
      this%transformation       = trim(transformation)
      this%es_type              = 'valence'
!
      records_in_memory           = .false.
!
!     Initialize convergence checker with default threshols
!
      this%convergence_checker = convergence_tool(energy_threshold   = 1.0d-3,   &
                                                    residual_threshold = 1.0d-3,   &
                                                    energy_convergence = .false.)
!
      call this%read_settings(records_in_memory)
!
      call this%print_es_settings()
!
      this%wf%n_singlet_states = this%n_singlet_states
!
      lindep_threshold = min(1.0d-11, 0.1d0*this%convergence_checker%residual_threshold)
      this%davidson = eigen_davidson_tool('cc_es_davidson',        &
                                            this%wf%n_es_amplitudes,      &
                                            this%n_singlet_states, &
                                            lindep_threshold,        &
                                            this%max_dim_red,      &
                                            records_in_memory)
!
   end function new_davidson_cc_es
!
!
   subroutine run_davidson_cc_es(this)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(davidson_cc_es), intent(inout) :: this
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
      call this%initialize_start_vector_tool()
      call this%initialize_projection_tool()
!
      call this%prepare_wf_for_excited_state()
!
      call this%davidson%initialize()
!
      call this%set_precondition_vector()
!
!     :: Set start vectors and handle restart ::
!
      call this%initialize_energies()
!
      call mem%alloc(c, this%wf%n_es_amplitudes)
!
      do trial = 1, this%n_singlet_states
!
         call this%set_initial_guesses(c, trial, trial)
         call this%davidson%set_trial(c, trial)
!
      enddo
!
      call mem%dealloc(c, this%wf%n_es_amplitudes)
!
!     :: Iterative loop ::
!
      call mem%alloc(converged, this%n_singlet_states)
      call mem%alloc(residual_lt_lindep, this%n_singlet_states)
!
      converged          = .false.
      residual_lt_lindep = .false.
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. this%max_iterations))
!
         iteration = iteration + 1
!
         write(timer_name, '(a,i0)') "Davidson: CC ES iteration ", iteration
         iteration_time = timings(trim(timer_name), pl="n")
         call iteration_time%turn_on()
!
!        Reduced space preparations
!
         if (this%davidson%red_dim_exceeds_max()) call this%davidson%set_trials_to_solutions()
!
         call this%davidson%update_reduced_dim()
!
         call this%davidson%orthonormalize_trial_vecs()
!
!        Print iteration information
!
         call output%printf('n', 'Iteration:               (i4)', &
                            ints=[iteration], fs='(/t3,a)')
         call output%printf('n', 'Reduced space dimension: (i4)', ints=[this%davidson%dim_red])
!
!        Transform new trial vectors and write them to file
!
         call mem%alloc(c, this%wf%n_es_amplitudes)
         call mem%alloc(residual, this%wf%n_es_amplitudes)
!
         do trial = this%davidson%first_new_trial(), this%davidson%last_new_trial()
!
            call this%davidson%get_trial(c, trial)
!
            call this%wf%construct_jacobian_transform(this%transformation, c, residual)
!
            if (this%projector%active) call this%projector%do_(residual)
!
            call this%davidson%set_transform(residual, trial)
!
         enddo ! Done transforming new trials
!
         call mem%dealloc(c, this%wf%n_es_amplitudes)
!
!        Solve problem in the reduced space
!
         call this%davidson%solve_reduced_problem()
!
!        Construct the full space solutions, test for convergence,
!        and use the residuals to construct the next trial vectors
!
         call mem%alloc(solution, this%wf%n_es_amplitudes)
!
         call output%printf('n', ' Root  Eigenvalue (Re)   Eigenvalue &
                            &(Im)    |residual|   Delta E (Re)', fs='(/t3,a)', ll=120)
         call output%print_separator('n', 72,'-')
!
         call construct_new_trial%turn_on()
!
         do n = 1, this%n_singlet_states
!
            call this%davidson%construct_solution(solution, n)
!
            call this%wf%save_excited_state(solution, n, n, this%transformation, this%energies(n))
!
            call this%davidson%construct_residual(residual, solution, n)
!
            if (this%projector%active) call this%projector%do_(residual) ! CVS projection,
                                                                             ! for instance
!
            residual_norm = get_l2_norm(residual, this%wf%n_es_amplitudes)
!
            converged(n) = this%convergence_checker%has_converged(residual_norm, &
                                 this%davidson%omega_re(n)-this%energies(n), iteration)
!
            residual_lt_lindep(n) = (residual_norm <= this%davidson%lindep_threshold)
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
            call output%printf('n', '(i4) (f16.12)  (f16.12)    (e11.4)  (e11.4) ',             &
                               ints=[n], ll=120,                                                &
                               reals=[this%davidson%omega_re(n),                              &
                                      this%davidson%omega_im(n),                              &
                                      residual_norm,                                            &
                                      abs(this%davidson%omega_re(n) - this%energies(n))])
!
         enddo ! Done constructing new trials from residuals
!
         call construct_new_trial%turn_off()
         call construct_new_trial%reset()
!
         call output%print_separator('n', 72,'-')
!
         call mem%dealloc(residual, this%wf%n_es_amplitudes)
         call mem%dealloc(solution, this%wf%n_es_amplitudes)
!
!        Update energies and save them
!
         this%energies = this%davidson%omega_re
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
         call this%wf%set_excitation_energies(this%energies, this%transformation)
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
         call mem%alloc(r, this%wf%n_es_amplitudes, this%n_singlet_states)
!
         do state = 1, this%n_singlet_states
!
            call this%davidson%construct_solution(r(:,state), state)
!
         enddo
!
         call this%print_summary(r)
!
         call mem%dealloc(r, this%wf%n_es_amplitudes, this%n_singlet_states)
!
      endif
!
      call this%davidson%cleanup()
!
      call mem%dealloc(converged, this%n_singlet_states)
      call mem%dealloc(residual_lt_lindep, this%n_singlet_states)
!
   end subroutine run_davidson_cc_es
!
!
   subroutine set_precondition_vector_davidson_cc_es(this)
!!
!!    Set precondition vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets precondition vector to orbital differences.
!!
      implicit none
!
      class(davidson_cc_es) :: this
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, this%wf%n_es_amplitudes)
      call this%wf%get_orbital_differences(preconditioner, this%wf%n_es_amplitudes)
      call this%davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, this%wf%n_es_amplitudes)
!
   end subroutine set_precondition_vector_davidson_cc_es
!
!
   subroutine read_settings_davidson_cc_es(this, records_in_memory)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      implicit none
!
      class(davidson_cc_es) :: this
      logical, intent(inout) :: records_in_memory
!
      call this%read_es_settings(records_in_memory)
      call this%read_davidson_settings()
!
   end subroutine read_settings_davidson_cc_es
!
!
   subroutine read_davidson_settings_davidson_cc_es(this)
!!
!!    Read Davidson settings
!!    Written by Andreas S. Skeidsvoll, Oct 2021
!!
      use global_in, only: input
!
      implicit none
!
      class(davidson_cc_es) :: this
!
      this%max_dim_red = max(100, 10*this%n_singlet_states)
!
      call input%get_keyword('max reduced dimension', &
                             'solver cc es',          &
                             this%max_dim_red)
!
   end subroutine read_davidson_settings_davidson_cc_es
!
!
end module davidson_cc_es_class
