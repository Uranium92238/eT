!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   implicit none
!
   type, extends(abstract_cc_es) :: davidson_cc_es
!
      integer :: max_dim_red
!
   contains
!     
      procedure, non_overridable :: run              => run_davidson_cc_es
!
      procedure, nopass :: set_precondition_vector   => set_precondition_vector_davidson_cc_es
!
      procedure :: read_settings                     => read_settings_davidson_cc_es
      procedure :: read_davidson_settings            => read_davidson_settings_davidson_cc_es
!
      procedure :: print_settings                    => print_settings_davidson_cc_es
!
      procedure :: set_start_vectors                 => set_start_vectors_davidson_cc_es
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
      implicit none
!
      type(davidson_cc_es) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      logical, intent(in) :: restart
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' excited state (' // trim(transformation) //')')
      call solver%timer%turn_on()
!
      solver%name_ = 'Davidson coupled cluster excited state solver'
      solver%tag   = 'Davidson'
!
      solver%description1 = 'A Davidson solver that calculates the lowest eigenvalues and &
               & the right or left eigenvectors of the Jacobian matrix, A. The eigenvalue &
               & problem is solved in a reduced space, the dimension of which is &
               & expanded until the convergence criteria are met.'
!
      solver%description2 = 'A complete description of the algorithm can be found in &
                                 & E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%n_singlet_states     = 0
      solver%max_iterations       = 100
      solver%eigenvalue_threshold = 1.0d-6
      solver%residual_threshold   = 1.0d-6
      solver%restart              = restart
      solver%max_dim_red          = 100 
      solver%transformation       = trim(transformation)
      solver%es_type              = 'valence'
      solver%records_in_memory    = .false.  
      solver%storage              = 'disk'
!
      call solver%read_settings()
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      wf%n_singlet_states = solver%n_singlet_states
!
      call solver%initialize_start_vector_tool(wf)
      call solver%initialize_projection_tool(wf)
!
      call solver%prepare_wf_for_excited_state(wf)
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
   end function new_davidson_cc_es
!
!
   subroutine print_settings_davidson_cc_es(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(davidson_cc_es) :: solver 
!
      call solver%print_es_settings()
!
      call output%printf('m', 'Max reduced space dimension:  (i11)', &
                         ints=[solver%max_dim_red], fs='(/t6,a)')
!
   end subroutine print_settings_davidson_cc_es
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
      logical, dimension(:), allocatable :: converged_eigenvalue
      logical, dimension(:), allocatable :: converged_residual
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
      type(eigen_davidson_tool) :: davidson 
!  
      real(dp) :: lindep_threshold
!
!     :: Preparations ::
!  
      lindep_threshold = min(solver%residual_threshold, 1.0d-11)
!
      davidson = eigen_davidson_tool('cc_es_davidson',                           &
                                       wf%n_es_amplitudes,                       &
                                       solver%n_singlet_states,                  &
                                       lindep_threshold,                         &
                                       solver%max_dim_red)
!
      call davidson%initialize_trials_and_transforms(solver%records_in_memory)
!
      call solver%set_precondition_vector(wf, davidson)
      call solver%set_start_vectors(wf, davidson)
!
!     :: Iterative loop :: 
!
      call mem%alloc(converged, solver%n_singlet_states)
      call mem%alloc(converged_eigenvalue, solver%n_singlet_states)
      call mem%alloc(converged_residual, solver%n_singlet_states)
!
      converged            = .false. 
      converged_eigenvalue = .false. 
      converged_residual   = .false. 
!
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1
         call davidson%iterate()
!
         call output%printf('n', 'Iteration:               (i4)', &
                            ints=[iteration], fs='(/t3,a)')
         call output%printf('n', 'Reduced space dimension: (i4)', ints=[davidson%dim_red])
!
!        Transform new trial vectors and write them to file
!
         call mem%alloc(c, wf%n_es_amplitudes)
!
         do trial = davidson%first_new_trial(), davidson%last_new_trial()
!
            call davidson%get_trial(c, trial)
!
            call wf%construct_Jacobian_transform(solver%transformation, c)
!
            if (solver%projector%active) call solver%projector%do_(c)
!
            call davidson%set_transform(c, trial)
!
         enddo ! Done transforming new trials 
!
         call mem%dealloc(c, wf%n_es_amplitudes)
!
!        Solve problem in the reduced space
!
         call davidson%solve_reduced_problem()
!
!        Construct the full space solutions, test for convergence, 
!        and use the residuals to construct the next trial vectors 
!
         call mem%alloc(residual, wf%n_es_amplitudes)
         call mem%alloc(solution, wf%n_es_amplitudes)
!
         call output%printf('n', ' Root  Eigenvalue (Re)   Eigenvalue &
                            &(Im)    |residual|   Delta E (Re)', fs='(/t3,a)', ll=120)
         call output%print_separator('n', 72,'-')
!
         converged_residual   = .true.
         converged_eigenvalue = .true.
!
         do n = 1, solver%n_singlet_states
!
            call davidson%construct_solution(solution, n) 
!
            call wf%save_excited_state(solution, n, solver%transformation) 
!
            call davidson%construct_residual(residual, solution, n) 
!
            if (solver%projector%active) call solver%projector%do_(residual) ! CVS projection, 
                                                                             ! for instance
!
            residual_norm = get_l2_norm(residual, wf%n_es_amplitudes)
!
            converged_residual(n)   = residual_norm <= solver%residual_threshold
!
            converged_eigenvalue(n) = dabs(davidson%omega_re(n) &
                                          - solver%energies(n)) <= solver%eigenvalue_threshold
!
            converged(n) = converged_residual(n) .and. converged_eigenvalue(n) .or. &
                           converged_residual(n) .and. iteration == 1
!
            if (.not. converged(n)) then 
!
               if (residual_norm <= lindep_threshold) then 
!
                  call output%warning_msg('Residual norm for root (i0) smaller than linear ' // &
                                          'dependence threshold, but energy and residual  '  // &
                                          'thresholds have not yet been met. No new trial ' // &
                                          'added for this root.', ints=[n])
!
               else
!
                  call davidson%construct_next_trial(residual, n)
!
               endif
!
            endif
!
            call output%printf('n', '(i4) (f16.12)  (f16.12)    (e11.4)  (e11.4) ',       &
                               ints=[n],                                                  &
                               reals=[davidson%omega_re(n),                               &
                                      davidson%omega_im(n),                               &
                                      residual_norm,                                      &
                                      dabs(davidson%omega_re(n) - solver%energies(n))],   &
                               ll=120)
!
         enddo ! Done constructing new trials from residuals
!
         call output%print_separator('n', 72,'-')
!
         call mem%dealloc(residual, wf%n_es_amplitudes)
         call mem%dealloc(solution, wf%n_es_amplitudes)
!
!        Update energies and save them
!
         solver%energies = davidson%omega_re
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies, solver%transformation)
!
!        Special case when residuals converge in first iteration, e.g. on restart
!        => Exit without testing energy convergence 
!
         if (all(converged_residual) .and. iteration == 1) then 
!
            call output%printf('m', 'Note: Residual(s) converged in first &
                               &iteration. ' // 'Energy convergence therefore &
                               &not tested in this calculation.')
!
         endif
!
      enddo ! End of iterative loop
!
!     :: Calculation summary :: 
!
      if (.not. all(converged)) then
!
         call output%error_msg("Did not converge in the max number of " // & 
                                 "iterations in run_davidson_cc_es.")
!
      else 
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
         call mem%alloc(r, wf%n_es_amplitudes, solver%n_singlet_states)
!
         do state = 1, solver%n_singlet_states
!
            call davidson%construct_solution(r(:,state), state) 
!
         enddo 
!
         call solver%print_summary(wf, r) 
!
         call mem%dealloc(r, wf%n_es_amplitudes, solver%n_singlet_states)
!
      endif
!
      call davidson%finalize_trials_and_transforms()
!
      call mem%dealloc(converged, solver%n_singlet_states)
      call mem%dealloc(converged_eigenvalue, solver%n_singlet_states)
      call mem%dealloc(converged_residual, solver%n_singlet_states)
!
   end subroutine run_davidson_cc_es
!
!
   subroutine set_start_vectors_davidson_cc_es(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cc_es) :: solver
!
      class(ccs) :: wf
!
      class(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: c
!
      integer :: trial, n_solutions_on_file
!
      if (solver%restart) then 
!
!        Read the solutions from file & set as initial trial vectors 
!
         call solver%determine_restart_transformation(wf) ! Read right or left?
         n_solutions_on_file = wf%get_n_excited_states_on_file(solver%restart_transformation)
!
         call output%printf('m', 'Requested restart - restarting (i0) (a0) &
                           &eigenvectors from file.', fs='(/t3,a)', &
                            ints=[n_solutions_on_file], &
                            chars=[solver%restart_transformation])
!
         call mem%alloc(c, wf%n_es_amplitudes)
!
         do trial = 1, n_solutions_on_file
!
            call wf%read_excited_state(c, trial, solver%restart_transformation)
            call davidson%set_trial(c, trial)
!
         enddo 
!
         call mem%dealloc(c, wf%n_es_amplitudes)
!
      else
!
         n_solutions_on_file = 0
!
      endif 
!
!     Sets whatever remains using the start vector tool 
!
      if (n_solutions_on_file .lt. solver%n_singlet_states) then 
!
         call mem%alloc(c, wf%n_es_amplitudes)
!
         do trial = n_solutions_on_file + 1, solver%n_singlet_states
!
            call solver%start_vectors%get(c, trial)
            call davidson%set_trial(c, trial)
!
         enddo 
!
         call mem%dealloc(c, wf%n_es_amplitudes)
!     
      endif
!
   end subroutine set_start_vectors_davidson_cc_es
!
!
   subroutine set_precondition_vector_davidson_cc_es(wf, davidson)
!!
!!    Set precondition vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets precondition vector to orbital differences.
!!
      implicit none
!
      class(eigen_davidson_tool) :: davidson
!
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(preconditioner, wf%n_es_amplitudes)
      call davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_es_amplitudes)
!
   end subroutine set_precondition_vector_davidson_cc_es
!
!
   subroutine read_settings_davidson_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_davidson_settings()
!
   end subroutine read_settings_davidson_cc_es
!
!
   subroutine read_davidson_settings_davidson_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_es) :: solver 
!
      call input%get_keyword_in_section('max reduced dimension', 'solver cc es', solver%max_dim_red)
!
   end subroutine read_davidson_settings_davidson_cc_es
!
!
end module davidson_cc_es_class
