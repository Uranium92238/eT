!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
!!    Davidson coupled cluster excited state solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
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
      type(eigen_davidson_tool) :: davidson 
!
   contains
!     
      procedure, non_overridable :: run              => run_davidson_cc_es
!
      procedure :: set_precondition_vector           => set_precondition_vector_davidson_cc_es
!
      procedure :: read_settings                     => read_settings_davidson_cc_es
      procedure :: read_davidson_settings            => read_davidson_settings_davidson_cc_es
!
      procedure :: print_settings                    => print_settings_davidson_cc_es
!
      procedure :: set_start_vectors                 => set_start_vectors_davidson_cc_es
      procedure :: transform_trial_vector            => transform_trial_vector_davidson_cc_es
!       
      procedure :: initialize_projection_tool        => initialize_projection_tool_davidson_cc_es
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
   function new_davidson_cc_es(transformation, wf) result(solver)
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
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' excited state (' // trim(transformation) //')')
      call solver%timer%turn_on()
!
      solver%name_ = 'Davidson coupled cluster excited state solver'
      solver%tag = 'Davidson'
      solver%author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
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
      solver%transformation       = 'right'
      solver%restart              = .false.
      solver%max_dim_red          = 100 
      solver%transformation       = trim(transformation)
      solver%es_type              = 'valence'
!
      call solver%read_settings()
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      wf%n_excited_states = solver%n_singlet_states
!
      solver%davidson = eigen_davidson_tool('cc_es_davidson', wf%n_es_amplitudes, solver%n_singlet_states, &
                                       solver%residual_threshold, solver%eigenvalue_threshold)
!
      call solver%initialize_start_vector_tool(wf)
      call solver%initialize_projection_tool(wf)
!
      call solver%prepare_wf_for_excited_state(wf)
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
      call output%printf('Max reduced space dimension: (i4)',  ints=[solver%max_dim_red], fs='(/t6,a)')
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
      logical :: converged
      logical :: converged_eigenvalue
      logical :: converged_residual
!
      integer :: iteration, trial, solution, state
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:), allocatable :: c_i
      real(dp), dimension(:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: r
!
      converged            = .false. 
      converged_eigenvalue = .false. 
      converged_residual   = .false. 
!
      iteration = 1
!
      call solver%set_precondition_vector(wf)
!
!     Construct first trial vectors
!
      call solver%set_start_vectors(wf)
!
!     Enter iterative loop
!
      do while (.not. converged .and. (iteration .le. solver%max_iterations))
!
         write(output%unit,'(/t3,a25,i4)') 'Iteration:               ', iteration
         write(output%unit,'(t3,a25,i4/)') 'Reduced space dimension: ', solver%davidson%dim_red
         flush(output%unit)
!
         write(output%unit,'(t3,a)') 'Root     Eigenvalue (Re)        Eigenvalue (Im)      Residual norm'
         write(output%unit,'(t3,a)') '-------------------------------------------------------------------'
!
         flush(output%unit)
!
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         do trial = solver%davidson%dim_red - solver%davidson%n_new_trials + 1, solver%davidson%dim_red
!
            call solver%davidson%read_trial(c_i, trial)
!
            call solver%transform_trial_vector(wf, c_i)
!
            call solver%davidson%projection(c_i)
!
            call solver%davidson%write_transform(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes)
!
!        Solve problem in reduced space
!
         call solver%davidson%construct_reduced_matrix()
         call solver%davidson%solve_reduced_problem()
!
!        Construct new trials and check if convergence criterion on residual is satisfied
!
         converged_residual = .true.
!
         solver%davidson%n_new_trials = 0
!
         call mem%alloc(X, wf%n_es_amplitudes)
!
         do solution = 1, solver%n_singlet_states
!
            call solver%davidson%construct_next_trial_vec(residual_norm, solution)
            call solver%davidson%construct_X(X, solution)
            call wf%save_excited_state(X, solution, solver%transformation)
!
            write(output%unit,'(t3,i2,5x,f16.12,7x,f16.12,11x,e11.4)') &
            solution, solver%davidson%omega_re(solution), solver%davidson%omega_im(solution), residual_norm
            flush(output%unit)
!
            if (residual_norm .gt. solver%residual_threshold) converged_residual = .false.
!
         enddo
!
         call mem%dealloc(X, wf%n_es_amplitudes)
!
         write(output%unit,'(t3,a)') '-------------------------------------------------------------------'
!
!        Check if convergence criterion on energy is satisfied
!
         converged_eigenvalue = .true.
!
         do solution = 1, solver%n_singlet_states
!
            if (abs(solver%davidson%omega_re(solution) - solver%energies(solution)) &
               .gt. solver%eigenvalue_threshold) converged_eigenvalue = .false.
!
         enddo
!
         if (solver%davidson%dim_red .ge. solver%max_dim_red) then
!
            call solver%davidson%set_trials_to_solutions()
!
         else
!
            solver%davidson%dim_red = solver%davidson%dim_red + solver%davidson%n_new_trials
!
         endif
!
!        Update energies and save them
!
         solver%energies = solver%davidson%omega_re
!
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies, solver%transformation)
!
!        Test for total convergence
!
         if (converged_residual) then 
!
!           Tests for convergence of energy or restart
!
            if (converged_eigenvalue) then
!
              converged = .true.
!
            elseif (iteration .eq. 1) then
!
               converged = .true.
               write(output%unit,'(/t3,a,/t3,a)') 'Note: residual(s) converged in first iteration.', &
                                                   'Energy convergence therefore not tested in this calculation.'
!
            endif
!
            flush(output%unit)
!
         endif ! Not yet converged
!
         iteration = iteration + 1       
!
      enddo
!
      if (converged) then
!
         write(output%unit,'(/t3,a, i3, a)') 'Convergence criterion met in ', iteration - 1, ' iterations!'
!
         call mem%alloc(r, wf%n_es_amplitudes, solver%n_singlet_states)
!
         do state = 1, solver%n_singlet_states
!
            call solver%davidson%construct_X(r(:,state), state)   
!
         enddo 
!
         call solver%print_summary(wf, r) 
!
         call mem%dealloc(r, wf%n_es_amplitudes, solver%n_singlet_states)
!
         write(output%unit,'(/t3,a)') '- Storing excited states to file.'
!
         call mem%alloc(X, wf%n_es_amplitudes)
!
         do solution = 1, solver%n_singlet_states
!
            call solver%davidson%construct_X(X, solution)
!
            call wf%save_excited_state(X, solution, solver%transformation)
!
         enddo
!
         call mem%dealloc(X, wf%n_es_amplitudes)
!
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies, solver%transformation)
!
      elseif (.not. converged ) then
!
         call output%error_msg("Did not converge in the max number of iterations.")
!
      endif
!
   end subroutine run_davidson_cc_es
!
!
   subroutine transform_trial_vector_davidson_cc_es(solver, wf, c_i)
!!
!!    Transform trial vector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Transforms the trial vector according to specified transformation routine.
!!
      class(davidson_cc_es), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c_i
!
      if (trim(solver%transformation) == 'right') then 
!
         call wf%jacobian_transform_trial_vector(c_i)
!
      elseif (trim(solver%transformation) == 'left') then 
!
         call wf%jacobian_transpose_transform_trial_vector(c_i)
!
      endif 
!
   end subroutine transform_trial_vector_davidson_cc_es
!
!
   subroutine set_start_vectors_davidson_cc_es(solver, wf)
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
      real(dp), dimension(:), allocatable :: c_i
!
      integer :: trial, n_solutions_on_file
!
      if (solver%restart) then 
!
!        Read the solutions from file & set as initial trial vectors 
!
         call wf%is_restart_safe('excited state')
!
         call solver%determine_restart_transformation(wf) ! Read right or left?
         n_solutions_on_file = wf%get_n_excited_states_on_file(solver%restart_transformation)
!
         call output%printf('Requested restart - there are (i0) (a0) eigenvectors on file.', &
                              ints=[n_solutions_on_file], chars=[solver%restart_transformation])
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         do trial = 1, n_solutions_on_file
!
            call wf%read_excited_state(c_i, trial, solver%restart_transformation)
            call solver%davidson%write_trial(c_i)
!
         enddo 
!
         call mem%dealloc(c_i, wf%n_es_amplitudes)
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
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         do trial = n_solutions_on_file + 1, solver%n_singlet_states
!
            call solver%start_vector_tool%get_vector(c_i, trial)
            call solver%davidson%write_trial(c_i)
!
         enddo 
!     
      endif
!
      call solver%davidson%orthonormalize_trial_vecs()
!
   end subroutine set_start_vectors_davidson_cc_es
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
!
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(preconditioner, wf%n_es_amplitudes)
      call solver%davidson%set_preconditioner(preconditioner)
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
   subroutine initialize_projection_tool_davidson_cc_es(solver, wf)
!!
!!    Initialize projection tool 
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none
!
      class(davidson_cc_es) :: solver 
!
      class(ccs) :: wf 
!
      if (trim(solver%es_type) == 'valence') then 
!
         solver%projection_tool = es_valence_projection_tool()
!
      elseif (trim(solver%es_type) == 'core') then 
!
         solver%projection_tool = es_cvs_projection_tool(wf, solver%davidson)
!
      elseif (trim(solver%es_type) == 'ionize') then 
!
         solver%projection_tool = es_ip_projection_tool(wf, solver%davidson)
!
      else 
!
         call output%error_msg('could not recognize excited state type in abstract_cc_es')
!
      endif
!
   end subroutine initialize_projection_tool_davidson_cc_es
!
!
end module davidson_cc_es_class
