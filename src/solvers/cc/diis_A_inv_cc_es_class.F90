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
module diis_A_inv_cc_es_class
!
!!
!!    Newton-Raphson coupled cluster excited state solver class module
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    Based on the NR solver for ground state (newton_raphson_cc_gs_class.F90)
!!    and the diis solver for excited states (diis_cc_es_class.F90)
!!
!!
!
   use kinds
   use file_class
   use ccs_class
   use diis_tool_class
   use abstract_cc_es_class, only: abstract_cc_es
!
   implicit none
!
   type, extends(abstract_cc_es) :: diis_A_inv_cc_es
!
      integer :: max_micro_iterations
!
   contains
!     
      procedure, non_overridable :: run            => run_diis_A_inv_cc_es
!
      procedure :: set_start_vectors               => set_start_vectors_diis_A_inv_cc_es
!
      procedure :: read_settings                   => read_settings_diis_A_inv_cc_es
      procedure :: read_nr_settings                => read_nr_settings_diis_A_inv_cc_es
!
      procedure :: print_settings                  => print_settings_diis_A_inv_cc_es
!
      procedure, private :: do_micro_iterations    => do_micro_iterations_diis_A_inv_cc_es
!
   end type diis_A_inv_cc_es
!
!
   interface diis_A_inv_cc_es 
!
      procedure :: new_diis_A_inv_cc_es
!
   end interface diis_A_inv_cc_es
!
!
contains
!
!
   function new_diis_A_inv_cc_es(transformation, wf) result(solver)
!!
!!    New NR CC ES 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(diis_A_inv_cc_es) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' excited state (' // trim(transformation) //')')
      call solver%timer%turn_on()
!
!     Set printables
!
      solver%name_  = 'Newton-Raphson coupled cluster excited state solver'
      solver%tag    = 'Newton-Raphson'
      solver%author = 'S. D. Folkestad and E. F. Kjønstad, 2019'
!
      solver%description1 = 'A DIIS solver for excited states ( left or &
                           &right eigenvectors of the Jacobian matrix, A) using & 
                           &A^-1 preconditioner. The eigenvalue &
                           &problem is solved by DIIS extrapolation of residuals for each &
                           &eigenvector until the convergence criteria are met.'
!
      solver%description2 = 'More on the DIIS algorithm can be found in &
                           &P. Pulay, Chemical Physics Letters, 73(2), 393-398 (1980).'
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
      solver%max_micro_iterations = 100
      solver%restart              = .false.
      solver%transformation       = trim(transformation)
!
      call solver%read_settings()
      call solver%print_settings()
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      call mem%alloc(solver%energies, solver%n_singlet_states)
      solver%energies = zero
!
      wf%n_excited_states = solver%n_singlet_states
!
   end function new_diis_A_inv_cc_es
!
!
   subroutine set_start_vectors_diis_A_inv_cc_es(solver, wf, R, orbital_differences)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(diis_A_inv_cc_es), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes, solver%n_singlet_states), intent(inout) :: R 
      real(dp), dimension(wf%n_es_amplitudes), intent(in)                             :: orbital_differences 
!
      real(dp), dimension(:), allocatable :: lowest_orbital_differences
!
      integer, dimension(:), allocatable :: lowest_orbital_differences_index
!
      integer :: state
!
      if (wf%bath_orbital) call output%error_msg('Bath orbitals can not be used in valence excitation calculation')
!
      call mem%alloc(lowest_orbital_differences, solver%n_singlet_states)
      call mem%alloc(lowest_orbital_differences_index, solver%n_singlet_states)
!
      call get_n_lowest(solver%n_singlet_states, wf%n_es_amplitudes, orbital_differences, &
                           lowest_orbital_differences, lowest_orbital_differences_index)
!
      call zero_array(R, (solver%n_singlet_states)*(wf%n_es_amplitudes))
!
      do state = 1, solver%n_singlet_states
!
         R(lowest_orbital_differences_index(state), state) = one
!
      enddo 
!
      call mem%dealloc(lowest_orbital_differences, solver%n_singlet_states)
      call mem%dealloc(lowest_orbital_differences_index, solver%n_singlet_states)      
!
   end subroutine set_start_vectors_diis_A_inv_cc_es
!
!
   subroutine read_settings_diis_A_inv_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_A_inv_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_nr_settings()
!
   end subroutine read_settings_diis_A_inv_cc_es
!
!
   subroutine read_nr_settings_diis_A_inv_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_A_inv_cc_es) :: solver 
!
      call input%get_keyword_in_section('max micro iterations', 'solver cc es', solver%max_micro_iterations)
!
   end subroutine read_nr_settings_diis_A_inv_cc_es
!
!
   subroutine print_settings_diis_A_inv_cc_es(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_A_inv_cc_es) :: solver 
!
      call solver%print_es_settings()
!
      call output%printf('Max micro iterations: (i3)',  ints=[solver%max_micro_iterations], fs='(/t6,a)')
!
   end subroutine print_settings_diis_A_inv_cc_es
!
!
!
   subroutine do_micro_iterations_diis_A_inv_cc_es(solver, wf, R, dX, micro_threshold, final_micro_iteration)
!!
!!    Do micro iterations 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      use linear_davidson_tool_class
!
      implicit none 
!
      class(diis_A_inv_cc_es), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: R 
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: dX  
!
      real(dp), intent(in) :: micro_threshold
!
      integer, intent(out) :: final_micro_iteration
!
      real(dp), dimension(:), allocatable :: eps, first_trial, c_i, rho_i
!
      real(dp) :: norm_trial, residual_norm
!
      type(linear_davidson_tool) :: davidson 
!
      integer :: micro_iteration
!
      logical :: converged_residual 
!
!     Initialize solver tool and set preconditioner 
!
      call mem%alloc(eps, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(eps, wf%n_es_amplitudes)
!
      davidson = linear_davidson_tool('cc_es_NR', wf%n_es_amplitudes, micro_threshold, -R)
!
      call davidson%set_preconditioner(eps)
      call mem%dealloc(eps, wf%n_es_amplitudes)
!
!     Set start vector / initial guess 
!
      call mem%alloc(first_trial, wf%n_es_amplitudes)
!
      call copy_and_scale(-one, R, first_trial, wf%n_es_amplitudes)
!
      call davidson%precondition(first_trial)
!
      norm_trial = get_l2_norm(first_trial, wf%n_es_amplitudes)
      call dscal(wf%n_es_amplitudes, one/norm_trial, first_trial, 1)
!
      call davidson%write_trial(first_trial, 'rewind')
      call mem%dealloc(first_trial, wf%n_es_amplitudes)
!
!     Enter iterative loop
!
   !   write(output%unit,'(/t6,a)') 'Micro-iter.  Residual norm'
   !   write(output%unit,'(t6,a)')  '--------------------------'
   !   flush(output%unit)
!
      micro_iteration = 1
      converged_residual = .false.
!
      do while (.not. converged_residual .and. (micro_iteration .le. solver%max_micro_iterations))
!
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, davidson%n_parameters)
         call mem%alloc(rho_i, davidson%n_parameters)
!
         call davidson%read_trial(c_i, davidson%dim_red)
         call dcopy(davidson%n_parameters, c_i, 1, rho_i, 1)
!
         if (solver%transformation == 'right') then
!
            call wf%jacobian_transform_trial_vector(rho_i)
!
         elseif (solver%transformation == 'left') then
!
            call wf%jacobian_transpose_transform_trial_vector(rho_i)
!
         endif
!
         if (micro_iteration == 1) then
!
            call davidson%write_transform(rho_i, 'rewind')
!
         else
!
            call davidson%write_transform(rho_i, 'append')
!
         endif
!
         call mem%dealloc(c_i, davidson%n_parameters)
         call mem%dealloc(rho_i, davidson%n_parameters)
!
!        Solve problem in reduced space
!
         call davidson%construct_reduced_matrix()
         call davidson%construct_reduced_gradient()
         call davidson%solve_reduced_problem()
!
!        Construct new trials and check if convergence criterion on residual is satisfied
!
         davidson%n_new_trials = 0
!
         call davidson%construct_next_trial_vec(residual_norm)
!
      !   write(output%unit,'(t6,2x,i3,8x,e11.4)') micro_iteration, residual_norm
      !   flush(output%unit)
!
         converged_residual = .true.
!
         if (residual_norm .gt. micro_threshold) converged_residual = .false.
!   
         davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
         if (.not. converged_residual) micro_iteration = micro_iteration + 1       
!
      enddo
!
   !   write(output%unit,'(t6,a/)')  '--------------------------'
   !   flush(output%unit)
!
      if (.not. converged_residual) then
!
         write(output%unit, '(/t6,a)')  'Warning: was not able to converge the equations in the given'
         write(output%unit, '(t6,a/)')  'number of maximum micro-iterations.'
         flush(output%unit)
!
      endif
!
      call davidson%construct_X(dX, 1)
      call davidson%cleanup()
!
      final_micro_iteration = micro_iteration
!
   end subroutine do_micro_iterations_diis_A_inv_cc_es
!
!
!
!
   subroutine run_diis_A_inv_cc_es(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_A_inv_cc_es) :: solver
!
      class(ccs) :: wf
!
      logical, dimension(:), allocatable :: converged
!
      logical, dimension(:), allocatable :: converged_eigenvalue
      logical, dimension(:), allocatable :: converged_residual
!
      real(dp), dimension(:), allocatable :: prev_energies 
      real(dp), dimension(:), allocatable :: residual_norms
!
      type(diis_tool), dimension(:), allocatable :: diis
!
      integer :: iteration, state, n_solutions_on_file, micro_iteration
!
      character(len=3) :: string_state
!
      real(dp) :: norm_X, micro_threshold
!
      real(dp), dimension(:), allocatable   :: eps, dX
      real(dp), dimension(:,:), allocatable :: X, R
!
      call solver%prepare_wf_for_excited_state(wf)
!
!     Initialize energies, residual norms, and convergence arrays 
!
      call mem%alloc(prev_energies, solver%n_singlet_states)
      call mem%alloc(residual_norms, solver%n_singlet_states)
!
      prev_energies     = zero 
      residual_norms    = zero 
!
      allocate(converged(solver%n_singlet_states))
      allocate(converged_residual(solver%n_singlet_states))
      allocate(converged_eigenvalue(solver%n_singlet_states))
!
      converged            = .false.
      converged_residual   = .false.
      converged_eigenvalue = .false.
!
!     Make nr tools array & initialize the individual nr tools 
!
      allocate(diis(solver%n_singlet_states))
!
      do state = 1, solver%n_singlet_states
!  
         write(string_state, '(i3.3)') state
         diis(state) = diis_tool('diis_A_inv_cc_es_' // string_state, wf%n_es_amplitudes, wf%n_es_amplitudes, 20)
!
      enddo 
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call mem%alloc(eps, wf%n_es_amplitudes)
      call wf%get_es_orbital_differences(eps, wf%n_es_amplitudes)
!
      call mem%alloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
!
      call solver%set_start_vectors(wf, X, eps) ! Use orbital differences (Koopman)
!
      if (solver%restart) then ! Overwrite all or some of the orbital differences 
!
         call wf%is_restart_safe('excited state')
         call solver%determine_restart_transformation(wf) ! Read right or left?
!
         n_solutions_on_file = wf%get_n_excited_states_on_file(solver%restart_transformation)
!
         call output%printf('Requested restart - there are (i0) (a0) eigenvectors on file.', &
                              ints=[n_solutions_on_file], chars=[solver%restart_transformation])
!
         do state = 1, n_solutions_on_file
!
            call wf%read_excited_state(X(:,state), state, solver%restart_transformation)
!
         enddo
!
         solver%energies = zero
         call wf%read_excitation_energies(n_solutions_on_file, solver%energies(1:n_solutions_on_file))
!
      endif
!
!     Enter iterative loop
!
      call mem%alloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
!
      iteration = 1
!
      do while (.not. all(converged) .and. (iteration .le. solver%max_iterations))   
!
         write(output%unit,'(/t3,a25,i4)') 'Iteration:               ', iteration
!
         write(output%unit,'(/t3,a)') 'Root  Eigenvalue (Re)  Residual norm   # micro iter '
         write(output%unit,'(t3,a)')  '----------------------------------------------------'
         flush(output%unit)
!
         do state = 1, solver%n_singlet_states
!
            if (.not. converged(state)) then 
!
!              Construct residual and energy and precondition the former 
!
               call wf%construct_excited_state_equation(X(:,state), R(:,state), solver%energies(state), &
                                                        solver%transformation)
!
               residual_norms(state) = get_l2_norm(R(:, state), wf%n_es_amplitudes)
!
               call mem%alloc(dX, wf%n_es_amplitudes)
!
               micro_threshold = 1.0d-2*residual_norms(state)
!
               call solver%do_micro_iterations(wf, R(:,state), dX, micro_threshold, micro_iteration) 
!
!              Update convergence logicals 
!
               converged_eigenvalue(state) = abs(solver%energies(state)-prev_energies(state)) &
                                                      .lt. solver%eigenvalue_threshold
               converged_residual(state)   = residual_norms(state) .lt. solver%residual_threshold
!
               converged(state) = converged_eigenvalue(state) .and. converged_residual(state)
!
!              Perform nr extrapolation to the optimal next guess for X,
!              then normalize it to avoid accumulating norm in X
!
               if (converged_residual(state) .and. iteration .eq. 1) then 
!
                  converged(state) = .true.
!
               endif
!
               if (.not. converged(state)) then
!
                  call daxpy(wf%n_es_amplitudes, one, dX, 1, X(1, state), 1)
!
                  call diis(state)%update(dX(:), X(:,state))
!
                  norm_X = get_l2_norm(X(:,state), wf%n_es_amplitudes)
!
                  call dscal(wf%n_es_amplitudes, one/norm_X, X(1,state), 1)
!
               endif 
!
               call mem%dealloc(dX, wf%n_es_amplitudes)
!
            endif 
!
            write(output%unit, '(i3,6x,f12.8,6x,e11.4,5x,i3)') state, solver%energies(state), residual_norms(state), micro_iteration
            flush(output%unit)
            flush(timing%unit)
!
         enddo
!
         if (.not. all(converged)) iteration = iteration + 1
!
!        Save excited states and excitation energies
!
         do state = 1, solver%n_singlet_states
!
            call wf%save_excited_state(X(:,state), state, solver%transformation)
!
         enddo 
!
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies, solver%transformation)
         prev_energies = solver%energies 
!
         write(output%unit,'(t3,a)')  '----------------------------------------------------'
!
      enddo 
!
      if (all(converged)) then 
!
         if (iteration .eq. 1) then 
!
            write(output%unit, '(/t3,a)')  'Note: residual of all states converged in first iteration.'
            write(output%unit, '(t3,a/)')  'Energy convergence has not been tested.'
!
         endif
!
         write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
         call solver%print_summary(wf, X) 
!
         write(output%unit, '(/t3,a)') '- Storing converged states to file.'       
!
         do state = 1, solver%n_singlet_states
!
            call wf%save_excited_state(X(:,state), state, solver%transformation)
!
         enddo 
!
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies, solver%transformation)
!
      endif 
!
      do state = 1, solver%n_singlet_states
!
         call diis(state)%cleanup()
!
      enddo
!
      call mem%dealloc(prev_energies, solver%n_singlet_states)
      call mem%dealloc(residual_norms, solver%n_singlet_states)
!
      deallocate(converged)
      deallocate(converged_residual)
      deallocate(converged_eigenvalue)
!
      call mem%dealloc(eps, wf%n_es_amplitudes)
      call mem%dealloc(X, wf%n_es_amplitudes, solver%n_singlet_states)
      call mem%dealloc(R, wf%n_es_amplitudes, solver%n_singlet_states)
!
   end subroutine run_diis_A_inv_cc_es
!
!
end module diis_A_inv_cc_es_class
