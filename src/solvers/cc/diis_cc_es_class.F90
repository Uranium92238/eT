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
module diis_cc_es_class
!
!!
!!    DIIS coupled cluster excited state solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
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
   type, extends(abstract_cc_es) :: diis_cc_es
!
      integer :: diis_dimension
!
   contains
!     
      procedure, non_overridable :: run            => run_diis_cc_es
!
      procedure :: set_start_vectors               => set_start_vectors_diis_cc_es
!
      procedure :: read_settings                   => read_settings_diis_cc_es
      procedure :: read_diis_settings              => read_diis_settings_diis_cc_es
!
      procedure :: print_settings                  => print_settings_diis_cc_es
!
   end type diis_cc_es
!
!
   interface diis_cc_es 
!
      procedure :: new_diis_cc_es
!
   end interface diis_cc_es
!
!
contains
!
!
   function new_diis_cc_es(transformation, wf) result(solver)
!!
!!    New DIIS CC ES 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(diis_cc_es) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' excited state (' // trim(transformation) //')')
      call solver%timer%turn_on()
!
!     Set printables
!
      solver%name_  = 'DIIS coupled cluster excited state solver'
      solver%tag    = 'DIIS'
      solver%author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      solver%description1 = 'A DIIS solver that solves for the lowest eigenvalues and &
                           & the right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                           & problem is solved by DIIS extrapolation of residuals for each &
                           & eigenvector until the convergence criteria are met.'
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
      solver%diis_dimension       = 20
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
   end function new_diis_cc_es
!
!
   subroutine print_settings_diis_cc_es(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_es) :: solver 
!
      call solver%print_es_settings()
!
      call output%printf('DIIS dimension: (i3)',  ints=[solver%diis_dimension], fs='(/t6,a)')
!
   end subroutine print_settings_diis_cc_es
!
!
   subroutine read_settings_diis_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_diis_settings()
!
   end subroutine read_settings_diis_cc_es
!
!
   subroutine read_diis_settings_diis_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_es) :: solver 
!
      call input%get_keyword_in_section('diis dimension', 'solver cc es', solver%diis_dimension)
!
   end subroutine read_diis_settings_diis_cc_es
!
!
   subroutine run_diis_cc_es(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es) :: solver
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
      integer :: iteration, state, amplitude, n_solutions_on_file
!
      character(len=3) :: string_state
!
      real(dp) :: norm_X
!
      real(dp), dimension(:), allocatable   :: eps
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
!     Make DIIS tools array & initialize the individual DIIS tools 
!
      allocate(diis(solver%n_singlet_states))
!
      do state = 1, solver%n_singlet_states
!  
         write(string_state, '(i3.3)') state
         diis(state) = diis_tool('diis_cc_es_' // string_state, wf%n_es_amplitudes, wf%n_es_amplitudes, solver%diis_dimension)
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
      iteration = 0
!
      do while (.not. all(converged) .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1   
!
         write(output%unit,'(/t3,a25,i4)') 'Iteration:               ', iteration
!
         write(output%unit,'(/t3,a)') 'Root     Eigenvalue (Re)     Residual norm    '
         write(output%unit,'(t3,a)')  '----------------------------------------------'
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
!$omp parallel do private(amplitude)
               do amplitude = 1, wf%n_es_amplitudes
!
                  R(amplitude, state) = -R(amplitude, state)/(eps(amplitude) - solver%energies(state))
!
               enddo
!$omp end parallel do 
!
!              Update convergence logicals 
!
               converged_eigenvalue(state) = abs(solver%energies(state)-prev_energies(state)) &
                                                      .lt. solver%eigenvalue_threshold
               converged_residual(state)   = residual_norms(state) .lt. solver%residual_threshold
!
               converged(state) = converged_eigenvalue(state) .and. converged_residual(state)
!
!              Perform DIIS extrapolation to the optimal next guess for X,
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
                  X(:,state) = X(:,state) + R(:,state)
!
                  call diis(state)%update(R(:,state), X(:,state))
!
                  norm_X = get_l2_norm(X(:,state), wf%n_es_amplitudes)
                  X(:,state) = X(:,state)/norm_X
!
               endif 
!
            endif 
!
            write(output%unit, '(i3,3x,f19.12,6x,e11.4)') state, solver%energies(state), residual_norms(state)
            flush(output%unit)
            flush(timing%unit)
!
         enddo
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
         write(output%unit,'(t3,a)')  '----------------------------------------------'     
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
      else 
!
         call output%error_msg("Did not converge in the max number of iterations.")
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
   end subroutine run_diis_cc_es
!
!
   subroutine set_start_vectors_diis_cc_es(solver, wf, R, orbital_differences)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(diis_cc_es), intent(in) :: solver 
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
      call mem%alloc(lowest_orbital_differences, solver%n_singlet_states)
      call mem%alloc(lowest_orbital_differences_index, solver%n_singlet_states)
!
      call get_n_lowest(solver%n_singlet_states, wf%n_es_amplitudes, orbital_differences, &
                           lowest_orbital_differences, lowest_orbital_differences_index)
!
      do state = 1, solver%n_singlet_states
!
         R(:,state) = zero
         R(lowest_orbital_differences_index(state), state) = one
!
      enddo 
!
      call mem%dealloc(lowest_orbital_differences, solver%n_singlet_states)
      call mem%dealloc(lowest_orbital_differences_index, solver%n_singlet_states)      
!
   end subroutine set_start_vectors_diis_cc_es
!
!
end module diis_cc_es_class
