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
   use file_class
   use ccs_class
   use eigen_davidson_tool_class
!
   implicit none
!
   type :: davidson_cc_es
!
      character(len=100) :: tag = 'Davidson coupled cluster excited state solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      character(len=500) :: description1 = 'A Davidson solver that calculates the lowest eigenvalues and &
                                           &the right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                                           &problem is solved in a reduced space, the dimension of which is &
                                           &expanded until the convergence criteria are met.'
!
      character(len=500) :: description2 = 'A complete description of the algorithm can be found in &
                                           &E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      integer :: max_iterations
!
      real(dp) :: eigenvalue_threshold  
      real(dp) :: residual_threshold  
!
      logical  :: restart = .false.
!
      integer :: n_singlet_states = 0
!
      character(len=40) :: transformation 
!
      real(dp), dimension(:,:), allocatable :: energies
!
      integer, dimension(:,:), allocatable :: start_vectors
!
   contains
!     
      procedure, non_overridable :: prepare          => prepare_davidson_cc_es
      procedure, non_overridable :: run              => run_davidson_cc_es
      procedure, non_overridable :: cleanup          => cleanup_davidson_cc_es
!
      procedure, nopass :: set_precondition_vector   => set_precondition_vector_davidson_cc_es
      procedure :: set_projection_vector             => set_projection_vector_davidson_cc_es
!
      procedure :: print_banner                      => print_banner_davidson_cc_es
!
      procedure :: read_settings                     => read_settings_davidson_cc_es
!
      procedure :: print_settings                    => print_settings_davidson_cc_es
      procedure :: print_summary                     => print_summary_davidson_cc_es
!
      procedure :: set_start_vectors                 => set_start_vectors_davidson_cc_es
      procedure :: transform_trial_vector            => transform_trial_vector_davidson_cc_es
!       
      procedure :: initialize_energies               => initialize_energies_davidson_cc_es
      procedure :: destruct_energies                 => destruct_energies_davidson_cc_es   
!
   end type davidson_cc_es
!
!
contains
!
!
<<<<<<< HEAD:src/solvers/cc/davidson_cc_es_solver_class.F90
   subroutine prepare_davidson_cc_es_solver(solver, transform)
=======
   subroutine prepare_davidson_cc_es(solver)
>>>>>>> feadc9e97efbaf4bd11d0ab9fa937d2aa6dd0263:src/solvers/cc/davidson_cc_es_class.F90
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
<<<<<<< HEAD:src/solvers/cc/davidson_cc_es_solver_class.F90
      class(davidson_cc_es_solver) :: solver
      character(len=*), optional :: transform
=======
      class(davidson_cc_es) :: solver
>>>>>>> feadc9e97efbaf4bd11d0ab9fa937d2aa6dd0263:src/solvers/cc/davidson_cc_es_class.F90
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
!
      call solver%read_settings()
!
      if (present(transform)) solver%transformation = transform
!
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      write(output%unit, '(/t3,a,a,a)') 'Solving for the ', trim(solver%transformation), ' eigenvectors.'
      flush(output%unit)
!
   end subroutine prepare_davidson_cc_es
!
!
   subroutine initialize_energies_davidson_cc_es(solver)
!!
!!    Initialize energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initialize excitation energies
!!
      implicit none
!
      class(davidson_cc_es) :: solver
!
      if (.not. allocated(solver%energies)) &
            call mem%alloc(solver%energies, solver%n_singlet_states, 1)
!
   end subroutine initialize_energies_davidson_cc_es
!
!
   subroutine destruct_energies_davidson_cc_es(solver)
!!
!!    Destruct energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Destruct excitation energies
!!
      implicit none
!
      class(davidson_cc_es) :: solver
!
      if (allocated(solver%energies)) &
            call mem%dealloc(solver%energies, solver%n_singlet_states, 1)
!
   end subroutine destruct_energies_davidson_cc_es
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
      write(output%unit, '(/t3,a)') '- Davidson CC excited state solver settings:'
!
      write(output%unit,'(/t6,a20,e9.2)') 'Energy threshold:   ', solver%eigenvalue_threshold
      write(output%unit,'(t6,a20,e9.2)')  'Residual threshold: ', solver%residual_threshold
      write(output%unit,'(/t6,a,i3,a)')   'Number of singlet states: ', solver%n_singlet_states
      write(output%unit, '(t6,a26,i3)')   'Max number of iterations: ', solver%max_iterations
      write(output%unit, '(t6,a20,a)')   'Transformation: ', solver%transformation
      flush(output%unit)
!
   end subroutine print_settings_davidson_cc_es
!
!
   subroutine print_summary_davidson_cc_es(solver, davidson, wf)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(davidson_cc_es), intent(in) :: solver 
      class(eigen_davidson_tool) :: davidson
!
      class(ccs), intent(in) :: wf 
!
      integer :: state 
!
      real(dp), dimension(:,:), allocatable :: r
!
      write(output%unit, '(/t3,a)') '- Excitation vector amplitudes:'
!
      call mem%alloc(r, wf%n_es_amplitudes, 1)
!
      do state = 1, solver%n_singlet_states
!
         write(output%unit, '(/t6,a21,i2)')    'Electronic state nr. ', state
!
         call davidson%read_solution(r, state)         
!
         write(output%unit, '(/t6,a30,f15.12)')  'Energy (Hartree):             ', davidson%get_eigenvalue(state)
         write(output%unit, '(t6,a30,f15.12)') 'Fraction singles (|r1|/|r|):  ', &
                        get_l2_norm(r(1:wf%n_t1,1),wf%n_t1)/get_l2_norm(r,wf%n_es_amplitudes)   
!
         call wf%print_dominant_x_amplitudes(r, 'r')
!
      enddo 
!
      call mem%dealloc(r, wf%n_es_amplitudes, 1)
!
      write(output%unit, '(/t3,a)') '- Electronic excitation energies:'
!
      write(output%unit, '(/t6,a)') '                                 Excitation energy            '
      write(output%unit, '(t6,a)')  '                     ------------------------------------------'
      write(output%unit, '(t6,a)')  'State                (Hartree)             (eV)                '
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
      do state = 1, solver%n_singlet_states
!
         write(output%unit, '(t6,i2,14x,f19.12,4x,f19.12)') state, davidson%get_eigenvalue(state), &
                                                           davidson%get_eigenvalue(state)*Hartree_to_eV
!
      enddo 
!
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
      write(output%unit, '(t6,a26,f11.8)') 'eV/Hartree (CODATA 2014): ', Hartree_to_eV
!
   end subroutine print_summary_davidson_cc_es
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
      type(eigen_davidson_tool) :: davidson
!
      integer :: iteration, trial, solution
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:), allocatable :: X
!
      call wf%prepare_for_jacobian()
!
      converged            = .false. 
      converged_eigenvalue = .false. 
      converged_residual   = .false. 
!
      iteration = 1
!
<<<<<<< HEAD:src/solvers/cc/davidson_cc_es_solver_class.F90
      call davidson%prepare(trim(wf%name_)//'_es_davidson_'//trim(solver%transformation), wf%n_es_amplitudes, &
                            solver%n_singlet_states, solver%residual_threshold, solver%eigenvalue_threshold)
=======
      call davidson%prepare('cc_es_davidson', wf%n_es_amplitudes, solver%n_singlet_states, &
                               solver%residual_threshold, solver%eigenvalue_threshold)
>>>>>>> feadc9e97efbaf4bd11d0ab9fa937d2aa6dd0263:src/solvers/cc/davidson_cc_es_class.F90
!
!     Construct first trial vectors
!
      call solver%set_start_vectors(wf, davidson)
!
      call solver%set_precondition_vector(wf, davidson)
      call solver%set_projection_vector(wf, davidson)
!
!     Enter iterative loop
!
      do while (.not. converged .and. (iteration .le. solver%max_iterations))
!
         write(output%unit,'(/t3,a25,i4)') 'Iteration:               ', iteration
         write(output%unit,'(t3,a25,i4/)') 'Reduced space dimension: ', davidson%dim_red
         flush(output%unit)
!
         write(output%unit,'(t3,a)') 'Root     Eigenvalue (Re)        Eigenvalue (Im)      Residual norm'
         write(output%unit,'(t3,a)') '-------------------------------------------------------------------'
!
         flush(output%unit)
!
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, wf%n_es_amplitudes, 1)
!
         do trial = davidson%dim_red - davidson%n_new_trials + 1, davidson%dim_red
!
            call davidson%read_trial(c_i, trial)
!
            call solver%transform_trial_vector(wf, c_i)
!
            call davidson%projection(c_i)
!
            call davidson%write_transform(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
!
!        Solve problem in reduced space
!
         call davidson%construct_reduced_matrix()
         call davidson%solve_reduced_problem()
!
!        Construct new trials and check if convergence criterion on residual is satisfied
!
         converged_residual = .true.
!
         davidson%n_new_trials = 0
!
         do solution = 1, solver%n_singlet_states
!
            call davidson%construct_next_trial_vec(residual_norm, solution)
!
            write(output%unit,'(t3,i2,5x,f16.12,7x,f16.12,11x,e11.4)') &
            solution, davidson%omega_re(solution, 1), davidson%omega_im(solution, 1), residual_norm
            flush(output%unit)
!
            if (residual_norm .gt. solver%residual_threshold) converged_residual = .false.
!
         enddo
!
         write(output%unit,'(t3,a)') '-------------------------------------------------------------------'
!
!        Check if convergence criterion on energy is satisfied
!
         converged_eigenvalue = .true.
!
         do solution = 1, solver%n_singlet_states
!
            if (abs(davidson%omega_re(solution, 1) - solver%energies(solution, 1)) &
               .gt. solver%eigenvalue_threshold) converged_eigenvalue = .false.
!
         enddo
!
         if (davidson%dim_red .ge. davidson%max_dim_red) then
!
            call davidson%set_trials_to_solutions()
!
         else
!
            davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
         endif
!
!        Update energies
!
         solver%energies = davidson%omega_re
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
         call solver%print_summary(davidson, wf)
!
         write(output%unit,'(/t3,a)') 'Storing excited states to file.'
!
         call mem%alloc(X, wf%n_es_amplitudes)
!
         do solution = 1, solver%n_singlet_states
!
            call davidson%read_solution(X, solution)
            call wf%save_excited_state(X, solution, solver%transformation)
!
         enddo
!
         call wf%save_excitation_energies(solver%n_singlet_states, solver%energies)
!
      elseif (.not. converged ) then
!
         write(output%unit,'(/t3,a)') 'Maximal number of iterations performed without reaching convergence!'
!
      endif
!
      call davidson%cleanup()
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
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: c_i
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
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: orbital_differences
      real(dp), dimension(:,:), allocatable :: lowest_orbital_differences
!
      integer, dimension(:,:), allocatable :: lowest_orbital_differences_index
!
      integer :: trial, n_solutions_on_file
!
      if (allocated(solver%start_vectors)) then
!
!        Initial trial vectors given on input
!
         call mem%alloc(c_i, wf%n_es_amplitudes, 1)
!
         do trial = 1, solver%n_singlet_states
!
            c_i = zero
            c_i(solver%start_vectors(trial, 1), 1) = one
!
            call davidson%write_trial(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
!
      else
!
         if (solver%restart) then 
!
!           Read the solutions from file & set as initial trial vectors 
!
            call wf%get_n_excited_states_on_file(solver%transformation, n_solutions_on_file)
!
            write(output%unit, '(/t3,a,i0,a)') 'Requested restart. There are ', n_solutions_on_file, &
                                                ' solutions on file.'
!
            call mem%alloc(c_i, wf%n_es_amplitudes, 1)
!
            do trial = 1, n_solutions_on_file
!
               call wf%restart_excited_state(c_i, trial, solver%transformation)
               call davidson%write_trial(c_i)
!
            enddo 
!
            call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
!
         else
!
            n_solutions_on_file = 0
!
         endif 
!
         if (n_solutions_on_file .lt. solver%n_singlet_states) then 
!
!           Compute the remaining start vectors using Koopman
!
            call mem%alloc(orbital_differences, wf%n_es_amplitudes, 1)
            call wf%get_es_orbital_differences(orbital_differences, wf%n_es_amplitudes)
!
            call mem%alloc(lowest_orbital_differences, solver%n_singlet_states, 1)
            call mem%alloc(lowest_orbital_differences_index, solver%n_singlet_states, 1)
!
            call get_n_lowest(solver%n_singlet_states, wf%n_es_amplitudes, orbital_differences, &
                           lowest_orbital_differences, lowest_orbital_differences_index)
!
            call mem%dealloc(lowest_orbital_differences, solver%n_singlet_states, 1)
            call mem%dealloc(orbital_differences, wf%n_es_amplitudes, 1)
!
            call mem%alloc(c_i, wf%n_es_amplitudes, 1)
!
            do trial = n_solutions_on_file + 1, solver%n_singlet_states
!
               c_i = zero
               c_i(lowest_orbital_differences_index(trial, 1), 1) = one
!
               call davidson%write_trial(c_i)
!
            enddo 
!
            call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
            call mem%dealloc(lowest_orbital_differences_index, solver%n_singlet_states, 1)
!
         endif
!
         call davidson%orthonormalize_trial_vecs()
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
!!    Sets precondition vector to orbital differences 
!!
      implicit none
!
!      class(davidson_cc_es) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_es_amplitudes, 1)
      call wf%get_es_orbital_differences(preconditioner, wf%n_es_amplitudes)
      call davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_es_amplitudes, 1)
!
   end subroutine set_precondition_vector_davidson_cc_es
!
!
   subroutine cleanup_davidson_cc_es(solver)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_es) :: solver
!
      write(output%unit, '(/t3,a,a,a)') 'Cleaning up ', trim(solver%tag), '.'
!
      call solver%destruct_energies()
      if (allocated(solver%start_vectors)) call mem%dealloc(solver%start_vectors, solver%n_singlet_states, 1)
!
   end subroutine cleanup_davidson_cc_es
!
!
   subroutine print_banner_davidson_cc_es(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(davidson_cc_es) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
      call long_string_print(solver%description2)
!
   end subroutine print_banner_davidson_cc_es
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
      integer :: n_specs, i, j, n_start_vecs
!
      character(len=100) :: line
!
      if (.not. requested_section('cc excited state')) then 
!
         call output%error_msg('number of excitations must be specified.')
!
      endif
!
      call move_to_section('cc excited state', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:19) == 'residual threshold:' ) then
!
            read(line(20:100), *) solver%residual_threshold
!
         elseif (line(1:17) == 'energy threshold:' ) then
!
            read(line(18:100), *) solver%eigenvalue_threshold
!
         elseif (line(1:15) == 'singlet states:' ) then
!
            read(line(16:100), *) solver%n_singlet_states
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         elseif (line(1:18) == 'right eigenvectors') then 
!
            solver%transformation = 'right'
!
         elseif (line(1:18) == 'left eigenvectors') then 
!
            solver%transformation = 'left'
!
         elseif (trim(line) == 'restart') then
!
            solver%restart = .true.
!
         elseif (line(1:14) == 'start vectors:') then
!
            line = line(15:100)
            line = remove_preceding_blanks(line)
            n_start_vecs = 0
!
            do j = 1, 86
!
               if (line(j:j) .ne. ' ') n_start_vecs = n_start_vecs + 1
!
            enddo
!
            call mem%alloc(solver%start_vectors, n_start_vecs, 1)
            read(line, *) solver%start_vectors
!
         endif
!
      enddo
!
      if (allocated(solver%start_vectors)) then
!
         if (n_start_vecs .ne. solver%n_singlet_states) then
!
            call output%error_msg('mismatch in number of start vectors and number of specified roots.')
!
         endif
!
      endif
!
   end subroutine read_settings_davidson_cc_es
!
!
   subroutine set_projection_vector_davidson_cc_es(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cc_es) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
!     Do nothing for regular excited states, but will be used in descendants
!     CVS and IP
!
      davidson%do_projection = .false.
!
      if (.false.) write(output%unit, *) wf%name_, solver%tag ! Hack to suppress unavoidable compiler warnings
!
   end subroutine set_projection_vector_davidson_cc_es
!
!
end module davidson_cc_es_class
