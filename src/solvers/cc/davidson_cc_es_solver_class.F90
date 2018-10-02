module davidson_cc_es_solver_class
!
!!
!!    Davidson coupled cluster excited state solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    
!
   use kinds
   use file_class
   use ccs_class
   use eigen_davidson_tool_class
!
   implicit none
!
   type :: davidson_cc_es_solver
!
      integer(i15) :: max_iterations = 100
!
      real(dp) :: eigenvalue_threshold = 1.0d-6
      real(dp) :: residual_threshold   = 1.0d-6
!
      logical      :: do_restart = .false.
!
      integer(i15) :: n_singlet_states = 0
!
      character(len=40) :: transformation 
!
      real(dp), dimension(:,:), allocatable :: energies
!
      integer(i15), dimension(:,:), allocatable :: start_vectors
!
      type(file) :: restart_file
!
   contains
!     
      procedure, non_overridable :: prepare  => prepare_davidson_cc_es_solver
      procedure, non_overridable :: run      => run_davidson_cc_es_solver
      procedure, non_overridable :: cleanup  => cleanup_davidson_cc_es_solver
!
      procedure :: print_banner              => print_banner_davidson_cc_es_solver
      procedure :: print_summary             => print_summary_davidson_cc_es_solver
!
      procedure :: read_settings             => read_settings_davidson_cc_es_solver
!
      procedure :: print_settings            => print_settings_davidson_cc_es_solver
!
      procedure :: set_start_vectors         => set_start_vectors_davidson_cc_es_solver
      procedure :: set_precondition_vector   => set_precondition_vector_davidson_cc_es_solver
      procedure :: set_projection_vector     => set_projection_vector_davidson_cc_es_solver
!
      procedure :: transform_trial_vector    => transform_trial_vector_davidson_cc_es_solver
!       
      procedure :: initialize_energies       => initialize_energies_davidson_cc_es_solver
      procedure :: destruct_energies         => destruct_energies_davidson_cc_es_solver   
!
      procedure :: restart                   => restart_davidson_cc_es_solver 
      procedure :: write_restart_file        => write_restart_file_cc_es_solver
!
   end type davidson_cc_es_solver
!
!
contains
!
!
   subroutine prepare_davidson_cc_es_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%max_iterations       = 100
      solver%eigenvalue_threshold = 1.0d-6
      solver%residual_threshold   = 1.0d-6
      solver%transformation       = 'right'
      solver%do_restart           = .false.
!
      call solver%read_settings()
!
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
      call solver%restart_file%init('davidson_cc_es_restart_info', 'sequential', 'formatted')
!
   end subroutine prepare_davidson_cc_es_solver
!
!
   subroutine write_restart_file_cc_es_solver(solver)
!!
!!    Write restart 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2018
!!
      implicit none 
!
      class(davidson_cc_es_solver) :: solver 
!
      call disk%open_file(solver%restart_file, 'write', 'rewind')
!
      write(solver%restart_file%unit, *) 'n_singlet_states'
      write(solver%restart_file%unit, *) solver%n_singlet_states
!
      call disk%close_file(solver%restart_file) 
!
   end subroutine write_restart_file_cc_es_solver
!
!
   subroutine restart_davidson_cc_es_solver(solver, davidson)
!!
!!    Restart 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2018
!!
      implicit none 
!
      class(davidson_cc_es_solver) :: solver 
!
      class(eigen_davidson_tool) :: davidson
!
      integer(i15) :: n_solutions_on_file 
!
!     Read in the number of solutions to restart from - according the restart file 
!
      call disk%open_file(solver%restart_file, 'read', 'rewind')
!
      n_solutions_on_file = 0
      read(solver%restart_file%unit, *) ! Empty read to skip banner
      read(solver%restart_file%unit, *) n_solutions_on_file
!
      call disk%close_file(solver%restart_file) 
!
!     Avoid reading too many solutions if fewer are requested than previously converged 
!
      if (n_solutions_on_file .gt. solver%n_singlet_states) then 
!
         n_solutions_on_file = solver%n_singlet_states
!
      endif 
!
!     Ask Davidson to restart - use the previous solutions as trial vectors 
!
      call davidson%restart_from_solutions(n_solutions_on_file)
!
!     For the remaining states, use orbital differences 
!
!     Todo... 
!
   end subroutine restart_davidson_cc_es_solver
!
!
   subroutine initialize_energies_davidson_cc_es_solver(solver)
!!
!!    Initialize energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initialize excitation energies
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      if (.not. allocated(solver%energies)) &
            call mem%alloc(solver%energies, solver%n_singlet_states, 1)
!
   end subroutine initialize_energies_davidson_cc_es_solver
!
!
   subroutine destruct_energies_davidson_cc_es_solver(solver)
!!
!!    Destruct energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Destruct excitation energies
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      if (allocated(solver%energies)) &
            call mem%dealloc(solver%energies, solver%n_singlet_states, 1)
!
   end subroutine destruct_energies_davidson_cc_es_solver
!
!
   subroutine print_settings_davidson_cc_es_solver(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(davidson_cc_es_solver) :: solver 
!
      write(output%unit, '(/t3,a)') '- Davidson CC excited state solver settings:'
!
      write(output%unit,'(/t6,a20,e9.2)') 'Energy threshold:   ', solver%eigenvalue_threshold
      write(output%unit,'(t6,a20,e9.2)')  'Residual threshold: ', solver%residual_threshold
      write(output%unit,'(/t6,a,i3,a)')   'Number of singlet states: ', solver%n_singlet_states
      write(output%unit, '(t6,a26,i3)')   'Max number of iterations: ', solver%max_iterations
      flush(output%unit)
!
   end subroutine print_settings_davidson_cc_es_solver
!
!
   subroutine run_davidson_cc_es_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      logical :: converged
      logical :: converged_eigenvalue
      logical :: converged_residual
!
      type(eigen_davidson_tool) :: davidson
!
      integer(i15) :: iteration, trial, solution, i
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      converged            = .false. 
      converged_eigenvalue = .false. 
      converged_residual   = .false. 
!
      iteration = 1
!
     call davidson%prepare(wf%name // '_es_davidson', wf%n_amplitudes, solver%n_singlet_states, &
                               solver%residual_threshold, solver%eigenvalue_threshold)
!
!     Construct first trial vectors
!
      if (solver%do_restart) then 
!
         call solver%restart(davidson)
!
      else
!
         call solver%set_start_vectors(wf, davidson)
!
      endif 
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
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, wf%n_amplitudes, 1)
!
         write(output%unit, *) 'n new trials:', davidson%n_new_trials
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
         call mem%dealloc(c_i, wf%n_amplitudes, 1)
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
            call davidson%construct_next_trial_vec(residual_norm, iteration, solution)
            write(output%unit,'(t3,i2,5x,f16.12,7x,f16.12,11x,e10.4)') &
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
         if (converged_residual) then ! Converged residual
!
!           Tests for convergence of energy or restart
!
            if (converged_eigenvalue) then
!
              converged = .true.
!
            elseif (iteration .eq. 1 .and. wf%name .eq. 'ccs') then
!
                  converged = .true.
                  write(output%unit,'(/t3,a,/t3,a)') 'Note: residual converged in first iteration.', &
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
      elseif (.not. converged .and. iteration == solver%max_iterations) then
!
         write(output%unit,'(/t3,a)') 'Maximal number of iterations performed without reaching convergence!'
!
      endif
!
      call davidson%cleanup()
      call solver%print_summary(wf)
!
   end subroutine run_davidson_cc_es_solver
!
!
   subroutine transform_trial_vector_davidson_cc_es_solver(solver, wf, c_i)
!!
!!    Transform trial vector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Transforms the trial vector according to specified transformation routine.
!!
      class(davidson_cc_es_solver), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c_i
!
      if (trim(solver%transformation) == 'right') then 
!
         call wf%jacobi_transform_trial_vector(c_i)
!
      elseif (trim(solver%transformation) == 'left') then 
!
         call wf%jacobi_transpose_transform_trial_vector(c_i)
!
      endif 
!
   end subroutine transform_trial_vector_davidson_cc_es_solver
!
!
   subroutine set_start_vectors_davidson_cc_es_solver(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: orbital_differences
      real(dp), dimension(:,:), allocatable :: lowest_orbital_differences
!
      integer(i15), dimension(:,:), allocatable :: lowest_orbital_differences_index
!
      integer(i15) :: trial
!
      if (allocated(solver%start_vectors)) then
!
!        Initial trial vectors given on input
!
         call mem%alloc(c_i, wf%n_amplitudes, 1)
!
         call davidson%rewind_trials()
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
         call mem%dealloc(c_i, wf%n_amplitudes, 1)
!
      else
!
!        Initial trial vectors given by Koopman
!
         call mem%alloc(orbital_differences, wf%n_amplitudes, 1)
         call wf%get_orbital_differences(orbital_differences)
!
         call mem%alloc(lowest_orbital_differences, solver%n_singlet_states, 1)
         call mem%alloc_int(lowest_orbital_differences_index, solver%n_singlet_states, 1)
!
         call get_n_lowest(solver%n_singlet_states, wf%n_amplitudes, orbital_differences, &
                           lowest_orbital_differences, lowest_orbital_differences_index)
!
         call mem%dealloc(lowest_orbital_differences, solver%n_singlet_states, 1)
         call mem%dealloc(orbital_differences, wf%n_amplitudes, 1)
!
         call mem%alloc(c_i, wf%n_amplitudes, 1)
!
         call davidson%rewind_trials()
!
         do trial = 1, solver%n_singlet_states
!
            c_i = zero
            c_i(lowest_orbital_differences_index(trial, 1), 1) = one
!
            call davidson%write_trial(c_i)
!
         enddo 
!
         call mem%dealloc(c_i, wf%n_amplitudes, 1)
         call mem%dealloc_int(lowest_orbital_differences_index, solver%n_singlet_states, 1)
!
      endif
!
   end subroutine set_start_vectors_davidson_cc_es_solver
!
!
   subroutine set_precondition_vector_davidson_cc_es_solver(solver, wf, davidson)
!!
!!    Set precondition vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets precondition vector to orbital differences 
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_amplitudes, 1)
      call wf%get_orbital_differences(preconditioner)
      call davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_amplitudes, 1)
!
   end subroutine set_precondition_vector_davidson_cc_es_solver
!
!
   subroutine set_projection_vector_davidson_cc_es_solver(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
!     Do nothing for regular excited states
!
   end subroutine set_projection_vector_davidson_cc_es_solver
!
!
   subroutine cleanup_davidson_cc_es_solver(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      call solver%write_restart_file() 
!
   end subroutine cleanup_davidson_cc_es_solver
!
!
   subroutine print_banner_davidson_cc_es_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(davidson_cc_es_solver) :: solver 
!
      write(output%unit, '(//t3,a)') ':: Davidson coupled cluster excited state solver'
      write(output%unit, '(t3,a)')   ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(/t3,a)')  'A Davidson solver that calculates the lowest eigenvalues and the'
      write(output%unit, '(t3,a)')   'right eigenvectors of the Jacobian matrix, A. The eigenvalue problem'
      write(output%unit, '(t3,a)')   'is solved in a reduced space, the dimension of which is expanded'
      write(output%unit, '(t3,a)')   'until the convergence criteria are met.'
!
      write(output%unit, '(/t3,a)')  'A complete description of the algorithm can be found in'
      write(output%unit, '(t3,a)')   'E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      flush(output%unit)
!
   end subroutine print_banner_davidson_cc_es_solver
!
!
   subroutine print_summary_davidson_cc_es_solver(solver, wf)
!!
!!    Print summary 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(davidson_cc_es_solver) :: solver 
!
      class(ccs) :: wf 
!
    !  call wf%print_wavefunction_summary()
    !  call solver%print_solver_summary()
!
   end subroutine print_summary_davidson_cc_es_solver
!
!
   subroutine read_settings_davidson_cc_es_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_es_solver) :: solver 
!
      integer(i15) :: n_specs, i, j, n_start_vecs
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
         elseif (line(1:25) == 'number of singlet states:' ) then
!
            read(line(26:100), *) solver%n_singlet_states
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
            solver%do_restart = .true.
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
            call mem%alloc_int(solver%start_vectors, n_start_vecs, 1)
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
   end subroutine read_settings_davidson_cc_es_solver
!
!
end module davidson_cc_es_solver_class
