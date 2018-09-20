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
      integer(i15) :: max_iterations = 50
!
      real(dp) :: eigenvalue_threshold = 1.0d-6
      real(dp) :: residual_threshold  = 1.0d-6
!
      logical :: restart = .false.
!
      integer(i15) :: n_singlet_states = 0
!
      real(dp), dimension(:,:), allocatable :: energies
!
      integer(i15), dimension(:,:), allocatable :: start_vectors
!
   contains
!     
      procedure :: prepare        => prepare_davidson_cc_es_solver
      procedure :: run            => run_davidson_cc_es_solver
      procedure :: cleanup        => cleanup_davidson_cc_es_solver
!
      procedure :: print_banner   => print_banner_davidson_cc_es_solver
      procedure :: print_summary  => print_summary_davidson_cc_es_solver
!
      procedure :: read_settings  => read_settings_davidson_cc_es_solver
!
      procedure :: print_settings => print_settings_davidson_cc_es_solver
!
       procedure :: set_start_vectors => set_start_vectors_davidson_cc_es_solver
!       
       procedure :: initialize_energies => initialize_energies_davidson_cc_es_solver
       procedure :: destruct_energies => destruct_energies_davidson_cc_es_solver
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
      solver%n_singlet_states = 4 ! FOR TESTING SHOULD BE SET ON INPUT
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%max_iterations = 50
!
      solver%eigenvalue_threshold = 1.0d-6
      solver%residual_threshold  = 1.0d-6
!
      solver%restart = .false.
!
      if (requested_section('excited state')) then
!
         call solver%read_settings()
!
      else
!
         call output%error_msg('number of excitations must be specified.')
!
      endif
!
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
   end subroutine prepare_davidson_cc_es_solver
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
      write(output%unit,'(/t6,a20,e9.2)') 'Energy threshold:   ', solver%eigenvalue_threshold
      write(output%unit,'(t6,a20,e9.2)') 'Residual threshold: ', solver%residual_threshold
      write(output%unit,'(/t6,a,i3,a)')  'Number of singlet states: ', solver%n_singlet_states
      write(output%unit, '(t6,a26,i3)') 'Max number of iterations: ', solver%max_iterations
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
      type(eigen_davidson_tool) :: davidson
!
      logical :: converged
      logical :: converged_eigenvalue
      logical :: converged_residual
!
      integer(i15) :: iteration, trial, solution, i
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:,:), allocatable :: preconditioner
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
      if (solver%restart) then
!
!        Restart
!
      else
!
         call solver%set_start_vectors(wf, davidson)
!
      endif
!
!     Prepare preconditioner
!
      call mem%alloc(preconditioner, wf%n_amplitudes, 1)
      call wf%get_orbital_differences(preconditioner)
      call davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_amplitudes, 1)
!
!     Enter iterative loop
!
      do while (.not. converged .and. (iteration .le. solver%max_iterations))
!
         write(output%unit,'(/t3,a,i3)') 'Iteration: ', iteration
         write(output%unit,'(t3,a,i4/)') 'Reduced space dimension: ', davidson%dim_red
         flush(output%unit)
!
         write(output%unit,'(t3,a)') 'Root     Eigenvalue (Re)        Eigenvalue (Im)      Residual norm'
         write(output%unit,'(t3,a)') '-------------------------------------------------------------------'
!
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, wf%n_amplitudes, 1)
!
         do trial = davidson%dim_red - davidson%n_new_trials + 1, davidson%dim_red
!
            call davidson%read_trial(c_i, trial)
!
            call wf%transform_trial_vector(c_i)
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
            elseif (iteration .eq. 1 .and. wf%name .eq. 'CCS') then
!
                  converged = .true.
                  write(output%unit,'(//t3,a,/t3,a)')'Note: residual converged in first iteration.', &
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
         write(output%unit,'(/t3,a, i3, a)') 'Convergence criterion met in ', iteration, ' iterations!'
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
   subroutine set_start_vectors_davidson_cc_es_solver(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
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
         c_i = zero
         c_i(solver%start_vectors(1, 1), 1) = one
!
         call davidson%write_trial(c_i, 'rewind')
!
         do trial = 2, solver%n_singlet_states
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
         c_i = zero
         c_i(lowest_orbital_differences_index(1, 1), 1) = one
!
         call davidson%write_trial(c_i, 'rewind')
!
         do trial = 2, solver%n_singlet_states
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
!     Nothing here yet 
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
      write(output%unit, '(t3,a/)')  ':: E. F. Kjønstad, S. D. Folkestad, 2018'

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
      integer(i15) :: n_specs, i
!
      character(len=100) :: line
!
      call move_to_section('excited state', n_specs)
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
         elseif (trim(line) == 'restart') then
!
            solver%restart = .true.
!
         endif
!
      enddo
!
   end subroutine read_settings_davidson_cc_es_solver
!
!
end module davidson_cc_es_solver_class
