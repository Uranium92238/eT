module diis_cc_es_solver_class
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
!
   implicit none
!
   type :: diis_cc_es_solver
!
      character(len=100) :: tag = 'DIIS coupled cluster excited state solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      character(len=500) :: description1 = 'A DIIS solver that solves for the lowest eigenvalues and &
                                           &the right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                                           &problem is solved by DIIS extrapolation of residuals for each &
                                           &eigenvector until the convergence criteria are met.'
!
      integer(i15) :: max_iterations
!
      real(dp) :: eigenvalue_threshold  
      real(dp) :: residual_threshold  
!
      integer(i15) :: n_singlet_states, diis_dimension
!
      real(dp), dimension(:,:), allocatable :: energies
      real(dp), dimension(:,:), allocatable :: residual_norms 
!
      character(len=40) :: transformation 
!
      integer(i15), dimension(:,:), allocatable :: start_vectors
!
   contains
!     
      procedure, non_overridable :: prepare        => prepare_diis_cc_es_solver
      procedure, non_overridable :: run            => run_diis_cc_es_solver
      procedure, non_overridable :: cleanup        => cleanup_diis_cc_es_solver
!
      procedure :: set_start_vectors               => set_start_vectors_diis_cc_es_solver
!
      procedure :: print_banner                    => print_banner_diis_cc_es_solver
!
      procedure :: read_settings                   => read_settings_diis_cc_es_solver
      procedure :: print_settings                  => print_settings_diis_cc_es_solver
!
   end type diis_cc_es_solver
!
!
contains
!
!
   subroutine prepare_diis_cc_es_solver(solver)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
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
      solver%diis_dimension       = 8
!
      call solver%read_settings()
      call solver%print_settings()
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
   end subroutine prepare_diis_cc_es_solver
!
!
   subroutine print_settings_diis_cc_es_solver(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_es_solver) :: solver 
!
      write(output%unit, '(/t3,a)') '- DIIS CC excited state solver settings:'
!
      write(output%unit,'(/t6,a20,e9.2)') 'Energy threshold:   ', solver%eigenvalue_threshold
      write(output%unit,'(t6,a20,e9.2)')  'Residual threshold: ', solver%residual_threshold
      write(output%unit,'(/t6,a,i3,a)')   'Number of singlet states: ', solver%n_singlet_states
      write(output%unit, '(t6,a26,i3)')   'Max number of iterations: ', solver%max_iterations
      flush(output%unit)
!
   end subroutine print_settings_diis_cc_es_solver
!
!
   subroutine read_settings_diis_cc_es_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_es_solver) :: solver 
!
      integer(i15) :: n_specs, i
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
         elseif (line(1:15) == 'diis dimension:' ) then
!
            read(line(16:100), *) solver%diis_dimension
!
         elseif (line(1:18) == 'right eigenvectors') then 
!
            solver%transformation = 'right'
!
         elseif (line(1:18) == 'left eigenvectors') then 
!
            solver%transformation = 'left'
!
         endif
!
      enddo
!
   end subroutine read_settings_diis_cc_es_solver
!
!
   subroutine cleanup_diis_cc_es_solver(solver)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
!
!     Nothing here yet...
!
   end subroutine cleanup_diis_cc_es_solver
!
!
   subroutine print_banner_diis_cc_es_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(diis_cc_es_solver) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a)')
!
   end subroutine print_banner_diis_cc_es_solver
!
!
   subroutine run_diis_cc_es_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      logical, dimension(:), allocatable :: converged
!
      logical, dimension(:), allocatable :: converged_eigenvalue
      logical, dimension(:), allocatable :: converged_residual
!
      real(dp), dimension(:), allocatable :: energies 
      real(dp), dimension(:), allocatable :: prev_energies 
!
      real(dp), dimension(:), allocatable :: residual_norms
!
      type(diis_tool), dimension(:), allocatable :: diis 
!
      integer(i15) :: iteration, trial, solution, state, amplitude
!
      character(len=3) :: string_state
!
      real(dp) :: norm_X
!
      real(dp), dimension(:,:), allocatable :: X, R, eps
!
!     Initialize energies, residual norms, and convergence arrays 
!
      allocate(energies(solver%n_singlet_states))
      allocate(prev_energies(solver%n_singlet_states))
      allocate(residual_norms(solver%n_singlet_states))
!
      energies       = zero 
      prev_energies  = zero 
      residual_norms = zero 
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
         call diis(state)%init('diis_cc_es_solver_' // string_state, wf%n_amplitudes, wf%n_amplitudes, solver%diis_dimension)
!
      enddo 
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call mem%alloc(eps, wf%n_amplitudes, 1)
      call wf%get_orbital_differences(eps)
!
      call mem%alloc(X, wf%n_amplitudes, solver%n_singlet_states)
      call solver%set_start_vectors(wf, X, energies, eps)
!
!     Enter iterative loop
!
      call mem%alloc(R, wf%n_amplitudes, solver%n_singlet_states)
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
               call wf%construct_excited_state_equation(X(:,state), R(:,state), energies(state))
!
!$omp parallel do private(amplitude)
               do amplitude = 1, wf%n_amplitudes
!
                  R(amplitude, state) = -R(amplitude, state)/(eps(amplitude, 1) - energies(state))
!
               enddo
!$omp end parallel do 
!
!              Update convergence logicals 
!
               residual_norms(state) = get_l2_norm(R(:, state), wf%n_amplitudes)
!
               converged_eigenvalue(state) = abs(energies(state)-prev_energies(state)) .lt. solver%eigenvalue_threshold
               converged_residual(state)   = residual_norms(state)                     .lt. solver%residual_threshold
!
               converged(state) = converged_eigenvalue(state) .and. converged_residual(state)
!
!              Perform DIIS extrapolation to the optimal next guess for X,
!              then normalize it to avoid accumulating norm in X
!
               X(:,state) = X(:,state) + R(:,state)
               call diis(state)%update(R(:,state), X(:,state))
!
               norm_X = get_l2_norm(X(:,state), wf%n_amplitudes)
               X(:,state) = X(:,state)/norm_X
!
            endif 
!
            write(output%unit, '(i3,3x,f19.12,6x,e10.4)') state, energies(state), residual_norms(state)
            flush(output%unit)
!
         enddo
!
         prev_energies = energies 
!
         write(output%unit,'(t3,a)')  '----------------------------------------------'     
!
      enddo 
!
      if (all(converged)) then 
!
         write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
         ! call solver%print_summary() ... make this when printing routines from new-eT have been merged 
!
      endif 
!
      deallocate(energies)
      deallocate(prev_energies)
      deallocate(residual_norms)
!
      deallocate(converged)
      deallocate(converged_residual)
      deallocate(converged_eigenvalue)
!
   end subroutine run_diis_cc_es_solver
!
!
   subroutine set_start_vectors_diis_cc_es_solver(solver, wf, R, energies, orbital_differences)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(diis_cc_es_solver), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, solver%n_singlet_states), intent(inout) :: R 
      real(dp), dimension(solver%n_singlet_states), intent(inout)                  :: energies 
      real(dp), dimension(wf%n_amplitudes, 1), intent(in)                          :: orbital_differences 
!
      real(dp), dimension(:,:), allocatable :: lowest_orbital_differences
!
      integer(i15), dimension(:,:), allocatable :: lowest_orbital_differences_index
!
      integer(i15) :: state
!
      call mem%alloc(lowest_orbital_differences, solver%n_singlet_states, 1)
      call mem%alloc(lowest_orbital_differences_index, solver%n_singlet_states, 1)
!
      call get_n_lowest(solver%n_singlet_states, wf%n_amplitudes, orbital_differences, &
                           lowest_orbital_differences, lowest_orbital_differences_index)
!
      do state = 1, solver%n_singlet_states
!
         R(:,state) = zero
         R(lowest_orbital_differences_index(state, 1), state) = one
         energies(state) = lowest_orbital_differences(state, 1)
!
      enddo 
!
      call mem%dealloc(lowest_orbital_differences, solver%n_singlet_states, 1)
      call mem%dealloc(lowest_orbital_differences_index, solver%n_singlet_states, 1)      
!
   end subroutine set_start_vectors_diis_cc_es_solver
!
!
end module diis_cc_es_solver_class
