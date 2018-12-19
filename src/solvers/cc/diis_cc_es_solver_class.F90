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
      integer(i15) :: n_singlet_states
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
      procedure :: initialize_energies             => initialize_energies_diis_cc_es_solver
      procedure :: destruct_energies               => destruct_energies_diis_cc_es_solver   
!
      procedure :: initialize_residual_norms       => initialize_residual_norms_diis_cc_es_solver
      procedure :: destruct_residual_norms         => destruct_residual_norms_diis_cc_es_solver  
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
!
      call solver%read_settings()
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
      call solver%initialize_residual_norms()
      solver%residual_norms = one 
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
   end subroutine prepare_diis_cc_es_solver
!
!
   subroutine initialize_energies_diis_cc_es_solver(solver)
!!
!!    Initialize energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
!
      if (.not. allocated(solver%energies)) &
            call mem%alloc(solver%energies, solver%n_singlet_states, 1)
!
   end subroutine initialize_energies_diis_cc_es_solver
!
!
   subroutine destruct_energies_diis_cc_es_solver(solver)
!!
!!    Destruct energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
!
      if (allocated(solver%energies)) &
            call mem%dealloc(solver%energies, solver%n_singlet_states, 1)
!
   end subroutine destruct_energies_diis_cc_es_solver
!
!
   subroutine initialize_residual_norms_diis_cc_es_solver(solver)
!!
!!    Initialize residual_norms
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
!
      if (.not. allocated(solver%residual_norms)) &
            call mem%alloc(solver%residual_norms, solver%n_singlet_states, 1)
!
   end subroutine initialize_residual_norms_diis_cc_es_solver
!
!
   subroutine destruct_residual_norms_diis_cc_es_solver(solver)
!!
!!    Destruct residual_norms
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_es_solver) :: solver
!
      if (allocated(solver%residual_norms)) &
            call mem%dealloc(solver%residual_norms, solver%n_singlet_states, 1)
!
   end subroutine destruct_residual_norms_diis_cc_es_solver
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
      logical :: converged
      logical :: converged_eigenvalue
      logical :: converged_residual
!
      type(diis_tool), dimension(:), allocatable :: diis 
!
      integer(i15) :: iteration, trial, solution, state, I
!
      character(len=3) :: string_state
!
      real(dp) :: ddot
!
      real(dp), dimension(:,:), allocatable :: X, R, eps
!
!     Create DIIS tools vector & initialize them 
!
      allocate(diis(solver%n_singlet_states))
!
      do state = 1, solver%n_singlet_states
!  
         write(string_state, '(i3.3)') state
         call diis(state)%init('diis_cc_es_solver_' // string_state, wf%n_amplitudes, wf%n_amplitudes, 8)
!
      enddo 
!
!     Make initial guess on the eigenvectors X = [X1 X2 X3 ...]
!
      call mem%alloc(eps, wf%n_amplitudes, 1)
      call wf%get_orbital_differences(eps)
!
      call mem%alloc(X, wf%n_amplitudes, solver%n_singlet_states)
      call solver%set_start_vectors(wf, X, eps)
!
!     Enter iterative loop
!
      call mem%alloc(R, wf%n_amplitudes, solver%n_singlet_states)
!
      converged = .false.
      iteration = 1
      do while (.not. converged .and. (iteration .le. solver%max_iterations))
!
         write(output%unit,'(/t3,a25,i4)') 'Iteration:               ', iteration
         flush(output%unit)
!
         write(output%unit,'(t3,a)') 'Root     Eigenvalue (Re)        Residual norm'
         write(output%unit,'(t3,a)') '----------------------------------------------'
!
!        Construct residual vector R = AX - (X^T AX / X^X) X = [R1 R2 R3 ...],
!        then precondition, print information and do DIIS update 
!
         do state = 1, solver%n_singlet_states
!
            call wf%construct_excited_state_equation(X(:,state), R(:,state), solver%energies(state, 1))
!
            do I = 1, wf%n_amplitudes
!
               R(I, state) = R(I, state)/(eps(I,1)-solver%energies(state, 1))
!
            enddo
!
            solver%residual_norms(state, 1) = sqrt(ddot(wf%n_amplitudes, R(1,state), 1, R(1,state), 1))
!
            write(output%unit, *) state, solver%energies(state, 1), solver%residual_norms(state, 1)
!
            X(:,state) = X(:,state) + R(:,state)
            call diis(state)%update(R(:,state), X(:,state))
!
         enddo
!
         iteration = iteration + 1        
!
      enddo 
!
   end subroutine run_diis_cc_es_solver
!
!
   subroutine set_start_vectors_diis_cc_es_solver(solver, wf, R, orbital_differences)
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
      real(dp), dimension(wf%n_amplitudes, 1), intent(in)                          :: orbital_differences 
!
      real(dp), dimension(:,:), allocatable :: lowest_orbital_differences
!
      integer(i15), dimension(:,:), allocatable :: lowest_orbital_differences_index
!
      integer(i15) :: state
!
      call mem%alloc(lowest_orbital_differences, solver%n_singlet_states, 1)
      call mem%alloc_int(lowest_orbital_differences_index, solver%n_singlet_states, 1)
!
      call get_n_lowest(solver%n_singlet_states, wf%n_amplitudes, orbital_differences, &
                           lowest_orbital_differences, lowest_orbital_differences_index)
!
      call mem%dealloc(lowest_orbital_differences, solver%n_singlet_states, 1)
!
      do state = 1, solver%n_singlet_states
!
         R(:,state) = zero
         R(lowest_orbital_differences_index(state, 1), state) = one
!
      enddo 
!
      call mem%dealloc_int(lowest_orbital_differences_index, solver%n_singlet_states, 1)      
!
   end subroutine set_start_vectors_diis_cc_es_solver
!
!
end module diis_cc_es_solver_class
