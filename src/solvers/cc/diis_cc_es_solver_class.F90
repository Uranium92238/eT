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
      procedure :: print_banner                    => print_banner_diis_cc_es_solver
!
      procedure :: read_settings                   => read_settings_diis_cc_es_solver
      procedure :: print_settings                  => print_settings_diis_cc_es_solver
!
      procedure :: initialize_energies             => initialize_energies_diis_cc_es_solver
      procedure :: destruct_energies               => destruct_energies_diis_cc_es_solver   
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
      type(diis_tool) :: diis 
!
      integer(i15) :: iteration, trial, solution
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      write(output%unit, *) 'Hello world!'
!
   end subroutine run_diis_cc_es_solver
!
!
end module diis_cc_es_solver_class
