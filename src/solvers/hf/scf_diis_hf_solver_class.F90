module scf_diis_hf_solver_class
!
!!
!!		Self-consistent field DIIS HF solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    A DIIS-accelerated Roothan-Hall self-consistent field solver. 
!!    In other words, it does a least-square fit to a zero gradient 
!!    using the previously recorded Fock matrices and associated 
!!    gradients. In each Roothan-Hall update, the fitted F matrix 
!!    is used instead of the one produced from the previously 
!!    obtained density matrix D. 
!!
!!    Supported wavefunctions: HF 
!!
!!    Note to developers: although steps have been taken to make 
!!    the solver general enough to allow for unrestricted HF theory,
!!    there are still explicit references to the AO Fock and density.
!!    These need to be replaced by AO-densities and Fock-matrices,
!!    where it might be reasonable to adopt a "number of densities"
!!    variable, such that we work with F = [F_a F_b], G = [G_a G_b] and 
!!    D = [D_a D_b] in the DIIS solver. Generalized and overwritten
!!    wavefunction routines can then control the handling of these
!!    arrays and their expected size (when constructing Fock
!!    cumulatively, for instance, where the previous density 
!!    is required). - Eirik, Sep 2018
!!    
!
   use kinds
   use diis_tool_class
   use file_class
   use hf_class
   use disk_manager_class
   use abstract_hf_solver_class
!
   implicit none
!
   type, extends(abstract_hf_solver) :: scf_diis_hf_solver
!
      integer(i15) :: diis_dimension
!
      logical :: converged 
!
      type(file) :: restart_file
      logical    :: do_restart 
!
   contains
!     
      procedure :: prepare                => prepare_scf_diis_hf_solver
      procedure :: run                    => run_scf_diis_hf_solver
      procedure :: cleanup                => cleanup_scf_diis_hf_solver
!
      procedure :: print_banner           => print_banner_scf_diis_hf_solver
      procedure :: print_summary          => print_summary_scf_diis_hf_solver
!
      procedure :: read_settings           => read_settings_scf_diis_hf_solver 
      procedure :: read_scf_diis_settings  => read_scf_diis_settings_scf_diis_hf_solver
!
      procedure :: print_scf_diis_settings => print_scf_diis_settings_scf_diis_hf_solver
      procedure :: write_restart_file      => write_restart_file_scf_diis_hf_solver
!
      procedure :: restart                 => restart_scf_diis_hf_solver 
!
   end type scf_diis_hf_solver
!
!
contains
!
!
   subroutine prepare_scf_diis_hf_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_hf_solver) :: solver
!
      class(hf) :: wf
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%do_restart         = .false.
      solver%diis_dimension     = 8
      solver%max_iterations     = 100
      solver%ao_density_guess   = 'SAD'
      solver%energy_threshold   = 1.0D-6
      solver%gradient_threshold = 1.0D-6
!
!     Read user's specified settings
!
      call solver%read_settings()
!
      call wf%set_screening_and_precision_thresholds(solver%gradient_threshold)
!
      write(output%unit, '(/t3,a/)') '- Hartree-Fock solver settings:'
!
      call solver%print_scf_diis_settings()
      call solver%print_hf_solver_settings()
!
!     Initialize the orbitals, density, and the Fock matrix (or matrices)
!
      call wf%initialize_fock()
      call wf%initialize_density()
      call wf%initialize_orbitals()
!
!     Prepare restart information file 
!
      call solver%restart_file%init('scf_diis_restart_info', 'sequential', 'formatted')
!
      if (solver%do_restart) then 
!
         call solver%restart(wf)
!
      else 
!
         write(output%unit, '(/t3,a,a,a)') '- Setting initial AO density to ', trim(solver%ao_density_guess), ':'
         call wf%set_initial_ao_density_guess(solver%ao_density_guess)
!
      endif
!
   end subroutine prepare_scf_diis_hf_solver
!
!
   subroutine print_scf_diis_settings_scf_diis_hf_solver(solver)
!!
!!    Print SCF-DIIS settings    
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(scf_diis_hf_solver) :: solver 
!
      write(output%unit, '(t6,a29,i2)') 'DIIS dimension:              ', solver%diis_dimension
!
   end subroutine print_scf_diis_settings_scf_diis_hf_solver
!
!
   subroutine run_scf_diis_hf_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_hf_solver) :: solver
!
      class(hf) :: wf
!
      type(diis_tool) :: diis_manager
!
      logical :: converged_energy
      logical :: converged_gradient
!
      real(dp) :: max_grad, energy, prev_energy, n_electrons
!
      real(dp) :: ddot
!
      integer(i15) :: iteration
!
      real(dp), dimension(:,:), allocatable :: D
      real(dp), dimension(:,:), allocatable :: F 
      real(dp), dimension(:,:), allocatable :: Po 
      real(dp), dimension(:,:), allocatable :: Pv 
      real(dp), dimension(:,:), allocatable :: G 
      real(dp), dimension(:,:), allocatable :: ao_fock 
      real(dp), dimension(:,:), allocatable :: h_wx 
      real(dp), dimension(:,:), allocatable :: prev_ao_density 
!
      integer(i15) :: n_s, i
!
      real(dp), dimension(:,:), allocatable     :: sp_eri_schwarz
      integer(i15), dimension(:,:), allocatable :: sp_eri_schwarz_list
!
!     :: Part I. Preparations. 
!
     ! write(output%unit, '(/t3,a)') ':: Running SCF-DIIS object'
!
!     Construct screening vectors for efficient Fock construction 
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     Initialize the DIIS manager object
!
      call diis_manager%init('hf_diis', (wf%n_ao)*(wf%n_ao + 1)/2, wf%n_ao**2, solver%diis_dimension)
!
!     Set the initial density guess and Fock matrix 
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(h_wx)
!
      call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!
      call wf%get_n_electrons_in_density(n_electrons)
!
      write(output%unit, '(/t6,a30,f17.12)') 'Energy of initial guess:      ', wf%energy
      write(output%unit, '(t6,a30,f17.12)')  'Number of electrons in guess: ', n_electrons
!
!     Do a Roothan-Hall update to ensure idempotentency of densities,
!     and use it to construct the first proper Fock matrix from which 
!     to begin cumulative construction 
!
      call wf%roothan_hall_update_orbitals() ! F => C
      call wf%update_ao_density()            ! C => D
!
      call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!
      call mem%alloc(ao_fock, wf%n_ao, wf%n_ao)          ! Holds Fock matrix temporarily
      call mem%alloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
      call mem%alloc(G, wf%n_ao, wf%n_ao)                ! Gradient 
      call mem%alloc(F, wf%n_ao*(wf%n_ao + 1)/2, 1)      ! Fock matrix packed 
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
!
      call wf%construct_projection_matrices(Po, Pv)
      call wf%construct_roothan_hall_gradient(G, Po, Pv)
!
      max_grad = get_abs_max(G, (wf%n_ao)**2)
!
      call packin(F, wf%ao_fock, wf%n_ao)
      call diis_manager%update(G, F)
!
!     Part II. Iterative SCF loop.
!
      solver%converged   = .false.
      converged_energy   = .false.
      converged_gradient = .false.
!
      prev_energy = zero
!
      write(output%unit, '(/t3,a)') 'Iteration    Energy (a.u.)        Max(grad.)    Delta E (a.u.)'
      write(output%unit, '(t3,a)')  '--------------------------------------------------------------'
!
      iteration = 1
!
      do while (.not. solver%converged .and. iteration .le. solver%max_iterations)         
!
!        Set energy and print information for current iteration
!
         energy = wf%energy
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e10.4,4x,e10.4)') iteration, wf%energy, &
                                          max_grad, abs(wf%energy-prev_energy)
         flush(output%unit)
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_gradient = max_grad                .lt. solver%gradient_threshold
!
         solver%converged = converged_gradient .and. converged_energy
!
         if (converged_gradient .and. iteration .eq. 1) solver%converged = .true.
!
         if (solver%converged) then
!
            write(output%unit, '(t3,a)')          '--------------------------------------------------------------'
            write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
!
            if (.not. converged_energy) then 
!
               write(output%unit, '(/t3,a,/t9,a)') 'Note: the gradient converged in the first iteration,', &
                                                          'so the energy convergence has not been tested!'
!
            endif
!
            call solver%print_summary(wf)
!
         else
!
            prev_energy     = wf%energy
            prev_ao_density = wf%ao_density
!
            call wf%roothan_hall_update_orbitals()     ! DIIS F => C
            call wf%update_ao_density()                ! C => D
!
            if (iteration .ne. 1) wf%ao_fock = ao_fock ! Restore F 
!
            call wf%update_fock_and_energy_cumulative(sp_eri_schwarz, sp_eri_schwarz_list, n_s, prev_ao_density, h_wx)
!
            call wf%construct_projection_matrices(Po, Pv)
            call wf%construct_roothan_hall_gradient(G, Po, Pv)
!
            max_grad = get_abs_max(G, (wf%n_ao)**2)
!
            call packin(F, wf%ao_fock, wf%n_ao)
            call diis_manager%update(G, F)
!
            ao_fock = wf%ao_fock  
            call wf%set_ao_fock(F)
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(G, wf%n_ao, wf%n_ao)           
      call mem%dealloc(F, wf%n_ao*(wf%n_ao + 1)/2, 1) 
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call diis_manager%finalize()
!
      if (.not. solver%converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         stop
!
      endif 
!
   end subroutine run_scf_diis_hf_solver
!
!
   subroutine cleanup_scf_diis_hf_solver(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_hf_solver) :: solver
!
      class(hf) :: wf
!
      logical :: do_mo_transformation
!
      integer(i15) :: i
!
!     Do a final Roothan-Hall step to transform the Fock matrix in the canonical MO basis 
!
      do_mo_transformation = .true.
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies, do_mo_transformation)
!
!     Save the orbitals to file & store restart information 
!
      call wf%save_orbital_coefficients()
!
      call solver%write_restart_file(wf)
!
   end subroutine cleanup_scf_diis_hf_solver
!
!
   subroutine write_restart_file_scf_diis_hf_solver(solver, wf)
!!
!!    Write restart file 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2018     
!!
!!    Stores information about the current calculation which may be 
!!    used, e.g., to spot inconsistencies on restarting from files 
!!    present in the folder.
!!
      implicit none 
!
      class(scf_diis_hf_solver), intent(inout) :: solver 
!
      class(hf), intent(in) :: wf 
!
      call disk%open_file(solver%restart_file, 'write', 'rewind')
!
      write(solver%restart_file%unit, *) 'n_ao n_mo energy_threshold gradient_threshold converged'
      write(solver%restart_file%unit, *) wf%n_ao
      write(solver%restart_file%unit, *) wf%n_mo
      write(solver%restart_file%unit, *) solver%energy_threshold
      write(solver%restart_file%unit, *) solver%gradient_threshold
      write(solver%restart_file%unit, *) solver%converged
!
      call disk%close_file(solver%restart_file) 
!
   end subroutine write_restart_file_scf_diis_hf_solver
!
!
   subroutine restart_scf_diis_hf_solver(solver, wf)
!!
!!    Restart 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(scf_diis_hf_solver), intent(in) :: solver 
!
      class(hf), intent(inout) :: wf 
!
      integer(i15) :: n_ao, n_mo 
!
!     Sanity checks 
!
      call disk%open_file(solver%restart_file, 'read', 'rewind')  
!
      read(solver%restart_file%unit, *) ! Empty read to skip banner 
!
      read(solver%restart_file%unit, *) n_ao 
      read(solver%restart_file%unit, *) n_mo
!
      call disk%close_file(solver%restart_file)
!
      if (n_ao .ne. wf%n_ao .or. n_mo .ne. wf%n_mo) call output%error_msg('Inconsistent dimensions on restart in SCF-DIIS.') 
!
!     Do restart 
!
      call wf%read_orbital_coefficients()
      call wf%update_ao_density()
!
   end subroutine restart_scf_diis_hf_solver
!
!
   subroutine print_banner_scf_diis_hf_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_diis_hf_solver) :: solver 
!
      write(output%unit, '(//t3,a)') ':: Self-consistent field DIIS Hartree-Fock solver'
      write(output%unit, '(t3,a/)')  ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(t3,a)')  'A DIIS-accelerated Roothan-Hall self-consistent field solver.'
      write(output%unit, '(t3,a)')  'In other words, a least-square fit toward a zero gradient vector' 
      write(output%unit, '(t3,a)')  'is performed using the previously recorded Fock matrices and the'
      write(output%unit, '(t3,a)')  'associated gradients. After each Roothan-Hall update of the density,'
      write(output%unit, '(t3,a)')  'a DIIS-fitted Fock matrix is used to get the next orbital coefficients,'
      write(output%unit, '(t3,a)')  'instead of the one produced directly from the AO density matrix.'

      flush(output%unit)
!
   end subroutine print_banner_scf_diis_hf_solver
!
!
   subroutine print_summary_scf_diis_hf_solver(solver, wf)
!!
!!    Print summary 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_diis_hf_solver) :: solver 
!
      class(hf) :: wf 
!
      call wf%print_wavefunction_summary()
!
   end subroutine print_summary_scf_diis_hf_solver
!
!
   subroutine read_settings_scf_diis_hf_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(scf_diis_hf_solver) :: solver 
!
      call solver%read_hf_solver_settings()
      call solver%read_scf_diis_settings()
!
   end subroutine read_settings_scf_diis_hf_solver
!
!
   subroutine read_scf_diis_settings_scf_diis_hf_solver(solver)
!!
!!    Read SCF DIIS settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads settings specific to the class. 
!!
      implicit none 
!
      class(scf_diis_hf_solver) :: solver 
!
      integer(i15) :: n_records, i 
!
      character(len=100) :: line, value 
!
      if (requested_section('hf')) then ! User has requested something 
!
         call move_to_section('hf', n_records)
!
         do i = 1, n_records
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:15) == 'diis dimension:') then
!
               value = line(16:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%diis_dimension
               cycle
!
            elseif (line(1:7) == 'restart') then
!
               solver%do_restart = .true.
               cycle
!
            endif 
!
         enddo
!
      endif 
!
   end subroutine read_scf_diis_settings_scf_diis_hf_solver
!
!
end module scf_diis_hf_solver_class
