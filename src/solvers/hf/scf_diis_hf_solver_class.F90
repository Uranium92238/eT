module scf_diis_hf_solver_class
!
!!
!!		Self-consistent field DIIS HF solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    A DIIS-accelerated Roothan-Hall self-consistent field solver. 
!!    In each Roothan-Hall update, a fitted F matrix is used instead 
!!    of the one produced from the previously obtained density. 
!!
!!    Supported wavefunctions: HF, UHF
!!    
!
   use kinds
   use diis_tool_class
   use timings_class
   use file_class
   use hf_class
   use disk_manager_class
   use abstract_hf_solver_class
!
   implicit none
!
   type, extends(abstract_hf_solver) :: scf_diis_hf_solver
!
      integer :: diis_dimension
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
      solver%tag = 'Self-consistent field DIIS Hartree-Fock solver'
      solver%author = 'E. F. Kjønstad and S, D. Folkestad, 2018'
      solver%description = 'A DIIS-accelerated Roothan-Hall self-consistent field solver. &
                                  &A least-square DIIS fit is performed on the previous Fock matrices and &
                                  &associated gradients. Following the Roothan-Hall update of the density, &
                                  &the DIIS-fitted Fock matrix is used to get the next orbital coefficients.' 
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
!     Read user's specified settings & set wavefunction screening based on them
!     (note that the screenings must be tighter for tighter gradient thresholds)
!     & print settings to output
!
      call solver%read_settings()
      call wf%set_screening_and_precision_thresholds(solver%gradient_threshold)
!
      write(output%unit, '(/t3,a/)') '- Hartree-Fock solver settings:'
!
      call solver%print_scf_diis_settings()
      call solver%print_hf_solver_settings()
      call wf%print_screening_settings()
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
      write(output%unit, '(t6,a29,i2)') 'DIIS dimension:               ', solver%diis_dimension
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
      integer :: iteration
!
      real(dp), dimension(:,:), allocatable :: F 
      real(dp), dimension(:,:), allocatable :: G 
      real(dp), dimension(:,:), allocatable :: ao_fock 
      real(dp), dimension(:,:), allocatable :: h_wx 
      real(dp), dimension(:,:), allocatable :: prev_ao_density 
!
      integer :: n_s
!
      real(dp), dimension(:,:), allocatable     :: sp_eri_schwarz
      integer, dimension(:,:), allocatable :: sp_eri_schwarz_list
!
      integer :: dim_gradient, dim_fock
!
      type(timings) :: iteration_timer, solver_timer 
!
!     :: Part I. Preparations. 
!
      call iteration_timer%init('SCF DIIS iteration time')
      call solver_timer%init('SCF DIIS solver time')
!
      call solver_timer%start()
!
!     Construct screening vectors for efficient Fock construction 
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     Initialize the DIIS manager object
!
      dim_fock     = ((wf%n_ao)*(wf%n_ao + 1)/2)*(wf%n_densities)
      dim_gradient = (wf%n_ao*(wf%n_ao - 1)/2)*(wf%n_densities)
!
      call diis_manager%init('hf_diis', dim_fock, dim_gradient, solver%diis_dimension)
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
      if (.not. solver%do_restart) then 
!         
         call wf%roothan_hall_update_orbitals() ! F => C
         call wf%update_ao_density()            ! C => D
!
         call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!
      endif
!
      call mem%alloc(ao_fock, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
      call mem%alloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
!
      call mem%alloc(G, wf%n_ao*(wf%n_ao - 1)/2, wf%n_densities)          
      call mem%alloc(F, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities) 
!
      call wf%get_packed_roothan_hall_gradient(G)
!
      max_grad = get_abs_max(G, dim_gradient)
!
      call wf%get_ao_fock(F)
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
         call iteration_timer%start()       
!
!        Set energy and print information for current iteration
!
         energy = wf%energy
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e11.4,4x,e11.4)') iteration, wf%energy, &
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
            call wf%print_wavefunction_summary()

!
         else
!
            prev_energy = wf%energy
            call wf%get_ao_density_sq(prev_ao_density)
!
            call wf%roothan_hall_update_orbitals()     ! DIIS F => C
            call wf%update_ao_density()                ! C => D
!
            if (iteration .ne. 1) call wf%set_ao_fock(ao_fock) ! Restore F 
!
            call wf%update_fock_and_energy_cumulative(sp_eri_schwarz, sp_eri_schwarz_list, n_s, prev_ao_density, h_wx)
!
            call wf%get_packed_roothan_hall_gradient(G)
!
            max_grad = get_abs_max(G, dim_gradient)
!
            call wf%get_ao_fock(F)
            call dcopy(dim_fock, F, 1, ao_fock, 1)
!
            call diis_manager%update(G, F)
            call wf%set_ao_fock(F)
!
         endif
!
         call iteration_timer%freeze()
         call iteration_timer%switch_off()       
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%dealloc(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call mem%dealloc(G, wf%n_ao*(wf%n_ao - 1)/2, wf%n_densities)          
      call mem%dealloc(F, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities) 
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(ao_fock, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
      call mem%dealloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
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
      call solver_timer%freeze()
      call solver_timer%switch_off()
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
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(solver%tag)
!
!     Do a final Roothan-Hall step to transform the Fock matrix in the canonical MO basis 
!
      do_mo_transformation = .true.
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies, do_mo_transformation)
!
!     Save the orbitals to file & store restart information 
!
      call wf%save_orbital_coefficients()
      call solver%write_restart_file(wf)
! ------- DEBUG
     flush(output%unit)
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
      write(solver%restart_file%unit, *) 'n_ao n_mo'
      write(solver%restart_file%unit, *) wf%n_ao
      write(solver%restart_file%unit, *) wf%n_mo
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
      integer :: n_ao, n_mo 
!
      write(output%unit, '(/t3,a)') '- Requested restart. Reading orbitals from file:'
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
      integer :: n_records, i 
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
