module scf_diis_solver_class
!
!!
!!		Self-consistent field DIIS solver class module
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
   use hf_solver_class
   use diis_tool_class
   use file_class
   use hf_class
   use disk_manager_class
!
   implicit none
!
   type, extends(hf_solver) :: scf_diis_solver
!
      integer(i15) :: diis_dimension = 8
!
   contains
!
      procedure :: run                    => run_scf_diis_solver
      procedure :: cleanup                => cleanup_scf_diis_solver
!
      procedure :: print_banner           => print_banner_scf_diis_solver
      procedure :: print_summary          => print_summary_scf_diis_solver
!
      procedure :: read_settings          => read_settings_scf_diis_solver 
      procedure :: read_scf_diis_settings => read_scf_diis_settings_scf_diis_solver
!
   end type scf_diis_solver
!
!
contains
!
!
   subroutine run_scf_diis_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
      type(diis_tool) :: diis_manager
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_residual
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
!     Print solver banner
!
      call solver%print_banner()
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
!     Initialize the orbitals, density, and the Fock matrix (or matrices),
!     and set the initial density guess and Fock matrix 
!
      call wf%initialize_fock()
      call wf%initialize_density()
      call wf%initialize_orbitals()
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(h_wx)
!
      write(output%unit, '(/t3,a,a,a)') 'Initial AO density is ', trim(solver%ao_density_guess), '.'
!
      call wf%set_initial_ao_density_guess(solver%ao_density_guess)
!
      call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s,  &
                                 h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                 solver%coulomb_precision)
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
      call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s,       &
                                 h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                 solver%coulomb_precision)  
!
      call mem%alloc(ao_fock, wf%n_ao, wf%n_ao)          ! Holds Fock matrix, when wf%ao_fock 
                                                         ! becomes the DIIS Fock matrix 
!
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
      iteration = 1
      converged = .false.
!
      converged_energy   = .false.
      converged_residual = .false.
!
      prev_energy = zero
!
      write(output%unit, '(/t3,a)') 'Iteration    Energy (a.u.)        Max(grad.)    Delta E (a.u.)'
      write(output%unit, '(t3,a)')  '--------------------------------------------------------------'
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)         
!
!        Set current energy
!
         energy = wf%energy
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e10.4,4x,e10.4)') iteration, wf%energy, &
                                          max_grad, abs(wf%energy-prev_energy)
         flush(output%unit)
!
!        Test for convergence:
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_residual = max_grad                .lt. solver%residual_threshold
!
         converged = converged_residual .and. converged_energy
!
         if (converged) then
!
            write(output%unit, '(t3,a)') '--------------------------------------------------------------'
            write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
!
            call solver%print_summary(wf)
!
         else
!
            prev_energy     = wf%energy
            prev_ao_density = wf%ao_density
!
            call wf%roothan_hall_update_orbitals() ! DIIS F => C
            call wf%update_ao_density()            ! C => D
!
            if (iteration .ne. 1) wf%ao_fock = ao_fock ! Restore non-DIIS F 
!
            call wf%update_fock_and_energy_cumulative(sp_eri_schwarz, sp_eri_schwarz_list,            &
                                                      n_s, prev_ao_density, h_wx, solver%coulomb_thr, &
                                                      solver%exchange_thr, solver%coulomb_precision)
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
      call wf%destruct_ao_overlap()
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call diis_manager%finalize()
!
      if (.not. converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         stop
!
      endif 
!
   end subroutine run_scf_diis_solver
!
!
   subroutine cleanup_scf_diis_solver(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
      logical :: do_mo_transformation
!
      integer(i15) :: i
!
!     When finished, we do a final Roothan-Hall step
!     to express the Fock matrix in the canonical MO basis 
!
    !  do_mo_transformation = .true.
    !  call wf%do_roothan_hall(do_mo_transformation)
!
   end subroutine cleanup_scf_diis_solver
!
!
   subroutine print_banner_scf_diis_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: Self-consistent field DIIS Hartree-Fock solver'
      write(output%unit, '(t3,a/)') ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(t3,a)')  'A DIIS-accelerated Roothan-Hall self-consistent field solver.'
      write(output%unit, '(t3,a)')  'In other words, a least-square fit toward a zero gradient vector' 
      write(output%unit, '(t3,a)')  'is performed using the previously recorded Fock matrices and the'
      write(output%unit, '(t3,a)')  'associated gradients. After each Roothan-Hall update of the density,'
      write(output%unit, '(t3,a)')  'a DIIS-fitted Fock matrix is used to get the next orbital coefficients,'
      write(output%unit, '(t3,a)')  'instead of the one produced directly from the AO density matrix.'

      flush(output%unit)
!
   end subroutine print_banner_scf_diis_solver
!
!
   subroutine print_summary_scf_diis_solver(solver, wf)
!!
!!    Print summary 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      class(hf) :: wf 
!
      write(output%unit, '(/t3,a,a,a)') 'Final ', trim(wf%name), ' energetics (a.u.):'
!
      write(output%unit, '(/t6,a26,f17.12)')  'Nuclear repulsion energy: ', wf%system%get_nuclear_repulsion()
      write(output%unit, '(t6,a26,f17.12)')   'Total electronic energy:  ', wf%energy
!
      call wf%print_orbital_energies()
!
   end subroutine print_summary_scf_diis_solver
!
!
   subroutine read_settings_scf_diis_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      call solver%read_hf_solver_settings()
      call solver%read_scf_diis_settings()
!
   end subroutine read_settings_scf_diis_solver
!
!
   subroutine read_scf_diis_settings_scf_diis_solver(solver)
!!
!!    Read SCF DIIS settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads settings specific to the class. 
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
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
            endif
!
         enddo
!
      endif 
!
   end subroutine read_scf_diis_settings_scf_diis_solver
!
!
end module scf_diis_solver_class
