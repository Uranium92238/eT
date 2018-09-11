module scf_hf_solver_class
!
!!
!!		Self-consistent field solver HF solver class module
!!		Written by Eirik F. Kjønstad, Sep 2018
!!
!!    A Roothan-Hall self-consistent field solver. In each iteration, 
!!    the Roothan-Hall equation (or equations for unrestricted HF theory)
!!    are solved to provide the next orbital coefficients. From the new 
!!    orbitals, a new density provides the next Fock matrix. The cycle 
!!    repeats until the solution is self-consistent (as measured by 
!!    the energy change).
!!
!!    Supported wavefunctions: HF, UHF    
!!
!
   use kinds
   use file_class
   use hf_class
   use disk_manager_class
   use abstract_hf_solver_class
!
   implicit none
!
   type, extends(abstract_hf_solver) :: scf_hf_solver
!
!     Nothing here yet, except variables in ancestor
!
   contains
!
      procedure :: prepare       => prepare_scf_hf_solver
      procedure :: run           => run_scf_hf_solver
      procedure :: cleanup       => cleanup_scf_hf_solver
!
      procedure :: print_banner  => print_banner_scf_hf_solver
      procedure :: print_summary => print_summary_scf_hf_solver
!
   end type scf_hf_solver
!
!
contains
!
!
   subroutine prepare_scf_hf_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_hf_solver) :: solver
!
      class(hf) :: wf
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
!     Initialize orbital coefficients, densities, and Fock matrices (plural for unrestricted methods)
!
      call wf%initialize_orbitals()
      call wf%initialize_density()
      call wf%initialize_fock()
!
   end subroutine prepare_scf_hf_solver
!
!
   subroutine run_scf_hf_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_hf_solver) :: solver
!
      class(hf) :: wf
!
      logical :: converged
      logical :: converged_energy
!
      real(dp) :: energy, prev_energy, n_electrons
!
      real(dp) :: ddot
!
      integer(i15) :: iteration
!
      real(dp), dimension(:,:), allocatable :: h_wx 
!
      integer(i15) :: n_s, i
!
      real(dp), dimension(:,:), allocatable     :: sp_eri_schwarz
      integer(i15), dimension(:,:), allocatable :: sp_eri_schwarz_list
!
!     :: Part I. Preparations
!
!     Construct ERI screening vector for efficient Fock construction 
!
      n_s = wf%system%n_s
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     Set initial AO density (or densities) guess
!
      write(output%unit, '(/t3,a,a,a)') 'Initial AO density is ', trim(solver%ao_density_guess), '.'
!
      call wf%set_initial_ao_density_guess(solver%ao_density_guess)
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
!     Update the orbitals and density to make sure the density is idempotent
!     (not the case for the standard atomic superposition density)
!
      call wf%roothan_hall_update_orbitals() ! F => C 
      call wf%update_ao_density()            ! C => D 
!
      call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!
!     :: Part II. Iterative SCF loop.
!
      iteration = 1
!
      converged        = .false.
      converged_energy = .false.
!
      prev_energy = zero
!
      write(output%unit, '(/t3,a)') 'Iteration    Energy (a.u.)        Delta E (a.u.)'
      write(output%unit, '(t3,a)')  '------------------------------------------------'
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)         
!
         energy = wf%energy
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e10.4)') iteration, energy, abs(energy-prev_energy)
         flush(output%unit)
!
!        Test for convergence:
!
         converged_energy = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged        = converged_energy
!
         if (converged) then
!
            write(output%unit, '(t3,a)') '------------------------------------------------'
            write(output%unit, '(/t3,a27,i3,a12)') 'Converged criterion met in ', iteration, ' iterations!'
!
            call solver%print_summary(wf)
!
         else
!
            call wf%roothan_hall_update_orbitals() ! F => C 
            call wf%update_ao_density()            ! C => D 
!
            prev_energy = wf%energy
            call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
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
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
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
   end subroutine run_scf_hf_solver
!
!
   subroutine cleanup_scf_hf_solver(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_hf_solver) :: solver
!
      class(hf) :: wf
!
!     Save AO density (or densities) to disk 
!
      call wf%save_ao_density()
!
!     Final deallocations of solver 
!     (note that we keep certain arrays in the wavefunction for later)
!
      call wf%destruct_ao_overlap()
      call wf%destruct_fock()
!
   end subroutine cleanup_scf_hf_solver
!
!
   subroutine print_banner_scf_hf_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_hf_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: Direct-integral Hartree-Fock self-consistent field solver'
      write(output%unit, '(t3,a/)') ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(t3,a)')  'A Roothan-Hall self-consistent field solver. In each iteration,' 
      write(output%unit, '(t3,a)')  'the Roothan-Hall equation (or equations for unrestricted HF theory)'
      write(output%unit, '(t3,a)')  'are solved to provide the next orbital coefficients. From the new'
      write(output%unit, '(t3,a)')  'orbitals, a new density provides the next Fock matrix. The cycle' 
      write(output%unit, '(t3,a)')  'repeats until the solution is self-consistent (as measured by' 
      write(output%unit, '(t3,a)')  'the energy change).'
!
      flush(output%unit)
!
   end subroutine print_banner_scf_hf_solver
!
!
   subroutine print_summary_scf_hf_solver(solver, wf)
!!
!!    Print summary 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_hf_solver) :: solver 
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
   end subroutine print_summary_scf_hf_solver
!
!
end module scf_hf_solver_class
