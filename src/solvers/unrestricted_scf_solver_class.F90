module unrestricted_scf_solver_class
!
!!
!!		Unrestricted self-consistent field solver class module
!!		Written by Eirik F. Kjønstad, Sep 2018
!!
!!    A Roothan-Hall self-consistent field solver. In each iteration, 
!!    the Pople–Nesbet–Berthier equations are solved,
!! 
!!       F^sigma C = S C e ----> C ------> D^sigma,   sigma = alpha, beta,
!!
!!    From the obtained C, the AO densities D^sigma are constructed, 
!!    and from them a new set of spin Fock matrices F^sigma are made.
!!    The cycle ends when the energy has converged to within a certain
!!    threshold. 
!!
!
   use kinds
   use file_class
   use uhf_class
   use disk_manager_class
   use hf_solver_class
!
   implicit none
!
   type, extends(hf_solver) :: unrestricted_scf_solver
!
!     Nothing here yet
!
   contains
!
      procedure :: run_2   => run_2_unrestricted_scf_solver
      procedure :: cleanup => cleanup_unrestricted_scf_solver
!
      procedure :: print_banner => print_banner_unrestricted_scf_solver
!
   end type unrestricted_scf_solver
!
!
contains
!
!
   subroutine run_2_unrestricted_scf_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(unrestricted_scf_solver) :: solver
!
      class(uhf) :: wf
!
      logical :: converged
      logical :: converged_energy
!
      real(dp) :: energy, prev_energy
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
!     Print solver banner
!
      call solver%print_banner()
!
!     :: Construct screening vectors, as well as a degeneracy vector, 
!     used to construct AO Fock efficiently
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     :: Get an initial AO density and initial MO coefficients
!     by solving the Roothan-Hall equations for the SAD density
!
      call wf%initialize_ao_fock()
      call wf%initialize_ao_density()
      call wf%initialize_mo_coefficients()
!
      call wf%initialize_ao_fock_a()
      call wf%initialize_ao_fock_b()
!
      call wf%initialize_ao_density_a()
      call wf%initialize_ao_density_b()
!
      call wf%initialize_orbital_coefficients_a()
      call wf%initialize_orbital_coefficients_b()
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call get_ao_h_xy(h_wx)
!
      call wf%set_ao_density_to_core_guess(h_wx)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha',   &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s,      &
                                       h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                       solver%coulomb_precision)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta',   &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s,      &
                                       h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                       solver%coulomb_precision)
!
      call wf%calculate_uhf_energy(h_wx)
!
      wf%ao_fock = wf%ao_fock_a 
      call wf%do_roothan_hall()
      wf%orbital_coefficients_a = wf%orbital_coefficients 
!
      wf%ao_fock = wf%ao_fock_b 
      call wf%do_roothan_hall()
      wf%orbital_coefficients_b = wf%orbital_coefficients 
!
      call wf%construct_ao_spin_density('alpha') ! D_alpha 
      call wf%construct_ao_spin_density('beta')  ! D_beta 
!
      call wf%form_ao_density() ! D_alpha + D_beta 
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha',   &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s,      &
                                       h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                       solver%coulomb_precision)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta',   &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s,      &
                                       h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                       solver%coulomb_precision)
!
      call wf%calculate_uhf_energy(h_wx)
!
      iteration = 1
      converged = .false.
!
      converged_energy   = .false.
!
      prev_energy = zero
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)        Delta E (a.u.)'
      write(output%unit, '(t3,a)') '------------------------------------------------'
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)         
!
         energy = wf%uhf_energy
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
            write(output%unit, '(/t3,a13,i3,a12/)') 'Converged in ', iteration, ' iterations!'
!
         else
!
!           Solve the Pople–Nesbet–Berthier coupled equations
!           to get new spin densities
!
            wf%ao_fock = wf%ao_fock_a 
            call wf%do_roothan_hall()
            wf%orbital_coefficients_a = wf%orbital_coefficients 
!
            wf%ao_fock = wf%ao_fock_b 
            call wf%do_roothan_hall()
            wf%orbital_coefficients_b = wf%orbital_coefficients 
!
            call wf%construct_ao_spin_density('alpha') ! D_alpha 
            call wf%construct_ao_spin_density('beta')  ! D_beta 
!
            call wf%form_ao_density() ! D_alpha + D_beta 
!
!           Construct the spin AO Fock matrices and compute the energy
!
            call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha',   &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s,      &
                                       h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                       solver%coulomb_precision)
!
            call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta',    &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s,      &
                                       h_wx, solver%coulomb_thr, solver%exchange_thr, &
                                       solver%coulomb_precision)
!
            prev_energy = wf%uhf_energy
            call wf%calculate_uhf_energy(h_wx)
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
      call wf%destruct_ao_overlap()
!
      call wf%destruct_ao_fock()
!
      call wf%destruct_ao_fock_a()
      call wf%destruct_ao_fock_b()
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
   end subroutine run_2_unrestricted_scf_solver
!
!
   subroutine cleanup_unrestricted_scf_solver(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(unrestricted_scf_solver) :: solver
!
      class(hf) :: wf
!
!     Nothing here yet, but: maybe write density if requested?
!
   end subroutine cleanup_unrestricted_scf_solver
!
!
   subroutine print_banner_unrestricted_scf_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(unrestricted_scf_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: Direct-integral unrestricted Hartree-Fock self-consistent field solver'
      write(output%unit, '(t3,a/)') ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(t3,a)') 'A Roothan-Hall self-consistent field solver. In each iteration,'
      write(output%unit, '(t3,a/)') 'the Pople–Nesbet–Berthier equations are solved,'
!
      write(output%unit, '(t3,a/)') '   F^sigma C = S C e ----> C ------> D^sigma,   sigma = alpha, beta,'
!
      write(output%unit, '(t3,a)') 'From the obtained C, the AO densities D^sigma are constructed,' 
      write(output%unit, '(t3,a)') 'and from them a new set of spin Fock matrices F^sigma are made.'
      write(output%unit, '(t3,a)') 'The cycle ends when the energy has converged to within a certain'
      write(output%unit, '(t3,a/)') 'threshold.'

      flush(output%unit)
!
   end subroutine print_banner_unrestricted_scf_solver
!
!
end module unrestricted_scf_solver_class
