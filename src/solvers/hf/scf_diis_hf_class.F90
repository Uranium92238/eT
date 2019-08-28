!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module scf_diis_hf_class
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
   type, extends(abstract_hf_solver) :: scf_diis_hf
!
      integer :: diis_dimension
!
      logical :: converged
!
      type(file) :: restart_file
      logical    :: restart
!
   contains
!
      procedure :: run                    => run_scf_diis_hf
      procedure :: cleanup                => cleanup_scf_diis_hf
!
      procedure :: read_settings           => read_settings_scf_diis_hf
      procedure :: read_scf_diis_settings  => read_scf_diis_settings_scf_diis_hf
!
      procedure :: print_scf_diis_settings => print_scf_diis_settings_scf_diis_hf
!
   end type scf_diis_hf
!
!
   interface scf_diis_hf 
!
      procedure :: new_scf_diis_hf
!
   end interface scf_diis_hf
!
!
contains
!
!
   function new_scf_diis_hf(wf, restart) result(solver)
!!
!!    New SCF DIIS 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(scf_diis_hf) :: solver
!
      class(hf) :: wf
!
      logical, intent(in) :: restart 
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
      solver%restart             = restart
      solver%diis_dimension      = 8
      solver%max_iterations      = 100
      solver%ao_density_guess    = 'SAD'
      solver%energy_threshold    = 1.0D-6
      solver%gradient_threshold  = 1.0D-6
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
      if (solver%restart) then
!
         write(output%unit, '(/t3,a)') '- Requested restart. Reading orbitals from file:'
!
         call wf%read_orbital_coefficients()
         call wf%update_ao_density()
         call wf%read_orbital_energies()
!
      else
!
         write(output%unit, '(/t3,a,a,a)') '- Setting initial AO density to ', trim(solver%ao_density_guess), ':'
         call wf%set_initial_ao_density_guess(solver%ao_density_guess)
!
      endif
!
   end function new_scf_diis_hf
!
!
   subroutine print_scf_diis_settings_scf_diis_hf(solver)
!!
!!    Print SCF-DIIS settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(scf_diis_hf) :: solver
!
      write(output%unit, '(t6,a29,i3)') 'DIIS dimension:               ', solver%diis_dimension
!
   end subroutine print_scf_diis_settings_scf_diis_hf
!
!
   subroutine run_scf_diis_hf(solver, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_hf) :: solver
!
      class(hf) :: wf
!
      type(diis_tool) :: diis
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
      integer :: dim_gradient, dim_fock
!
      type(timings) :: iteration_timer, solver_timer
!
!     :: Part I. Preparations.
!
      iteration_timer = new_timer('SCF DIIS iteration time')
      solver_timer = new_timer('SCF DIIS solver time')
!
      call solver_timer%turn_on()
!
!     Construct screening vectors for efficient Fock construction
!
      n_s = wf%system%get_n_shells()
!
      if (.not. allocated(wf%sp_eri_schwarz)) call mem%alloc(wf%sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      if (.not. allocated(wf%sp_eri_schwarz_list)) call mem%alloc(wf%sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz()
!
!     Initialize the DIIS manager object
!
      dim_fock     = ((wf%n_ao)*(wf%n_ao + 1)/2)*(wf%n_densities)
      dim_gradient = (wf%n_ao*(wf%n_ao - 1)/2)*(wf%n_densities)
!
      diis = diis_tool('hf_diis', dim_fock, dim_gradient, solver%diis_dimension)
!
!     Set the initial density guess and Fock matrix
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(h_wx)
!
      call wf%update_fock_and_energy(h_wx)
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
      if (.not. solver%restart) then
!
         call wf%roothan_hall_update_orbitals() ! F => C
         call wf%update_ao_density()            ! C => D
!
         call wf%update_fock_and_energy(h_wx)
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
      call diis%update(G, F)
!
!     Part II. Iterative SCF loop.
!
      solver%converged   = .false.
      converged_energy   = .false.
      converged_gradient = .false.
!
      prev_energy = zero
!
      write(output%unit, '(/t3,a)') 'Iteration       Energy (a.u.)      Max(grad.)    Delta E (a.u.)'
      write(output%unit, '(t3,a)')  '---------------------------------------------------------------'
!
      iteration = 1
!
      do while (.not. solver%converged .and. iteration .le. solver%max_iterations)
!
         call iteration_timer%turn_on()
!
!        Set energy and print information for current iteration
!
         energy = wf%energy
!
         write(output%unit, '(t3,i4,9x,f17.12,4x,e11.4,4x,e11.4)') iteration, wf%energy, &
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
            write(output%unit, '(t3,a)')          '---------------------------------------------------------------'
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
            prev_energy = wf%energy
            call wf%get_ao_density_sq(prev_ao_density)
!
            call wf%roothan_hall_update_orbitals()     ! DIIS F => C
!
            call wf%save_orbital_coefficients()
            call wf%save_orbital_energies()
!
            call wf%update_ao_density()                ! C => D
!
            if (iteration .ne. 1) call wf%set_ao_fock(ao_fock) ! Restore F
!
!           calling a wrapper for cumulative or no_cumulative depending on options
!
            call wf%update_fock_and_energy(h_wx,prev_ao_density)
!
            call wf%get_packed_roothan_hall_gradient(G)
!
            max_grad = get_abs_max(G, dim_gradient)
!
            call wf%get_ao_fock(F)
            call dcopy(dim_fock, F, 1, ao_fock, 1)
!
            call diis%update(G, F)
            call wf%set_ao_fock(F)
!
         endif
!
         call iteration_timer%turn_off()
         call iteration_timer%reset()
!
         iteration = iteration + 1
!
      enddo
!
   !   call mem%dealloc(wf%sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
   !   call mem%dealloc(wf%sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call mem%dealloc(G, wf%n_ao*(wf%n_ao - 1)/2, wf%n_densities)
      call mem%dealloc(F, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(ao_fock, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
      call mem%dealloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
!
      call diis%cleanup()
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
      call solver_timer%turn_off()
!
   end subroutine run_scf_diis_hf
!
!
   subroutine cleanup_scf_diis_hf(solver, wf)
!!
!! 	Cleanup
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_hf) :: solver
!
      class(hf) :: wf
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(solver%tag)
!
!     MO transform the AO Fock matrix
!
      call wf%initialize_mo_fock()
      call wf%construct_mo_fock()
!
!     Save the orbitals to file & store restart information
!
      call wf%save_orbital_coefficients()
      call wf%save_orbital_energies()
!
!     Save AO density (or densities) to disk 
!
      call wf%save_ao_density()
!
   end subroutine cleanup_scf_diis_hf
!
!
   subroutine read_settings_scf_diis_hf(solver)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      implicit none
!
      class(scf_diis_hf) :: solver
!
      call solver%read_hf_solver_settings()
      call solver%read_scf_diis_settings()
!
   end subroutine read_settings_scf_diis_hf
!
!
   subroutine read_scf_diis_settings_scf_diis_hf(solver)
!!
!!    Read SCF DIIS settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Reads settings specific to the class.
!!
      implicit none
!
      class(scf_diis_hf) :: solver
!
      call input%get_keyword_in_section('diis dimension', 'solver hf', solver%diis_dimension)
!
   end subroutine read_scf_diis_settings_scf_diis_hf
!
!
end module scf_diis_hf_class
