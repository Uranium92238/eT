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
module scf_hf_class
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
   type, extends(abstract_hf_solver) :: scf_hf
!
      character(len=400) :: warning
!
   contains
!
      procedure :: prepare       => prepare_scf_hf
      procedure :: run           => run_scf_hf
      procedure :: cleanup       => cleanup_scf_hf
!
      procedure :: print_banner  => print_banner_scf_hf
!
   end type scf_hf
!
!
contains
!
!
   subroutine prepare_scf_hf(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_hf) :: solver
!
      class(hf) :: wf
!
!     Print solver banner
!
      solver%tag = 'Self-consistent field solver'
      solver%author = 'E. F. Kjønstad and S, D. Folkestad, 2018'
!
      solver%description = 'A Roothan-Hall self-consistent field solver. In each iteration, &
                                  &the Roothan-Hall equation (or equations for unrestricted HF theory) &
                                  &are solved to provide the next orbital coefficients. From the new &
                                  &orbitals, a new density provides the next Fock matrix. The cycle &
                                  &repeats until the solution is self-consistent (as measured by &
                                  &the energy change).' 
!
      solver%warning = 'Warning: We recommend to use the SCF-DIIS algorithm instead, which &
                              &supports a gradient threshold and typically converges much faster. &
                              &Use only when absolutely necessary!'
!
      call solver%print_banner()
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call wf%set_screening_and_precision_thresholds(solver%gradient_threshold)
      call wf%print_screening_settings()
!
!     Initialize orbital coefficients, densities, and Fock matrices (plural for unrestricted methods)
!
      call wf%initialize_orbitals()
      call wf%initialize_density()
      call wf%initialize_fock()
!
!     Set initial AO density guess
!
      write(output%unit, '(/t3,a,a,a)') '- Setting initial AO density to ', trim(solver%ao_density_guess), ':'
      call wf%set_initial_ao_density_guess(solver%ao_density_guess)
!
   end subroutine prepare_scf_hf
!
!
   subroutine run_scf_hf(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_hf) :: solver
!
      class(hf) :: wf
!
      logical :: converged
      logical :: converged_energy
!
      real(dp) :: energy, prev_energy, n_electrons
!
      integer :: iteration
!
      real(dp), dimension(:,:), allocatable :: h_wx 
!
      integer :: n_s
!
      real(dp), dimension(:,:), allocatable     :: sp_eri_schwarz
      integer, dimension(:,:), allocatable :: sp_eri_schwarz_list
!
!     :: Part I. Preparations
!
!     Construct ERI screening vector for efficient Fock construction 
!
      n_s = wf%system%n_s
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
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
         write(output%unit, '(t3,i3,10x,f17.12,4x,e11.4)') iteration, energy, abs(energy-prev_energy)
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
            call wf%print_wavefunction_summary()
            flush(output%unit)
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
      call mem%dealloc(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
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
   end subroutine run_scf_hf
!
!
   subroutine cleanup_scf_hf(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_hf) :: solver
!
      class(hf) :: wf
!
      logical :: do_mo_transformation
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(solver%tag)
      flush(output%unit)
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
!     Do a final Roothan-Hall step to transform the Fock matrix in the canonical MO basis 
!
      do_mo_transformation = .true.
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies, do_mo_transformation)
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
   end subroutine cleanup_scf_hf
!
!
   subroutine print_banner_scf_hf(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_hf) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%warning,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
      call long_string_print(solver%description)
!
   end subroutine print_banner_scf_hf
!
!
end module scf_hf_class
