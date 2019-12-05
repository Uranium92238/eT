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
   use abstract_hf_solver_class
!
   use memory_manager_class, only : mem
   use hf_class, only : hf
   use array_utilities, only: get_abs_max
!
   implicit none
!
   type, extends(abstract_hf_solver) :: scf_hf
!
      character(len=400) :: warning
!
      logical :: restart
!
   contains
!
      procedure :: run           => run_scf_hf
!
      procedure :: print_banner  => print_banner_scf_hf
      procedure :: prepare       => prepare_scf_hf
!
   end type scf_hf
!
!
   interface scf_hf
!
      procedure :: new_scf_hf
      procedure :: new_scf_hf_from_parameters
!
   end interface scf_hf 
!
!
contains
!
!
   function new_scf_hf(wf, restart) result(solver)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(scf_hf) :: solver
!
      class(hf) :: wf
!
      logical, intent(in) :: restart
! 
!     Set standards 
!
      solver%max_iterations      = 100
      solver%ao_density_guess    = 'SAD'
      solver%energy_threshold    = 1.0D-6
      solver%gradient_threshold  = 1.0D-6
      solver%restart             = restart
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call solver%prepare(wf)
!
   end function new_scf_hf
!
!
   function new_scf_hf_from_parameters(wf, restart,            &
                                        max_iterations,        &
                                        ao_density_guess,      &
                                        energy_threshold,      &
                                        gradient_threshold) result(solver)
!!
!!    New SCF from parameters
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      type(scf_hf) :: solver
!
      class(hf) :: wf
!
      logical,            intent(in) :: restart
      integer,            intent(in) :: max_iterations
      character(len=200), intent(in) :: ao_density_guess
      real(dp),           intent(in) :: energy_threshold
      real(dp),           intent(in) :: gradient_threshold
!
!     Set settings from parameters
!
      solver%max_iterations      = max_iterations
      solver%ao_density_guess    = ao_density_guess
      solver%energy_threshold    = energy_threshold
      solver%gradient_threshold  = gradient_threshold
      solver%restart             = restart
!
      call solver%prepare(wf)
!
   end function new_scf_hf_from_parameters
!
!
   subroutine prepare_scf_hf(solver, wf)
!!
!!    Prepare SCF
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
      call wf%set_screening_and_precision_thresholds(solver%gradient_threshold)
!
      call output%printf('- Hartree-Fock solver settings:',fs='(/t3,a)', pl='minimal')
!
      call solver%print_hf_solver_settings()
      call wf%print_screening_settings()
!
!     Initialize orbital coefficients, densities, and Fock matrices (plural for unrestricted methods)
!
      call wf%initialize_orbitals()
      call wf%initialize_density()
      call wf%initialize_fock()
!
!     The initial (idempotent) density is obtained either from one Roothan-Hall step, or
!     by reading orbital coefficients (restart)
!
      if (solver%restart) then
!
         call output%printf('- Requested restart. Reading orbitals from file',fs='(/t3,a)', pl='minimal')
!
         call wf%is_restart_safe('ground state')
         call wf%read_for_scf_restart()
!
      else 
!
         call output%printf('- Setting initial AO density to '//trim(solver%ao_density_guess),&
            fs='(/t3,a)',pl='minimal')
!
         call wf%write_scf_restart()
         call wf%set_initial_ao_density_guess(solver%ao_density_guess)
         call wf%prepare_for_roothan_hall()
!
      endif
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
      logical :: converged_energy, converged_gradient 
!
      real(dp) :: energy, prev_energy, max_grad 
!
      integer :: iteration
!
      real(dp), dimension(:,:), allocatable :: h_wx
!
      if (wf%n_ao == 1) then 
!
         call solver%run_single_ao(wf)
         call solver%print_summary(wf)
         return
!
      endif 
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(h_wx)
!
!     Construct first Fock from initial idempotent density 
!
      call wf%update_fock_and_energy(h_wx)
!
!     :: Part I. Iterative SCF loop.
!
      iteration = 1
!
      converged            = .false.
      converged_energy     = .false.
      converged_gradient   = .false.
!
      prev_energy = zero
!
      call output%printf('Iteration       Energy (a.u.)      Max(grad.)    Delta E (a.u.)',fs='(/t3,a)',pl='normal') 
      call output%print_separator('n', 63,'-')
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)
!
         energy = wf%energy
!
         max_grad = wf%get_max_roothan_hall_gradient()
!
!        Print current iteration information
!
         call output%printf('(i4)  (f25.12)    (e11.4)    (e11.4)', &
               ints=[iteration], &
               reals=[wf%energy, max_grad, abs(wf%energy-prev_energy)], &
               fs='(t3,a)', pl='normal')
!
!        Test for convergence:
!
         converged_energy     = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_gradient   = max_grad                .lt. solver%gradient_threshold
!
         converged            = converged_energy .and. converged_gradient
!
         if (converged_gradient .and. iteration .eq. 1) converged = .true.
!
         if (converged) then
!
            call output%print_separator('n', 63,'-')
            call output%printf('Convergence criterion met in (i0) iterations!', ints=[iteration], fs='(/t3,a)', pl='normal') 
!
            call solver%print_summary(wf)
!
         else
!
            call wf%roothan_hall_update_orbitals() ! F => C
            call wf%update_ao_density()            ! C => D
!
            prev_energy = wf%energy
            call wf%update_fock_and_energy(h_wx)
!
         endif
!
         call wf%save_orbital_coefficients()
         call wf%save_orbital_energies()
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      if (.not. converged) then
!
          call output%print_separator('n', 63,'-')
!
         call output%error_msg('Was not able to converge the equations in the given number of maximum iterations.')
!
      endif
!
   end subroutine run_scf_hf
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
      call output%printf(':: (a0)', pl='normal', fs='(//t3,a)',  chars=[trim(solver%tag)])
      call output%printf(':: (a0)', pl='normal', fs='(t3,a)',    chars=[trim(solver%author)])
      call output%printf('(a0)',    pl='normal', ffs='(/t3,a)',  chars=[trim(solver%warning)])
      call output%printf('(a0)',    pl='normal', ffs='(/t3,a)',  chars=[trim(solver%description)])
!
   end subroutine print_banner_scf_hf
!
!
end module scf_hf_class
