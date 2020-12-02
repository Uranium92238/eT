!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module mo_scf_diis_class
!
!!
!!    MO Self-consistent field DIIS HF solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Adapted from the scf_diis_hf_class
!!
!!    A DIIS-accelerated Roothan-Hall self-consistent field solver.
!!    In each Roothan-Hall update, a fitted F matrix is used instead
!!    of the one produced from the previously obtained density.
!!
!
   use abstract_hf_solver_class
   use kinds
   use diis_tool_class,          only : diis_tool
   use timings_class,            only : timings
   use sequential_file_class,    only : sequential_file
   use hf_class,                 only : hf
   use memory_manager_class,     only : mem
   use convergence_tool_class,   only : convergence_tool
!
   implicit none
!
   type, extends(abstract_hf_solver) :: mo_scf_diis
!
      integer :: diis_dimension
!
      logical :: converged
!
      logical :: restart
!
      class(diis_tool), allocatable :: diis
!
      character(len=200) :: storage
      logical :: records_in_memory
!
      logical  :: cumulative
      real(dp) :: cumulative_threshold
!
   contains
!
      procedure :: prepare                      => prepare_mo_scf_diis
      procedure :: run                          => run_mo_scf_diis
!
      procedure :: read_settings                => read_settings_mo_scf_diis
      procedure :: read_mo_scf_diis_settings    => read_mo_scf_diis_settings_mo_scf_diis
!
      procedure :: print_settings               => print_settings_mo_scf_diis
      procedure :: print_mo_scf_diis_settings   => print_mo_scf_diis_settings_mo_scf_diis
!
      procedure, private :: update_diis_history => update_diis_history_mo_scf_diis
!
   end type mo_scf_diis
!
!
   interface mo_scf_diis 
!
      procedure :: new_mo_scf_diis
!
   end interface mo_scf_diis
!
!
contains
!
!
   function new_mo_scf_diis(wf, restart) result(solver)
!!
!!    New MO SCF DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(mo_scf_diis) :: solver
!
      class(hf) :: wf
!
      logical, intent(in) :: restart
!
!     Set standard settings
!
      solver%restart                = restart
      solver%diis_dimension         = 8
      solver%max_iterations         = 100
      solver%ao_density_guess       = 'SAD'
      solver%storage                = 'memory'
      solver%cumulative             = .false.
      solver%cumulative_threshold   = 1.0d0
!
!     Initialize convergence checker with default threshols
!
      solver%convergence_checker = convergence_tool(energy_threshold   = 1.0d-7,   &
                                                    residual_threshold = 1.0d-7,   &
                                                    energy_convergence = .false.)
!
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call solver%prepare(wf)
!
   end function new_mo_scf_diis
!
!
   subroutine prepare_mo_scf_diis(solver, wf)
!!
!!    Prepare SCF-DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(mo_scf_diis) :: solver
!
      class(hf) :: wf
!
      solver%name_         = 'MO Self-consistent field DIIS Hartree-Fock solver'
      solver%tag           = 'MO SCF DIIS'
!
      solver%description   = 'A DIIS-accelerated Roothan-Hall self-consistent field solver. &
                              &A least-square DIIS fit is performed on the previous Fock matrices and &
                              &associated gradients. Following the Roothan-Hall update of the density, &
                              &the DIIS-fitted Fock matrix is used to get the next orbital coefficients.'
!
      call solver%print_banner()
!
!     Set wavefunction screenings.
!     Note that the screenings must be tighter for tighter gradient thresholds)
!     and print settings to output
!
      call wf%set_screening_and_precision_thresholds(solver%convergence_checker%residual_threshold)
!
      call solver%print_settings(wf)
!
!     Initialize the orbitals, density, and the Fock matrix (or matrices)
!
      call wf%initialize_fock()
      call wf%initialize_density()
      call wf%initialize_orbitals()
!
      if (solver%restart) then
!
         call output%printf('m', '- Requested restart. Reading orbitals from file', &
                            fs='(/t3,a)')
!
         call wf%is_restart_safe
         call wf%read_for_scf_restart_mo()
!
      else
!
         call output%printf('m', '- Setting initial AO density to (a0)', &
                            chars=[solver%ao_density_guess], fs='(/t3,a)')
!
         call wf%write_scf_restart()
         call wf%set_initial_ao_density_guess(solver%ao_density_guess)
         call wf%prepare_for_roothan_hall_mo() ! NOTE: this routine does wavefunction 
!                                                      specific preparations for RH, 
!                                                      which changes in MLHF!
!
      endif
!
!     Determine whether to store records in memory or on file
!
      if (trim(solver%storage) == 'memory') then 
!
         solver%records_in_memory = .true.
!
      elseif (trim(solver%storage) == 'disk') then 
!
         solver%records_in_memory = .false.
!
      else 
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                 trim(solver%storage))
!
      endif 
!
   end subroutine prepare_mo_scf_diis
!
!
   subroutine run_mo_scf_diis(solver, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
      use array_utilities, only : get_abs_max
!
      implicit none
!
      class(mo_scf_diis) :: solver
!
      class(hf) :: wf
!
      real(dp) :: max_grad, energy, prev_energy
!
      integer :: iteration
!
      real(dp), dimension(:,:), allocatable  :: F
      real(dp), dimension(:,:), allocatable  :: G
      real(dp), dimension(:,:), allocatable  :: prev_ao_density
!
      integer :: dim_gradient, dim_fock
!
      type(timings) :: iteration_timer, solver_timer
!
!     :: Part I. Preparations.
!
      iteration_timer = timings('MO SCF DIIS iteration time', pl='normal')
      solver_timer = timings('MO SCF DIIS solver time', pl='minimal')
!
      call solver_timer%turn_on()
!
      call wf%update_fock_and_energy_mo()
!
!     Initialize the DIIS manager object
!
      dim_fock     = (wf%n_mo)*(wf%n_mo)
      dim_gradient = (wf%n_o)*(wf%n_v)
!
      solver%diis = diis_tool('mo_hf_diis',           &
                              dim_fock,               &
                              dim_gradient,           &
                              solver%diis_dimension,  &
                              accumulate=.false.,     &
                              erase_history=.true.)
!
      call solver%diis%initialize_storers(solver%records_in_memory)
!
      call mem%alloc(G, wf%n_v, wf%n_o)
      call mem%alloc(F, wf%n_mo, wf%n_mo)
      call mem%alloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
!
      call wf%get_roothan_hall_mo_gradient(G)
!
      max_grad = get_abs_max(G, dim_gradient)
!
      call wf%get_mo_fock(F)           ! Sets F to wf%mo_fock
      call solver%diis%update(G, F)    ! First diis step -> Update F
!
!     Part II. Iterative SCF loop.
!
      solver%converged = .false.
!
      prev_energy = zero
!
      call output%printf('n', 'Iteration       Energy (a.u.)      Max(grad.)    &
                         &Delta E (a.u.)', fs='(/t3,a)')
      call output%print_separator('n', 63,'-')
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
         call output%printf('n', '(i4)  (f25.12)    (e11.4)    (e11.4)', &
                            ints=[iteration], reals=[wf%energy, max_grad, &
                            abs(wf%energy-prev_energy)], fs='(t3,a)')
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         solver%converged = solver%convergence_checker%has_converged(max_grad, energy-prev_energy, iteration)
!
         if (solver%converged) then
!
            call output%print_separator('n', 63,'-')
!
            call output%printf('n', 'Convergence criterion met in (i0) iterations!', &
                               ints=[iteration], fs='(t3,a)')
!
         else
!
            prev_energy = wf%energy
!
!           Switch to cumulative Fock construction?
!
            if (.not. solver%cumulative .and. &
                max_grad .lt. solver%cumulative_threshold) then 
!
               solver%cumulative = .true.
!
               call output%printf('v', 'Switching to Fock construction using &
                                  &density differences.', fs='(t3,a)')
!
            endif
!
            if (solver%cumulative) call wf%get_ao_density_sq(prev_ao_density)
!
!           Update Fock, coefficients, density and gradient
!
            call wf%roothan_hall_update_orbitals_mo()  ! DIIS F => C
            call wf%update_ao_density()                ! C => D
!
!           Construct updated Fock matrix from the density 
!
            if (solver%cumulative) then 
!
               call wf%update_fock_and_energy_mo(prev_ao_density) ! MO F
!
            else 
!
               call wf%update_fock_and_energy_mo() ! MO F
!
            endif
!      
            call wf%get_roothan_hall_mo_gradient(G)    ! MO G
            max_grad = get_abs_max(G, dim_gradient)
!
!           Update DIIS history
!
            call solver%update_diis_history(wf, iteration)
!
!           DIIS step
!
            call wf%get_mo_fock(F)              ! Sets F to wf%mo_fock
!
            call solver%diis%update(G, F)       ! Diis update -> F
            call wf%set_mo_fock(F)              ! Set wf%mo_fock to F
!
         endif
!
         call iteration_timer%turn_off()
         call iteration_timer%reset()
!
         iteration = iteration + 1
!
!        Save the orbitals to file
!
         call wf%save_orbital_coefficients()
         call wf%save_orbital_energies()
!
      enddo
!
      call mem%dealloc(G, wf%n_v, wf%n_o)
      call mem%dealloc(F, wf%n_mo, wf%n_mo)
      call mem%dealloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
!
!     Initialize engine (make final deallocations, and other stuff)
!
      if (.not. solver%converged) then
!
         call output%print_separator('n', 63,'-')
         call output%error_msg('Was not able to converge the equations in the given number of maximum iterations.')
!
      endif
!
      call wf%flip_final_orbitals()
!
      call solver%diis%finalize_storers()
!
      call solver_timer%turn_off()
!
   end subroutine run_mo_scf_diis
!
!
   subroutine print_mo_scf_diis_settings_mo_scf_diis(solver)
!!
!!    Print SCF-DIIS settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(mo_scf_diis) :: solver
!
      call output%printf('m', 'DIIS dimension:               (i11)', &
                         ints=[solver%diis_dimension], fs='(/t6,a)')
!
   end subroutine print_mo_scf_diis_settings_mo_scf_diis
!
!
   subroutine read_settings_mo_scf_diis(solver)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      implicit none
!
      class(mo_scf_diis) :: solver
!
      call solver%read_hf_solver_settings()
      call solver%read_mo_scf_diis_settings()
!
   end subroutine read_settings_mo_scf_diis
!
!
   subroutine read_mo_scf_diis_settings_mo_scf_diis(solver)
!!
!!    Read SCF DIIS settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Reads settings specific to the class.
!!
      implicit none
!
      class(mo_scf_diis) :: solver
!
      call input%get_keyword_in_section('diis dimension', 'solver scf', solver%diis_dimension)
      call input%get_keyword_in_section('storage', 'solver scf', solver%storage)
!
   end subroutine read_mo_scf_diis_settings_mo_scf_diis
!
!
   subroutine update_diis_history_mo_scf_diis(solver, wf, iteration)
!!
!!    Update DIIS history 
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Updates the DIIS history, such that it is 
!!    is in the current MO basis.
!!
!!    This is done by transforming the MO quantities
!!    of the DIIS history (previous Fock matrices and gradients)
!!
!!       F_new = W^T F W,
!!
!!       g_new = W^T h W,
!!
      implicit none
!
      class(mo_scf_diis) :: solver
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: iteration
!
      real(dp), dimension(:,:), allocatable  :: fock_i, X_i, Y_i
      real(dp), dimension(:,:), allocatable  :: gradient_i
!
      integer :: history_dim, i, current_diis_dim
!     
      call mem%alloc(gradient_i, wf%n_v, wf%n_o)
      call mem%alloc(fock_i, wf%n_mo, wf%n_mo)
      call mem%alloc(X_i, wf%n_mo, wf%n_mo)
      call mem%alloc(Y_i, wf%n_v, wf%n_o)
!
      current_diis_dim = solver%diis%get_dim_G()
!
      if (current_diis_dim .lt. iteration) then
!
         history_dim = current_diis_dim
!
      else
!
         history_dim = current_diis_dim - 1
!
      endif
!
      do i = 1, history_dim
!
!        Get diis history in previous mo basis
!
         call solver%diis%read_x(fock_i, i)
         call solver%diis%read_e(gradient_i, i)
!
!        Transform to current mo basis
!
         call dgemm('T', 'N',          &
                     wf%n_mo,          &
                     wf%n_mo,          &
                     wf%n_mo,          &
                     one,              &
                     wf%W_mo_update,   &
                     wf%n_mo,          &
                     fock_i,           &
                     wf%n_mo,          &
                     zero,             &
                     X_i,              &
                     wf%n_mo)
!
         call dgemm('N', 'N',          &
                     wf%n_mo,          &
                     wf%n_mo,          &
                     wf%n_mo,          &
                     one,              &
                     X_i,              &
                     wf%n_mo,          &
                     wf%W_mo_update,   &
                     wf%n_mo,          &
                     zero,             &
                     fock_i,           &
                     wf%n_mo)
!
         call dgemm('T', 'N',                                  &
                     wf%n_v,                                   &
                     wf%n_o,                                   &
                     wf%n_v,                                   &
                     one,                                      &
                     wf%W_mo_update(wf%n_o + 1, wf%n_o + 1),   &
                     wf%n_mo,                                  &
                     gradient_i,                               &
                     wf%n_v,                                   &
                     zero,                                     &
                     Y_i,                                      &
                     wf%n_v)
!
         call dgemm('N', 'N',          &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_o,           &
                     one,              &
                     Y_i,              &
                     wf%n_v,           &
                     wf%W_mo_update,   &
                     wf%n_mo,          &
                     zero,             &
                     gradient_i,       &
                     wf%n_v)
!
!        Set diis history in current mo basis
!
         call solver%diis%write_e(gradient_i, i)
         call solver%diis%write_x(fock_i, i)
!
      enddo
!
      call mem%dealloc(gradient_i, wf%n_v, wf%n_o)
      call mem%dealloc(fock_i, wf%n_mo, wf%n_mo)
      call mem%dealloc(X_i, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y_i, wf%n_v, wf%n_o)
!
   end subroutine update_diis_history_mo_scf_diis
!
!
   subroutine print_settings_mo_scf_diis(solver, wf)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad, Dec 2019
!!
      implicit none
!
      class(mo_scf_diis), intent(in) :: solver
!
      class(hf), intent(in) :: wf
!
      call output%printf('m', '- Hartree-Fock solver settings:',fs='(/t3,a)')
!
      call solver%print_mo_scf_diis_settings()
      call solver%convergence_checker%print_settings()
      call wf%print_screening_settings()
!
   end subroutine print_settings_mo_scf_diis
!
!
end module mo_scf_diis_class
