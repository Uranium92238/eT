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
   use abstract_hf_solver_class
!
   use diis_tool_class,          only : diis_tool
   use timings_class,            only : timings
   use hf_class,                 only : hf
   use memory_manager_class,     only : mem
   use array_utilities,          only : get_abs_max
   use convergence_tool_class,   only : convergence_tool
!
   implicit none
!
   type, extends(abstract_hf_solver) :: scf_diis_hf
!
      integer :: diis_dimension
!
      logical  :: cumulative
      real(dp) :: cumulative_threshold
!
      logical :: converged
      logical :: restart
!
      logical :: crop ! Standard DIIS if false; CROP variant of DIIS if true
!
      logical :: records_in_memory
      character(len=200) :: storage 
!
   contains
!
      procedure, private :: prepare       => prepare_scf_diis_hf
      procedure :: run                    => run_scf_diis_hf
!
      procedure :: read_settings           => read_settings_scf_diis_hf
      procedure :: read_scf_diis_settings  => read_scf_diis_settings_scf_diis_hf
!
      procedure :: print_settings            => print_settings_scf_diis_hf
      procedure :: print_scf_diis_settings   => print_scf_diis_settings_scf_diis_hf
!
   end type scf_diis_hf
!
!
   interface scf_diis_hf 
!
      procedure :: new_scf_diis_hf
      procedure :: new_scf_diis_hf_from_parameters
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
!     Set standard settings
!
      solver%restart                = restart
      solver%diis_dimension         = 8
      solver%max_iterations         = 100
      solver%ao_density_guess       = 'SAD'
      solver%storage                = 'memory'
      solver%cumulative             = .false.
      solver%cumulative_threshold   = 1.0d0
      solver%crop                   = .false.
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
   end function new_scf_diis_hf
!
!
   function new_scf_diis_hf_from_parameters(wf, restart,                &
                                                diis_dimension,         &
                                                max_iterations,         &
                                                ao_density_guess,       &
                                                energy_threshold,       &
                                                residual_threshold,     &
                                                storage,                &
                                                cumulative_threshold,   &
                                                crop,                   &
                                                energy_convergence) result(solver)
!!
!!    New SCF DIIS from parameters
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      type(scf_diis_hf) :: solver
!
      class(hf) :: wf
!
      logical,            intent(in) :: restart 
      integer,            intent(in) :: diis_dimension
      integer,            intent(in) :: max_iterations
      character(len=200), intent(in) :: ao_density_guess
      real(dp),           intent(in) :: energy_threshold
      real(dp),           intent(in) :: residual_threshold
      character(len=200), intent(in) :: storage 
      real(dp),           intent(in) :: cumulative_threshold
      logical,            intent(in) :: crop 
      logical,            intent(in) :: energy_convergence
!
!     Set settings from parameters
!
      solver%restart                = restart
      solver%diis_dimension         = diis_dimension
      solver%max_iterations         = max_iterations
      solver%ao_density_guess       = ao_density_guess
      solver%storage                = storage
      solver%cumulative_threshold   = cumulative_threshold
      solver%crop                   = crop
!
      solver%convergence_checker = convergence_tool(energy_threshold, residual_threshold, energy_convergence)
!
      call solver%prepare(wf)
!
   end function new_scf_diis_hf_from_parameters
!
!
   subroutine prepare_scf_diis_hf(solver, wf)
!!
!!    Prepare SCF-DIIS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_hf) :: solver
!
      class(hf) :: wf
!
      solver%timer = timings('SCF DIIS solver time', pl='minimal')
      call solver%timer%turn_on()
!
      solver%name_       = 'Self-consistent field DIIS Hartree-Fock solver'
      solver%tag         = 'SCF DIIS'
      solver%description = 'A DIIS-accelerated Roothan-Hall self-consistent field solver. &
                           &A least-square DIIS fit is performed on the previous Fock matrices and &
                           &associated gradients. Following the Roothan-Hall update of the density, &
                           &the DIIS-fitted Fock matrix is used to get the next orbital coefficients.'
!
      call solver%print_banner()
!
!     Set wavefunction screenings
!     (note that the screenings must be tighter for tighter gradient thresholds)
!     & print settings to output
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
         call wf%read_for_scf_restart()
!
      else
!
         call output%printf('m', '- Setting initial AO density to (a0)', &
                            chars=[solver%ao_density_guess], fs='(/t3,a)')
!
         call wf%write_scf_restart()
         call wf%set_initial_ao_density_guess(solver%ao_density_guess)
         call wf%prepare_for_roothan_hall()
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
   end subroutine prepare_scf_diis_hf
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
      call output%printf('m', 'DIIS dimension:               (i11)', &
                         ints=[solver%diis_dimension], fs='(/t6,a)')
!
      call output%printf('m', 'Cumulative Fock threshold:    (e11.2)', &
                         reals=[solver%cumulative_threshold], fs='(t6,a)')
!
      if (solver%crop) then 
!
         call output%printf('m', 'Enabled CROP in the DIIS algorithm.', fs='(/t6,a)')
!
      endif
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
      real(dp) :: max_grad, energy, prev_energy
!
      integer :: iteration
!
      real(dp), dimension(:,:), allocatable :: F
      real(dp), dimension(:,:), allocatable :: G
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: prev_ao_density
!
      integer :: dim_gradient, dim_fock
!
      type(timings), allocatable :: iteration_timer
!
      if (wf%n_ao == 1) then 
!
         call solver%run_single_ao(wf)
         return
!
      endif 
!
!     :: Part I. Preparations.
!
      iteration_timer = timings('SCF DIIS iteration time', pl='normal')
!
!     Initialize the DIIS manager object
!
      dim_fock     = ((wf%n_ao)*(wf%n_ao + 1)/2)*(wf%n_densities)
      dim_gradient = (wf%n_mo*(wf%n_mo - 1)/2)*(wf%n_densities)
!
      diis = diis_tool('hf_diis',                           &
                        dim_fock,                           &
                        dim_gradient,                       &
                        dimension_=solver%diis_dimension,   &
                        crop=solver%crop)
!
      call diis%initialize_storers(solver%records_in_memory)
!
!     Set the initial density guess and Fock matrix
!
      call wf%update_fock_and_energy()
!
      call mem%alloc(ao_fock, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
      call mem%alloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
!
      call mem%alloc(G, wf%n_mo*(wf%n_mo - 1)/2, wf%n_densities)
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
      solver%converged = .false.
      prev_energy = zero
!
      call output%printf('n', 'Iteration       Energy (a.u.)      Max(grad.)    &
                         &Delta E (a.u.)', fs='(/t3,a)')
      call output%print_separator('n', 63, '-')
!
      iteration = 0
!
      do while (.not. solver%converged .and. iteration .le. solver%max_iterations)
!
         iteration = iteration + 1
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
         solver%converged = solver%convergence_checker%has_converged(max_grad, wf%energy-prev_energy, iteration)
!
         if (solver%converged) then
!
            call output%print_separator('n', 63, '-', fs='(t3,a)')
!
            call output%printf('n', 'Convergence criterion met in (i0) iterations!', &
                               ints=[iteration], fs='(t3,a)')
!
         else
!
!           Switch to cumulative Fock construction?
!
            if (.not. solver%cumulative .and. &
                  max_grad .lt. solver%cumulative_threshold) then 
!
               solver%cumulative = .true.
               call output%printf('v', 'Switching to Fock construction using &
                                  &density differences.', fs='(t3,a)')
!
            endif
!
            if (solver%cumulative) call wf%get_ao_density_sq(prev_ao_density)
!            
!           Update the orbitals and density by solving the Roothan-Hall problem
!
            prev_energy = wf%energy
!
            call wf%roothan_hall_update_orbitals()     ! DIIS F => C
            call wf%update_ao_density()                ! C => D
!
!           Restore Fock to non-DIIS-extrapolated Fock matrix 
!
            if (iteration .gt. 1) call wf%set_ao_fock(ao_fock) 
!
!           Construct updated Fock matrix from the density 
!
            if (solver%cumulative) then 
!
               call wf%update_fock_and_energy(prev_ao_density)
!
            else 
!
               call wf%update_fock_and_energy()
!
            endif
!
!           Construct current gradient and the max norm of it 
!
            call wf%get_packed_roothan_hall_gradient(G)
            max_grad = get_abs_max(G, dim_gradient)
!
!           Keep a copy of non-extrapolated Fock matrix 
!
            call wf%get_ao_fock(F)
            call dcopy(dim_fock, F, 1, ao_fock, 1)
!
!           Perform the DIIS extrapolation to get the next Fock matrix for Roothan-Hall
!
            call diis%update(G, F)
            call wf%set_ao_fock(F)
!
         endif
!
         call wf%save_orbital_coefficients()
         call wf%save_orbital_energies()
!
         call iteration_timer%turn_off()
         call iteration_timer%reset()
!
      enddo
!
      call mem%dealloc(G, wf%n_mo*(wf%n_mo - 1)/2, wf%n_densities)
      call mem%dealloc(F, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
!
      call mem%dealloc(ao_fock, wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities)
      call mem%dealloc(prev_ao_density, wf%n_ao**2, wf%n_densities)
!
      if (.not. solver%converged) then
!
         call output%print_separator('n', 63, '-',fs='(t3,a)')
         call output%error_msg('Was not able to converge the equations       &
                               &in the given number of maximum iterations.')
!
      endif
!
      call wf%flip_final_orbitals()
!
      call diis%finalize_storers()
      call solver%timer%turn_off()
!
   end subroutine run_scf_diis_hf
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
      call input%get_keyword_in_section('diis dimension', 'solver scf', solver%diis_dimension)
      call input%get_keyword_in_section('storage', 'solver scf', solver%storage)
!
      call input%get_keyword_in_section('cumulative fock threshold', &
                                        'solver scf', solver%cumulative_threshold)
!
      if (input%requested_keyword_in_section('crop', 'solver scf')) then 
!
         solver%crop = .true.
!
      endif
!
   end subroutine read_scf_diis_settings_scf_diis_hf
!
!
   subroutine print_settings_scf_diis_hf(solver, wf)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad, Dec 2019
!!
      implicit none
!
      class(scf_diis_hf), intent(in) :: solver
!
      class(hf), intent(in) :: wf
!
      call output%printf('m', '- Hartree-Fock solver settings:',fs='(/t3,a)')
!
      call solver%print_scf_diis_settings()
      call solver%convergence_checker%print_settings()
      call wf%print_screening_settings()
!
   end subroutine print_settings_scf_diis_hf
!
!
end module scf_diis_hf_class
