!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module hf_class
!
!!
!!    Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Wavefunction class for restricted (closed-shell) Hartree-Fock theory. 
!!
!
   use wavefunction_class
!
   use reordering
!
   use timings_class,         only : timings
   use array_utilities,       only : get_abs_max, sandwich
   use array_utilities,       only : full_cholesky_decomposition_system
   use array_utilities,       only : get_n_highest
   use array_utilities,       only : identity_array
   use output_file_class,     only : output_file
   use sequential_file_class, only : sequential_file
   use array_utilities,       only : zero_array
   use interval_class,        only : interval
   use omp_lib
!
   use stream_file_class, only: stream_file
!
   implicit none
!
!  Hartree-Fock hf
!
   type, extends(wavefunction) :: hf
!
      real(dp), dimension(:,:), allocatable :: ao_density
      real(dp), dimension(:,:,:), allocatable :: previous_ao_density
!
      integer :: n_densities ! For RHF, there is one AO density. 
                             ! For UHF, there are two AO densities (alpha, beta).
!
      real(dp) :: coulomb_threshold  = 1.0D-12   ! Screening threshold (Fock, Coulomb)
      real(dp) :: exchange_threshold = 1.0D-10   ! Screening threshold (Fock, exchange)
      real(dp) :: integral_cutoff    = 1.0D-12   ! Default: sqrt(epsilon) 
!
      real(dp) :: cumulative_fock_threshold = 1.0d0
!
      real(dp), dimension(:,:), allocatable :: W_mo_update ! Eigenvectors for 
                                                           ! Roothan-Hall in MO basis
!
!     Frozen orbital variables. Frozen orbitals are typically frozen core or frozen HF orbitals.
!
      integer :: n_frozen_hf_o
      integer :: n_frozen_hf_orbitals
      integer :: n_frozen_core_orbitals
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_frozen_hf
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_fc
!
      type(stream_file)     :: orbital_file
!
      logical :: frozen_core
      logical :: frozen_hf_mos
!
      logical :: cumulative_fock
!
      logical :: plot_active_density
!
      integer :: gradient_dimension
      type(output_file) :: mo_information_file
!
   contains
!
      procedure :: print_banner                                => print_banner_hf
      procedure :: print_summary                               => print_summary_hf
      procedure :: print_energy                                => print_energy_hf
!
!     Read, save of orbital energies and coefficients 
!
      procedure :: read_orbitals                               => read_orbitals_hf
      procedure :: save_orbitals                               => save_orbitals_hf
!
!     Read, write restart information, and check for safe restart 
!
      procedure :: read_for_scf_restart                        => read_for_scf_restart_hf
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                                     => cleanup_hf
      procedure :: prepare_for_cc                              => prepare_for_cc_hf
!
      procedure :: read_settings                               => read_settings_hf
      procedure :: read_hf_settings                            => read_hf_settings_hf
!
!     AO Fock and energy related routines
!
      procedure :: read_frozen_orbitals_settings &
                => read_frozen_orbitals_settings_hf
!
      procedure :: construct_ao_G                              => construct_ao_G_hf
      procedure :: construct_ao_G_thread_terms                 => construct_ao_G_thread_terms_hf
      procedure :: construct_ao_G_thread_terms_mo_screened     => construct_ao_G_thread_terms_mo_screened_hf
      procedure :: construct_ao_G_1der                         => construct_ao_G_1der_hf
!
      procedure :: construct_coulomb_ao_G                      => construct_coulomb_ao_G_hf
      procedure :: construct_exchange_ao_G                     => construct_exchange_ao_G_hf
!
      procedure :: set_ao_fock                                 => set_ao_fock_hf
      procedure :: get_ao_fock                                 => get_ao_fock_hf
      procedure :: calculate_hf_energy_from_fock               => calculate_hf_energy_from_fock_hf
      procedure :: calculate_hf_energy_from_G                  => calculate_hf_energy_from_G_hf
      procedure :: initialize_fock                             => initialize_fock_hf
      procedure :: destruct_fock                               => destruct_fock_hf
      procedure :: update_G_non_cumulative                     => update_G_non_cumulative_hf
      procedure :: update_G_cumulative                         => update_G_cumulative_hf
      procedure :: update_fock_and_energy                      => update_fock_and_energy_hf
      procedure :: initialize_frozen_CCT                       => initialize_frozen_CCT_hf
      procedure :: destruct_frozen_CCT                         => destruct_frozen_CCT_hf
!
!     AO Density related routines
!
      procedure :: construct_ao_density                        => construct_ao_density_hf
      procedure :: decompose_ao_density                        => decompose_ao_density_hf
      procedure :: get_ao_density                              => get_ao_density_hf
      procedure :: set_ao_density                              => set_ao_density_hf
      procedure :: initialize_density                          => initialize_density_hf
      procedure :: update_ao_density                           => update_ao_density_hf
      procedure :: save_ao_density                             => save_ao_density_hf
      procedure :: get_ao_density_sq                           => get_ao_density_sq_hf
      procedure :: set_initial_ao_density_guess                => set_initial_ao_density_guess_hf
      procedure :: set_ao_density_to_core_guess                => set_ao_density_to_core_guess_hf
      procedure :: get_n_electrons_in_density                  => get_n_electrons_in_density_hf
      procedure :: construct_shp_density_schwarz               => construct_shp_density_schwarz_hf
!
!     MO orbital related routines
!
      procedure :: initialize_orbitals                         => initialize_orbitals_hf
      procedure :: print_orbitals_and_energies                 => print_orbitals_and_energies_hf
      procedure :: save_orbital_info                           => save_orbital_info_hf
      procedure :: print_orbitals_from_coefficients            => print_orbitals_from_coefficients_hf
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density                       => initialize_ao_density_hf
      procedure :: destruct_ao_density                         => destruct_ao_density_hf
!
!     Gradient and Hessian related routines
!
      procedure :: construct_projection_matrices               => construct_projection_matrices_hf
!
      procedure :: construct_roothan_hall_gradient             => construct_roothan_hall_gradient_hf
      procedure :: get_packed_roothan_hall_gradient            => get_packed_roothan_hall_gradient_hf
!
      procedure :: construct_molecular_gradient                => construct_molecular_gradient_hf
!
!     Integral related routines
!
      procedure :: set_n_mo                                    => set_n_mo_hf
!
      procedure :: set_screening_and_precision_thresholds      => set_screening_and_precision_thresholds_hf
      procedure :: print_screening_settings                    => print_screening_settings_hf
!
      procedure :: prepare_for_roothan_hall                    => prepare_for_roothan_hall_hf
      procedure :: prepare                                     => prepare_hf
!
!     Properties
!
      procedure :: calculate_expectation_value                 => calculate_expectation_value_hf
!
!     Frozen core
!
      procedure :: prepare_mos                                 => prepare_mos_hf
      procedure :: prepare_frozen_fock_terms                   => prepare_frozen_fock_terms_hf
!
      procedure :: remove_core_orbitals                        => remove_core_orbitals_hf
      procedure :: remove_frozen_hf_orbitals                   => remove_frozen_hf_orbitals_hf
!
      procedure :: get_full_idempotent_density                 => get_full_idempotent_density_hf
!
      procedure :: construct_mo_fock_fc_term                   => construct_mo_fock_fc_term_hf
      procedure :: construct_mo_fock_frozen_hf_term            => construct_mo_fock_frozen_hf_term_hf
!
      procedure :: initialize_orbital_coefficients_fc          => initialize_orbital_coefficients_fc_hf
      procedure :: destruct_orbital_coefficients_fc            => destruct_orbital_coefficients_fc_hf
!
      procedure :: initialize_orbital_coefficients_frozen_hf   => initialize_orbital_coefficients_frozen_hf_hf
      procedure :: destruct_orbital_coefficients_frozen_hf     => destruct_orbital_coefficients_frozen_hf_hf
!
      procedure :: diagonalize_fock_frozen_hf_orbitals         => diagonalize_fock_frozen_hf_orbitals_hf
      procedure :: get_n_active_hf_atoms                       => get_n_active_hf_atoms_hf
!
      procedure :: flip_final_orbitals                         => flip_final_orbitals_hf
!
      procedure :: prepare_for_scf => prepare_for_scf_hf
      procedure :: get_energy => get_energy_hf
!
      procedure :: get_F => get_F_hf
      procedure :: get_gradient => get_gradient_hf
      procedure :: set_C_and_e => set_C_and_e_hf
      procedure :: can_do_cumulative_fock => can_do_cumulative_fock_hf
      procedure :: set_orbital_coefficients_from_reduced_ao_C &
                => set_orbital_coefficients_from_reduced_ao_C_hf
!
      procedure :: set_orbital_coefficients_from_orthonormal_ao_C &
                => set_orbital_coefficients_from_orthonormal_ao_C_hf
!
      procedure :: construct_initial_idempotent_density &
                => construct_initial_idempotent_density_hf
!
      procedure, non_overridable :: ao_to_orthonormal_ao_transformation => ao_to_orthonormal_ao_transformation_hf
!
      procedure :: ao_to_reduced_ao_transformation &
                => ao_to_reduced_ao_transformation_hf
!
      procedure :: is_restart_possible => is_restart_possible_hf
!
      procedure :: calculate_frozen_dipole_moment              => calculate_frozen_dipole_moment_hf
      procedure :: calculate_frozen_quadrupole_moment          => calculate_frozen_quadrupole_moment_hf
!
      procedure :: control_gradient_convergence                => control_gradient_convergence_hf
!
   end type hf
!
   interface
!
      include "frozen_orbital_hf_interface.F90"
      include "ao_fock_hf_interface.F90"
      include "set_get_hf_interface.F90"
      include "file_handling_hf_interface.F90"
      include "initialize_destruct_hf_interface.F90"
!
   end interface
!
!
   interface hf
!
      procedure :: new_hf
!
   end interface hf
!
!
contains
!
!
   function new_hf() result(wf)
!!
!!    New HF
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(hf) :: wf
!
      wf%name_ = 'rhf'
!
      wf%cumulative_fock_threshold  = 1.0d0
      wf%cumulative_fock            = .false.
!
      call wf%read_settings()
      call wf%print_banner()
!
   end function new_hf
!
!
   subroutine read_for_scf_restart_hf(wf)
!!
!!    Read for SCF restart
!!    Written by Sarai D. Folkestad and Linda Goletto, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%read_orbitals()
      call wf%update_ao_density()
!
   end subroutine read_for_scf_restart_hf
!
!
   subroutine print_energy_hf(wf)
!!
!!    Print energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified by Tommaso Giovannini, March 2019
!!
!!    Prints information related to the wavefunction, most of which is meaningful
!!    only for a properly converged wavefunction. Should be overwritten in descendants
!!    if more or less or other information is present.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp) :: homo_lumo_gap, nuclear_repulsion
!
      if (wf%n_v > 0) then
!
         homo_lumo_gap = wf%orbital_energies(wf%n_o + 1) - wf%orbital_energies(wf%n_o)
         call output%printf('m', 'HOMO-LUMO gap:             (f19.12)', &
                            reals=[homo_lumo_gap], fs='(/t6,a)')
!
      endif
!
      nuclear_repulsion = wf%get_nuclear_repulsion()    
!
      call output%printf('m', 'Nuclear repulsion energy:  (f19.12)', &
                         reals=[nuclear_repulsion], fs='(t6,a)')

!
      call output%printf('m', 'Electronic energy:         (f19.12)', &
                         reals=[wf%energy - nuclear_repulsion], fs='(t6,a)')
      call output%printf('m', 'Total energy:              (f19.12)', &
                         reals=[wf%energy], fs='(t6,a)')
!
      if (wf%embedded) call wf%embedding%print_energy(wf%ao, wf%ao_density)
!
   end subroutine print_energy_hf
!
!
   subroutine read_settings_hf(wf)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Designed to be overwritten by descendants with
!!    descendant specific settings (see e.g. UHF).
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%read_hf_settings()
      call wf%read_frozen_orbitals_settings()
!
   end subroutine read_settings_hf
!
!
   subroutine read_hf_settings_hf(wf)
!!
!!    Read HF settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call input%get_keyword('coulomb threshold',              &
                                        'solver scf',          &
                                        wf%coulomb_threshold)
!
      call input%get_keyword('exchange threshold',             &
                                        'solver scf',          &
                                        wf%exchange_threshold)
!
      call input%get_keyword('integral cutoff',                &
                                        'solver scf',          &
                                        wf%integral_cutoff)      
!
      call input%get_keyword('cumulative fock threshold',      &
                                        'solver scf',          &
                                        wf%cumulative_fock_threshold)  
!
   end subroutine read_hf_settings_hf
!
!
   subroutine print_summary_hf(wf, print_mo_info)
!!
!!    Print Summary
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none 
!
      class(hf), intent(inout) :: wf
!
      logical, intent(in) :: print_mo_info
!
      call output%printf('m', '- Summary of '// &
                         &trim(convert_to_uppercase(wf%name_))// ' wavefunction &
                         &energetics (a.u.):', fs='(/t3,a)')
!
      call wf%print_energy()
!
      if (print_mo_info) call wf%save_orbital_info()
!
   end subroutine print_summary_hf
!
!
   subroutine save_orbital_info_hf(wf)
!!
!!    Save orbital info
!!    Written by Alexander C. Paul, Nov 2020
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      wf%mo_information_file = output_file('eT.mo_information.out')
      call wf%mo_information_file%open_('rewind')
!
      call wf%print_orbitals_and_energies(wf%orbital_energies,     &
                                          wf%orbital_coefficients, &
                                          '- Molecular orbital')
!
      call wf%mo_information_file%close_()
!
   end subroutine save_orbital_info_hf
!
!
   subroutine print_orbitals_and_energies_hf(wf, mo_energies, mo_coefficients, prefix)
!!
!!    Print orbitals and energies
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!    Modified by Alexander C. Paul to print all MOs to file, Dec 2019
!!    Modified by ALexander C. Paul prints MO energies to file as well, Oct 2020
!!
!!    Prints the orbital energies and the orbitals with atom & orbital information given.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo), intent(in) :: mo_energies
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(in) :: mo_coefficients
!
      character(len=*), intent(in) :: prefix
      character(len=200) :: header
!
      write(header, '(a,a)') trim(prefix), ' energies'
!
      call wf%mo_information_file%print_vector('m', header, wf%n_mo, mo_energies, &
                                               fs='(f16.12)', columns=4)
!
      write(header, '(a,a)') trim(prefix), ' coefficients'
      call wf%mo_information_file%printf('m', header, fs='(//t3,a)')
!
      call wf%print_orbitals_from_coefficients(mo_coefficients)
!
   end subroutine print_orbitals_and_energies_hf
!
!
   subroutine print_orbitals_from_coefficients_hf(wf, orbital_coefficients)
!!
!!    Print orbitals from coefficients
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!
!!    Modified by Alexander C. Paul, Dec 2019. Printing of l and m_l.
!!    Modified by Eirik F. Kjønstad, Sep 2020. Printing now delegated to AO tool.
!!
!!    Prints the orbitals from coefficients with atom & orbital information given.
!!
!!       orbital_coefficients: (n_ao, n_mo) array containing orbital coefficients
!!       the_file:             output file where the coefficients are to be printed
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(hf),                             intent(in) :: wf
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(in) :: orbital_coefficients 
!
      integer, parameter :: n_entries = 5
!
      integer :: mo_offset, first_mo, last_mo
!
!     Print at most n_entries (5) MOs at a time 
!
      do mo_offset = 1, wf%n_mo, n_entries
!
         first_mo = mo_offset
         last_mo  = min(mo_offset + n_entries - 1, wf%n_mo)
!
!        Print the current set of MO vectors 
!
         call wf%ao%print_ao_vectors(C        = orbital_coefficients(:, first_mo : last_mo),&
                                     out_file = wf%mo_information_file,                     &
                                     m        = last_mo - first_mo + 1,                     &
                                     offset   = first_mo - 1)
!
      enddo
!
   end subroutine print_orbitals_from_coefficients_hf
!
!
   subroutine set_initial_ao_density_guess_hf(wf, guess)
!!
!!    Set initial AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Sets initial AO density (or densities) to the
!!    appropriate initial guess requested by the solver.
!!
      implicit none
!
      class(hf) :: wf
!
      character(len=*) :: guess
!
      real(dp), dimension(:,:), allocatable :: h_wx
!
      if (trim(guess) == 'sad' .or. trim(guess) == 'SAD') then
!
         call wf%ao%get_sad_guess(wf%ao_density)
!
      elseif (trim(guess) == 'core' .or. trim(guess) == 'CORE') then
!
         call mem%alloc(h_wx, wf%ao%n, wf%ao%n)
         call wf%ao%get_oei('hamiltonian', h_wx)
!
         call wf%set_ao_density_to_core_guess(h_wx)
!
         call mem%dealloc(h_wx, wf%ao%n, wf%ao%n)
!
      else
!
         call output%error_msg('Guess AO density ' // trim(guess) // ' is currently not supported.')
!
      endif
!
      call zero_array(wf%previous_ao_density, wf%ao%n**2*wf%n_densities)
!
   end subroutine set_initial_ao_density_guess_hf
!
!
   subroutine update_ao_density_hf(wf)
!!
!!    Update AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Updates the AO density (or densities, if unrestricted) based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(hf) :: wf
!
      call dcopy(wf%ao%n**2, wf%ao_density, 1, wf%previous_ao_density, 1)
!
      call wf%construct_ao_density()
!
   end subroutine update_ao_density_hf
!
!
   subroutine cleanup_hf(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
!     Deallocations
!
      call wf%save_ao_density()
!
      call wf%destruct_orbital_energies()
      call wf%destruct_orbital_coefficients()
      call wf%destruct_fock()
      call wf%destruct_ao_density()
!
      call wf%destruct_mo_fock()
!
      call wf%destruct_mo_fock_frozen()
      call wf%destruct_orbital_coefficients_fc()
      call wf%destruct_orbital_coefficients_frozen_hf()
      call wf%destruct_frozen_CCT()
!
      deallocate(wf%ao)
      if (wf%embedded) deallocate(wf%embedding)
!
   end subroutine cleanup_hf
!
!
   subroutine prepare_for_cc_hf(wf)
!!
!!    Prepare for CC
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Prepares frozen Fock terms and sets the HF energy
!!
      implicit none
!
      class(hf) :: wf
!
      wf%hf_energy = wf%energy
!
!     Check if there are frozen Fock contributions; if not, return
!
      if ((.not. wf%frozen_core) .and. &
          (.not. wf%frozen_hf_MOs) .and. &
          (.not. wf%embedded)) then
!
         call output%printf('v', 'No frozen fock contributions!', fs='(/t3,a)')
         wf%exists_frozen_fock_terms = .false.
         return
!
      endif
!
      wf%exists_frozen_fock_terms = .true.
!
      call wf%initialize_frozen_CCT()
      call zero_array(wf%frozen_CCT, wf%ao%n**2)
!
!     Change the MOs if frozen core or frozen hf is requested
!
      call wf%prepare_mos()
!
!     Prepare frozen Fock terms from frozen core and frozen HF
!
      call wf%calculate_frozen_dipole_moment()
      call wf%calculate_frozen_quadrupole_moment()
!
      call wf%prepare_frozen_fock_terms()
!
   end subroutine prepare_for_cc_hf
!
!
   subroutine construct_ao_density_hf(wf)
!!
!!    Construct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Routine which calculates D_αβ = 2 sum_i C_αi C_βi, where C are the MO coefficients.
!!
      implicit none
!
      class(hf) :: wf
!
     call dgemm('N', 'T',                   &
                 wf%ao%n,                   &
                 wf%ao%n,                   &
                 wf%n_o,                    &
                 two,                       &
                 wf%orbital_coefficients,   &
                 wf%ao%n,                   &
                 wf%orbital_coefficients,   &
                 wf%ao%n,                   &
                 zero,                      &
                 wf%ao_density,             &
                 wf%ao%n)
!
   end subroutine construct_ao_density_hf
!
!
   real(dp) function calculate_hf_energy_from_G_hf(wf, half_GD_wx, h_wx) result(hf_energy)
!!
!!    Calculate HF energy from G(D)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates the Hartree-Fock energy,
!!
!!       E = Tr(h D) + 1/4 * Tr(D G(D)),
!!
!!    where D is the AO density and
!!
!!       G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ
!!
!!    The traces are calculated as dot products (since B is symmetric here):
!!
!!       Tr(AB) = sum_x (AB)_xx = sum_xy A_xy B_yx = sum_xy A_xy B_xy.
!!
!!    NOTE: the nuclear repulsion is not included
!!
!!    Modified by Sarai D. Folkestad, Oct 2019
!!
!!    Changed from subroutine to function.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: half_GD_wx
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: h_wx
!
      hf_energy = zero
!
      hf_energy = hf_energy + ddot((wf%ao%n)**2, h_wx, 1, wf%ao_density, 1)
      hf_energy = hf_energy + half*ddot((wf%ao%n)**2, wf%ao_density, 1, half_GD_wx, 1)

   end function calculate_hf_energy_from_G_hf
!
!
   function calculate_hf_energy_from_fock_hf(wf, F_wx, h_wx) result(hf_energy)
!!
!!    Calculate HF energy from Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates the Hartree-Fock energy,
!!
!!       E = Tr(h D) + 1/4 * Tr(D G(D)) = 1/2 * Tr (F D) + 1/2 * Tr(h D),
!!
!!    where D is the AO density and
!!
!!       G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ
!!
!!    The traces are calculated as dot products (since B is symmetric here):
!!
!!       Tr(AB) = sum_x (AB)_xx = sum_xy A_xy B_yx = sum_xy A_xy B_xy = A-dot-B.
!!
!!    Adds the embedding contribution if there is embedding
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: F_wx
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: h_wx
!
      real(dp) :: hf_energy
!
      hf_energy = wf%get_nuclear_repulsion()
!
      hf_energy = hf_energy + half*ddot((wf%ao%n)**2, h_wx, 1, wf%ao_density, 1)
      hf_energy = hf_energy + half*ddot((wf%ao%n)**2, wf%ao_density, 1, F_wx, 1)
!
      if (wf%embedded) hf_energy = hf_energy + wf%embedding%get_energy(wf%ao, wf%ao_density)
!
   end function calculate_hf_energy_from_fock_hf
!
!
   subroutine decompose_ao_density_hf(wf)
!!
!!    Decompose AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Does a Cholesky decomposition of the AO density matrix,
!!
!!       D^AO_wx = sum_J L_w,J L_J,x^T,
!!
!!    and sets the MO coefficients accordingly.
!!
!!    Modified by Rolf H. Myhre and Alexander C. Paul, Oct 2019
!!    included sanity check
!!
      implicit none
!
      class(hf) :: wf
!
      integer, dimension(:), allocatable :: used_diag
!
      integer :: rank
!
      call mem%alloc(used_diag, wf%ao%n)
!
      call dscal(wf%ao%n**2, half, wf%ao_density, 1)
!
      call full_cholesky_decomposition_system(wf%ao_density,            &
                                              wf%orbital_coefficients,  &
                                              wf%ao%n,                  &
                                              rank,                     &
                                              1.0d-12,                  &
                                              used_diag)
!
      call dscal(wf%ao%n**2, two, wf%ao_density, 1)
!
      if (any(used_diag .gt. wf%ao%n) .or. any(used_diag .le. 0)) then
!
         call output%printf('m', 'Something went wrong when decomposing the AO density.', &
                            fs='(/t3,a)')
!
         call output%error_msg('Trying to access elements outside of an array.')
!
      end if
!
      call mem%dealloc(used_diag, wf%ao%n)
!
   end subroutine decompose_ao_density_hf
!
!
   subroutine construct_projection_matrices_hf(wf, Po, Pv, D)
!!
!!    Construct projection matrices
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs Po = 1/2 D S
!!               Pv = 1 - Po
!!
!!    D: AO density matrix, S: AO overlap matrix. Both are assumed
!!    to be allocated and properly set.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: Po
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: Pv
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
!
      type(timings), allocatable :: timer 
!
      timer = timings('Construction of projection matrices (Po, Pv)', 'verbose')
      call timer%turn_on()
!
!     Po = 1/2 D S
!
      call dgemm('N','N',        &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  half,          &
                  D,             &
                  wf%ao%n,       &
                  wf%ao%s,       &
                  wf%ao%n,       &
                  zero,          &
                  Po,            &
                  wf%ao%n)
!
!     Pv = I - Po
!
      call identity_array(Pv, wf%ao%n)
!
      call daxpy(wf%ao%n**2, -one, Po, 1, Pv, 1)
!
      call timer%turn_off()
!
   end subroutine construct_projection_matrices_hf
!
!
   subroutine get_packed_roothan_hall_gradient_hf(wf, G)
!!
!!    Get packed Roothan-Hall gradient
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Constructs and returns the RH gradient as an anti-symmetric packed array
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo*(wf%n_mo - 1)/2, wf%n_densities), intent(inout) :: G
!
      real(dp), dimension(:,:), allocatable :: G_sq, Po, Pv
!
      call mem%alloc(Po, wf%ao%n, wf%ao%n)
      call mem%alloc(Pv, wf%ao%n, wf%ao%n)
!
      call wf%construct_projection_matrices(Po, Pv, wf%ao_density)
!
      call mem%alloc(G_sq, wf%n_mo, wf%n_mo)
!
      call wf%construct_roothan_hall_gradient(G_sq, Po, Pv, wf%ao_fock)
!
      call mem%dealloc(Po, wf%ao%n, wf%ao%n)
      call mem%dealloc(Pv, wf%ao%n, wf%ao%n)
!
      call packin_anti(G(:,1), G_sq, wf%n_mo)
!
      call mem%dealloc(G_sq, wf%n_mo, wf%n_mo)
!
   end subroutine get_packed_roothan_hall_gradient_hf
!
!
   subroutine construct_roothan_hall_gradient_hf(wf, G, Po, Pv, F)
!!
!!    Construct Roothan-Hall gradient
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       G = Fov - Fvo = Po^T F Pv - Pv^T F Po,
!!
!!    where Po = 1/2 D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: Po
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: Pv
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: F
!
      real(dp), dimension(wf%n_mo, wf%n_mo) :: G
!
      real(dp), dimension(:, :), allocatable :: tmp, G_ao 
!
      type(timings), allocatable :: timer 
!
      timer = timings('Construct RH gradient', 'verbose')
      call timer%turn_on()
!
!     G_ao = Fov = Po^T F Pv
!
      call mem%alloc(G_ao, wf%ao%n, wf%ao%n)
!
      call dcopy(wf%ao%n**2, F, 1, G_ao, 1)
!
      call sandwich(G_ao, Po, Pv, wf%ao%n)
!
!     tmp = Fvo = Pv^T F Po
!
      call mem%alloc(tmp, wf%ao%n, wf%ao%n)
!
      call dcopy(wf%ao%n**2, F, 1, tmp, 1)
!
      call sandwich(tmp, Pv, Po, wf%ao%n)
!
!     G_ao = G_ao - tmp = Fov - Fvo
!
      call daxpy(wf%ao%n**2, -one, tmp, 1, G_ao, 1)
!
      call mem%dealloc(tmp, wf%ao%n, wf%ao%n)
!
!     Transform G_ao to OAO pivot basis: G = P^T G_ao P 
!
      call wf%ao%orthonormal_ao_pivot_basis_transformation(G_ao, G)
!
      call mem%dealloc(G_ao, wf%ao%n, wf%ao%n)
!
      call timer%turn_off()
!
   end subroutine construct_roothan_hall_gradient_hf
!
!
!
   function get_n_electrons_in_density_hf(wf) result(n)
!!
!!    Get number of electrons in density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the number of electrons in the current AO density, a useful check on 
!!    the sensibility of the AO density. It is by default printed by some solvers.
!!
!!    Calculated as 
!!
!!       n = Tr(D S) =  sum_wx D_wx S_xw = sum_wx D_wx S_wx 
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp) :: n
!
      real(dp) :: ddot
!
      n = ddot(wf%ao%n**2, wf%ao_density, 1, wf%ao%s, 1)
!
   end function get_n_electrons_in_density_hf
!
!
   subroutine set_ao_density_to_core_guess_hf(wf, h_wx)
!!
!!    Set AO density to core guess
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Solves the Roothan-Hall equation ignoring the two-
!!    electron AO density dependent part of the Fock matrix,
!!    giving an initial density from the resulting orbital
!!    coefficients.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: h_wx
!
      call dcopy(wf%ao%n**2, h_wx, 1, wf%ao_fock, 1)
!
      call wf%construct_initial_idempotent_density()
!
   end subroutine set_ao_density_to_core_guess_hf
!
!
   subroutine print_screening_settings_hf(wf)
!!
!!    Print screening settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call output%printf('n', '- Screening and integral thresholds:', &
                        fs='(/t3,a)')
!
      call output%printf('n', 'Coulomb screening threshold:  (e11.4)', &
                         reals=[wf%coulomb_threshold], fs='(/t6,a)')
!
      call output%printf('n', 'Exchange screening threshold: (e11.4)', &
                         reals=[wf%exchange_threshold], fs='(t6,a)')
!
      call output%printf('n', 'Integral cutoff:              (e11.4)', &
                         reals=[wf%integral_cutoff], fs='(t6,a)')
!
      call output%printf('n', 'Cumulative Fock threshold:    (e11.4)', &
                         reals=[wf%cumulative_fock_threshold], fs='(t6,a)')
!
   end subroutine print_screening_settings_hf
!
!
   subroutine set_n_mo_hf(wf)
!!
!!    Set number of molecular orbitals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      wf%n_mo = wf%ao%get_n_orthonormal_ao()
!
      wf%n_o = (wf%ao%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      call output%printf('m', '- Molecular orbital details:', &
                         fs='(/t3,a)')
!
      call output%printf('m', 'Number of occupied orbitals:  (i8)', ints=[wf%n_o], fs='(/t6,a)')
      call output%printf('m', 'Number of virtual orbitals:   (i8)', ints=[wf%n_v], fs='(t6,a)')
      call output%printf('m', 'Number of molecular orbitals: (i8)', ints=[wf%n_mo], fs='(t6,a)')
!
   end subroutine set_n_mo_hf
!
!
   subroutine construct_molecular_gradient_hf(wf, E_qk)
!!
!!    Contruct molecular gradient
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
!!    Constructs the molecular gradient,
!!
!!       E^x = Tr[D h^x] + (1/2)Tr[D G^x(D)] - (1/2)Tr[D F D S^x] + h_nuc^x.
!!
!!    Here, x denotes the energy in the x direction. In the code,
!!    x = (q,k), where q denotes the component (x,y, or z) and k
!!    denotes the atom (k = 1,2,3,...,n_atoms).
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(3, wf%n_atomic_centers), intent(inout) :: E_qk ! Molecular gradient
!
      real(dp), dimension(:,:,:,:), allocatable :: h_wxqk, h_wxqk_nuc
      real(dp), dimension(:,:,:,:), allocatable :: G_wxqk
      real(dp), dimension(:,:,:,:), allocatable :: s_wxqk
!
      real(dp), dimension(:,:), allocatable :: DFD, FD, E_kq
!
      real(dp), dimension(3, wf%n_atomic_centers) :: TrDh_qk, TrDG_qk, TrDFDS_qk
!
      real(dp) :: ddot
!
      integer :: k, q
!
      type(timings) :: s_timer, h_timer, G_timer, G_timer_sym, non_integral_timer, timer 
!
!     Initialize timers
!
      timer = timings('HF gradient', pl='normal')
!
      s_timer = timings('HF gradient - 1st derivative of S', pl='verbose')
!
      h_timer = timings('HF gradient - 1st derivative of h', pl='verbose')
!
      G_timer = timings('HF gradient - 1st derivative of G(D) - integrals', pl='verbose')
!
      G_timer_sym = timings('HF gradient - 1st derivative of G(D) - symmetrization', pl='verbose')
!
      non_integral_timer = timings('HF gradient - non-integral-time', pl='verbose')
!
      call timer%turn_on()
!
!     Construct h_nuc^x, and the AO integral derivatives, h^x, S^x, and G^x(D)
!
      E_qk = wf%get_nuclear_repulsion_1der() ! 1der = first derivative
!
      call mem%alloc(h_wxqk, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
      call mem%alloc(G_wxqk, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
      call mem%alloc(s_wxqk, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
!
      call s_timer%turn_on()
!
      call wf%ao%get_oei_1der('overlap', s_wxqk)
!
      call s_timer%turn_off()
!
      call G_timer%turn_on()
!
      call wf%construct_ao_G_1der(G_wxqk, wf%ao_density)
!
      call G_timer%turn_off()
      call G_timer_sym%turn_on()
!
      do k = 1, wf%n_atomic_centers

         do q = 1, 3

           call symmetric_sum(G_wxqk(:,:,q,k), wf%ao%n)
           call dscal(wf%ao%n**2, half, G_wxqk(:,:,q,k), 1)

         enddo

      enddo
!
      call G_timer_sym%turn_off()
!
      call h_timer%turn_on()
!
      call wf%ao%get_oei_1der('kinetic', h_wxqk)
!
      call mem%alloc(h_wxqk_nuc, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)

      call wf%ao%get_oei_1der('nuclear', h_wxqk_nuc)

      call daxpy(wf%ao%n**2 * 3 * wf%n_atomic_centers, one, h_wxqk_nuc, 1, h_wxqk, 1)

      call mem%dealloc(h_wxqk_nuc, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
!
      call h_timer%turn_off()
!
      call non_integral_timer%turn_on()
!
!     Construct D F D
!
      call mem%alloc(FD, wf%ao%n, wf%ao%n)
!
      call dgemm('N','N',        &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  one,           &
                  wf%ao_fock,    &
                  wf%ao%n,       &
                  wf%ao_density, &
                  wf%ao%n,       &
                  zero,          &
                  FD,            &
                  wf%ao%n)
!
      call mem%alloc(DFD, wf%ao%n, wf%ao%n)
!
      call dgemm('N','N',        &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  one,           &
                  wf%ao_density, &
                  wf%ao%n,       &
                  FD,            &
                  wf%ao%n,       &
                  zero,          &
                  DFD,           &
                  wf%ao%n)
!
      call mem%dealloc(FD, wf%ao%n, wf%ao%n)
!
!     Perform the traces, adding the contributions to the gradient
!
      do k = 1, wf%n_atomic_centers
         do q = 1, 3
!
            TrDh_qk(q,k)   = ddot(wf%ao%n**2, wf%ao_density, 1, h_wxqk(:,:,q,k), 1)
            TrDG_qk(q,k)   = ddot(wf%ao%n**2, wf%ao_density, 1, G_wxqk(:,:,q,k), 1)
            TrDFDS_qk(q,k) = ddot(wf%ao%n**2, DFD, 1, s_wxqk(:,:,q,k), 1)
!
            E_qk(q,k) = E_qk(q,k)            &
                      + TrDh_qk(q,k)         &
                      + half*TrDG_qk(q,k)    &
                      - half*TrDFDS_qk(q,k)
!
         enddo
      enddo
!
      call mem%alloc(E_kq, wf%n_atomic_centers, 3)
!
      call sort_12_to_21(E_qk, E_kq, 3, wf%n_atomic_centers)
!
      call output%print_matrix('n', 'Molecular gradient (Hartree/bohr)', E_kq, wf%n_atomic_centers, 3)
!
      call mem%dealloc(E_kq, wf%n_atomic_centers, 3)
!
      call mem%dealloc(DFD, wf%ao%n, wf%ao%n)
!
      call non_integral_timer%turn_off()
!
      call mem%dealloc(h_wxqk, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
      call mem%dealloc(G_wxqk, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
      call mem%dealloc(s_wxqk, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers)
!
      call timer%turn_off()
!
   end subroutine construct_molecular_gradient_hf
!
!
   subroutine prepare_for_roothan_hall_hf(wf)
!!
!!    Prepare for Roothan-Hall
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs the necessary preparations needed to solve
!!    of the Roothan-Hall equation in the iterative cycle
!!    construct Fock - calculate energy - Roothan-Hall-update orbitals -
!!    update the density:
!!
!!    - constructs the ao Fock matrix and
!!      performs a Roothan-Hall step to get the
!!      initial idempotent density;
!!    - prints the number of electrons and the energy
!!      of the initial guess.
!!
!!    NOTE: this routine is overwritten for MLHF!
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: n_electrons
!
      call wf%update_fock_and_energy(cumulative = .false.)
!
      n_electrons = wf%get_n_electrons_in_density()
!
      call output%printf('m', 'Energy of initial guess:      (f25.12)', &
                         reals=[wf%energy], fs='(/t6, a)')
      call output%printf('m', 'Number of electrons in guess: (f25.12)', &
                         reals=[n_electrons], fs='(t6, a)')
!
      call wf%construct_initial_idempotent_density()
!
   end subroutine prepare_for_roothan_hall_hf
!
!
   subroutine prepare_hf(wf, centers, embedding)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initializes files, writes the restart file used for consistency checks
!!    and constructs screening vectors
!!
      use atomic_center_class, only: atomic_center
!
      implicit none
!
      class(hf) :: wf
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers 
!
      logical, intent(in), optional :: embedding 
!
      wf%orbital_file = stream_file('orbital_coefficients')
!
      call wf%prepare_ao_tool(centers)
      call wf%prepare_embedding(embedding)
      if (wf%embedded) call wf%embedding%print_description
!
      wf%n_densities = 1
!
      call wf%set_n_mo()
!
      wf%gradient_dimension = wf%n_mo*(wf%n_mo - 1)/2
!
   end subroutine prepare_hf
!
!
   function calculate_expectation_value_hf(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Calculate the expectation value of a one-electron
!!    operator Â
!!
!!    Modified by Anders Hutcheson and Linda Goletto, Oct 2019
!!
!!    The expectation values is calculated as:
!!
!!       < A > = < HF | Â | HF > = sum_pq A_pq D_pq
!!
!!    where A_pq are the integrals in the AO basis
!!    and D_pq is the a one-electron density matrix in the AO basis
!!
      implicit none
!  
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: A
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: density
!
      real(dp) :: expectation_value
!
      real(dp) :: ddot
!
      expectation_value = ddot(wf%ao%n**2, A, 1, density, 1)
!
   end function calculate_expectation_value_hf
!
!
   subroutine print_banner_hf(wf)
!!
!!    Print banner
!!    Sarai D. Folkestad, Dec 2019
!!
      use string_utilities, only : convert_to_uppercase
!
      implicit none
!
      class(hf) :: wf
!
      character(len=200) :: name_
!
      name_ = trim(convert_to_uppercase(wf%name_)) // ' wavefunction'
!
      call output%printf('m', ':: (a0)', chars=[name_], fs='(//t3,a)')
      call output%print_separator('minimal', len_trim(name_) + 6, '=')
!
   end subroutine print_banner_hf
!
!
   subroutine prepare_frozen_fock_terms_hf(wf)
!!
!!    Prepare frozen Fock contributions
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    This routine prepares the frozen Fock contributions
!!    to coupled cluster.
!!
!!    Possible contributions to frozen fock:
!!
!!       - Frozen core
!!
!!       - Frozen HF orbitals
!!
!!       - Embedding QM/MM or PCM fock
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: mo_mm_fock
      real(dp), dimension(:,:), allocatable :: mo_fc_fock
      real(dp), dimension(:,:), allocatable :: mo_frozen_hf_fock
!
      call wf%initialize_mo_fock_frozen()
      call zero_array(wf%mo_fock_frozen, wf%n_mo**2)
!
!     Contribution from frozen core orbitals
!
      if (wf%frozen_core) then
!
         call mem%alloc(mo_fc_fock, wf%n_mo, wf%n_mo)
!
         call wf%construct_mo_fock_fc_term(mo_fc_fock)
         call daxpy(wf%n_mo**2, one, mo_fc_fock, 1, wf%mo_fock_frozen, 1)
!
         call mem%dealloc(mo_fc_fock, wf%n_mo, wf%n_mo)
!
      endif
!
!     Contributions from frozen HF orbitals
!
      if (wf%frozen_hf_mos) then
!
         call mem%alloc(mo_frozen_hf_fock, wf%n_mo, wf%n_mo)
!
         call wf%construct_mo_fock_frozen_hf_term(mo_frozen_hf_fock)
         call daxpy(wf%n_mo**2, one, mo_frozen_hf_fock, 1, wf%mo_fock_frozen, 1)
!
         call mem%dealloc(mo_frozen_hf_fock, wf%n_mo, wf%n_mo)
!
      endif
!
!     Contribution from embedding
!
      if (wf%embedded) then
!
        call mem%alloc(mo_mm_fock, wf%n_mo, wf%n_mo)
!
        call wf%mo_transform(wf%ao%v, mo_mm_fock)
!
        call daxpy(wf%n_mo**2, one, mo_mm_fock, 1, wf%mo_fock_frozen, 1)
!
        call mem%dealloc(mo_mm_fock, wf%n_mo, wf%n_mo)
!
     endif
!
   end subroutine prepare_frozen_fock_terms_hf
!
!
   subroutine read_frozen_orbitals_settings_hf(wf)
!!
!!    Read frozen orbitals 
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Reads the frozen orbitals section of the input
!!
!!     - Frozen core 
!!
!!    - Frozen hf orbitals
!!
!!    - plot active density
!!
!!    This routine is used at HF level to prepare mos and 
!!    frozen fock contributions.
!!
!!    This routine is read at cc level to figure out if there
!!    should be frozen fock contributions
!!
      implicit none
!
      class(hf) :: wf
!
      wf%frozen_core    = .false.
      wf%frozen_hf_mos  = .false.
!
      wf%plot_active_density = .false.
!
      wf%frozen_core = input%is_keyword_present('core', 'frozen orbitals')
      wf%frozen_hf_mos = input%is_keyword_present('hf', 'frozen orbitals')
!
      wf%plot_active_density = input%is_keyword_present('plot hf active density', &
            'visualization')
!
      if (wf%plot_active_density .and. .not.  (wf%frozen_core .or. wf%frozen_hf_mos)) &
         call output%warning_msg('no active density for CC to plot in HF, no plots produced.')
!
   end subroutine read_frozen_orbitals_settings_hf
!
!
   subroutine flip_final_orbitals_hf(wf)
!!
!!    Flip final orbitals
!!    Written by Sarai D. Folkestad, Apr 2020
!!
!!    Ensures that orbitals will have the same
!!    signs each time HF converged.
!!
!!    Finds the first element of each MO with magnitude larger than
!!    0.01 and flips it if it is < 0
!!
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      integer :: p, w
!
      do p = 1, wf%n_mo
!
         do w = 1, wf%ao%n
!
            if (abs(wf%orbital_coefficients(w,p)) .gt. 1.0d-2) then
!  
!              Found the first significant contribution to MO p
!
               if (wf%orbital_coefficients(w,p) .lt. 0.0d0) then
!
!                 First contribution is negative, so we flip MO p
!
                  call dscal(wf%ao%n, -one, wf%orbital_coefficients(1,p), 1)
!
               endif

               exit
!
            endif
!
         enddo
!
      enddo
!
      call wf%save_orbitals()
!
   end subroutine flip_final_orbitals_hf
!
!
   subroutine prepare_for_scf_hf(wf, restart, skip, ao_density_guess, &
      gradient_threshold, dim_, n_densities, gradient_dimension)
!!
!!    Prepare for SCF
!!    Written by Sarai D. Folkestad
!!
!!    1. Start guess
!!
!!    2. Sets dimensions
!!
!!    3. Prints screenings
!!
      implicit none
!
      class(hf)                     :: wf
!
      logical, intent(in)           :: restart   
      logical, intent(in)           :: skip   
      character(len=*), intent(in)  :: ao_density_guess 
      real(dp)                      :: gradient_threshold
      integer, intent(out)          :: dim_
      integer, intent(out)          :: n_densities
      integer, intent(out)          :: gradient_dimension
!
      wf%cumulative_fock = .false. 
!
      call wf%set_screening_and_precision_thresholds(gradient_threshold)
!
      call wf%initialize_orbitals()
      call wf%initialize_density()
      call wf%initialize_fock()
!
      call zero_array(wf%ao_density, wf%ao%n**2)
!
      if (restart .or. skip) then
!
         call output%printf('m', '- Requested restart. Reading orbitals from file', &
                            fs='(/t3,a)')
!
         call wf%read_for_scf_restart()
!
         if (skip) then

            if (.not. wf%control_gradient_convergence(gradient_threshold)) &
               call output%error_msg('cannot skip scf when gradient has not converged.')
         endif
!
      else 
!
         call output%printf('m', '- Setting initial AO density to &
                            &'//trim(ao_density_guess), fs='(/t3,a)')
!
         call wf%set_initial_ao_density_guess(ao_density_guess)
         call wf%prepare_for_roothan_hall()
!
      endif 
!
      dim_ = wf%n_mo
      n_densities = wf%n_densities
      gradient_dimension = wf%gradient_dimension
!
      call wf%print_screening_settings()
!
   end subroutine prepare_for_scf_hf
!
!
   function get_energy_hf(wf) result(E)
!!
!!    Get energy
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(hf) :: wf
      real(dp)  :: E
!
      E = wf%energy
!
   end function get_energy_hf
!
!
   subroutine set_C_and_e_hf(wf, C, e)
!!
!!    Set C and e
!!    Written by Sarai D. Folkestad, 2020
!! 
!!    Sets the orbital coefficients from the orbital 
!!    coefficients in the orthonormal AO basis (C)
!!
!!    Sets the orbital energies (e)
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, wf%n_densities), intent(in)  :: C
      real(dp), dimension(wf%n_mo, wf%n_densities), intent(in)           :: e
!
!     Transform back the solutions to original basis, C = P L^-T (L^-1 P^T C) = P L^-T C'
!
      call wf%set_orbital_coefficients_from_orthonormal_ao_C(C(:,:,1), wf%orbital_coefficients)
!
      call dcopy(wf%n_mo, e, 1, wf%orbital_energies, 1)
!
      call wf%update_ao_density()! C => D
!
      call wf%save_orbitals()
!
   end subroutine set_C_and_e_hf
!
!
   subroutine get_F_hf(wf, F_packed)
!!
!!    Get F
!!    Written by Sarai D. Folkestad
!!
!!    Constructs the Fock matrix in the orthonormal AO basis
!!    and returns it packed
!!
      use reordering, only: packin 
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_mo*(wf%n_mo + 1)/2*wf%n_densities), intent(out) :: F_packed
!
      real(dp), dimension(:,:), allocatable :: F
!
      call wf%update_fock_and_energy(wf%cumulative_fock)
!
      call mem%alloc(F, wf%n_mo, wf%n_mo)
!
      call wf%ao_to_orthonormal_ao_transformation(F, wf%ao_fock)
      call packin(F_packed, F, wf%n_mo)
!
      call mem%dealloc(F, wf%n_mo, wf%n_mo)
!
   end subroutine get_F_hf
!
!
   subroutine get_gradient_hf(wf, G)
!!
!!    Get gradient
!!    Written by Sarai D. Folkestad
!!
!!    Returns the packed gradient
!!
!!    If the gradient norm is sufficiently small, 
!!    'cumulative_fock' is enabled
!!
      class(hf) :: wf
!
      real(dp), dimension(wf%gradient_dimension), intent(out) :: G
!
      call wf%get_packed_roothan_hall_gradient(G)
!
      wf%cumulative_fock = wf%can_do_cumulative_fock(G)
!
   end subroutine get_gradient_hf
!
!
   subroutine construct_initial_idempotent_density_hf(wf)
!!
!!    Construct initial idempotent density
!!    Written by Sarai D. Folkestad
!!
!
      use array_utilities, only: generalized_diagonalization_symmetric
!
      implicit none
!
      class(hf), intent(inout)               :: wf
!
      real(dp), dimension(:,:), allocatable :: F, S
!
      call mem%alloc(F, wf%n_mo, wf%n_mo)
      call mem%alloc(S, wf%n_mo, wf%n_mo)
!
      call wf%ao_to_reduced_ao_transformation(F, wf%ao_fock)
      call wf%ao%get_reduced_ao_metric(S) ! RAO metric
!
      call generalized_diagonalization_symmetric(F, S, wf%n_mo, wf%orbital_energies)
      call wf%set_orbital_coefficients_from_reduced_ao_C(F, wf%orbital_coefficients)
!
      call wf%save_orbitals()
!
      call mem%dealloc(F, wf%n_mo, wf%n_mo)
      call mem%dealloc(S, wf%n_mo, wf%n_mo)
!
      call wf%update_ao_density()
!
   end subroutine construct_initial_idempotent_density_hf
!
!
   subroutine set_orbital_coefficients_from_reduced_ao_C_hf(wf, C, orbital_coefficients)
!!
!!    Set orbital coefficients from reduced AO C
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2020 
!! 
!!    Sets the orbital coefficients from the orbital coefficients
!!    in the reduced ao basis 
!!
!!       orbital_coefficients = P C
!!
!!    where
!!
!!       P^T S P = L L^T
!!
!!    defines the Cholesky decomposition of the 
!!    atomic orbital overlap matrix S.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)   :: C
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(out)  :: orbital_coefficients
!
      call dgemm('N','N',                 &
                  wf%ao%n,                &
                  wf%n_mo,                &
                  wf%n_mo,                &
                  one,                    &
                  wf%ao%P,                &
                  wf%ao%n,                &
                  C,                      & 
                  wf%n_mo,                &
                  zero,                   &
                  orbital_coefficients,   &
                  wf%ao%n)
!
   end subroutine set_orbital_coefficients_from_reduced_ao_C_hf
!
!
   subroutine set_orbital_coefficients_from_orthonormal_ao_C_hf(wf, C, orbital_coefficients)
!!
!!    Set orbital coefficients from orthonormal AO C
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2020 
!! 
!!    Sets the orbital coefficients from the orbital coefficients
!!    in the reduced ao basis 
!!
!!       orbital_coefficients = P L^-T C
!!
!!    where
!!
!!       P^T S P = L L^T
!!
!!    defines the Cholesky decomposition of the 
!!    atomic orbital overlap matrix S.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)  :: C
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(out)  :: orbital_coefficients
!
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: Y
!
      integer :: info
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call dcopy(wf%n_mo**2, wf%ao%L, 1, X, 1)
!
      call dtrtri('l','n', wf%n_mo, X, wf%n_mo, info)
!
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      call dgemm('T','N',  &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  X,       &
                  wf%n_mo, &
                  C,       &
                  wf%n_mo, &
                  zero,    &
                  Y,       &
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
!
      call dgemm('N','N',                 &
                  wf%ao%n,                &
                  wf%n_mo,                &
                  wf%n_mo,                &
                  one,                    &
                  wf%ao%P,                &
                  wf%ao%n,                &
                  Y,                      &
                  wf%n_mo,                &
                  zero,                   &
                  orbital_coefficients,   &
                  wf%ao%n)
!
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
!
   end subroutine set_orbital_coefficients_from_orthonormal_ao_C_hf
!
!
   subroutine ao_to_reduced_ao_transformation_hf(wf, X_reduced_ao, X_ao)
!!
!!    AO to RAO transformation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2020
!!
!!    Transforms a matrix from the AO basis to the
!!    reduced AO (RAO) basis
!!
!!       X_reduced_ao = P^T X_ao P, 
!!
!!    where P is a projection onto the linearly idependent AO basis
!!    and is obtained from the Cholesky decomposition  of 
!!    the atomic orbital overlap matrix S:
!!
!!       P^T S P = L L^T
!!
      implicit none
!
      class(hf) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out)  :: X_reduced_ao
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in)   :: X_ao
!
      real(dp), dimension(:,:), allocatable :: XP
!
      call mem%alloc(XP, wf%ao%n, wf%n_mo)
!
      call dgemm('N','N',    &
                  wf%ao%n,   &
                  wf%n_mo,   &
                  wf%ao%n,   &
                  one,       &
                  X_ao,      &
                  wf%ao%n,   &
                  wf%ao%P,   &
                  wf%ao%n,   &
                  zero,      &
                  XP,        &
                  wf%ao%n)
!
      call dgemm('T','N',    &
                  wf%n_mo,   &
                  wf%n_mo,   &
                  wf%ao%n,   &
                  one,       &
                  wf%ao%P,   &
                  wf%ao%n,   &
                  XP,        &
                  wf%ao%n,   &
                  zero,      &
                  X_reduced_ao,     & 
                  wf%n_mo)
!
      call mem%dealloc(XP, wf%ao%n, wf%n_mo)      
!
   end subroutine ao_to_reduced_ao_transformation_hf
!
!
   subroutine ao_to_orthonormal_ao_transformation_hf(wf, X_orthonormal_ao, X_ao)
!!
!!    AO to OAO transformation
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Transforms a matrix from the AO basis to the
!!    orthonormal AO (OAO) basis 
!!
!!       X_orthonormal_ao = L^-1 P^T X_ao P L^-T, 
!!
!!    where
!!
!!       P^T S P = L L^T
!!
!!    defines the Cholesky decomposition of the 
!!    atomic orbital overlap matrix S.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: X_orthonormal_ao
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in)  :: X_ao
!
      real(dp), dimension(:,:), allocatable :: X, L_inv
      real(dp), dimension(:,:), allocatable :: X_reduced_ao
!
      integer :: info
!
      call mem%alloc(L_inv, wf%n_mo, wf%n_mo)
      call dcopy(wf%n_mo**2, wf%ao%L, 1, L_inv, 1)
!
      call dtrtri('l','n', wf%n_mo, L_inv, wf%n_mo, info)
!
      if (info .ne. 0) call output%error_msg('problem with inversion of Cholesky factor of S')
!
!     Construct reduced AO Fock matrix, F' = P^T F P
!
      call mem%alloc(X_reduced_ao, wf%n_mo, wf%n_mo)
!
      call wf%ao_to_reduced_ao_transformation(X_reduced_ao, X_ao)
!
!     Construct orthogonal AO Fock matrix, X_orthonormal_ao = L^-1 X_reduced_ao L^-T
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'N',    &
                  wf%n_mo,    &
                  wf%n_mo,    &
                  wf%n_mo,    &
                  one,        &
                  L_inv,      &
                  wf%n_mo,    &
                  X_reduced_ao,      &
                  wf%n_mo,    &
                  zero,       &
                  X,          &
                  wf%n_mo)
!
      call mem%dealloc(X_reduced_ao, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'T',    &
                  wf%n_mo,    &
                  wf%n_mo,    &
                  wf%n_mo,    &
                  one,        &
                  X,          &
                  wf%n_mo,    &
                  L_inv,      &
                  wf%n_mo,    &
                  zero,       &
                  X_orthonormal_ao,      &
                  wf%n_mo)
!
      call mem%dealloc(L_inv, wf%n_mo, wf%n_mo)
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
!
   end subroutine ao_to_orthonormal_ao_transformation_hf
!
!
   function is_restart_possible_hf(wf) result(restart_is_possible)
!!
!!    Is restart possible
!!    Written by Alexander C. Paul, Okt 2020
!!
!!    Checks if restart files exist
!!
      implicit none
!
      class(hf) :: wf
!
      logical :: restart_is_possible
!
      restart_is_possible = (wf%orbital_file%exists())
!
   end function is_restart_possible_hf
!
!
   function control_gradient_convergence_hf(wf, threshold) result(converged)
!!
!!    Control gradient convergence
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Checks if gradient has converged to the given threshold.
!!
!!    Prints energy and gradient.
!!
      implicit none
!
      class(hf),  intent(inout)  :: wf
      real(dp),   intent(in)     :: threshold
!
      logical :: converged
      real(dp) :: max_gradient
!
      real(dp), dimension(:,:), allocatable :: G
!
      if (wf%n_mo == 1) then
!
         converged = .true.
         return
!
      endif
!
      converged = .false.
!
      call wf%update_fock_and_energy(cumulative = .false.)
!
      call mem%alloc(G, wf%n_mo*(wf%n_mo-1)/2, wf%n_densities)
!
      call wf%get_packed_roothan_hall_gradient(G)
      max_gradient = get_abs_max(G, wf%n_mo*(wf%n_mo-1)/2*wf%n_densities)
!
      call mem%dealloc(G, wf%n_mo*(wf%n_mo-1)/2, wf%n_densities)
!
      if (max_gradient .lt. threshold) converged = .true.
!
   end function control_gradient_convergence_hf
!
!
   pure function can_do_cumulative_fock_hf(wf, G) result(cumulative)
!!
!!    Do cumulative Fock
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp), dimension(wf%gradient_dimension), intent(in) :: G
!
      logical :: cumulative
!
      cumulative = .false.
      if (get_abs_max(G, wf%gradient_dimension) .lt. wf%cumulative_fock_threshold) &
         cumulative = .true.
!
   end function can_do_cumulative_fock_hf
!
!
end module hf_class
