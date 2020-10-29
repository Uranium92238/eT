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
   use timings_class,   only : timings
   use array_utilities, only : get_abs_max, sandwich
   use array_utilities, only : full_cholesky_decomposition_system
   use array_utilities, only : get_n_highest
   use array_utilities, only : identity_array
!
   use libint_initialization, only : set_coulomb_precision_c
!
   use omp_lib
!
   implicit none
!
!  Hartree-Fock hf
!
   type, extends(wavefunction) :: hf
!
      real(dp), dimension(:,:), allocatable :: ao_density
!
      integer :: n_densities ! For RHF, there is one AO density. 
                             ! For UHF, there are two AO densities (alpha, beta).
!
      real(dp), dimension(:,:), allocatable :: ao_h
!
      real(dp) :: coulomb_threshold  = 1.0D-12   ! Screening threshold (Fock, Coulomb)
      real(dp) :: exchange_threshold = 1.0D-10   ! Screening threshold (Fock, exchange)
      real(dp) :: libint_epsilon     = 1.0D-20   ! Libint electron repulsion integral 
                                                 ! precision given approximately by sqrt(epsilon)
      real(dp) :: integral_cutoff    = 1.0D-10   ! Default: sqrt(epsilon) 
!
      real(dp), dimension(:,:), allocatable :: shp_eri_schwarz      ! Screening for (wx | yz)
      integer,  dimension(:,:), allocatable :: shp_eri_schwarz_list ! Indices for screening 
!
      real(dp), dimension(:,:), allocatable :: W_mo_update ! Eigenvectors for 
                                                           ! Roothan-Hall in MO basis
!
      real(dp), dimension(:,:), allocatable :: ao_overlap   ! S_wx = < w | x >
!
      real(dp) :: linear_dep_threshold = 1.0D-6 ! threshold for decomposition 
                                                ! of AO overlap: S = (P L) (P L)^T
!
      real(dp), dimension(:,:), allocatable :: cholesky_ao_overlap      ! L in S = (P L) (P L)^T
      real(dp), dimension(:,:), allocatable :: pivot_matrix_ao_overlap  ! P in S = (P L) (P L)^T
!
!     QM/MM quantities 
!
      real(dp) :: electrostatic_energy_qmmm
      real(dp) :: energy_qmmm
      real(dp) :: energy_scf_qmmm
!
      real(dp) :: electrostatic_energy_qmpcm
      real(dp) :: energy_scf_qmpcm
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
      type(sequential_file) :: orbital_coefficients_file
      type(sequential_file) :: orbital_energies_file
!
      logical :: frozen_core
      logical :: frozen_hf_mos
!
!     Restart files 
!
      type(sequential_file) :: restart_file
      type(sequential_file) :: orbital_information_file
!
      logical :: plot_active_density
!
   contains
!
      procedure :: print_banner                                => print_banner_hf
!
      procedure :: print_summary                               => print_summary_hf
!
!     Read, save of orbital energies and coefficients 
!
      procedure :: read_orbital_coefficients                   => read_orbital_coefficients_hf
      procedure :: save_orbital_coefficients                   => save_orbital_coefficients_hf
      procedure :: read_orbital_energies                       => read_orbital_energies_hf
      procedure :: save_orbital_energies                       => save_orbital_energies_hf
!
!     Read, write restart information, and check for safe restart 
!
      procedure :: read_for_scf_restart                        => read_for_scf_restart_hf
      procedure :: read_for_scf_restart_mo                     => read_for_scf_restart_mo_hf
      procedure :: is_restart_safe                             => is_restart_safe_hf
      procedure :: write_scf_restart                           => write_scf_restart_hf
      procedure :: write_orbital_information                   => write_orbital_information_hf
      procedure :: is_restart_possible                         => is_restart_possible_hf
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                                     => cleanup_hf
      procedure :: prepare_for_cc                              => prepare_for_cc_hf
!
      procedure :: read_settings                               => read_settings_hf
      procedure :: read_hf_settings                            => read_hf_settings_hf
      procedure :: construct_ao_overlap                        => construct_ao_overlap_hf
      procedure :: decompose_ao_overlap                        => decompose_ao_overlap_hf
      procedure :: print_energy                                => print_energy_hf
      procedure :: print_energy_mm                             => print_energy_mm_hf
      procedure :: print_energy_pcm                            => print_energy_pcm_hf
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
      procedure :: construct_mo_fock                           => construct_mo_fock_hf
      procedure :: set_ao_fock                                 => set_ao_fock_hf
      procedure :: get_ao_fock                                 => get_ao_fock_hf
      procedure :: get_fock_ov                                 => get_fock_ov_hf
      procedure :: calculate_hf_energy_from_fock               => calculate_hf_energy_from_fock_hf
      procedure :: calculate_mm_energy_terms                   => calculate_mm_energy_terms_hf
      procedure :: calculate_pcm_energy_terms                  => calculate_pcm_energy_terms_hf
      procedure :: calculate_hf_energy_from_G                  => calculate_hf_energy_from_G_hf
      procedure :: initialize_fock                             => initialize_fock_hf
      procedure :: destruct_fock                               => destruct_fock_hf
      procedure :: update_fock_and_energy_non_cumulative       => update_fock_and_energy_non_cumulative_hf
      procedure :: update_fock_and_energy_cumulative           => update_fock_and_energy_cumulative_hf
      procedure :: update_fock_and_energy                      => update_fock_and_energy_hf
      procedure :: add_mm_fock_terms                           => add_mm_fock_terms_hf
      procedure :: add_pcm_fock_term                           => add_pcm_fock_term_hf
      procedure :: initialize_ao_h                             => initialize_ao_h_hf
      procedure :: destruct_ao_h                               => destruct_ao_h_hf
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
      procedure :: set_ao_density_to_sad                       => set_ao_density_to_sad_hf
      procedure :: set_ao_density_to_core_guess                => set_ao_density_to_core_guess_hf
      procedure :: get_n_electrons_in_density                  => get_n_electrons_in_density_hf
      procedure :: construct_shp_density_schwarz               => construct_shp_density_schwarz_hf
!
!     MO orbital related routines
!
      procedure :: do_roothan_hall                             => do_roothan_hall_hf
      procedure :: initialize_orbitals                         => initialize_orbitals_hf
      procedure :: roothan_hall_update_orbitals                => roothan_hall_update_orbitals_hf
      procedure :: print_orbital_energies                      => print_orbital_energies_hf
      procedure :: print_orbitals                              => print_orbitals_hf
      procedure :: print_orbitals_from_coefficients            => print_orbitals_from_coefficients_hf
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density                       => initialize_ao_density_hf
      procedure :: initialize_ao_overlap                       => initialize_ao_overlap_hf
      procedure :: initialize_pivot_matrix_ao_overlap          => initialize_pivot_matrix_ao_overlap_hf
      procedure :: initialize_cholesky_ao_overlap              => initialize_cholesky_ao_overlap_hf
!
      procedure :: destruct_ao_density                         => destruct_ao_density_hf
      procedure :: destruct_ao_overlap                         => destruct_ao_overlap_hf
      procedure :: destruct_pivot_matrix_ao_overlap            => destruct_pivot_matrix_ao_overlap_hf
      procedure :: destruct_cholesky_ao_overlap                => destruct_cholesky_ao_overlap_hf
!
!     Gradient and Hessian related routines
!
      procedure :: construct_projection_matrices               => construct_projection_matrices_hf
!
      procedure :: construct_roothan_hall_gradient             => construct_roothan_hall_gradient_hf
      procedure :: get_packed_roothan_hall_gradient            => get_packed_roothan_hall_gradient_hf
      procedure :: get_max_roothan_hall_gradient               => get_max_roothan_hall_gradient_hf
!
      procedure :: construct_molecular_gradient                => construct_molecular_gradient_hf
!
!     Integral related routines
!
      procedure :: initialize_shp_eri_schwarz                  => initialize_shp_eri_schwarz_hf
      procedure :: destruct_shp_eri_schwarz                    => destruct_shp_eri_schwarz_hf
!
      procedure :: initialize_shp_eri_schwarz_list             => initialize_shp_eri_schwarz_list_hf
      procedure :: destruct_shp_eri_schwarz_list               => destruct_shp_eri_schwarz_list_hf
!
      procedure :: construct_shp_eri_schwarz                   => construct_shp_eri_schwarz_hf
      procedure :: get_n_sig_eri_shp                           => get_n_sig_eri_shp_hf
!
      procedure :: set_n_mo                                    => set_n_mo_hf
!
      procedure :: set_screening_and_precision_thresholds      => set_screening_and_precision_thresholds_hf
      procedure :: print_screening_settings                    => print_screening_settings_hf
!
      procedure :: prepare_for_roothan_hall                    => prepare_for_roothan_hall_hf
      procedure :: prepare                                     => prepare_hf
      procedure :: prepare_qmmm                                => prepare_qmmm_hf
!
!     Zop
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
!     MO-SCF routines
!
      procedure :: get_max_roothan_hall_mo_gradient            => get_max_roothan_hall_mo_gradient_hf
!  
      procedure :: update_fock_and_energy_mo                   => update_fock_and_energy_mo_hf
      procedure :: get_roothan_hall_mo_gradient                => get_roothan_hall_mo_gradient_hf
!  
      procedure :: do_roothan_hall_mo                          => do_roothan_hall_mo_hf
!  
      procedure :: initialize_W_mo_update                      => initialize_W_mo_update_hf
      procedure :: destruct_W_mo_update                        => destruct_W_mo_update_hf
!  
      procedure :: roothan_hall_update_orbitals_mo             => roothan_hall_update_orbitals_mo_hf
      procedure :: prepare_for_roothan_hall_mo                 => prepare_for_roothan_hall_mo_hf
!
      procedure :: flip_final_orbitals                         => flip_final_orbitals_hf
!
   end type hf
!
   interface
!
      include "frozen_orbital_hf_interface.F90"
      include "mo_hf_interface.F90"
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
   function new_hf(system) result(wf)
!!
!!    New HF
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(hf) :: wf
!
      class(molecular_system), target, intent(in) :: system
!
      wf%name_ = 'rhf'
!
      wf%system => system
!
      call wf%read_settings()
      call wf%print_banner()
      call wf%prepare()
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
      call wf%read_orbital_coefficients()
      call wf%update_ao_density()
      call wf%read_orbital_energies()
!
   end subroutine read_for_scf_restart_hf
!
!
   subroutine is_restart_safe_hf(wf)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(hf) :: wf
!
      integer :: n_ao, n_densities, n_electrons
!
      call wf%restart_file%open_('read', 'rewind')
!
      call wf%restart_file%read_(n_ao)
      call wf%restart_file%read_(n_densities)
      call wf%restart_file%read_(n_electrons)
!
      call wf%restart_file%close_
!
      if (n_ao .ne. wf%n_ao) then
         call output%error_msg('Attempted to restart HF with an inconsistent number ' // &
                               'of atomic orbitals.')
      endif
!
      if (n_densities .ne. wf%n_densities) then
         call output%error_msg('Attempted to restart HF with an inconsistent number ' // &
                               'of atomic densities (likely a HF/UHF inconsistency).')
      endif
!
      if (n_electrons .ne. wf%system%get_n_electrons()) then
         call output%error_msg('Attempted to restart HF with an inconsistent number ' // &
                               'of electrons.')
      endif
!
   end subroutine is_restart_safe_hf
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
      restart_is_possible = (wf%restart_file%exists() &
                       .and. wf%orbital_coefficients_file%exists() &
                       .and. wf%orbital_energies_file%exists())
!
   end function is_restart_possible_hf
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
      nuclear_repulsion = wf%system%get_total_nuclear_repulsion()
!
      call output%printf('m', 'Nuclear repulsion energy:  (f19.12)', &
                         reals=[nuclear_repulsion], fs='(t6,a)')
      call output%printf('m', 'Electronic energy:         (f19.12)', &
                         reals=[wf%energy - nuclear_repulsion], fs='(t6,a)')
      call output%printf('m', 'Total energy:              (f19.12)', &
                         reals=[wf%energy], fs='(t6,a)')
!
      if(wf%system%mm_calculation) call wf%print_energy_mm()
      if(wf%system%pcm_calculation) call wf%print_energy_pcm()
!
   end subroutine print_energy_hf
!
!
   subroutine print_energy_mm_hf(wf)
!!
!!    Print hf summary for QM/MM calculations
!!    Written by Tommaso Giovannini, March 2019
!!
!!    Modified by Sarai D. Folkestad, Oct 2019
!!      
!!    Moved the calculation of some non-polerizable QM/MM energy terms to 
!!    calculate_mm_energy_terms_hf
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call output%printf('m', '- Summary of QM/MM energetics:', fs='(/t3,a)')
      call output%printf('m', 'a.u.             eV     kcal/mol', fs='(t42,a)')
      call output%printf('m', 'QM/MM SCF Contribution: (f22.12)', &
                         reals=[wf%energy_scf_qmmm], fs='(t6,a)')
      call output%printf('m', 'QM/MM Electrostatic Energy:(f19.12)(f12.5) (f9.3)', &
                         reals=[wf%electrostatic_energy_qmmm, &
                         wf%electrostatic_energy_qmmm*Hartree_to_eV, &
                         wf%electrostatic_energy_qmmm*Hartree_to_kcalmol], fs='(t6,a)')
!
   end subroutine print_energy_mm_hf
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
      call input%get_keyword_in_section('coulomb threshold',   &
                                        'solver scf',          &
                                        wf%coulomb_threshold)
!
      call input%get_keyword_in_section('exchange threshold',  &
                                        'solver scf',          &
                                        wf%exchange_threshold)
!
      call input%get_keyword_in_section('integral precision',  &
                                        'solver scf',          &
                                        wf%libint_epsilon)
!
      call input%get_keyword_in_section('integral cutoff',     &
                                        'solver scf',          &
                                        wf%integral_cutoff)      
!
   end subroutine read_hf_settings_hf
!
!
   subroutine print_orbital_energies_hf(wf)
!!
!!    Print orbital energies
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints the current orbital energies to output
!!    in a hopefully readable way.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      call output%print_vector('normal', '- Molecular orbital energies', wf%n_mo, wf%orbital_energies, &
                              fs='(f16.12)', columns=4)
!
   end subroutine print_orbital_energies_hf
!
!
   subroutine print_summary_hf(wf)
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
      call output%printf('m', '- Summary of '// &
                         &trim(convert_to_uppercase(wf%name_))// ' wavefunction &
                         &energetics (a.u.):', fs='(/t3,a)')
!
      call wf%print_energy()
      call wf%print_orbital_energies()
!
   end subroutine print_summary_hf
!
!
   subroutine print_orbitals_hf(wf)
!!
!!    Print orbitals
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!    Modified by Alexander C. Paul to print all MOs to file, Dec 2019
!!
!!    Prints the orbitals with atom & orbital information given.
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(hf), intent(in) :: wf
!
      type(output_file) :: mo_coefficient_file
!
      mo_coefficient_file = output_file('mo_coefficients.out')
      call mo_coefficient_file%open_()
!
      call wf%print_orbitals_from_coefficients(wf%orbital_coefficients, &
                                               mo_coefficient_file)
!
      call mo_coefficient_file%close_()
!
   end subroutine print_orbitals_hf
!
!
   subroutine print_orbitals_from_coefficients_hf(wf, orbital_coefficients, the_file)
!!
!!    Print orbitals from coefficients
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!    Modified by Alexander C. Paul to print l and m_l, Dec 2019
!!
!!    Prints the orbitals from coefficients with atom & orbital information given.
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(hf),                             intent(in) :: wf
      type(output_file),                     intent(in) :: the_file
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(in) :: orbital_coefficients 
!
      integer, parameter :: n_entries  = 5
!
      integer :: mo_offset, first_mo, last_mo, mo, n_mo
      integer :: atom, shell, ao, l, index_in_shell
!
      logical :: adv
!
      character(len=2)  :: symbol
      character(len=8)  :: ang_mom
      character(len=50) :: n_format_string
      character(len=7)  :: format_string = '(f10.6)'
!
!     Print n_entries columns with wf%n_mo rows
!
      do mo_offset = 1, wf%n_mo, n_entries
!
         first_mo = mo_offset
         last_mo  = min(mo_offset + n_entries - 1, wf%n_mo)
!
         n_mo = last_mo - first_mo + 1
!
!        Print header: AO Atom l ml 1 2 3 4 ..
!
         call the_file%printf('n', '  AO    Atom  l m_l', fs='(/t3,a)', adv=.false.)
!
         do mo = first_mo, last_mo
!
            adv = (mo == last_mo)
!
            call the_file%printf('n', '(i4)', adv=adv, ints=[mo], fs='(8x,a)')
!
         enddo
!
         call the_file%print_separator(pl='normal', n=83, symbol='-')
!
!        Print content: AO, Element symbol, orbital and orbital coefficients
!
         do atom = 1, wf%system%n_atoms
!
            symbol  = trim( wf%system%atoms(atom)%symbol )
!
            do shell = 1, wf%system%atoms(atom)%n_shells
!
               l = wf%system%atoms(atom)%shells(shell)%l
!
               do ao = wf%system%atoms(atom)%shells(shell)%first, &
                       wf%system%atoms(atom)%shells(shell)%last
!
                  index_in_shell = ao - wf%system%atoms(atom)%shells(shell)%first + 1
!
                  call wf%system%atoms(atom)%shells(shell)%get_angular_momentum_label( &
                                 l, index_in_shell, ang_mom, wf%system%atoms(atom)%cartesian)
!
!                 Setup the string for printf to print the right number of reals
!
                  n_format_string = repeat('  ' // format_string, n_mo)
!
                  call the_file%printf('n', '(i4) (i4) (b2)  '// ang_mom // n_format_string, &
                                      reals=[orbital_coefficients(ao, first_mo:last_mo)], &
                                      ints=[ao, atom], chars=[symbol], ll=85)
!
               enddo
!
            enddo
!
         enddo
!
         call the_file%print_separator(pl='normal', n=83, symbol='-')
!
      enddo
!
   end subroutine print_orbitals_from_coefficients_hf
!
!
   subroutine mo_transform_hf(wf, X_wx, Y_pq)
!!
!!    MO transform
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Performs MO transformation of X and saves the result in Y:
!!
!!       Y_pq = sum_wx C_wp X_wx C_xq
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: X_wx
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Y_pq
!
      real(dp), dimension(:,:), allocatable :: Z_wq ! = sum_x X_wx C_xq
!
      call mem%alloc(Z_wq, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                 &
                  wf%n_ao,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  X_wx,                    &
                  wf%n_ao,                 &
                  wf%orbital_coefficients, & ! C_xq
                  wf%n_ao,                 &
                  zero,                    &
                  Z_wq,                    &
                  wf%n_ao)
!
      call dgemm('T', 'N',                 &
                  wf%n_mo,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  wf%orbital_coefficients, & ! C_wp
                  wf%n_ao,                 &
                  Z_wq,                    &
                  wf%n_ao,                 &
                  zero,                    &
                  Y_pq,                    &
                  wf%n_mo)
!
      call mem%dealloc(Z_wq, wf%n_ao, wf%n_mo)
!
   end subroutine mo_transform_hf
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
         call wf%set_ao_density_to_sad()
!
      elseif (trim(guess) == 'core' .or. trim(guess) == 'CORE') then
!
         call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
         call wf%get_ao_h_wx(h_wx)
!
         call wf%set_ao_density_to_core_guess(h_wx)
!
         call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      else
!
         call output%error_msg('Guess AO density ' // trim(guess) // ' is currently not supported.')
!
      endif
!
   end subroutine set_initial_ao_density_guess_hf
!
!
   subroutine roothan_hall_update_orbitals_hf(wf)
!!
!!    Roothan-Hall update of orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    This routine guides the construction of new orbital coefficients
!!    from the current AO Fock matrix (or matrices if the hf
!!    is unrestricted).
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies)
!
   end subroutine roothan_hall_update_orbitals_hf
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
      call wf%destruct_ao_overlap()
      call wf%destruct_fock()
      call wf%destruct_ao_density()
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
      call wf%destruct_shp_eri_schwarz()
      call wf%destruct_shp_eri_schwarz_list()
      call wf%destruct_ao_h()
!
      call wf%destruct_W_mo_update()
      call wf%destruct_mo_fock()
!
      call wf%destruct_mm_matrices()
      call wf%destruct_pcm_matrices()
      call wf%destruct_mo_fock_frozen()
      call wf%destruct_orbital_coefficients_fc()
      call wf%destruct_orbital_coefficients_frozen_hf()
      call wf%destruct_frozen_CCT()
!
   end subroutine cleanup_hf
!
!
   subroutine prepare_for_cc_hf(wf)
!!
!!    Prepare for CC
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Prepares frozen fock terms, 
!!    and places energy in hf_energy
!!
      implicit none
!
      class(hf) :: wf
!
      wf%hf_energy = wf%energy
!
!     Check if there are frozen Fock contributions, if not
!     return
!
      if ((.not. wf%frozen_core) .and. &
          (.not. wf%frozen_hf_MOs) .and. &
          (.not. wf%system%pcm_calculation) .and. &  
          (.not. wf%system%mm_calculation)) then
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
      call zero_array(wf%frozen_CCT, wf%n_ao**2)
!
!     Change the MOs if frozen core or frozen hf 
!     is requested
!
      call wf%prepare_mos()
!
!     Prepare frozen Fock terms from frozen core 
!     and frozen HF
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
!!    Routine which calculates D_αβ = sum_i C_αi C_βi,
!!    where C are the MO coefficients.
!!
!!    Density is packed
!!
      implicit none
!
      class(hf) :: wf
!
     call dgemm('N', 'T',                   &
                 wf%n_ao,                   &
                 wf%n_ao,                   &
                 wf%n_o,                    &
                 two,                       &
                 wf%orbital_coefficients,   &
                 wf%n_ao,                   &
                 wf%orbital_coefficients,   &
                 wf%n_ao,                   &
                 zero,                      &
                 wf%ao_density,             &
                 wf%n_ao)
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: half_GD_wx
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      hf_energy = zero
!
      hf_energy = hf_energy + ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      hf_energy = hf_energy + half*ddot((wf%n_ao)**2, wf%ao_density, 1, half_GD_wx, 1)

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
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: F_wx
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp) :: hf_energy
!
      hf_energy = wf%system%get_nuclear_repulsion()
!
      hf_energy = hf_energy + half*ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      hf_energy = hf_energy + half*ddot((wf%n_ao)**2, wf%ao_density, 1, F_wx, 1)
!
   end function calculate_hf_energy_from_fock_hf
!
!
   subroutine construct_ao_overlap_hf(wf)
!!
!!    Construct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%get_ao_s_wx(wf%ao_overlap)
!
   end subroutine construct_ao_overlap_hf
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
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      integer :: rank, j
!
      call mem%alloc(used_diag, wf%n_ao)
!
      call dscal(wf%n_ao**2, half, wf%ao_density, 1)
!
      call full_cholesky_decomposition_system(wf%ao_density, wf%orbital_coefficients, &
                                              wf%n_ao, rank, 1.0d-12, used_diag)
!
      call dscal(wf%n_ao**2, two, wf%ao_density, 1)
!
!     Make permutation matrix P
!
      call mem%alloc(perm_matrix, wf%n_ao, wf%n_ao)
!
      if (any(used_diag .gt. wf%n_ao) .or. any(used_diag .le. 0)) then
!
         call output%printf('m', 'Something went wrong when decomposing the AO density.', &
                            fs='(/t3,a)')
!
         call output%error_msg('Trying to access elements outside of an array.')
!
      end if
!
      call zero_array(perm_matrix, wf%n_ao**2)
!
!$omp parallel do private(j)
      do j = 1, wf%n_ao
!
         perm_matrix(used_diag(j), j) = one
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(used_diag, wf%n_ao)
      call mem%dealloc(perm_matrix, wf%n_ao, wf%n_ao)
!
   end subroutine decompose_ao_density_hf
!
!
   subroutine decompose_ao_overlap_hf(wf)
!!
!!    Decompose AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs a Cholesky decomposition of the AO overlap matrix S,
!!    to within a given threshold:
!!
!!       P^T S P = L L^T
!!
!!    The routine allocates and sets P and L, the 'permutation matrix'
!!    and the 'cholesky ao overlap'. Moreover, it sets the number of
!!    linearly independent AOs, wf%n_mo, which is sometimes less than
!!    wf%n_ao. From P and L, we can transform equations to the linearly
!!    independent basis and back.
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
      real(dp), dimension(:, :), allocatable :: L
!
      integer :: i, j
!
      type(timings), allocatable :: timer
!
      timer = timings('Cholesky decomposition of AO overlap', 'normal')
      call timer%turn_on()
!
      call mem%alloc(used_diag, wf%n_ao)
!
!$omp parallel do private(j)
      do j = 1, wf%n_ao 
!
         used_diag(j) = 0
!
      enddo
!$omp end parallel do
!
      call mem%alloc(L, wf%n_ao, wf%n_ao) ! Full Cholesky vector
!
      call zero_array(L, wf%n_ao**2)
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, wf%n_mo, &
                                              wf%linear_dep_threshold, used_diag)
!
      call wf%initialize_cholesky_ao_overlap()
!
!$omp parallel do private(j, i)
      do j = 1, wf%n_mo
         do i = 1, wf%n_mo 
!
            wf%cholesky_ao_overlap(i, j) = L(i, j)
!
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(L, wf%n_ao, wf%n_ao)
!
!     Make permutation matrix P
!
      call wf%initialize_pivot_matrix_ao_overlap()
!
      call zero_array(wf%pivot_matrix_ao_overlap, wf%n_ao*wf%n_mo)
!
      if (wf%n_mo .gt. wf%n_ao .or. wf%n_mo .le. 0 .or. &
          any(used_diag .gt. wf%n_ao) .or. any(used_diag .le. 0)) then
!
         call output%printf('m', 'Something went wrong when decomposing the AO overlap.', &
                            fs='(/t3,a)')
!
         call output%printf('m', 'Did you compile with wrong type of integers &
                            &in setup? For example system native BLAS  with &
                            &default 64-bit integers.', ffs='(/t3,a)')
!
         call output%printf('m', 'If that is the case, use setup with --int32 &
                            &or install MKL.')
!
         call output%error_msg('Failed to decompose AO overlap.')
!
      end if
!
!$omp parallel do private(j)
      do j = 1, wf%n_mo
!
         wf%pivot_matrix_ao_overlap(used_diag(j), j) = one
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(used_diag, wf%n_ao)
!
      call timer%turn_off()
!
   end subroutine decompose_ao_overlap_hf
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
      real(dp), dimension(wf%n_ao, wf%n_ao) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao) :: Pv
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      integer :: w, x
!
      type(timings), allocatable :: timer 
!
      timer = timings('Construction of projection matrices (Po, Pv)', 'verbose')
      call timer%turn_on()
!
!     Po = 1/2 D S
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  half,          &
                  D,             &
                  wf%n_ao,       &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  zero,          &
                  Po,            &
                  wf%n_ao)
!
!     Pv = I
!
      call zero_array(Pv, wf%n_ao**2)
!
!$omp parallel do private(x)
      do x = 1, wf%n_ao
!
         Pv(x, x) = one
!
      enddo
!$omp end parallel do
!
!     Pv = I - Po
!
!$omp parallel do private(x, w)
      do x = 1, wf%n_ao
         do w = 1, wf%n_ao
!
            Pv(w, x) = Pv(w, x) - Po(w, x)
!
         enddo
      enddo
!$omp end parallel do
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
!!    Constructs and returns the gradient as an
!!    anti-symmetrized packed vector.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo*(wf%n_mo - 1)/2, wf%n_densities), intent(inout) :: G
!
      real(dp), dimension(:,:), allocatable :: G_sq, Po, Pv
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
!
      call wf%construct_projection_matrices(Po, Pv, wf%ao_density)
!
      call mem%alloc(G_sq, wf%n_mo, wf%n_mo)
!
      call wf%construct_roothan_hall_gradient(G_sq, Po, Pv, wf%ao_fock)
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      call packin_anti(G(:,1), G_sq, wf%n_mo)
!
      call mem%dealloc(G_sq, wf%n_mo, wf%n_mo)
!
   end subroutine get_packed_roothan_hall_gradient_hf
!
!
   function get_max_roothan_hall_gradient_hf(wf) result(max_gradient)
!!
!!    Get max Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Constructs and returns the absolute maximum
!!    of the HF gradient
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp) :: max_gradient
!
      real(dp), dimension(:), allocatable :: G
!
      call mem%alloc(G, wf%n_ao*(wf%n_ao - 1)/2*wf%n_densities)
!
      call wf%get_packed_roothan_hall_gradient(G)
!
      max_gradient = get_abs_max(G, wf%n_ao*(wf%n_ao - 1)/2*wf%n_densities)
!
      call mem%dealloc(G, wf%n_ao*(wf%n_ao - 1)/2*wf%n_densities)
!
   end function get_max_roothan_hall_gradient_hf
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: F
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
!     Construct tmp = Fov = Po^T F Pv and set G = tmp = Fov
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
      call mem%alloc(G_ao, wf%n_ao, wf%n_ao)
!
      call dcopy(wf%n_ao**2, F, 1, tmp, 1)
      call sandwich(tmp, Po, Pv, wf%n_ao)
!
      call dcopy(wf%n_ao**2, tmp, 1, G_ao, 1)
!
!     Construct tmp = Fvo = Pv^T F Po and set H = H - tmp = Fov - Fvo
!
      call dcopy(wf%n_ao**2, F, 1, tmp, 1)
      call sandwich(tmp, Pv, Po, wf%n_ao)
!
      call daxpy(wf%n_ao**2, -one, tmp, 1, G_ao, 1)
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
!     Transform G to linearly independent AO basis 
!
      call mem%alloc(tmp, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                      &
                  wf%n_ao,                      &
                  wf%n_mo,                      &
                  wf%n_ao,                      &
                  one,                          &
                  G_ao,                         &
                  wf%n_ao,                      &
                  wf%pivot_matrix_ao_overlap,   & ! P 
                  wf%n_ao,                      &
                  zero,                         &
                  tmp,                          &
                  wf%n_ao)
!
      call dgemm('T', 'N',                      &
                  wf%n_mo,                      &
                  wf%n_mo,                      &
                  wf%n_ao,                      &
                  one,                          &
                  wf%pivot_matrix_ao_overlap,   & ! P 
                  wf%n_ao,                      &
                  tmp,                          & ! G P  
                  wf%n_ao,                      &
                  zero,                         &
                  G,                            &
                  wf%n_mo)
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_mo)
      call mem%dealloc(G_ao, wf%n_ao, wf%n_ao)
!
      call timer%turn_off()
!
   end subroutine construct_roothan_hall_gradient_hf
!
!
   subroutine do_roothan_hall_hf(wf, F, C, e)
!!
!!    Do Roothan-Hall
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the equation F C = S C e for the orbital coefficients C.
!!    More precisely, it solves the equation in a linearly independent
!!    subspace,
!!
!!       P^T F P (P^T C) = P^T S P P^T C e = L L^T (P^T C) e,
!!
!!    which differs from the AO basis if linear dependence is present
!!    to within some given threshold. The components of the equation is
!!    given by the Cholesky decomposition
!!
!!       P^T S P = L L^T,
!!
!!    where P is referred to as the 'permutation matrix' and L the 'cholesky
!!    ao overlap' (note that these are member variables of the solver which
!!    must be set by a call to solver%decompose_ao_overlap(wf)). The number
!!    of linearly independent orbitals is wf%n_mo, whereas the full number
!!    is wf%n_ao.
!!
!!    Default is to not transform the Fock matrix to the MO basis. If
!!    do_mo_transformation is passed and is set true, the MO Fock matrix
!!    is initialized and transformed to the MO basis:
!!
!!       F_mo = C'^T P^T F P C'
!!
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)     :: F
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(inout)  :: C
      real(dp), dimension(wf%n_mo), intent(inout)           :: e
!
      real(dp), dimension(:), allocatable   :: work
      real(dp), dimension(:,:), allocatable :: metric
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: FP
!
      integer :: info
!
      type(timings), allocatable :: timer 
!
      timer = timings('Solve FC = SCe', 'verbose')
      call timer%turn_on()
!
      call mem%alloc(metric, wf%n_mo, wf%n_mo)
!
      call dgemm('N','T',                 &
                  wf%n_mo,                &
                  wf%n_mo,                &
                  wf%n_mo,                &
                  one,                    &
                  wf%cholesky_ao_overlap, &
                  wf%n_mo,                &
                  wf%cholesky_ao_overlap, &
                  wf%n_mo,                &
                  zero,                   &
                  metric,                 & ! metric = L L^T
                  wf%n_mo)
!
!     Allocate reduced space matrices
!
      call mem%alloc(ao_fock, wf%n_mo, wf%n_mo)
!
!     Construct reduced space Fock matrix, F' = P^T F P,
!     which is to be diagonalized over the metric L L^T
!
      call mem%alloc(FP, wf%n_ao, wf%n_mo)
!
      call dgemm('N','N',                       &
                  wf%n_ao,                      &
                  wf%n_mo,                      &
                  wf%n_ao,                      &
                  one,                          &
                  F,                            &
                  wf%n_ao,                      &
                  wf%pivot_matrix_ao_overlap,   &
                  wf%n_ao,                      &
                  zero,                         &
                  FP,                           &
                  wf%n_ao)
!
      call dgemm('T','N',                       &
                  wf%n_mo,                      &
                  wf%n_mo,                      &
                  wf%n_ao,                      &
                  one,                          &
                  wf%pivot_matrix_ao_overlap,   &
                  wf%n_ao,                      &
                  FP,                           &
                  wf%n_ao,                      &
                  zero,                         &
                  ao_fock,                      & ! F' = P^T F P
                  wf%n_mo)
!
      call mem%dealloc(FP, wf%n_ao, wf%n_mo)
!
!     Solve F'C' = L L^T C' e
!
      info = 0
!
      call mem%alloc(work, 4*wf%n_mo)
      call zero_array(work, 4*wf%n_mo)
!
      call dsygv(1, 'V', 'L',       &
                  wf%n_mo,          &
                  ao_fock,          & ! ao_fock on entry, orbital coefficients on exit
                  wf%n_mo,          &
                  metric,           &
                  wf%n_mo,          &
                  e,                &
                  work,             &
                  4*(wf%n_mo),      &
                  info)
!
      call mem%dealloc(metric, wf%n_mo, wf%n_mo)
      call mem%dealloc(work, 4*wf%n_mo)
!
      if (info .ne. 0) call output%error_msg('Could not solve Roothan-Hall equations.')
!
!     Transform back the solutions to original basis, C = P (P^T C) = P C'
!
      call dgemm('N','N',                       &
                  wf%n_ao,                      &
                  wf%n_mo,                      &
                  wf%n_mo,                      &
                  one,                          &
                  wf%pivot_matrix_ao_overlap,   &
                  wf%n_ao,                      &
                  ao_fock,                      & ! orbital coefficients
                  wf%n_mo,                      &
                  zero,                         &
                  C,                            &
                  wf%n_ao)
!
      call mem%dealloc(ao_fock, wf%n_mo, wf%n_mo)
!
      call timer%turn_off()
!
   end subroutine do_roothan_hall_hf
!
!
   subroutine set_ao_density_to_sad_hf(wf)
!!
!!    Set AO density to SAD
!!    Written by Eirik F. Kjønstad, Aug-Sep 2018
!!
!!    The routine uses the pre-computed library of atomic densities to
!!    assemble the superposition of atomic densities (SAD) start guess.
!!    The library is assumed to be in the location given by the environ-
!!    ment variable "ET_SAD_DIR". Note that the library is located in the
!!    directory
!!
!!       eT/src/molecular_system/sad/ => set ET_SAD_DIR accordingly
!!
!!    The SAD guess is based on ground state UHF calculations, where
!!    the valence electrons are evenly spread out in the degenerate HOMO
!!    orbitals (to ensure a spherically symmetric AO density, which in
!!    turn ensures a rotationally invariant SAD guess for the molecule).
!!
!!    The algorithm is based on the procedure implemented in Psi4, although
!!    we perform spherical averaging of the densities in the UHF calculations
!!    in order to get a unique rotationally invariant SAD density. This is
!!    especially important for us, since we will use the initial idempotent
!!    density that results from a single Roothan-Hall step in multilevel HF
!!    calculations. The atomic calculations are performed with the ground
!!    state multiplicities, as listed in Griffiths, David J.,
!!    "Introduction to Quantum Mechanics." (1995).
!!
      implicit none
!
      class(hf) :: wf
!
      integer :: I, n_ao_on_atom, first_ao_on_atom, last_ao_on_atom, n_s_on_atom
!
      real(dp), dimension(:,:), allocatable :: atomic_density
!
      call zero_array(wf%ao_density, wf%n_ao**2)
!
      do I = 1, wf%system%n_atoms
!
         n_ao_on_atom     = wf%system%atoms(I)%n_ao
         n_s_on_atom      = wf%system%atoms(I)%n_shells
!
         first_ao_on_atom = wf%system%atoms(I)%shells(1)%first
         last_ao_on_atom  = wf%system%atoms(I)%shells(n_s_on_atom)%last
!
         call mem%alloc(atomic_density, n_ao_on_atom, n_ao_on_atom)
!
         call wf%system%atoms(I)%read_atomic_density(atomic_density)
!
         wf%ao_density(first_ao_on_atom : last_ao_on_atom, &
                       first_ao_on_atom : last_ao_on_atom) = atomic_density(1 : n_ao_on_atom, &
                                                                            1 : n_ao_on_atom)
!
         call mem%dealloc(atomic_density, n_ao_on_atom, n_ao_on_atom)
!
      enddo
!
   end subroutine set_ao_density_to_sad_hf
!
!
   subroutine get_n_electrons_in_density_hf(wf, n_electrons)
!!
!!    Get number of electrons in density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the number of electrons in the current
!!    AO density, a useful check on the sensibility of the
!!    initial AO density guess. It is by default printed
!!    by some solvers.
!!
      implicit none
!
      class(hf), intent(in)   :: wf
      real(dp), intent(inout) :: n_electrons
!
      integer :: ao
!
      real(dp), dimension(:,:), allocatable :: DS
!
      call mem%alloc(DS, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  zero,          &
                  DS,            &
                  wf%n_ao)
!
      n_electrons = zero
!
!$omp parallel do private(ao) reduction(+:n_electrons)
      do ao = 1, wf%n_ao
!
         n_electrons = n_electrons + DS(ao, ao)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(DS, wf%n_ao, wf%n_ao)
!
   end subroutine get_n_electrons_in_density_hf
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      call dcopy(wf%n_ao**2, h_wx, 1, wf%ao_fock, 1)
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies)
      call wf%construct_ao_density()
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
      call output%printf('n', 'Coulomb screening threshold:  (e11.4)', &
                         reals=[wf%coulomb_threshold], fs='(/t6,a)')
!
      call output%printf('n', 'Exchange screening threshold: (e11.4)', &
                         reals=[wf%exchange_threshold], fs='(t6,a)')
!
      call output%printf('n', 'Fock precision:               (e11.4)', &
                         reals=[wf%libint_epsilon], fs='(t6,a)')
!
      call output%printf('n', 'Integral cutoff:              (e11.4)', &
                         reals=[wf%integral_cutoff], fs='(t6,a)')
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
      call output%printf('n', '- Cholesky decomposition of AO overlap to get &
                         &linearly independent orbitals:', ffs='(/t3,a)')
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap()
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      call output%printf('m', '- Orbital details:', &
                         fs='(/t3,a)')
!
      call output%printf('m', 'Number of occupied orbitals:  (i8)', ints=[wf%n_o], fs='(/t6,a)')
      call output%printf('m', 'Number of virtual orbitals:   (i8)', ints=[wf%n_v], fs='(t6,a)')
      call output%printf('m', 'Number of molecular orbitals: (i8)', ints=[wf%n_mo], fs='(t6,a)')
      call output%printf('m', 'Number of atomic orbitals:    (i8)', ints=[wf%n_ao], fs='(t6,a)')
!
      if (wf%n_mo .lt. wf%n_ao) &
         call output%printf('m', '- Removed (i0) AOs due to linear dependencies', &
            ints=[wf%n_ao - wf%n_mo], fs='(/t3,a)')
!
   end subroutine set_n_mo_hf
!
!
   subroutine calculate_mm_energy_terms_hf(wf) 
!!
!!    Calculate MM energy contributions
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Calculates QM/MM energy contributions.
!!    Adds the QM/MM contribution to the HF energy.
!!
!!    Modified by Sarai D. Folkestad, Oct 2019
!!
!!    Moved some energy calculation from the print_energy_mm_hf routine
!!    and added it to the non-polarizable part of this routine.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(:), allocatable :: nopol_potential
!
      wf%electrostatic_energy_qmmm    = zero
!
!     energy terms due to non-polarizable QM/MM
!
      if(wf%system%mm%forcefield.eq.'non-polarizable') then
!
!        change the nuclear repulsion
!
         wf%energy = wf%energy + wf%system%get_nuclear_repulsion_mm()
!
         wf%energy = wf%energy + one/two*ddot((wf%n_ao)**2, &
                                       half*wf%nopol_h_wx, 1,wf%ao_density, 1)
!
         wf%energy_scf_qmmm  = + one/two*ddot((wf%n_ao)**2, one*wf%nopol_h_wx, 1, wf%ao_density, 1)
!
        call mem%alloc(nopol_potential, wf%system%mm%n_atoms)
!
        call wf%construct_ao_electrostatics(0, 1, 'prop', wf%system%mm%n_atoms, wf%system%mm%coordinates, &
                     property_points=nopol_potential, ao_density=wf%ao_density)
!
        wf%electrostatic_energy_qmmm = zero
!
        wf%electrostatic_energy_qmmm = wf%electrostatic_energy_qmmm - &
                                          dot_product(wf%system%mm%charge,nopol_potential)
!
        wf%energy_qmmm = wf%electrostatic_energy_qmmm
!
        call mem%dealloc(nopol_potential, wf%system%mm%n_atoms)
!
      endif
!
!     energy terms due to polarizable QM/FQ
!
      if(wf%system%mm%forcefield .eq. 'fq') then
!
         wf%energy = wf%energy - one/four*ddot((wf%n_ao)**2, wf%ao_density, &
                                          1, wf%pol_emb_fock, 1) &
                               - one/two*dot_product(wf%system%mm%pol_emb_lhs,wf%system%mm%pol_emb_rhs)
!
         wf%energy_scf_qmmm  = - one/two*dot_product(wf%system%mm%pol_emb_lhs,wf%system%mm%pol_emb_rhs)
!
         wf%system%mm%pol_emb_rhs(1:wf%system%mm%n_atoms) = wf%system%mm%pol_emb_rhs(1:wf%system%mm%n_atoms) &
                                                            + wf%system%mm%chi
         wf%electrostatic_energy_qmmm = - dot_product(wf%system%mm%pol_emb_lhs,wf%system%mm%pol_emb_rhs)
!
      endif
!
   end subroutine calculate_mm_energy_terms_hf
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
      real(dp), dimension(3, wf%system%n_atoms), intent(inout) :: E_qk ! Molecular gradient
!
      real(dp), dimension(:,:,:,:), allocatable :: h_wxqk
      real(dp), dimension(:,:,:,:), allocatable :: G_wxqk
      real(dp), dimension(:,:,:,:), allocatable :: s_wxqk
!
      real(dp), dimension(:,:), allocatable :: DFD, FD
!
      real(dp), dimension(3, wf%system%n_atoms) :: TrDh_qk, TrDG_qk, TrDFDS_qk
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
      E_qk = wf%system%get_nuclear_repulsion_1der() ! E_qk = h_nuc_qk
!
      call mem%alloc(h_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%alloc(G_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%alloc(s_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
!
      call s_timer%turn_on()
!
      call wf%get_ao_s_wx_1der(s_wxqk)
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
      do k = 1, wf%system%n_atoms

         do q = 1, 3

           call symmetric_sum(G_wxqk(:,:,q,k), wf%n_ao)
           call dscal(wf%n_ao**2, half, G_wxqk(:,:,q,k), 1)

         enddo

      enddo
!
      call G_timer_sym%turn_off()
!
      call h_timer%turn_on()
!
      call wf%get_ao_h_wx_1der(h_wxqk)
!
      call h_timer%turn_off()
!
      call non_integral_timer%turn_on()
!
!     Construct D F D
!
      call mem%alloc(FD, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_fock,    &
                  wf%n_ao,       &
                  wf%ao_density, &
                  wf%n_ao,       &
                  zero,          &
                  FD,            &
                  wf%n_ao)
!
      call mem%alloc(DFD, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  FD,            &
                  wf%n_ao,       &
                  zero,          &
                  DFD,           &
                  wf%n_ao)
!
      call mem%dealloc(FD, wf%n_ao, wf%n_ao)
!
!     Perform the traces, adding the contributions to the gradient
!
      do k = 1, wf%system%n_atoms
         do q = 1, 3
!
            TrDh_qk(q,k)   = ddot(wf%n_ao**2, wf%ao_density, 1, h_wxqk(:,:,q,k), 1)
            TrDG_qk(q,k)   = ddot(wf%n_ao**2, wf%ao_density, 1, G_wxqk(:,:,q,k), 1)
            TrDFDS_qk(q,k) = ddot(wf%n_ao**2, DFD, 1, s_wxqk(:,:,q,k), 1)
!
            E_qk(q,k) = E_qk(q,k)            &
                      + TrDh_qk(q,k)         &
                      + half*TrDG_qk(q,k)    &
                      - half*TrDFDS_qk(q,k)
!
         enddo
      enddo
!
      call output%printf('m', 'Molecular gradient (Hartree/bohr):', fs='(/t6,a/)')
!
      do k = 1, wf%system%n_atoms
!
         call output%printf('m', '(a2)(f19.12)(f19.12)(f19.12)', &
                            chars=[wf%system%atoms(k)%symbol], reals=[E_qk(1, &
                            k), E_qk(2, k), E_qk(3, k)], fs='(t6,a)')
!
      enddo
!
      call mem%dealloc(DFD, wf%n_ao, wf%n_ao)
!
      call non_integral_timer%turn_off()
!
      call mem%dealloc(h_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%dealloc(G_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%dealloc(s_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
!
      call timer%turn_off()
!
   end subroutine construct_molecular_gradient_hf
!
!
   subroutine calculate_pcm_energy_terms_hf(wf)
!!
!!    Calculate HF energy from F for QM/MM methods
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Adds the QM/MM contribution to the HF energy. 
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      wf%electrostatic_energy_qmpcm = zero
!
      wf%energy = wf%energy - one/four*ddot((wf%n_ao)**2, wf%ao_density, 1, wf%pcm_fock, 1) &
                            - one/two*dot_product(wf%system%pcm%charges,wf%system%pcm%pcm_rhs)
!
      wf%energy_scf_qmpcm  = - one/two*dot_product(wf%system%pcm%charges,wf%system%pcm%pcm_rhs)
!
      wf%electrostatic_energy_qmpcm = two * wf%energy_scf_qmpcm
!
!
   end subroutine calculate_pcm_energy_terms_hf
!
!
   subroutine print_energy_pcm_hf(wf)
!!
!!    Print wavefunction summary for QM/MM calculations
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!      
      call output%printf('m', '- Summary of QM/PCM energetics:', fs='(/t3,a)')
      call output%printf('m', 'a.u.             eV     kcal/mol', fs='(t42,a)')
      call output%printf('m', 'QM/PCM SCF Contribution: (f21.12)', &
                         reals=[wf%energy_scf_qmpcm], fs='(t6,a)')
      call output%printf('m', 'QM/PCM Electrostatic Energy:(f18.12)(f12.5) (f9.3)', &
                         reals=[wf%electrostatic_energy_qmpcm, &
                         wf%electrostatic_energy_qmpcm*Hartree_to_eV, &
                         wf%electrostatic_energy_qmpcm*Hartree_to_kcalmol], fs='(t6,a)')
!
!
   end subroutine print_energy_pcm_hf
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
!!    NOTE: this routine is overwritten
!!          for MLHF!
!!
      implicit none
!
      class(hf) :: wf
!

      real(dp) :: n_electrons
!
      call wf%update_fock_and_energy()
!
      call wf%get_n_electrons_in_density(n_electrons)
!
      call output%printf('m', 'Energy of initial guess:      (f25.12)', &
                         reals=[wf%energy], fs='(/t6, a)')
      call output%printf('m', 'Number of electrons in guess: (f25.12)', &
                         reals=[n_electrons], fs='(t6, a)')
!
!     Update the orbitals and density to make sure the density is idempotent
!     (not the case for the standard atomic superposition density)
!
      call wf%roothan_hall_update_orbitals() ! F => C
      call wf%update_ao_density() ! C => D
!
   end subroutine prepare_for_roothan_hall_hf
!
!
   subroutine prepare_hf(wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initializes files, writes the restart file used for consistency checks
!!    and constructs screening vectors
!!
      implicit none
!
      class(hf) :: wf
!
      wf%n_ao        = wf%system%get_n_aos()
      wf%n_densities = 1
!
      call wf%set_n_mo()
!
      wf%orbital_coefficients_file = sequential_file('orbital_coefficients')
      wf%orbital_energies_file = sequential_file('orbital_energies')
!
      wf%restart_file = sequential_file('scf_restart_file')
!
      call wf%initialize_shp_eri_schwarz()
      call wf%initialize_shp_eri_schwarz_list()
!
      call wf%construct_shp_eri_schwarz()
!
      if (wf%system%mm_calculation) call wf%prepare_qmmm()
      if (wf%system%pcm_calculation) call wf%initialize_pcm_matrices()
!
      call wf%initialize_ao_h()
      call wf%get_ao_h_wx(wf%ao_h)
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: A
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: density
!
      real(dp) :: expectation_value
!
      real(dp) :: ddot
!
      expectation_value = ddot(wf%n_ao**2, A, 1, density, 1)
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
      if (wf%system%mm_calculation) then
!
         call output%printf('m', 'This is a QM/MM calculation', fs='(/t3,a)')
         call wf%system%mm%print_description()
!
      endif
!
      if (wf%system%pcm_calculation)  then
!
         call output%printf('m', 'This is a PCM calculation', fs='(/t3,a)')
         call wf%system%pcm%print_description_and_settings()
!
      endif
!
   end subroutine print_banner_hf
!
!
   subroutine prepare_qmmm_hf(wf)
!!
!!    Prepare QM/MM 
!!    Written by Tommaso Giovannini
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%initialize_mm_matrices()
!
      if (wf%system%mm%forcefield .eq. 'non-polarizable') then
!
         call wf%construct_ao_electrostatics(0, 0, 'fock', wf%system%mm%n_atoms, &
                                             wf%system%mm%coordinates,           &
                                             elec_fock = wf%nopol_h_wx,          &
                                             charges = wf%system%mm%charge)
!         
         call output%print_matrix('debug', 'Electrostatic Embedding h:', &
                                  wf%nopol_h_wx, wf%n_ao, wf%n_ao)
      endif
!
   end subroutine prepare_qmmm_hf
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
!!       - MM fock (QM/MM)
!!
!!       - PCM fock
!!
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: mo_pcm_fock
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
!     MM
!
      if (wf%system%mm_calculation) then
!
         call mem%alloc(mo_mm_fock, wf%n_mo, wf%n_mo)
!
         if (wf%system%mm%forcefield .eq. 'non-polarizable') then
!
            call wf%mo_transform(wf%nopol_h_wx, mo_mm_fock)
!
         elseif (wf%system%mm%forcefield .eq. 'fq') then
!
            call wf%mo_transform(wf%pol_emb_fock, mo_mm_fock)
!
         endif
!
         call daxpy(wf%n_mo**2, half, mo_mm_fock, 1, wf%mo_fock_frozen, 1)
!
         call mem%dealloc(mo_mm_fock, wf%n_mo, wf%n_mo)
!
      endif
!
!     PCM
!
      if (wf%system%pcm_calculation) then
!
         call mem%alloc(mo_pcm_fock, wf%n_mo, wf%n_mo)
!
         call wf%mo_transform(wf%pcm_fock, mo_pcm_fock)
!
         call daxpy(wf%n_mo**2, half, mo_pcm_fock, 1, wf%mo_fock_frozen, 1)
!
         call mem%dealloc(mo_pcm_fock, wf%n_mo, wf%n_mo)
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
      wf%frozen_core = input%requested_keyword_in_section('core', 'frozen orbitals')
      wf%frozen_hf_mos = input%requested_keyword_in_section('hf', 'frozen orbitals')
!
      wf%plot_active_density = input%requested_keyword_in_section('plot hf active density', &
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
         do w = 1, wf%n_ao
!
            if (abs(wf%orbital_coefficients(w,p)) .gt. 1.0d-2) then
!  
!              Found the first significant contribution to MO p
!
               if (wf%orbital_coefficients(w,p) .lt. 0.0d0) then
!
!                 First contribution is negative, so we flip MO p
!
                  call dscal(wf%n_ao, -one, wf%orbital_coefficients(1,p), 1)
!
               endif
!
               exit
!
            endif
!
         enddo
!
      enddo
!
   end subroutine flip_final_orbitals_hf
!
!
end module hf_class
