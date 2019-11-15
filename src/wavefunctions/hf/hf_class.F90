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
module hf_class
!
!!
!!    Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
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
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: mo_fock
      real(dp), dimension(:,:), allocatable :: W_mo_update ! Eigenvectors for Roothan-Hall in MO basis
!
      real(dp), dimension(:,:), allocatable :: ao_overlap
      real(dp), dimension(:,:), allocatable :: cholesky_ao_overlap
      real(dp), dimension(:,:), allocatable :: pivot_matrix_ao_overlap
!
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz
      integer,  dimension(:,:), allocatable :: sp_eri_schwarz_list
!
!     Declarations for QM/MM
!
      real(dp) :: electrostatic_energy_qmmm
      real(dp) :: energy_qmmm
      real(dp) :: energy_scf_qmmm
!
      real(dp) :: electrostatic_energy_qmpcm
      real(dp) :: energy_scf_qmpcm
!      
      real(dp) :: linear_dep_threshold = 1.0D-6
!
      real(dp) :: coulomb_threshold    = 1.0D-12   ! screening threshold
      real(dp) :: exchange_threshold   = 1.0D-10   ! screening threshold
      real(dp) :: libint_epsilon       = 1.0D-20   ! ε for libint, integral precision given
                                                   ! approximately by sqrt(ε)
!
      type(sequential_file) :: restart_file
      type(sequential_file) :: orbital_information_file
!
      integer :: n_densities
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
   contains
!
!
      procedure :: read_orbital_coefficients                   => read_orbital_coefficients_hf
      procedure :: save_orbital_coefficients                   => save_orbital_coefficients_hf
      procedure :: read_orbital_energies                       => read_orbital_energies_hf
      procedure :: save_orbital_energies                       => save_orbital_energies_hf
!
      procedure :: read_for_scf_restart                        => read_for_scf_restart_hf
      procedure :: is_restart_safe                             => is_restart_safe_hf
      procedure :: write_scf_restart                           => write_scf_restart_hf
      procedure :: write_orbital_information                   => write_orbital_information_hf
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                                     => cleanup_hf
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
      procedure :: construct_ao_fock                           => construct_ao_fock_hf
!
      procedure :: construct_ao_G                              => construct_ao_G_hf
      procedure, private :: construct_ao_G_thread_terms        => construct_ao_G_thread_terms_hf
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
      procedure :: update_fock_and_energy_no_cumulative        => update_fock_and_energy_no_cumulative_hf
      procedure :: update_fock_and_energy_cumulative           => update_fock_and_energy_cumulative_hf
      procedure :: update_fock_and_energy                      => update_fock_and_energy_hf
      procedure :: update_fock_mm                              => update_fock_mm_hf
      procedure :: update_fock_pcm                             => update_fock_pcm_hf
!
!     AO Density related routines
!
      procedure :: construct_ao_density                        => construct_ao_density_hf
      procedure :: rotate_ao_density                           => rotate_ao_density_hf
      procedure :: purify_ao_density                           => purify_ao_density_hf
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
      procedure :: construct_sp_density_schwarz                => construct_sp_density_schwarz_hf
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
      procedure :: initialize_ao_fock                          => initialize_ao_fock_hf
      procedure :: initialize_mo_fock                          => initialize_mo_fock_hf
      procedure :: initialize_ao_overlap                       => initialize_ao_overlap_hf
      procedure :: initialize_pivot_matrix_ao_overlap          => initialize_pivot_matrix_ao_overlap_hf
      procedure :: initialize_cholesky_ao_overlap              => initialize_cholesky_ao_overlap_hf
!
      procedure :: destruct_ao_density                         => destruct_ao_density_hf
      procedure :: destruct_ao_fock                            => destruct_ao_fock_hf
      procedure :: destruct_mo_fock                            => destruct_mo_fock_hf
      procedure :: destruct_ao_overlap                         => destruct_ao_overlap_hf
      procedure :: destruct_pivot_matrix_ao_overlap            => destruct_pivot_matrix_ao_overlap_hf
      procedure :: destruct_cholesky_ao_overlap                => destruct_cholesky_ao_overlap_hf
!
!     Gradient and Hessian related routines
!
      procedure :: construct_projection_matrices               => construct_projection_matrices_hf
      procedure :: project_redundant_rotations                 => project_redundant_rotations_hf
!
      procedure :: construct_roothan_hall_hessian              => construct_roothan_hall_hessian_hf
      procedure :: construct_roothan_hall_gradient             => construct_roothan_hall_gradient_hf
      procedure :: get_packed_roothan_hall_gradient            => get_packed_roothan_hall_gradient_hf
      procedure :: get_max_roothan_hall_gradient               => get_max_roothan_hall_gradient_hf
!
      procedure :: construct_molecular_gradient                => construct_molecular_gradient_hf
!
!     Integral related routines
!
      procedure :: initialize_sp_eri_schwarz                   => initialize_sp_eri_schwarz_hf
      procedure :: destruct_sp_eri_schwarz                     => destruct_sp_eri_schwarz_hf
!
      procedure :: initialize_sp_eri_schwarz_list              => initialize_sp_eri_schwarz_list_hf
      procedure :: destruct_sp_eri_schwarz_list                => destruct_sp_eri_schwarz_list_hf
!
      procedure :: construct_sp_eri_schwarz                    => construct_sp_eri_schwarz_hf
      procedure :: get_n_sig_eri_sp                            => get_n_sig_eri_sp_hf
!
      procedure :: set_n_mo                                    => set_n_mo_hf
!
      procedure :: set_screening_and_precision_thresholds      => set_screening_and_precision_thresholds_hf
      procedure :: print_screening_settings                    => print_screening_settings_hf
!
      procedure :: prepare_for_roothan_hall                    => prepare_for_roothan_hall_hf
      procedure :: prepare                                     => prepare_hf
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
      procedure :: construct_mo_fock_fc_term                   => construct_mo_fock_fc_term_hf
      procedure :: construct_mo_fock_frozen_hf_term            => construct_mo_fock_frozen_hf_term_hf
!
      procedure :: initialize_orbital_coefficients_fc          => initialize_orbital_coefficients_fc_hf
      procedure :: destruct_orbital_coefficients_fc            => destruct_orbital_coefficients_fc_hf
!
      procedure :: initialize_orbital_coefficients_frozen_hf   => initialize_orbital_coefficients_frozen_hf_hf
      procedure :: destruct_orbital_coefficients_frozen_hf     => destruct_orbital_coefficients_frozen_hf_hf
!
!     MO scf routines
!
      procedure :: get_max_roothan_hall_mo_gradient         => get_max_roothan_hall_mo_gradient_hf
!
      procedure :: update_fock_and_energy_mo                => update_fock_and_energy_mo_hf
      procedure :: get_roothan_hall_mo_gradient             => get_roothan_hall_mo_gradient_hf
!
      procedure :: get_mo_fock                              => get_mo_fock_hf
      procedure :: set_mo_fock                              => set_mo_fock_hf
!
      procedure :: do_roothan_hall_mo                       => do_roothan_hall_mo_hf
!
      procedure :: initialize_W_mo_update                   => initialize_W_mo_update_hf
      procedure :: destruct_W_mo_update                     => destruct_W_mo_update_hf
!
      procedure :: roothan_hall_update_orbitals_mo          => roothan_hall_update_orbitals_mo_hf
      procedure :: prepare_for_roothan_hall_mo              => prepare_for_roothan_hall_mo_hf
!
   end type hf
!
   interface
!
      include "frozen_orbital_hf_interface.F90"
      include "mo_hf_interface.F90"
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
      wf%name_ = 'hf'
!
      wf%system => system
!
      call wf%read_settings()
!
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
!     Allocate active mo specific arrays
!     and construct them
!
      call wf%initialize_W_mo_update()
      call wf%initialize_mo_fock()
!
      call identity_array(wf%W_mo_update, wf%n_mo)
!
   end subroutine read_for_scf_restart_hf
!
!
   subroutine is_restart_safe_hf(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(hf) :: wf
!
      character(len=*), intent(in) :: task
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
         call output%error_msg('attempted to restart HF with an inconsistent number ' // &
                               'of atomic orbitals for task ' // trim(task))
      endif
!
      if (n_densities .ne. wf%n_densities) then
         call output%error_msg('attempted to restart HF with an inconsistent number ' // &
                               'of atomic densities (likely a HF/UHF inconsistency) for task ' // &
                               trim(task))
      endif
!
      if (n_electrons .ne. wf%system%get_n_electrons()) then
         call output%error_msg('attempted to restart HF with an inconsistent number ' // &
                               'of electrons for task ' // trim(task))
      endif
!
   end subroutine is_restart_safe_hf
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
         call output%printf('HOMO-LUMO gap:             (f19.12)', pl='m', fs='(/t6,a)', reals=[homo_lumo_gap])
!
      endif
!
      nuclear_repulsion = wf%system%get_total_nuclear_repulsion()
!
      call output%printf('Nuclear repulsion energy:  (f19.12)', pl='m', fs='(t6,a)',  reals=[nuclear_repulsion])
      call output%printf('Electronic energy:         (f19.12)', pl='m', fs='(t6,a)',  reals=[wf%energy - nuclear_repulsion])
      call output%printf('Total energy:              (f19.12)', pl='m', fs='(t6,a)',  reals=[wf%energy])
!
      if(wf%system%mm_calculation) call wf%print_energy_mm()
      if(wf%system%pcm_calculation) call wf%print_energy_pcm()
!
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
      call output%printf('- Summary of QM/MM energetics:', pl='m', fs='(/t3,a)')
      call output%printf('a.u.             eV     kcal/mol', pl='m', fs='(t42,a)')
      call output%printf('QM/MM SCF Contribution: (f22.12)', pl='m', &
                          reals=[wf%energy_scf_qmmm], fs='(t6,a)')
      call output%printf('QM/MM Electrostatic Energy:(f19.12)(f12.5) (f9.3)', pl='m', &
                          reals=[wf%electrostatic_energy_qmmm,               &
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
      call input%get_keyword_in_section('coulomb threshold', 'solver scf', wf%coulomb_threshold)
      call input%get_keyword_in_section('exchange threshold', 'solver scf', wf%exchange_threshold)
      call input%get_keyword_in_section('integral precision', 'solver scf', wf%libint_epsilon)
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
      call output%print_vector('normal', '- Molecular orbital energies', wf%n_ao, wf%orbital_energies, &
                              fs='(f16.12)', columns=4)
!
   end subroutine print_orbital_energies_hf
!
!
   subroutine print_orbitals_hf(wf, orbital_list)
!!
!!    Print orbitals
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!
!!    Prints the orbitals with atom & orbital information given.
!!
!!    orbital_list: A list of integers determining which MO to print
!!                     [1,3,4] -> print MO 1, 3 and 4
!!
      implicit none
!
      class(hf),             intent(in)           :: wf
      integer, dimension(:), intent(in), optional :: orbital_list
!
      integer, dimension(:), allocatable :: orbital_list_local
!
      integer :: n_orbitals
      integer :: mo, first_mo, last_mo
!
!     Set orbital_list
!
      if (.not. present(orbital_list)) then
!
!        Default: orbital_list = [n_o - 9, ..., n_o + 10]
!
         first_mo = max(1, wf%n_o - 9)
         last_mo = min(wf%n_mo, wf%n_o + 10)
!
         n_orbitals = last_mo - first_mo + 1
!
         call mem%alloc(orbital_list_local, n_orbitals)
!
         do mo = 1, n_orbitals
!
            orbital_list_local(mo) = first_mo + mo - 1
!
         enddo
!
      else
!
!        Input: orbital_list
!
         n_orbitals = size(orbital_list)
!
         call mem%alloc(orbital_list_local, n_orbitals)
!
         orbital_list_local = orbital_list
!
      endif
!
!     Print MO and coefficients
!
      call output%printf('- Printing molecular orbitals ((i0) from total (i0))', pl='normal', fs='(/t3,a)', &
                        ints=[n_orbitals, wf%n_mo])
!
      call wf%print_orbitals_from_coefficients(orbital_list_local, wf%orbital_coefficients)
!
      call mem%dealloc(orbital_list_local, n_orbitals)
!
   end subroutine print_orbitals_hf
!
!
   subroutine print_orbitals_from_coefficients_hf(wf, orbital_list, orbital_coefficients)
!!
!!    Print orbitals from coefficients
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!
!!    Prints the orbitals from coefficients with atom & orbital information given.
!!
      implicit none
!
      class(hf),                             intent(in) :: wf
      integer, dimension(:),                 intent(in) :: orbital_list
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(in) :: orbital_coefficients 
!
      integer, parameter :: n_entries  = 5
!
      character(len=1), dimension(6), parameter :: angular_momentums = ['s', 'p', 'd', 'f', 'g', 'h']
!
      integer :: n_orbitals, mo_offset, first_mo, last_mo, mo
      integer :: atom, shell, ao, l
!
      logical :: adv
!
      character(len=2) :: symbol, ang_mom
      real(dp) :: orb_coeff
!
!     Sanity check
!
      if (any(orbital_list < 1) .or. any(orbital_list > wf%n_mo)) then
!
         call output%error_msg('Tried to print non-existent orbital')
!
      endif
!
      n_orbitals = size(orbital_list)
!
!     Print `n_entries` columns with `wf%n_mo` rows
!
      do mo_offset = 1, n_orbitals, n_entries
!
         first_mo = mo_offset
         last_mo  = min(mo_offset + n_entries - 1, n_orbitals)
!
!        Print header: AO Atom 1 2 3 4 ..
!
         call output%printf('  AO    Atom', pl='normal', fs='(/t3,a)', adv=.false.)
!
         do mo = first_mo, last_mo
!
            adv = (mo == last_mo)
!
            call output%printf('(i4)', pl='normal', fs='(8x,a)', adv=adv, &
                              ints=[orbital_list(mo)] )
!
         enddo
!
         call output%print_separator(pl='normal', n=76, symbol='-')
!
!        Print content: AO, Atom and orbital coefficients
!
         do atom = 1, wf%system%n_atoms
!
            do shell = 1, wf%system%atoms(atom)%n_shells
!
               l = wf%system%atoms(atom)%shells(shell)%l
!
               do ao = wf%system%atoms(atom)%shells(shell)%first, wf%system%atoms(atom)%shells(shell)%last
!
                  symbol  = trim( wf%system%atoms(atom)%symbol )
                  ang_mom = trim( angular_momentums(l+1) )
!
                  call output%printf('(i4) (i4) (a2)((a1))', pl='normal', adv=.false., &
                                    ints=[ao, atom], chars=[symbol, ang_mom])
!
                  do mo = first_mo, last_mo
!
                     adv = (mo == last_mo)
                     orb_coeff = orbital_coefficients(ao, orbital_list(mo))
!
                     call output%printf('(f10.6)', pl='normal', fs='(2x,a)', adv=adv, &
                                       reals=[orb_coeff] )
!
                  enddo
!
               enddo
!
            enddo
!
         enddo
!
         call output%print_separator(pl='normal', n=76, symbol='-')
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
   subroutine initialize_orbitals_hf(wf)
!!
!!    Initialize orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the arrays associated with the orbital
!!    coefficients. In spin-unrestricted hfs, this
!!    will include alpha and beta coefficients, though these
!!    are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      wf%orbital_coefficients = zero
      wf%orbital_energies     = zero
!
   end subroutine initialize_orbitals_hf
!
!
   subroutine initialize_density_hf(wf)
!!
!!    Initialize density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO density (or densities).
!!    In spin-unrestricted hfs, this alpha and beta densities,
!!    though these are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_ao_density()
!
      wf%ao_density = zero
!
   end subroutine initialize_density_hf
!
!
   subroutine initialize_fock_hf(wf)
!!
!!    Initialize Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO Fock matrix (or matrices).
!!    In spin-unrestricted hfs, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_ao_fock()
!
      wf%ao_fock = zero
!
   end subroutine initialize_fock_hf
!
!
   subroutine destruct_fock_hf(wf)
!!
!!    Destruct Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO Fock matrix (or matrices).
!!    In spin-unrestricted hfs, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%destruct_ao_fock()
!
   end subroutine destruct_fock_hf
!
!
   subroutine update_fock_and_energy_no_cumulative_hf(wf, h_wx)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified by Tommaso Giovannini, May 2019 for QM/MM
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: h_wx
!
      real(dp), dimension(:, :), allocatable            :: h_wx_eff
!
      call mem%alloc(h_wx_eff, wf%n_ao,wf%n_ao)
!
      h_wx_eff = h_wx
!
      if(wf%system%mm_calculation .and. wf%system%mm%forcefield .eq. 'non-polarizable') &
         call wf%update_h_wx_mm(h_wx_eff)
!
      call wf%construct_ao_fock(wf%ao_density, wf%ao_fock, h_wx_eff)
!
      call mem%dealloc(h_wx_eff, wf%n_ao,wf%n_ao)
!
      if(wf%system%mm_calculation .and. wf%system%mm%forcefield .ne. 'non-polarizable') &
         call wf%update_fock_mm()
!         
      if(wf%system%pcm_calculation) call wf%update_fock_pcm()
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, h_wx)
!
      if(wf%system%mm_calculation) call wf%calculate_mm_energy_terms()
      if(wf%system%pcm_calculation) call wf%calculate_pcm_energy_terms()
!
   end subroutine update_fock_and_energy_no_cumulative_hf
!
!
   subroutine update_fock_and_energy_cumulative_hf(wf, prev_ao_density, h_wx)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed. By cumulatively
!!    we mean using the density change to build the Fock matrix
!!    in the iterative loop.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in) :: prev_ao_density
!
      logical :: cumulative
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
!
      cumulative = .true.
      call wf%construct_ao_fock(wf%ao_density, wf%ao_fock, h_wx, cumulative)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, h_wx)
!
      if(wf%system%mm_calculation) call wf%calculate_mm_energy_terms()
!
   end subroutine update_fock_and_energy_cumulative_hf
!
!
   subroutine update_fock_and_energy_hf(wf, h_wx, prev_ao_density)
!!
!!    Wrapper for Update Fock and energy
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Call either cumulative or no_cumulative updating depending on options
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
      if (.not.present(prev_ao_density)) then
!
          call wf%update_fock_and_energy_no_cumulative(h_wx)
!
      else
!
         if(.not.wf%system%mm_calculation.and..not.wf%system%pcm_calculation) then
!
            call wf%update_fock_and_energy_cumulative(prev_ao_density, h_wx)
!
         else
!
            if(wf%system%mm%forcefield.eq.'non-polarizable') then
!
               call wf%update_fock_and_energy_cumulative(prev_ao_density, h_wx)
!
            else
!
               call wf%update_fock_and_energy_no_cumulative(h_wx)
!
            endif
!
         endif
!
      endif
!
!
   end subroutine update_fock_and_energy_hf
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
   subroutine save_ao_density_hf(wf)
!!
!!    Save AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Save the AO density based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(hf) :: wf
!
      type(sequential_file) :: ao_density_file
!
      ao_density_file = sequential_file('ao_density')
      call ao_density_file%open_('write', 'rewind')
!
      call ao_density_file%write_(wf%ao_density, wf%n_ao*wf%n_ao)
!
      call ao_density_file%close_
!
   end subroutine save_ao_density_hf
!
!
   subroutine get_ao_density_sq_hf(wf, D)
!!
!!    Get AO density squared
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Returns the unpacked AO density matrix D
!!    (or density matrices in descendants, see overwriting routines)
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao**2,wf%n_densities), intent(inout) :: D
!
      call dcopy(wf%n_ao**2, wf%ao_density, 1, D, 1)
!
   end subroutine get_ao_density_sq_hf
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
      call wf%prepare_mos()
      call wf%prepare_frozen_fock_terms()
!
!     Save orbital information in orbital_information_file for CC
!
      call wf%write_orbital_information()
!
!     Save MM and PCM matrices for CC
!
      if (wf%system%mm_calculation) call wf%write_mm_matrices()
      if (wf%system%pcm_calculation) call wf%write_pcm_matrices()
!
!     Deallocations
!
      call wf%destruct_orbital_energies()
      call wf%destruct_orbital_coefficients()
      call wf%destruct_ao_overlap()
      call wf%destruct_fock()
      call wf%destruct_ao_density()
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
      call wf%destruct_sp_eri_schwarz()
      call wf%destruct_sp_eri_schwarz_list()
!
      call wf%destruct_W_mo_update()
      call wf%destruct_mo_fock()
!
      call wf%destruct_mm_matrices()
      call wf%destruct_pcm_matrices()
!
   end subroutine cleanup_hf
!
!
   subroutine initialize_ao_density_hf(wf)
!!
!!    Initialize AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_density)) call mem%alloc(wf%ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_hf
!
!
   subroutine initialize_ao_fock_hf(wf)
!!
!!    Initialize AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_fock)) call mem%alloc(wf%ao_fock, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_fock_hf
!
!
   subroutine initialize_mo_fock_hf(wf)
!!
!!    Initialize MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%mo_fock)) call mem%alloc(wf%mo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_mo_fock_hf
!
!
   subroutine initialize_ao_overlap_hf(wf)
!!
!!    Initialize AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_overlap)) call mem%alloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_overlap_hf
!
!
   subroutine initialize_pivot_matrix_ao_overlap_hf(wf)
!!
!!    Initialize pivot matrix AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%pivot_matrix_ao_overlap)) call mem%alloc(wf%pivot_matrix_ao_overlap, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_pivot_matrix_ao_overlap_hf
!
!
   subroutine initialize_cholesky_ao_overlap_hf(wf)
!!
!!    Initialize cholesky vectors AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%cholesky_ao_overlap)) call mem%alloc(wf%cholesky_ao_overlap, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_cholesky_ao_overlap_hf
!
!
   subroutine initialize_sp_eri_schwarz_hf(wf)
!!
!!    Initialize shell pair eri schwarz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%sp_eri_schwarz)) &
         call mem%alloc(wf%sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2, 2)
!
   end subroutine initialize_sp_eri_schwarz_hf
!
!
   subroutine destruct_sp_eri_schwarz_hf(wf)
!!
!!    Destruct shell pair eri schwarz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%sp_eri_schwarz)) &
         call mem%dealloc(wf%sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2, 2)
!
   end subroutine destruct_sp_eri_schwarz_hf
!
!
   subroutine initialize_sp_eri_schwarz_list_hf(wf)
!!
!!    Initialize shell pair eri schwarz list
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%sp_eri_schwarz_list)) &
         call mem%alloc(wf%sp_eri_schwarz_list,wf%system%n_s*(wf%system%n_s + 1)/2, 3)
!
   end subroutine initialize_sp_eri_schwarz_list_hf
!
!
   subroutine destruct_sp_eri_schwarz_list_hf(wf)
!!
!!    Destruct shell pair eri schwarz list
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%sp_eri_schwarz_list)) &
         call mem%dealloc(wf%sp_eri_schwarz_list, wf%system%n_s*(wf%system%n_s + 1)/2, 3)
!
   end subroutine destruct_sp_eri_schwarz_list_hf
!
!
   subroutine destruct_ao_overlap_hf(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_overlap)) call mem%dealloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_overlap_hf
!
!
   subroutine destruct_ao_density_hf(wf)
!!
!!    Destruct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_density)) call mem%dealloc(wf%ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_hf
!
!
   subroutine destruct_ao_fock_hf(wf)
!!
!!    Destruct AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_fock)) call mem%dealloc(wf%ao_fock, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_hf
!
!
   subroutine destruct_mo_fock_hf(wf)
!!
!!    Destruct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%mo_fock)) call mem%dealloc(wf%mo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_mo_fock_hf
!
!
   subroutine destruct_pivot_matrix_ao_overlap_hf(wf)
!!
!!    Destruct pivot matrix AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%pivot_matrix_ao_overlap)) call mem%dealloc(wf%pivot_matrix_ao_overlap, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_pivot_matrix_ao_overlap_hf
!
!
   subroutine destruct_cholesky_ao_overlap_hf(wf)
!!
!!    Initialize cholesky vectors AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%cholesky_ao_overlap)) call mem%dealloc(wf%cholesky_ao_overlap, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_cholesky_ao_overlap_hf
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
   subroutine construct_sp_eri_schwarz_hf(wf)
!!
!!    Construct shell-pair electronic-repulsion-integral Schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of g_wxwx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively.
!!
      implicit none
!
      class(hf) :: wf
!
      integer, dimension(:),  allocatable :: sp_eri_schwarz_index_list
      real(dp), dimension(:), allocatable :: sorted_sp_eri_schwarz
!
!     Local variables
!
      integer :: s1, s2, s1s2
!
      real(dp) :: maximum
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      type(interval) :: A_interval, B_interval
!
!     Set the maximum element in each shell pair
!
      call set_coulomb_precision_c(1.0d-50)
!
!$omp parallel do private(s1, s2, s1s2, A_interval, B_interval, g, maximum) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
            s1s2 = (max(s1,s2)*(max(s1,s2)-3)/2) + s1 + s2
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call wf%system%construct_ao_g_wxyz(g, s1, s2, s1, s2)
!
            maximum = get_abs_max(g, ((A_interval%length)*(B_interval%length))**2)
!
            wf%sp_eri_schwarz(s1s2, 1) = sqrt(maximum)
!
            wf%sp_eri_schwarz_list(s1s2, 1) = s1
            wf%sp_eri_schwarz_list(s1s2, 2) = s2
!
         enddo
      enddo
!$omp end parallel do
!
!     Sort the sp_eri_schwarz vector and use the resulting index list
!     to resort the sp_eri_schwarz_list matrix
!
      call mem%alloc(sp_eri_schwarz_index_list, wf%system%n_s*(wf%system%n_s + 1)/2)
      sp_eri_schwarz_index_list = 0
!
      call mem%alloc(sorted_sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
      sorted_sp_eri_schwarz = wf%sp_eri_schwarz(:, 1)
!
      call get_n_highest(wf%system%n_s*(wf%system%n_s + 1)/2, wf%system%n_s*(wf%system%n_s + 1)/2, wf%sp_eri_schwarz,&
                      sorted_sp_eri_schwarz, sp_eri_schwarz_index_list)
!
      wf%sp_eri_schwarz(:, 2) = wf%sp_eri_schwarz(:, 1)
      wf%sp_eri_schwarz(:, 1) = sorted_sp_eri_schwarz
      call mem%dealloc(sorted_sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
      wf%sp_eri_schwarz_list(:,3) = sp_eri_schwarz_index_list
      call mem%dealloc(sp_eri_schwarz_index_list, wf%system%n_s*(wf%system%n_s + 1)/2)
!
   end subroutine construct_sp_eri_schwarz_hf
!
!
   subroutine construct_sp_density_schwarz_hf(wf, sp_density_schwarz, D)
!!
!!    Construct shell-pair density schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of D_wx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s) :: sp_density_schwarz
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), dimension(:,:), allocatable, target :: D_red
      real(dp), dimension(:,:), contiguous, pointer :: D_red_p => null()
!
      type(interval) :: A_interval, B_interval
!
      integer :: s1, s2
      integer :: n_threads = 1, thread = 0
!
      real(dp) :: maximum
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(D_red,wf%system%max_shell_size**2,n_threads)
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D_red_p, maximum, thread) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
!$          thread = omp_get_thread_num()
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            D_red_p(1:A_interval%length,1:B_interval%length) => D_red(1:A_interval%length*B_interval%length,thread+1)
!
            D_red_p = D(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            maximum = get_abs_max(D_red_p, (A_interval%length)*(B_interval%length))
!
            nullify(D_red_p)
!
            sp_density_schwarz(s1, s2) = maximum
            sp_density_schwarz(s2, s1) = maximum
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_red,wf%system%max_shell_size**2,n_threads)
!
   end subroutine construct_sp_density_schwarz_hf
!
!
   subroutine get_n_sig_eri_sp_hf(wf, n_sig_sp)
!!
!!    Get number of significant ERI shell-pairs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the number of significant shell pairs. The threshold
!!    determines how small the largest element of g_wxwx in a shell
!!    pair AB (w in A, x in B) to be ignored completely in the Fock
!!    construction loop.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(inout) :: n_sig_sp
!
      integer :: s1s2
!
      n_sig_sp = 0
!
      do s1s2 = 1, wf%system%n_s*(wf%system%n_s + 1)/2
!
         if (wf%sp_eri_schwarz(s1s2, 1)*wf%sp_eri_schwarz(1, 1) .lt. sqrt(wf%libint_epsilon)) then
!
            exit
!
         else
!
            n_sig_sp = n_sig_sp + 1
!
         endif
!
      enddo
!
   end subroutine get_n_sig_eri_sp_hf
!
!
   subroutine construct_ao_fock_hf(wf, D, ao_fock, h_wx, cumulative)
!!
!!    Construct AO Fock matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates
!!
!!       F_αβ = h_αβ + sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where D is the AO density. This routine is integral direct, and
!!    it calculates the Hartree-Fock energy by default.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: D
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: ao_fock
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      logical, intent(in), optional :: cumulative
!
      logical :: local_cumulative
!
      real(dp), dimension(:,:), allocatable :: G
!
      type(timings) :: ao_fock_timer
!
      ao_fock_timer = timings('AO Fock construction')
      call ao_fock_timer%turn_on()
!
!     Set whether to accumulate into Fock (density differences)
!     or to construct the entire Fock matrix
!
      local_cumulative = .false.
      if (present(cumulative)) then
!
         if (cumulative) then
!
            local_cumulative = .true.
!
         else
!
            local_cumulative = .false.
!
         endif
!
      endif
!
!     Construct the two electron part of the Fock matrix (G),
!     and add the contribution to the Fock matrix
!
      call mem%alloc(G, wf%n_ao, wf%n_ao)
      call wf%construct_ao_G(D, G)
!
      if (.not. local_cumulative) call zero_array(ao_fock, wf%n_ao**2)
!
      call daxpy(wf%n_ao**2, one, G, 1, ao_fock, 1)
      call mem%dealloc(G, wf%n_ao, wf%n_ao)
!
!     Add the one-electron contribution F =+ h
!
      if (.not. local_cumulative) call daxpy(wf%n_ao**2, one, h_wx, 1, ao_fock, 1)
!
      call ao_fock_timer%turn_off()
!
   end subroutine construct_ao_fock_hf
!
!
!
   subroutine construct_ao_G_hf(wf, D, G)
!!
!!    Construct AO G(D)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Modified by Sarai D. Folkestad, Oct 2019 for G construction
!!    only
!!
!!    Calculates
!!
!!       G(D)_αβ = sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where D is the inactive AO Density.
!!
      class(hf)   :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: D
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: G
!
      integer :: thread = 0, n_threads = 1
!
      real(dp), dimension(:,:), allocatable :: sp_Density_schwarz
      real(dp), dimension(:,:,:), allocatable :: G_thread
!
      integer :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      type(timings) :: G_timer
!
      G_timer = timings('G construction time')
      call G_timer%turn_on()
!
!     Construct the density screening vector and the maximum element in the density
!
      call mem%alloc(sp_Density_schwarz, wf%system%n_s, wf%system%n_s)
!
      call wf%construct_sp_Density_schwarz(sp_Density_schwarz, D)
      max_D_schwarz = get_abs_max(sp_Density_schwarz, wf%system%n_s**2)
!
!     Compute number of significant ERI shell pairs (the G construction
!     only loops over these shell pairs) and the maximum element
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
      max_eri_schwarz = get_abs_max(wf%sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
!     Construct the two electron part of the Fock matrix, using the screening vectors
!     and parallelizing over available threads (each gets its own copy of the Fock matrix)
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(G_thread, wf%n_ao, wf%n_ao, n_threads) ! [G(thread 1) G(thread 2) ...]
      call zero_array(G_thread, (wf%n_ao**2)*n_threads)
!
      call wf%construct_ao_G_thread_terms(G_thread, D, n_threads, max_D_schwarz, max_eri_schwarz,      &
                                         sp_density_schwarz, n_sig_sp,                 &
                                         wf%coulomb_threshold, wf%exchange_threshold,  &
                                         wf%libint_epsilon, wf%system%shell_limits)
!
      call mem%dealloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix,
!     and symmetrize the result
!
      call zero_array(G, (wf%n_ao**2))
!
      do thread = 1, n_threads
!
         call daxpy(wf%n_ao**2, one, G_thread(1, 1, thread), 1, G, 1)
!
      enddo
!
      call mem%dealloc(G_thread, wf%n_ao, wf%n_ao, n_threads)
!
      call symmetric_sum(G, wf%n_ao)
      call dscal(wf%n_ao**2, half, G, 1)
!
      call G_timer%turn_off()
!
   end subroutine construct_ao_G_hf
!
!
   subroutine construct_ao_G_thread_terms_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,    &
                                          sp_density_schwarz, n_sig_sp, coulomb_thr, &
                                          exchange_thr, precision_thr, shells)
!!
!!    Construct AO G
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ =+ sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ (= G(D)_αβ),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get G(D)_αβ.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads, n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                             &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34,     &
!$omp w, x, y, z, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset,         &
!$omp temp4, temp5, temp6, temp7, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,        &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                                &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4), &
                           sp_density_schwarz_s1s2)
!
               temp8 = max(sp_density_schwarz_s3s2, &
                           sp_density_schwarz_s3s1, &
                           sp_density_schwarz(s4,s2), &
                           sp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr .and. temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4,         &
                  precision_thr/max(temp7,temp8), thread, skip, shells(s1)%length, shells(s2)%length, &
                  shells(s3)%length, shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*(shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*(shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           temp3 = eighth*temp*d3
                           temp4 = eighth*temp*d4
                           temp5 = eighth*temp*d5
                           temp6 = eighth*temp*d6
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
!
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
!
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call set_coulomb_precision_c(wf%libint_epsilon)
!
   end subroutine construct_ao_G_thread_terms_hf
!
!
   subroutine construct_coulomb_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,     &
                                                   sp_density_schwarz, &
                                                   n_sig_sp, coulomb_thr, precision_thr, shells)
!!
!!    AO Fock Coulomb construction loop
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the Coulomb two-electron part of the Fock matrix,
!!
!!       F_αβ = F_αβ + sum_γδ g_αβγδ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get the Coulomb part of G(D)_αβ.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads,  n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: d1, d2, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp7, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp1, temp2, d1, d2, thread, thread_offset,                            &
!$omp temp7, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,                                &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                            &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4), sp_density_schwarz_s1s2)
!
               if (temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp7, thread, skip, shells(s1)%length, shells(s2)%length,    &
                  shells(s3)%length, shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*(shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*(shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_coulomb_ao_G_hf
!
!
   subroutine construct_exchange_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz, &
                                          sp_density_schwarz, &
                                           n_sig_sp, exchange_thr, precision_thr, shells)
!!
!!    AO Fock exchange construction loop
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ = F_αβ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get the exchange part of G(D)_αβ.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads,  n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp3, temp4, temp5, temp6, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp3, d3, d4, d5, d6, thread, thread_offset,                           &
!$omp temp4, temp5, temp6, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,           &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                            &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. exchange_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. exchange_thr) cycle ! Screened out shell pair
!
               temp8 = max(sp_density_schwarz_s3s2,   &
                           sp_density_schwarz_s3s1,   &
                           sp_density_schwarz(s4,s2), &
                           sp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp8, thread, skip, shells(s1)%length, shells(s2)%length,    &
                  shells(s3)%length, shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*(shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*(shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp3 = eighth*temp*d3
                           temp4 = eighth*temp*d4
                           temp5 = eighth*temp*d5
                           temp6 = eighth*temp*d6
!
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_exchange_ao_G_hf
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
   subroutine construct_ao_G_1der_hf(wf, G_ao, D_ao)
!!
!!    Construct AO G 1der
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβqk = sum_γδ g_αβγδqk D_γδ - 1/2 * sum_γδ g_αδγβqk D_γδ (= G(D)_αβqk),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: G_ao
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_ao
!
      real(dp), dimension((wf%system%max_shell_size**4)*3*4), target :: g_ABCDqk
      real(dp), dimension(:,:,:,:,:,:), pointer, contiguous :: g_ABCDqk_p
!
      integer :: A, B, C, D, D_max, w, x, y, z, n_sig_sp, AB, AB_packed
      integer :: w_red, x_red, y_red, z_red, tot_dim, k, q, n_threads, thread
!
      real(dp) :: d1, d2, d3, d4, d5, d6
!
      integer, dimension(4) :: atoms
!
      real(dp) :: deg, deg_AB, deg_CD, deg_AB_CD
!
      real(dp), dimension(3,4) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
!
      real(dp), dimension(:,:,:,:,:), allocatable :: G_ao_t
!
!$    n_threads = omp_get_max_threads()
      call mem%alloc(G_ao_t, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms, n_threads)
      G_ao_t = zero
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
!
!$omp parallel do private(A, B, C, D, D_max, atoms, deg, deg_CD, deg_AB, deg_AB_CD, g_ABCDqk, g_ABCDqk_p, &
!$omp w, x, y, z, w_red, x_red, y_red, z_red, temp, temp1, temp2, temp3, temp4, temp5, temp6, &
!$omp d1, d2, d3, d4, d5, d6, thread, q, k, tot_dim, AB, AB_packed) schedule(dynamic)
      do AB = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
!
         if (wf%sp_eri_schwarz(AB, 1)*wf%sp_eri_schwarz(1, 1) < wf%coulomb_threshold) cycle
         AB_packed = wf%sp_eri_schwarz_list(AB, 3)
!
         A = wf%sp_eri_schwarz_list(AB_packed, 1)
         B = wf%sp_eri_schwarz_list(AB_packed, 2)
!
         atoms(1) = wf%system%shell2atom(A)
         atoms(2) = wf%system%shell2atom(B)
!
         deg_AB = real(2-B/A, kind=dp)
!
            do C = 1, A
!
               D_max = (C/A)*B + (1-C/A)*C
               atoms(3) = wf%system%shell2atom(C)
!
               do D = 1, D_max
!
                  deg_CD    = real(2-D/C, kind=dp)
                  deg_AB_CD = min(1-C/A+2-min(D/B,B/D), 2)

                  deg = deg_AB*deg_CD*deg_AB_CD ! Shell degeneracy
!
                  atoms(4) = wf%system%shell2atom(D)
!
                  call wf%system%construct_ao_g_wxyz_1der(g_ABCDqk, A, B, C, D)
!
                  tot_dim = (wf%system%shell_limits(A)%length)*(wf%system%shell_limits(B)%length)&
                              *(wf%system%shell_limits(C)%length)*(wf%system%shell_limits(D)%length)&
                              *3*4
!
                  g_ABCDqk(1:tot_dim) = deg*g_ABCDqk(1:tot_dim)
!
                  g_ABCDqk_p(1 : wf%system%shell_limits(A)%length, 1 : wf%system%shell_limits(B)%length, &
                             1 : wf%system%shell_limits(C)%length, 1 : wf%system%shell_limits(D)%length, &
                             1 : 3, 1 : 4) => g_ABCDqk(1 : tot_dim)
!
                  do z = wf%system%shell_limits(D)%first, wf%system%shell_limits(D)%last
!
                     z_red = z - wf%system%shell_limits(D)%first + 1
!
                     do y = wf%system%shell_limits(C)%first, wf%system%shell_limits(C)%last
!
                        y_red = y - wf%system%shell_limits(C)%first + 1
!
                        d1 = D_ao(y, z)
!
                        do x = wf%system%shell_limits(B)%first, wf%system%shell_limits(B)%last
!
                           x_red = x - wf%system%shell_limits(B)%first + 1
!
                           d3 = D_ao(x, y)
                           d5 = D_ao(x, z)
!
                           do w = wf%system%shell_limits(A)%first, wf%system%shell_limits(A)%last
!
                              d2 = D_ao(w, x)
                              d4 = D_ao(w, y)
                              d6 = D_ao(w, z)
!
                              w_red = w - wf%system%shell_limits(A)%first + 1
!
                              do k = 1, 4
                                 do q = 1, 3
                                    temp(q, k) = g_ABCDqk_p(w_red, x_red, y_red, z_red, q, k)
                                 enddo
                              enddo
!
                              temp1 = half*temp*d1
                              temp2 = half*temp*d2
!
                              temp3 = eighth*temp*d3
                              temp4 = eighth*temp*d4
                              temp5 = eighth*temp*d5
                              temp6 = eighth*temp*d6
!
                              do k = 1, 4
                                 do q = 1, 3
!
                                    G_ao_t(w, x, q, atoms(k), thread+1) = G_ao_t(w, x, q, atoms(k), thread+1) + temp1(q, k)
                                    G_ao_t(y, x, q, atoms(k), thread+1) = G_ao_t(y, x, q, atoms(k), thread+1) - temp6(q, k)
                                    G_ao_t(y, z, q, atoms(k), thread+1) = G_ao_t(y, z, q, atoms(k), thread+1) + temp2(q, k)
                                    G_ao_t(w, z, q, atoms(k), thread+1) = G_ao_t(w, z, q, atoms(k), thread+1) - temp3(q, k)
                                    G_ao_t(x, z, q, atoms(k), thread+1) = G_ao_t(x, z, q, atoms(k), thread+1) - temp4(q, k)
                                    G_ao_t(w, y, q, atoms(k), thread+1) = G_ao_t(w, y, q, atoms(k), thread+1) - temp5(q, k)
!
                                 enddo
                              enddo
!
                           enddo
                        enddo
                     enddo
                  enddo
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
      G_ao = zero
!
      do thread = 1, n_threads
!
         call daxpy(3*wf%system%n_atoms*wf%n_ao**2, one, G_ao_t(1,1,1,1,thread), 1, G_ao, 1)
!
      enddo
!
      call mem%dealloc(G_ao_t, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms, n_threads)
!
   end subroutine construct_ao_G_1der_hf
!
!
   subroutine set_ao_density_hf(wf, D)
!!
!!    Set AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO density from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:) :: D ! Packed
!
      call squareup(D, wf%ao_density, wf%n_ao)
!
   end subroutine set_ao_density_hf
!
!
   subroutine set_ao_fock_hf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
      call squareup(F(:,1), wf%ao_fock, wf%n_ao)
!
   end subroutine set_ao_fock_hf
!
!
   subroutine get_ao_fock_hf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock from input
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
      call packin(F(:,1), wf%ao_fock, wf%n_ao)
!
   end subroutine get_ao_fock_hf
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
   subroutine construct_mo_fock_hf(wf)
!!
!!    Construct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the MO Fock matrix F_pq using the current AO
!!    Fock and the orbital coefficients.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%ao_fock,                &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  X,                         & ! X = F^ao C
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  wf%mo_fock,                & ! F = C^T F^ao C
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_ao, wf%n_mo)
!
   end subroutine construct_mo_fock_hf
!
!
   subroutine get_fock_ov_hf(wf, F)
!!
!!    Get HF equations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the occupied-virtual block of the Fock MO matrix,
!!    and returns the result in the array F.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v)   :: F ! F_ia
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_v)
!
      call dgemm('N', 'N',                                  &
                  wf%n_ao,                                  &
                  wf%n_v,                                   &
                  wf%n_ao,                                  &
                  one,                                      &
                  wf%ao_fock,                               &
                  wf%n_ao,                                  &
                  wf%orbital_coefficients(1, wf%n_o + 1),   &
                  wf%n_ao,                                  &
                  zero,                                     &
                  X,                                        &
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_o,                    &
                  wf%n_v,                    &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  F,                         &
                  wf%n_o)
!
      call mem%dealloc(X, wf%n_ao, wf%n_v)
!
   end subroutine get_fock_ov_hf
!
!
   subroutine get_ao_density_hf(wf, D)
!!
!!    Get AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Packs the AO density into D.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:) :: D
!
      call packin(D(:,1), wf%ao_density, wf%n_ao)
!
   end subroutine get_ao_density_hf
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
      wf%ao_density = half*wf%ao_density
      call full_cholesky_decomposition_system(wf%ao_density, wf%orbital_coefficients, &
                                              wf%n_ao, rank, 1.0d-12, used_diag)
      wf%ao_density = two*wf%ao_density
!
!     Make permutation matrix P
!
      call mem%alloc(perm_matrix, wf%n_ao, wf%n_ao)
!
      if (any(used_diag .gt. wf%n_ao) .or. any(used_diag .le. 0)) then
!
         call output%printf('Something went wrong when decomposing the AO density.', &
                             pl='m', fs='(/t3,a)')
!
         call output%error_msg('Trying to access elements outside of an array.')
!
      end if
!
      perm_matrix = zero
!
      do j = 1, wf%n_ao
!
         perm_matrix(used_diag(j), j) = one
!
      enddo
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
      integer :: j
!
      call mem%alloc(used_diag, wf%n_ao)
      used_diag = 0
!
      call mem%alloc(L, wf%n_ao, wf%n_ao) ! Full Cholesky vector
      L = zero
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, wf%n_mo, &
                                              wf%linear_dep_threshold, used_diag)
!
      call wf%initialize_cholesky_ao_overlap()
!
      wf%cholesky_ao_overlap(:,:) = L(1:wf%n_mo, 1:wf%n_mo)
!
      call mem%dealloc(L, wf%n_ao, wf%n_ao)
!
!     Make permutation matrix P
!
      call wf%initialize_pivot_matrix_ao_overlap()
!
      wf%pivot_matrix_ao_overlap = zero
!
      if (wf%n_mo .gt. wf%n_ao .or. wf%n_mo .le. 0 .or. &
          any(used_diag .gt. wf%n_ao) .or. any(used_diag .le. 0)) then
!
         call output%printf('Something went wrong when decomposing the AO overlap.', &
                            pl='m', fs='(/t3,a)')
!
         call output%printf('Did you compile with wrong type of integers in setup? &
                            &For example system native BLAS  &
                            &with default 64-bit integers.', &
                            pl='m', ffs='(/t3,a)')
!
         call output%printf('If that is the case, use setup with --int32 or install MKL.', &
                            pl='m')
!
         call output%error_msg('Failed to decompose AO overlap.')
!
      end if
!
      do j = 1, wf%n_mo
!
         wf%pivot_matrix_ao_overlap(used_diag(j), j) = one
!
      enddo
!
      call mem%dealloc(used_diag, wf%n_ao)
!
   end subroutine decompose_ao_overlap_hf
!
!
   subroutine rotate_ao_density_hf(wf, X)
!!
!!    Rotate AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs an update of the AO density according to a first-order
!!    truncation of the BCH expansion:
!!
!!       D^AO <- exp(-X S) D^AO exp(S X) ~ D^AO + [D^AO, X]_S = D_AO + C.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: X
!
      real(dp), dimension(:,:), allocatable :: M
      real(dp), dimension(:,:), allocatable :: C
!
!     Construct C = [D^AO, X]_S = D^AO S X - X S D^AO
!
      call mem%alloc(M, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  X,             &
                  wf%n_ao,       &
                  zero,          &
                  M,             & ! M = S X
                  wf%n_ao)
!
      call mem%alloc(C, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  M,             &
                  wf%n_ao,       &
                  zero,          &
                  C,             & ! C = D^AO M = D^AO S X
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  wf%ao_density, &
                  wf%n_ao,       &
                  zero,          &
                  M,             & ! M = S D^AO
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  -one,          &
                  X,             &
                  wf%n_ao,       &
                  M,             &
                  wf%n_ao,       &
                  one,           &
                  C,             & ! C = C - X M = C - X S D^AO = D^AO S X - X S D^AO
                  wf%n_ao)
!
      wf%ao_density = wf%ao_density + C
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
      call mem%dealloc(C, wf%n_ao, wf%n_ao)
!
   end subroutine rotate_ao_density_hf
!
!
   subroutine purify_ao_density_hf(wf, threshold)
!!
!!    Purify AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Purifies a non-idempotent AO density matrix - typically arising from
!!    the non-exact rotation of the density - by the following fixed point algorithm:
!!
!!       D^AO <- 3/4 D^AO S D^AO - 1/2 D^AO S D^AO S D^AO
!!
!!    To check whether the final result is consistent, it is possible to verify that
!!    1/2 D^AO S is idempotent (which was done during debug of routine).
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(:,:), allocatable :: M ! Arrays to temporarily hold matrix products
      real(dp), dimension(:,:), allocatable :: N ! Arrays to temporarily hold matrix products
!
      real(dp), dimension(:,:), allocatable :: prev_ao_density ! Holds previous density matrix
!
      real(dp) :: ddot, error
!
      logical :: pure = .false.
!
      integer :: iteration
      integer, parameter :: max_iterations = 50
!
      iteration = 1
!
      call mem%alloc(M, wf%n_ao, wf%n_ao)
      call mem%alloc(N, wf%n_ao, wf%n_ao)
      call mem%alloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
      pure = .false.
      error = zero
!
      do while (.not. pure .and. iteration .le. max_iterations)
!
         prev_ao_density = wf%ao_density
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_overlap, &
                     wf%n_ao,       &
                     wf%ao_density, &
                     wf%n_ao,       &
                     zero,          &
                     M,             & ! M = S D^AO
                     wf%n_ao)
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_density, &
                     wf%n_ao,       &
                     M,             &
                     wf%n_ao,       &
                     zero,          &
                     N,             & ! N = D^AO M = D^AO S D^AO
                     wf%n_ao)
!
         wf%ao_density = (three/two)*N ! D^AO = 3 D^AO S D^AO
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_overlap, &
                     wf%n_ao,       &
                     N,             &
                     wf%n_ao,       &
                     zero,          &
                     M,             & ! M = S N = S D^AO S D^AO
                     wf%n_ao)
!
         call dgemm('N','N',          &
                     wf%n_ao,         &
                     wf%n_ao,         &
                     wf%n_ao,         &
                     one,             &
                     prev_ao_density, &
                     wf%n_ao,         &
                     M,               &
                     wf%n_ao,         &
                     zero,            &
                     N,               & ! N = D^AO M = D^AO S D^AO S D^AO
                     wf%n_ao)
!
         wf%ao_density = wf%ao_density - (one/two)*N
!
         M = wf%ao_density - prev_ao_density
!
         error = sqrt(ddot((wf%n_ao)**2, M, 1, M, 1))
!
         if (error .lt. threshold) then
!
            pure = .true.
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      if (.not. pure) then
!
         call output%printf('Error: could not purify AO density. Final error: ', &
                             pl='m', reals=[error], fs='(/t3,a)')
         stop
!
      endif
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
      call mem%dealloc(N, wf%n_ao, wf%n_ao)
      call mem%dealloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine purify_ao_density_hf
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
      Pv = zero
!
      do x = 1, wf%n_ao
!
         Pv(x, x) = one
!
      enddo
!
!     Pv = I - Po
!
      do x = 1, wf%n_ao
         do w = 1, wf%n_ao
!
            Pv(w, x) = Pv(w, x) - Po(w, x)
!
         enddo
      enddo
!
   end subroutine construct_projection_matrices_hf
!
!
   subroutine project_redundant_rotations_hf(wf, X, Po, Pv)
!!
!!    Project redundant rotations
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Here, X is an antisymmetric rotations matrix on entering the routine,
!!    where some parameters are redundant for rotating the AO density. On exit,
!!    the redundant parameters have been projected out of X.
!!
!!    To achieve this, we set
!!
!!       X <- Po X Pv^T + Pv X Po^T,
!!
!!    where Po = 1/2 D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: X
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct
!
!        tmp = X Pv^T   =>   tmp^T = Pv X^T = - Pv X
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'T', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_ao, &
                  Pv,      &
                  wf%n_ao, &
                  zero,    &
                  tmp,     &
                  wf%n_ao)
!
!     X = Po X Pv^T = Po tmp
!
      call dgemm('N', 'N', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  Po,      &
                  wf%n_ao, &
                  tmp,     &
                  wf%n_ao, &
                  zero,    &
                  X,       &
                  wf%n_ao)
!
!     X = X + Pv X Po^T = X - (-Pv X) Po^T = X - tmp^T Po^T
!
      call dgemm('T', 'T', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  -one,    &
                  tmp,     &
                  wf%n_ao, &
                  Po,      &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_ao)
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine project_redundant_rotations_hf
!
!
   subroutine construct_roothan_hall_hessian_hf(wf, H, Po, Pv)
!!
!!    Construct Roothan-Hall Hessian
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall Hessian,
!!
!!       H = Fvv - Foo = Pv^T F Pv - Po^T F Po,
!!
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: H
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fvv = Pv^T F Pv and set H = tmp = Fvv
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = wf%ao_fock
      call sandwich(tmp, Pv, Pv, wf%n_ao)
!
      H = tmp
!
!     Construct tmp = Foo = Po^T F Po and set H = H - tmp = Fvv - Foo
!
      tmp = wf%ao_fock
      call sandwich(tmp, Po, Po, wf%n_ao)
!
      H = H - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine construct_roothan_hall_hessian_hf
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
      real(dp), dimension(wf%n_ao*(wf%n_ao - 1)/2, wf%n_densities), intent(inout) :: G
!
      real(dp), dimension(:,:), allocatable :: G_sq, Po, Pv
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
!
      Po = zero
      Pv = zero
!
      call wf%construct_projection_matrices(Po, Pv, wf%ao_density)
!
      call mem%alloc(G_sq, wf%n_ao, wf%n_ao)
      G_sq = zero
!
      call wf%construct_roothan_hall_gradient(G_sq, Po, Pv, wf%ao_fock)
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      G = zero
      call packin_anti(G(:,1), G_sq, wf%n_ao)
!
      call mem%dealloc(G_sq, wf%n_ao, wf%n_ao)
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
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
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
      real(dp), dimension(wf%n_ao, wf%n_ao) :: G
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fov = Po^T F Pv and set G = tmp = Fov
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = F
      call sandwich(tmp, Po, Pv, wf%n_ao)
!
      G = tmp
!
!     Construct tmp = Fvo = Pv^T F Po and set H = H - tmp = Fov - Fvo
!
      tmp = F
      call sandwich(tmp, Pv, Po, wf%n_ao)
!
      G = G - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
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
      real(dp), dimension(:,:), allocatable :: prev_C
!
      real(dp) :: ddot, orbital_dotprod
!
      integer :: info, p
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
      work = zero
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
      if (info .ne. 0)  call output%error_msg('Error: could not solve Roothan-Hall equations.')
!
!     Transform back the solutions to original basis, C = P (P^T C) = P C'
!
      call mem%alloc(prev_C, wf%n_ao, wf%n_mo)
      prev_C = C
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
!     Test for orbitals that were approximately sign-flipped by dsygv,
!     resetting them if this is the case (this changes the orbitals
!     minimally, thus preserving the sweetness of whatever CC guess
!     is on file).
!
      do p = 1, wf%n_mo
!
         orbital_dotprod = ddot(wf%n_ao, prev_C(1, p), 1, C(1, p), 1)/ddot(wf%n_ao, C(1, p), 1, C(1, p), 1)
!
         if (orbital_dotprod .lt. zero) then
!
            C(:,p) = -C(:,p)
!
         endif
!
      enddo
!
      call mem%dealloc(prev_C, wf%n_ao, wf%n_mo)
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
      wf%ao_density = zero
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
         wf%ao_density(first_ao_on_atom:last_ao_on_atom, first_ao_on_atom:last_ao_on_atom) &
                                             = atomic_density(1:n_ao_on_atom, 1:n_ao_on_atom)
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
      do ao = 1, wf%n_ao
!
         n_electrons = n_electrons + DS(ao, ao)
!
      enddo
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
      wf%ao_fock = h_wx
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
      call output%printf('Coulomb screening threshold:  (e11.4)', pl='n', &
                          reals=[wf%coulomb_threshold], fs='(/t6,a)')
!
      call output%printf('Exchange screening threshold: (e11.4)', pl='n', &
                          reals=[wf%exchange_threshold], fs='(t6,a)')
!
      call output%printf('ERI integral precision:       (e11.4)', pl='n', &
                          reals=[wf%libint_epsilon], fs='(t6,a)')
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
      call output%printf('- Cholesky decomposition of AO overlap to get linearly independent orbitals:', &
                         pl='n', fs='(/t3,a)', ll=100)
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap()
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      call output%printf('Number of occupied orbitals:  (i8)', ints=[wf%n_o], fs='(/t6,a)', pl='minimal')
      call output%printf('Number of virtual orbitals:   (i8)', ints=[wf%n_v], fs='(t6,a)', pl='minimal')
      call output%printf('Number of molecular orbitals: (i8)', ints=[wf%n_mo], fs='(t6,a)', pl='minimal')
      call output%printf('Number of atomic orbitals:    (i8)', ints=[wf%n_ao], fs='(t6,a)', pl='minimal')
!
      if (wf%n_mo .lt. wf%n_ao) &
         call output%printf('Removed (i0) AOs due to linear dep.', ints=[wf%n_ao - wf%n_mo], fs='(/t6,a)', pl='minimal')
!
   end subroutine set_n_mo_hf
!
!
   subroutine set_screening_and_precision_thresholds_hf(wf, gradient_threshold)
!!
!!    Set screening and precision thresholds
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the screening thresholds for Coulomb and exchange
!!    integrals given the convergence threshold for the gradient
!!
!!       coulomb_threshold  = gradient_threshold * 1.0d-3
!!       exchange_threshold = gradient_threshold * 1.0d-3
!!
!!       libint_epsilon = (gradient_threshold * 1.0d-3)**2
!!
!!    unless stricter thresholds are already set on input or by default.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), intent(in) :: gradient_threshold
!
      if (wf%coulomb_threshold .gt. gradient_threshold*1.0d-3)  wf%coulomb_threshold  = gradient_threshold*1.0d-3
      if (wf%exchange_threshold .gt. gradient_threshold*1.0d-3) wf%exchange_threshold = gradient_threshold*1.0d-3
!
      wf%coulomb_threshold  = min(wf%coulomb_threshold,  1.0d-12)
      wf%exchange_threshold = min(wf%exchange_threshold, 1.0d-10)
!
      if (wf%libint_epsilon .gt. (wf%coulomb_threshold)**2) wf%libint_epsilon = (wf%coulomb_threshold)**2
!
   end subroutine set_screening_and_precision_thresholds_hf
!
!
   subroutine update_fock_mm_hf(wf)
!!
!!    Update Fock with polarizable QM/MM terms
!!    For now: QM/FQ model (see mm_class and output file)
!!    Written by Tommaso Giovannini, July 2019 for QM/MM
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:), allocatable                   :: potential_points
      integer :: i
!
      if(wf%system%mm%forcefield.eq.'fq') then
!
         if(.not.allocated(potential_points)) call mem%alloc(potential_points, wf%system%mm%n_atoms)
!
         call zero_array(wf%pol_emb_fock,wf%n_ao*wf%n_ao)
         call zero_array(wf%system%mm%pol_emb_rhs,wf%system%mm%n_variables)
!
!        electrostatic potential contracted with density : \sum_i V_mu(D_mu)(r_i)
!
         call wf%construct_ao_electrostatics(0,1,'prop',wf%system%mm%n_atoms,wf%system%mm%coordinates, &
                                             property_points=potential_points,ao_density=wf%ao_density) 
!
!        rhs for fq: -chi - V(D)
!
         wf%system%mm%pol_emb_rhs(1:wf%system%mm%n_atoms) = -wf%system%mm%chi + potential_points
!
!        solve q=D^-1 (-chi - V(D))
!
         call dgemm('N', 'N',                   &
                     wf%system%mm%n_variables,  &
                     1,                         &
                     wf%system%mm%n_variables,  &
                     one,                       &
                     wf%system%mm%fq_matrix,    &
                     wf%system%mm%n_variables,  &
                     wf%system%mm%pol_emb_rhs,  &
                     wf%system%mm%n_variables,  &
                     zero,                      &
                     wf%system%mm%pol_emb_lhs,  &
                     wf%system%mm%n_variables)
!
!
         call output%print_separator('verbose', 67, fs='(/t3,a)')
!
         call output%printf('Atom          FQ LHS             FQ RHS        QM Potential@FQs', &
                            pl='v', fs='(t6,a)')
!
         do i = 1, wf%system%mm%n_atoms
!
            call output%printf('(i4)      (e13.6)      (e13.6)      (e13.6)', pl='v', &
                               fs='(t6,a)', ints=[i], reals=[wf%system%mm%pol_emb_lhs(i), &
                               wf%system%mm%pol_emb_rhs(i), potential_points(i)])
!
         enddo
!
         call output%print_separator('verbose', 67)
!
!
!        put FQ charges into charge (I am discrading langrangian multipliers)
!
         wf%system%mm%charge = wf%system%mm%pol_emb_lhs(1:wf%system%mm%n_atoms)
!
!        Fock creation: F_munu = \sum_i q_i V_munu(r_i)
!
         call wf%construct_ao_electrostatics(0,0,'fock',wf%system%mm%n_atoms,wf%system%mm%coordinates, & 
                                             elec_fock=wf%pol_emb_fock,charges=wf%system%mm%charge) 
!
         wf%ao_fock = wf%ao_fock + half * wf%pol_emb_fock
!
!
         call output%print_matrix('debug', 'QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'FQ Fock', wf%pol_emb_fock, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM/FQ Fock', wf%ao_fock, wf%n_ao, wf%n_ao)
!
         call mem%dealloc(potential_points, wf%system%mm%n_atoms)
!
      else
!
         call output%error_msg('The only available polarizable force field is fq')
!
      endif
!
   end subroutine update_fock_mm_hf
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
      type(timings) :: s_timer, h_timer, G_timer, G_timer_sym, non_integral_timer
!
!     Initialize timers
!
      s_timer = timings('HF gradient - 1st derivative-integrals of S')
      h_timer = timings('HF gradient - 1st derivative-integrals of h')
      G_timer = timings('HF gradient - 1st derivative-integrals of G(D) - integrals')
      G_timer_sym = timings('HF gradient - 1st derivative-integrals of G(D) - symmetrization')
      non_integral_timer = timings('HF gradient - non-integral-time')
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
      s_wxqk = zero
      call wf%get_ao_s_wx_1der(s_wxqk)
!
      call s_timer%turn_off()
!
      call G_timer%turn_on()
!
      G_wxqk = zero
      call wf%construct_ao_G_1der(G_wxqk, wf%ao_density)
!
      call G_timer%turn_off()
      call G_timer_sym%turn_on()
!
      do k = 1, wf%system%n_atoms

         do q = 1, 3

           call symmetric_sum(G_wxqk(:,:,q,k), wf%n_ao)
           G_wxqk(:,:,q,k) = half*G_wxqk(:,:,q,k)

         enddo

      enddo
!
      call G_timer_sym%turn_off()
!
      call h_timer%turn_on()
!
      h_wxqk = zero
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
            TrDh_qk(q,k) = ddot(wf%n_ao**2, wf%ao_density, 1, h_wxqk(:,:,q,k), 1)
            TrDG_qk(q,k) = ddot(wf%n_ao**2, wf%ao_density, 1, G_wxqk(:,:,q,k), 1)
            TrDFDS_qk(q,k) = ddot(wf%n_ao**2, DFD, 1, s_wxqk(:,:,q,k), 1)
!
            E_qk(q,k) = E_qk(q,k)   &
              + TrDh_qk(q,k)        &
              + half*TrDG_qk(q,k)   &
              - half*TrDFDS_qk(q,k)
!
         enddo
      enddo
!
      call output%printf('Molecular gradient (Hartree/bohr):', fs='(/t6,a/)', pl='m')
!
      do k = 1, wf%system%n_atoms
!
         call output%printf('(a2)(f19.12)(f19.12)(f19.12)',         &
                           chars=[wf%system%atoms(k)%symbol],       &
                           reals=[E_qk(1,k), E_qk(2,k), E_qk(3,k)], &
                           fs='(t6,a)', pl='m')
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
   end subroutine construct_molecular_gradient_hf
!
!
   subroutine update_fock_pcm_hf(wf)
!!
!!    Update Fock PCM 
!!    Written by Tommaso Giovannini, Oct 2019 
!!
!!    The QM Fock is updated with the contributions coming 
!!    from the PCM:
!!       q*V_αβ
!!
!!    Done by interfacing to PCMSolver
!!
      use pcmsolver
      implicit none
!
      class(hf) :: wf
!
      character(kind=c_char, len=*), parameter :: mep_lbl = 'NucMEP'
      character(kind=c_char, len=*), parameter :: asc_lbl = 'NucASC'
      integer :: i
!     
! 
      if(.not.allocated(wf%system%pcm%pcm_rhs))   call mem%alloc(wf%system%pcm%pcm_rhs, wf%system%pcm%n_tesserae)
!
      call zero_array(wf%pcm_fock,wf%n_ao*wf%n_ao)
      call zero_array(wf%system%pcm%pcm_rhs,wf%system%pcm%n_tesserae)
!      
!     electrostatic potential contracted with density : \sum_i V_mu(D_mu)(r_i)
!
      call wf%construct_ao_electrostatics(0,1,'prop',wf%system%pcm%n_tesserae,wf%system%pcm%grid_coord*bohr_to_angstrom, &
                                          property_points=wf%system%pcm%pcm_rhs,ao_density=wf%ao_density) 
!      
!     solve q=D^-1 (V(D)) 
! 
      call pcmsolver_set_surface_function(wf%system%pcm%pcm_context, int(wf%system%pcm%n_tesserae, kind=c_int), &
                                          -wf%system%pcm%pcm_rhs,pcmsolver_fstring_to_carray(mep_lbl))
!                                          
      call pcmsolver_compute_asc(wf%system%pcm%pcm_context, &
                                 pcmsolver_fstring_to_carray(mep_lbl), &
                                 pcmsolver_fstring_to_carray(asc_lbl), &
                                 irrep=0_c_int)
!                                 
      call pcmsolver_get_surface_function(wf%system%pcm%pcm_context, int(wf%system%pcm%n_tesserae, kind=c_int), &
                                          wf%system%pcm%charges,pcmsolver_fstring_to_carray(asc_lbl))
!
      call output%print_separator('verbose', 67, fs='(/t3,a)')
!
      call output%printf('Atom         PCM ASC            PCM RHS', &
                         pl='v', fs='(t6,a)')
!
      do i = 1, wf%system%pcm%n_tesserae
!
         call output%printf('(i4)      (e13.6)      (e13.6)', pl='v', &
                            fs='(t6,a)', ints=[i], reals=[wf%system%pcm%charges(i), &
                            wf%system%pcm%pcm_rhs(i)])
!
      enddo
!
      call output%print_separator('verbose', 67)
!
!     Fock creation: F_munu = \sum_i q_i V_munu(r_i)
!
      call wf%construct_ao_electrostatics(0,0,'fock',wf%system%pcm%n_tesserae,wf%system%pcm%grid_coord*bohr_to_angstrom, & 
                                          elec_fock=wf%pcm_fock,charges=wf%system%pcm%charges) 
!
      wf%ao_fock = wf%ao_fock + half*wf%pcm_fock
!
      call output%print_matrix('debug', 'QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'PCM Fock', wf%pcm_fock, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM/PCM Fock', wf%ao_fock, wf%n_ao, wf%n_ao)
!
!
   end subroutine update_fock_pcm_hf
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
      call output%printf('- Summary of QM/PCM energetics:', pl='m', fs='(/t3,a)')
      call output%printf('a.u.             eV     kcal/mol', pl='m', fs='(t42,a)')
      call output%printf('QM/PCM SCF Contribution: (f21.12)', pl='m', &
                          reals=[wf%energy_scf_qmpcm], fs='(t6,a)')
      call output%printf('QM/PCM Electrostatic Energy:(f18.12)(f12.5) (f9.3)', pl='m', &
                          reals=[wf%electrostatic_energy_qmpcm,               &
                                 wf%electrostatic_energy_qmpcm*Hartree_to_eV, &
                                 wf%electrostatic_energy_qmpcm*Hartree_to_kcalmol], fs='(t6,a)')
!
!
   end subroutine print_energy_pcm_hf!
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
      real(dp), dimension(:,:), allocatable :: h_wx

      real(dp) :: n_electrons
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(h_wx)
!
      call wf%update_fock_and_energy(h_wx)
!
      call wf%get_n_electrons_in_density(n_electrons)
!
      call output%printf('Energy of initial guess:      (f25.12)', reals=[wf%energy], fs='(/t6, a)',pl='minimal')
      call output%printf('Number of electrons in guess: (f25.12)', reals=[n_electrons], fs='(t6, a)',pl='minimal')
!
!     Update the orbitals and density to make sure the density is idempotent
!     (not the case for the standard atomic superposition density)
!
      call wf%roothan_hall_update_orbitals() ! F => C
      call wf%update_ao_density()            ! C => D
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
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
      call wf%initialize_sp_eri_schwarz()
      call wf%initialize_sp_eri_schwarz_list()
!
      call wf%construct_sp_eri_schwarz()
!
      if(wf%system%mm_calculation) call wf%initialize_mm_matrices()
      if(wf%system%pcm_calculation) call wf%initialize_pcm_matrices()
!
   end subroutine prepare_hf
!
!
   subroutine write_scf_restart_hf(wf)
!!
!!    Write HF restart file
!!    Written by Linda Goletto, Oct 2019
!!
!!    Writes a file used for consistency checks when restarting
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%restart_file%open_('write', 'rewind')
!
      call wf%restart_file%write_(wf%n_ao)
      call wf%restart_file%write_(wf%n_densities)
      call wf%restart_file%write_(wf%system%get_n_electrons())
!
      call wf%restart_file%close_
!
   end subroutine write_scf_restart_hf
!
!
   subroutine write_orbital_information_hf(wf)
!!
!!    Write HF information file
!!    Written by Linda Goletto, Oct 2019
!!
!!    Writes orbital information
!!    Used for CC and MP2 and for restarting MLHF
!!
      implicit none
!
      class(hf) :: wf
!
      type(sequential_file) :: CC_orbitals_file, CC_orbital_energies_file
!
      wf%orbital_information_file = sequential_file('orbital_information')
      call wf%orbital_information_file%open_('write', 'rewind')
!
      call wf%orbital_information_file%write_(wf%n_o)
      call wf%orbital_information_file%write_(wf%n_v)
      call wf%orbital_information_file%write_(wf%n_ao)
      call wf%orbital_information_file%write_(wf%n_mo)
      call wf%orbital_information_file%write_(wf%energy)
!
      call wf%orbital_information_file%close_
!
      CC_orbitals_file = sequential_file('cc_orbital_coefficients')
      call CC_orbitals_file%open_('write', 'rewind')
!
      call CC_orbitals_file%write_(wf%orbital_coefficients, wf%n_ao*wf%n_mo)
!
      call CC_orbitals_file%close_('keep')
!
      CC_orbital_energies_file = sequential_file('cc_orbital_energies')
      call CC_orbital_energies_file%open_('write', 'rewind')
!
      call CC_orbital_energies_file%write_(wf%orbital_energies, wf%n_mo)
!
      call CC_orbital_energies_file%close_('keep')
!
   end subroutine write_orbital_information_hf
!
!
   subroutine save_orbital_coefficients_hf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%orbital_coefficients_file%open_('write', 'rewind')
!
      call wf%orbital_coefficients_file%write_(wf%orbital_coefficients, wf%n_ao*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
!
   end subroutine save_orbital_coefficients_hf
!
!
   subroutine read_orbital_coefficients_hf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%orbital_coefficients_file%open_('read', 'rewind')
!
      call wf%orbital_coefficients_file%read_(wf%orbital_coefficients, wf%n_ao*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
!
   end subroutine read_orbital_coefficients_hf
!
!
   subroutine save_orbital_energies_hf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%orbital_energies_file%open_('write', 'rewind')
!
      call wf%orbital_energies_file%write_(wf%orbital_energies, wf%n_mo)
!
      call wf%orbital_energies_file%close_
!
   end subroutine save_orbital_energies_hf
!
!
   subroutine read_orbital_energies_hf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%orbital_energies_file%open_('read', 'rewind')
!
      call wf%orbital_energies_file%read_(wf%orbital_energies, wf%n_mo)
!
      call wf%orbital_energies_file%close_
!
   end subroutine read_orbital_energies_hf
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
end module hf_class
