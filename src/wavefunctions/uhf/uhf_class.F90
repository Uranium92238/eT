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
module uhf_class
!
!!
!!    Unrestricted Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!
   use hf_class
!
   use omp_lib
!
   implicit none
!
!
   type, extends(hf) :: uhf
!
      integer :: n_alpha ! Number of alpha electrons 
      integer :: n_beta  ! Number of beta electrons 
!
      real(dp), dimension(:,:), allocatable :: ao_density_a ! Alpha density 
      real(dp), dimension(:,:), allocatable :: ao_density_b ! Beta density 
!
      real(dp), dimension(:,:), allocatable :: ao_fock_a    ! Alpha Fock 
      real(dp), dimension(:,:), allocatable :: ao_fock_b    ! Beta Fock 
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_a  ! Alpha orbital coeffs.
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_b  ! Beta orbital coeffs.
!
      real(dp), dimension(:), allocatable :: orbital_energies_a ! Alpha orbital energies 
      real(dp), dimension(:), allocatable :: orbital_energies_b ! Beta orbital energies
!
      logical :: fractional_uniform_valence = .false. ! If true, smears electrons for degenerate 
                                                      ! HOMOs. Used to generate spherically 
                                                      ! symmetrical SAD guess. Gives non-UHF 
                                                      ! energy and should be false if the 
                                                      ! user wants standard UHF numbers.
!
   contains
!
!     Preparation routines
!
      procedure :: prepare                               => prepare_uhf
      procedure :: determine_n_alpha_and_n_beta          => determine_n_alpha_and_n_beta_uhf
      procedure :: read_settings                         => read_settings_uhf
      procedure :: read_uhf_settings                     => read_uhf_settings_uhf
!
!     AO Fock and energy related routines
!
      procedure :: initialize_fock                       => initialize_fock_uhf
      procedure :: construct_ao_spin_fock                => construct_ao_spin_fock_uhf
      procedure :: calculate_uhf_energy                  => calculate_uhf_energy_uhf
      procedure :: update_fock_and_energy_non_cumulative => update_fock_and_energy_non_cumulative_uhf
      procedure :: update_fock_and_energy_cumulative     => update_fock_and_energy_cumulative_uhf
      procedure :: update_fock_and_energy                => update_fock_and_energy_uhf
      procedure :: update_fock_mm                        => update_fock_mm_uhf
      procedure :: update_fock_pcm                       => update_fock_pcm_uhf
      procedure :: set_ao_fock                           => set_ao_fock_uhf
      procedure :: get_ao_fock                           => get_ao_fock_uhf
      procedure :: print_energy                          => print_energy_uhf
!
!     AO Density related routines
!
      procedure :: initialize_density                    => initialize_density_uhf
      procedure :: set_initial_ao_density_guess          => set_initial_ao_density_guess_uhf
      procedure :: save_ao_density                       => save_ao_density_uhf
      procedure :: update_ao_density                     => update_ao_density_uhf
      procedure :: form_ao_density                       => form_ao_density_uhf
      procedure :: construct_ao_spin_density             => construct_ao_spin_density_uhf
      procedure :: set_ao_density_to_core_guess          => set_ao_density_to_core_guess_uhf
      procedure :: get_homo_degeneracy                   => get_homo_degeneracy_uhf
      procedure :: get_ao_density_sq                     => get_ao_density_sq_uhf
!
      procedure :: construct_mo_fock                     => construct_mo_fock_uhf
!
!     MO orbital related routines
!

      procedure :: initialize_orbitals                   => initialize_orbitals_uhf
      procedure :: roothan_hall_update_orbitals          => roothan_hall_update_orbitals_uhf
      procedure :: save_orbital_info                     => save_orbital_info_uhf
      procedure :: save_orbital_coefficients             => save_orbital_coefficients_uhf
      procedure :: read_orbital_coefficients             => read_orbital_coefficients_uhf
      procedure :: save_orbital_energies                 => save_orbital_energies_uhf
      procedure :: print_spin                            => print_spin_uhf
      procedure :: print_summary                         => print_summary_uhf
      procedure :: read_orbital_energies                 => read_orbital_energies_uhf
!
!     Roothan-Hall gradient 
!
      procedure :: get_packed_roothan_hall_gradient      => get_packed_roothan_hall_gradient_uhf
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density_a               => initialize_ao_density_a_uhf
      procedure :: initialize_ao_density_b               => initialize_ao_density_b_uhf
!
      procedure :: initialize_ao_fock_a                  => initialize_ao_fock_a_uhf
      procedure :: initialize_ao_fock_b                  => initialize_ao_fock_b_uhf
!
      procedure :: initialize_orbital_coefficients_a     => initialize_orbital_coefficients_a_uhf
      procedure :: initialize_orbital_coefficients_b     => initialize_orbital_coefficients_b_uhf
!
      procedure :: initialize_orbital_energies_a         => initialize_orbital_energies_a_uhf
      procedure :: initialize_orbital_energies_b         => initialize_orbital_energies_b_uhf
!
      procedure :: destruct_ao_density                   => destruct_ao_density_uhf
      procedure :: destruct_ao_density_a                 => destruct_ao_density_a_uhf
      procedure :: destruct_ao_density_b                 => destruct_ao_density_b_uhf
!
      procedure :: destruct_fock                         => destruct_fock_uhf
      procedure :: destruct_ao_fock_a                    => destruct_ao_fock_a_uhf
      procedure :: destruct_ao_fock_b                    => destruct_ao_fock_b_uhf
!
      procedure :: destruct_orbital_coefficients         => destruct_orbital_coefficients_uhf
      procedure :: destruct_orbital_coefficients_a       => destruct_orbital_coefficients_a_uhf
      procedure :: destruct_orbital_coefficients_b       => destruct_orbital_coefficients_b_uhf
!
      procedure :: destruct_orbital_energies             => destruct_orbital_energies_uhf
      procedure :: destruct_orbital_energies_a           => destruct_orbital_energies_a_uhf
      procedure :: destruct_orbital_energies_b           => destruct_orbital_energies_b_uhf
!
      procedure :: set_n_mo                              => set_n_mo_uhf
!
      procedure :: cleanup                               => cleanup_uhf 
!
      procedure :: get_spin_contamination               => get_spin_contamination_uhf
      procedure :: get_exact_s2                         => get_exact_s2_uhf
!
   end type uhf
!
!
   interface
!
      include "ao_fock_uhf_interface.F90"
      include "set_get_uhf_interface.F90"
      include "file_handling_uhf_interface.F90"
      include "initialize_destruct_uhf_interface.F90"
!
   end interface
!
!
   interface uhf 
!
      procedure :: new_uhf 
      procedure :: new_uhf_from_parameters
!
   end interface uhf 
!
!
contains 
!
!
   function new_uhf(system) result(wf)
!!
!!    New UHF
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(uhf) :: wf
!
      class(molecular_system), target, intent(in) :: system
!
      wf%system => system
!
      wf%name_ = 'uhf'
!
      call wf%read_settings()
!
      call wf%print_banner()
!
      call wf%prepare()
!
   end function new_uhf
!
!
   function new_uhf_from_parameters(system, fractional_uniform_valence) result(wf)
!!
!!    New UHF from parameters
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      type(uhf) :: wf
!
      class(molecular_system), target, intent(in) :: system
      logical,                         intent(in) :: fractional_uniform_valence
!
      wf%system => system
!
      wf%name_ = 'uhf'
!
      wf%fractional_uniform_valence = fractional_uniform_valence
!
      call wf%prepare()
!
   end function new_uhf_from_parameters
!
!
   subroutine prepare_uhf(wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initializes files
!!    and constructs screening vectors
!!
      implicit none
!
      class(uhf) :: wf
!
      wf%n_ao = wf%system%get_n_aos()
!
      call wf%set_n_mo()
!
      if (wf%fractional_uniform_valence) then
!
         call output%printf('m', 'Requested fractional uniform valence. Valence &
                            &electrons will be  distributed evenly in the &
                            &highest molecular orbitals (if plural).', ffs='(/t3,a)')
!
      endif
!
      wf%orbital_coefficients_file = sequential_file('orbital_coefficients')
      wf%orbital_energies_file = sequential_file('orbital_energies')
      wf%restart_file = sequential_file('scf_restart_file')
!
      call wf%initialize_shp_eri_schwarz() 
      call wf%initialize_shp_eri_schwarz_list()
!
      call wf%construct_shp_eri_schwarz()
!
      if (wf%system%mm_calculation) call wf%prepare_qmmm()
      if(wf%system%pcm_calculation) call wf%initialize_pcm_matrices()
!
      wf%frozen_core = .false.
      wf%frozen_hf_mos = .false.
!
      call wf%initialize_ao_h()
      call wf%get_ao_h_wx(wf%ao_h)
!
   end subroutine prepare_uhf
!
!
   subroutine set_initial_ao_density_guess_uhf(wf, guess)
!!
!!    Set initial AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Sets initial AO density (or densities) to the
!!    appropriate initial guess requested by the solver.
!!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(uhf) :: wf
!
      character(len=*) :: guess
!
      real(dp) :: alpha_prefactor, beta_prefactor 
!
      if (trim(guess) == 'sad' .or. trim(guess) == 'SAD') then
!
         call wf%set_ao_density_to_sad()
!
         alpha_prefactor = real(wf%n_alpha, kind=dp)/real(wf%n_alpha + wf%n_beta, kind=dp)
         beta_prefactor  = real(wf%n_beta,  kind=dp)/real(wf%n_alpha + wf%n_beta, kind=dp)
!
         call copy_and_scale(alpha_prefactor, wf%ao_density, wf%ao_density_a, wf%n_ao**2)
         call copy_and_scale(beta_prefactor,  wf%ao_density, wf%ao_density_b, wf%n_ao**2)
!
      elseif (trim(guess) == 'core' .or. trim(guess) == 'CORE') then
!
         call wf%set_ao_density_to_core_guess(wf%ao_h)
!
      else
!
         call output%error_msg('Guess AO density ' // trim(guess) // ' is currently not supported.')
!
      endif
!
   end subroutine set_initial_ao_density_guess_uhf
!
!
   subroutine get_packed_roothan_hall_gradient_uhf(wf, G)
!!
!!    Get packed Roothan-Hall gradient
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Constructs and returns the gradient as an
!!    anti-symmetrized packed vector. For UHF,
!!    both the alpha and beta gradients are
!!    returned as follows: [G_a G_b]
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo*(wf%n_mo - 1)/2, wf%n_densities), intent(inout) :: G
!
      real(dp), dimension(:,:), allocatable :: G_sq
      real(dp), dimension(:), allocatable :: G_pck
      real(dp), dimension(:,:), allocatable :: Po, Pv
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
      call mem%alloc(G_pck, wf%n_mo*(wf%n_mo - 1)/2)
      call mem%alloc(G_sq, wf%n_mo, wf%n_mo)
!
!     Alpha gradient
!
      call wf%construct_projection_matrices(Po, Pv, wf%ao_density_a)
      call wf%construct_roothan_hall_gradient(G_sq, Po, Pv, wf%ao_fock_a)
      call packin_anti(G_pck, G_sq, wf%n_ao)
      call dcopy(wf%n_ao*(wf%n_ao - 1)/2, G_pck, 1, G, 1)
!
!     Beta gradient
!
      call wf%construct_projection_matrices(Po, Pv, wf%ao_density_b)
      call wf%construct_roothan_hall_gradient(G_sq, Po, Pv, wf%ao_fock_b)
      call packin_anti(G_pck, G_sq, wf%n_ao)
      call dcopy(wf%n_ao*(wf%n_ao - 1)/2, G_pck, 1, G(1, 2), 1)
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
      call mem%dealloc(G_pck, wf%n_mo*(wf%n_mo - 1)/2)
      call mem%dealloc(G_sq, wf%n_mo, wf%n_mo)
!
   end subroutine get_packed_roothan_hall_gradient_uhf
!
!
   subroutine read_settings_uhf(wf)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Designed to be overwritten by descendants.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%read_hf_settings()
      call wf%read_uhf_settings()
!
   end subroutine read_settings_uhf
!
!
   subroutine read_uhf_settings_uhf(wf)
!!
!!    Read UHF settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Reads settings specific to the wavefunction.
!!
      implicit none
!
      class(uhf) :: wf
!
      wf%fractional_uniform_valence = &
            input%requested_keyword_in_section('fractional uniform valence', 'hf')
!
   end subroutine read_uhf_settings_uhf
!
!
   subroutine roothan_hall_update_orbitals_uhf(wf)
!!
!!    Roothan-Hall update of orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    This routine guides the construction of new orbital coefficients
!!    from the current AO Fock matrix (or matrices if the wavefunction
!!    is unrestricted).
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%do_roothan_hall(wf%ao_fock_a, wf%orbital_coefficients_a, wf%orbital_energies_a)
      call wf%do_roothan_hall(wf%ao_fock_b, wf%orbital_coefficients_b, wf%orbital_energies_b)
!
   end subroutine roothan_hall_update_orbitals_uhf
!
!
   subroutine update_ao_density_uhf(wf)
!!
!!    Update AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Updates the AO density (or densities, if unrestricted) based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%construct_ao_spin_density('alpha') ! Make D_alpha
      call wf%construct_ao_spin_density('beta')  ! Make D_beta
!
      call wf%form_ao_density()                  ! Set D = D_alpha + D_beta
!
   end subroutine update_ao_density_uhf
!
!
   subroutine set_ao_density_to_core_guess_uhf(wf, h_wx)
!!
!!    Set AO density to core guess
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Solves the Roothan-Hall equation ignoring the two-
!!    electron AO density dependent part of the Fock matrix.
!!
!!    Based on the orbital coefficients, the routine constructs
!!    the associated AO spin densities.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      call dcopy(wf%n_ao**2, h_wx, 1, wf%ao_fock, 1)
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies)
!
      call dcopy(wf%n_ao*wf%n_mo, wf%orbital_coefficients, 1, wf%orbital_coefficients_a, 1)
      call dcopy(wf%n_ao*wf%n_mo, wf%orbital_coefficients, 1, wf%orbital_coefficients_b, 1)
!
      call dcopy(wf%n_mo, wf%orbital_energies, 1, wf%orbital_energies_a, 1)
      call dcopy(wf%n_mo, wf%orbital_energies, 1, wf%orbital_energies_b, 1)
!
      call wf%construct_ao_spin_density('alpha')
      call wf%construct_ao_spin_density('beta')
!
      call wf%form_ao_density()
!
   end subroutine set_ao_density_to_core_guess_uhf
!
!
   subroutine determine_n_alpha_and_n_beta_uhf(wf)
!!
!!    Determine the number of alpha and beta electrons
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Determines the number of alpha and beta elctrons
!!    from the provided multiplicity and the number of
!!    electrons in the system. We assume that n_alpha
!!    is greater than or equal to n_beta without loss
!!    of generality.
!!
      implicit none
!
      class(uhf) :: wf
!
      integer :: n_alpha_m_n_beta, n_alpha_p_n_beta
!
      n_alpha_m_n_beta = wf%system%multiplicity - 1 ! 2 S + 1 - 1 = 2 S
      n_alpha_p_n_beta = wf%system%n_electrons
!
      wf%n_alpha = (n_alpha_p_n_beta + n_alpha_m_n_beta)/2
      wf%n_beta  = wf%system%n_electrons - wf%n_alpha
!
      if ((wf%n_alpha + wf%n_beta) .ne. wf%system%n_electrons) then
!
         call output%error_msg('Given multiplicity and number of electrons is inconsistent.')
!
      endif
!
   end subroutine determine_n_alpha_and_n_beta_uhf
!
!
   subroutine construct_ao_spin_density_uhf(wf, sigma)
!!
!!    Construct AO spin density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Constructs either the alpha or beta density,
!!    depending on the value of 'spin' (='alpha' or 'beta').
!!    The routine assumes that the associated orbital coefficients
!!    are allocated and set. Typically, the construction of these
!!    densities will follow a Roothan-Hall step.
!!
!!    We employ the normalization so that D = D_alpha + D_beta is
!!    equal to the restricted HF density if the alpha and beta
!!    densities are equal.
!!
!!    If fractional_uniform_valence is set to true, the routine
!!    will distribute the HOMO electrons (evenly, in fractions) among
!!    the degenerate HOMO orbitals. Note that this restricts the UHF
!!    densities to be spherically symmetric if the zeroth iteration
!!    density possesses this symmetry.
!!
      implicit none
!
      class(uhf) :: wf
!
      character(len=*), intent(in) :: sigma
!
      integer :: n_homo_electrons
      integer :: n_homo_orbitals
      integer :: homo_first, homo_last
      real(dp)     :: electrons_to_fill
      integer :: alpha, beta, i
!
      if (trim(sigma) == 'alpha') then
!
         if (.not. wf%fractional_uniform_valence) then ! Standard
!
            call zero_array(wf%ao_density_a, wf%n_ao**2)
            if (wf%n_alpha .eq. 0) return
!
            call dgemm('N', 'T',                   &
                        wf%n_ao,                   &
                        wf%n_ao,                   &
                        wf%n_alpha,                &
                        one,                       &
                        wf%orbital_coefficients_a, &
                        wf%n_ao,                   &
                        wf%orbital_coefficients_a, &
                        wf%n_ao,                   &
                        zero,                      &
                        wf%ao_density_a,           &
                        wf%n_ao)
!
         else ! Smear out HOMO electrons, if degenerate
!
            call wf%get_homo_degeneracy(wf%orbital_energies_a, homo_first, homo_last, &
                                          n_homo_orbitals, n_homo_electrons, wf%n_alpha)
!
            call zero_array(wf%ao_density_a, wf%n_ao**2)
            if (wf%n_alpha .eq. 0) return
!
            electrons_to_fill = real(n_homo_electrons, kind=dp)/real(n_homo_orbitals, kind=dp)
!
            do alpha = 1, wf%n_ao
               do beta = 1, wf%n_ao
!
                  do i = 1, homo_first - 1
!
                     wf%ao_density_a(alpha,beta) = wf%ao_density_a(alpha,beta) &
                        + wf%orbital_coefficients_a(alpha,i)*wf%orbital_coefficients_a(beta,i)
!
                  enddo
!
                  do i = homo_first, homo_last
!
                     wf%ao_density_a(alpha,beta) = wf%ao_density_a(alpha,beta)  &
                        + electrons_to_fill*wf%orbital_coefficients_a(alpha,i)*wf%orbital_coefficients_a(beta,i)
!
                  enddo
!
               enddo
            enddo
!
         endif
!
      elseif (trim(sigma) == 'beta') then
!
         if (.not. wf%fractional_uniform_valence) then ! Standard
!
            call zero_array(wf%ao_density_b, wf%n_ao**2)
            if (wf%n_beta .eq. 0) return
!
            call dgemm('N', 'T',                   &
                        wf%n_ao,                   &
                        wf%n_ao,                   &
                        wf%n_beta,                 &
                        one,                       &
                        wf%orbital_coefficients_b, &
                        wf%n_ao,                   &
                        wf%orbital_coefficients_b, &
                        wf%n_ao,                   &
                        zero,                      &
                        wf%ao_density_b,           &
                        wf%n_ao)
!
         else ! Smear out HOMO electrons, if degenerate
!
            call wf%get_homo_degeneracy(wf%orbital_energies_b, homo_first, homo_last, &
                                          n_homo_orbitals, n_homo_electrons, wf%n_beta)
!
            call zero_array(wf%ao_density_b, wf%n_ao**2)
            if (wf%n_beta .eq. 0) return
!
            electrons_to_fill = real(n_homo_electrons, kind=dp)/real(n_homo_orbitals, kind=dp)
!
            do alpha = 1, wf%n_ao
               do beta = 1, wf%n_ao
!
                  do i = 1, homo_first - 1
!
                     wf%ao_density_b(alpha, beta) = wf%ao_density_b(alpha, beta) &
                        + wf%orbital_coefficients_b(alpha, i)*wf%orbital_coefficients_b(beta, i)
!
                  enddo
!
                  do i = homo_first, homo_last
!
                     wf%ao_density_b(alpha, beta) = wf%ao_density_b(alpha, beta) &
                        + electrons_to_fill*wf%orbital_coefficients_b(alpha, i)*wf%orbital_coefficients_b(beta, i)
!
                  enddo
!
               enddo
            enddo
!
         endif
!
      else ! Not good: sigma is neither alpha nor beta
!
         call output%error_msg('Did not recognize spin variable in construct_ao_spin_density:' // trim(sigma))
!
      endif
!
   end subroutine construct_ao_spin_density_uhf
!
!
   subroutine get_homo_degeneracy_uhf(wf, energies, homo_first, homo_last, n_homo_orbitals, n_homo_electrons, n_electrons)
!!
!!    Get HOMO degeneracy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    When filling electrons in orbitals (one in each), this routine determines
!!    the degeneracy of the highest molecular orbital and the number of electrons
!!    that should be distributed to them.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_mo), intent(in) :: energies
!
      integer, intent(inout) :: n_homo_orbitals, n_homo_electrons, homo_first, homo_last
!
      integer, intent(in) :: n_electrons
!
      integer :: homo, n_below, n_above, I
!
      real(dp), parameter :: threshold = 1.0D-12
!
!     Set standard non-degenerate values
!
      homo             = n_electrons
      homo_first       = homo
      homo_last        = homo
      n_homo_orbitals  = 1
      n_homo_electrons = 1
!
      if (n_electrons .eq. 0) return
!
      n_below = 0
      do I = 1, homo - 1
!
         if (abs(energies(homo) - energies(I)) .le. threshold) then
!
            n_below = n_below + 1
!
         endif
!
      enddo
!
      n_above = 0
      do I = homo + 1, wf%n_mo
!
         if (abs(energies(homo) - energies(I)) .le. threshold) then
!
            n_above = n_above + 1
!
         endif
!
      enddo
!
      n_homo_orbitals  = n_below + n_above + 1
      n_homo_electrons = n_below + 1
      homo_first = homo - n_below
      homo_last  = homo + n_above
!
   end subroutine get_homo_degeneracy_uhf
!
!
   subroutine form_ao_density_uhf(wf)
!!
!!    Form AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Given that the spin densities have been constructed,
!!    this routine sets
!!
!!       D = D_alpha + D_beta
!!
!!    Note that it is not constructed from the orbital
!!    coefficients and should be preceded by calls to
!!    "construct_ao_spin_density" to make sure that the
!!    spin densities are up-to-date.
!!
      implicit none
!
      class(uhf) :: wf
!
      call dcopy(wf%n_ao**2, wf%ao_density_a, 1, wf%ao_density, 1)
      call daxpy(wf%n_ao**2, one, wf%ao_density_b, 1, wf%ao_density, 1)
!
   end subroutine form_ao_density_uhf
!
!
   subroutine calculate_uhf_energy_uhf(wf, h_wx)
!!
!!    Calculate UHF energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the unrestricted Hartree-Fock energy,
!!
!!       E = E_alpha + E_beta = 1/2 Tr (D^alpha h) + 1/2 Tr (D^alpha F)
!!                            + 1/2 Tr (D^beta h)  + 1/2 Tr (D^beta F)
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: h_wx
!
      real(dp) :: ddot
!
      wf%energy = wf%system%get_nuclear_repulsion()
!
      if (wf%system%mm_calculation) call wf%calculate_mm_energy_terms()
      if (wf%system%pcm_calculation) call wf%calculate_pcm_energy_terms()
!      
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_a, 1, h_wx, 1)
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_a, 1, wf%ao_fock_a, 1)
!
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_b, 1, h_wx, 1)
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_b, 1, wf%ao_fock_b, 1)
!
   end subroutine calculate_uhf_energy_uhf
!
!
   subroutine set_n_mo_uhf(wf)
!!
!!    Set number of molecular orbitals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      call output%printf('n', '- Cholesky decomposition of AO overlap to get &
                         &linearly independent orbitals:', ll=100, fs='(/t3,a)')
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap()
!
      wf%n_densities = 2 ! [D_alpha, D_beta]
!
      call wf%determine_n_alpha_and_n_beta()
!
!
      call output%printf('m', '- Orbital details:', &
                         fs='(/t3,a)')
!
      call output%printf('m', 'Number of alpha electrons:        (i8)', &
                         ints=[wf%n_alpha], fs='(/t6,a)')
      call output%printf('m', 'Number of beta electrons:         (i8)', &
                         ints=[wf%n_beta], fs='(t6,a)')
!
      call output%printf('m', 'Number of virtual alpha orbitals: (i8)', &
                         ints=[wf%n_mo - wf%n_alpha], fs='(t6,a)')
      call output%printf('m', 'Number of virtual beta orbitals:  (i8)', &
                         ints=[wf%n_mo - wf%n_beta], fs='(t6,a)')
!
      call output%printf('m', 'Number of molecular orbitals:     (i8)', &
                         ints=[wf%n_mo], fs='(t6,a)')
      call output%printf('m', 'Number of atomic orbitals:        (i8)', &
                         ints=[wf%n_ao], fs='(t6,a)')
!
      if (wf%n_mo .lt. wf%n_ao) &
         call output%printf('m', 'Removed (i0) AOs due to linear dep.', &
                            ints=[wf%n_ao - wf%n_mo], fs='(/t6,a)')
!
   end subroutine set_n_mo_uhf
!
!
   subroutine print_energy_uhf(wf)
!!
!!    Print energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints information related to the wavefunction,
!!    most of which is meaningful only for a properly
!!    converged wavefunction. Should be overwritten in
!!    descendants if more or less or other information
!!    is present.
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      real(dp) :: homo_lumo_gap_a
      real(dp) :: homo_lumo_gap_b
      real(dp) :: nuclear_repulsion
!
      if (wf%n_alpha > 0 .and. wf%n_alpha < wf%n_mo) then 
!
         homo_lumo_gap_a = wf%orbital_energies_a(wf%n_alpha + 1) - wf%orbital_energies_a(wf%n_alpha)
         call output%printf('m', 'HOMO-LUMO gap (alpha):     (f19.12)', &
                            reals=[homo_lumo_gap_a], fs='(/t6,a)')
!
      endif 
!
      if (wf%n_beta > 0 .and. wf%n_beta < wf%n_mo) then 
!
         homo_lumo_gap_b = wf%orbital_energies_b(wf%n_beta + 1) - wf%orbital_energies_b(wf%n_beta)
         call output%printf('m', 'HOMO-LUMO gap (beta):      (f19.12)', &
                            reals=[homo_lumo_gap_b], fs='(t6,a)')
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
   end subroutine print_energy_uhf
!
!
   subroutine save_orbital_info_uhf(wf)
!!
!!    Make orbital info file
!!    Written by Alexander C. Paul Nov 2020
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      type(output_file) :: mo_information_file
!
      mo_information_file = output_file('mo_information.out')
      call mo_information_file%open_('rewind')
!
      call wf%print_orbitals_and_energies(mo_information_file, wf%orbital_energies_a, &
                                          wf%orbital_coefficients_a, &
                                          '- Alpha molecular orbital')
!
      call mo_information_file%print_separator('m', 83, '=', fs='(//t3,a/)')
!
      call wf%print_orbitals_and_energies(mo_information_file, wf%orbital_energies_b, &
                                          wf%orbital_coefficients_b, &
                                          '- Beta molecular orbital')
!
      call mo_information_file%close_()
!
   end subroutine save_orbital_info_uhf
!
!
   subroutine cleanup_uhf(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%save_ao_density()
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
      call wf%destruct_shp_eri_schwarz()
      call wf%destruct_shp_eri_schwarz_list()
      call wf%destruct_ao_h()
!
      call wf%destruct_W_mo_update()
      call wf%destruct_mo_fock()
!
      call wf%destruct_mm_matrices()
      call wf%destruct_pcm_matrices()
!
   end subroutine cleanup_uhf
!
!
   function get_spin_contamination_uhf(wf) result(spin_contamination)
!!
!!    Get spin contamination
!!    Written by Sarai D. Folkestad, 2020
!!
!!    S^2 - S_z(S_z + 1) = n_beta - sum_ij |<phi_i^alpha|phi_j^beta>|^2
!!                       = n_beta - Tr([D^alpha S][D^beta S])
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp) :: spin_contamination
!
      real(dp) :: ddot
!
      real(dp), dimension(:,:), allocatable :: DS_alpha, DS_beta
!
      spin_contamination = real(wf%n_beta,dp)
!
      call mem%alloc(DS_alpha, wf%n_ao, wf%n_ao)
      call mem%alloc(DS_beta, wf%n_ao, wf%n_ao)
!
      call dgemm('n', 'n',          &
                  wf%n_ao,          &
                  wf%n_ao,          &
                  wf%n_ao,          &
                  one,              &
                  wf%ao_density_a,  &
                  wf%n_ao,          &
                  wf%ao_overlap,    &
                  wf%n_ao,          &
                  zero,             &
                  DS_alpha,         &
                  wf%n_ao)
!
      call dgemm('t', 't',          &
                  wf%n_ao,          &
                  wf%n_ao,          &
                  wf%n_ao,          &
                  one,              &
                  wf%ao_overlap,    &
                  wf%n_ao,          &
                  wf%ao_density_b,  &
                  wf%n_ao,          &
                  zero,             &
                  DS_beta,          &
                  wf%n_ao)
!
      spin_contamination = spin_contamination &
            - ddot(wf%n_ao**2, DS_alpha, 1, DS_beta, 1)
!
      call mem%dealloc(DS_alpha, wf%n_ao, wf%n_ao)
      call mem%dealloc(DS_beta, wf%n_ao, wf%n_ao)
!
   end function get_spin_contamination_uhf
!
!
   pure function get_exact_s2_uhf(wf) result(spin)
!!
!!    Get exact S^2
!!    Written by Sarai D. Folkestad, 2020
!!
!!    S^2 = S_z(S_z + 1) = (n_alpha - n_beta)/2*((n_alpha - n_beta)/2 + 1)
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp) :: spin, sz
!
      spin = half*real(wf%n_alpha - wf%n_beta,dp)*(half*real(wf%n_alpha - wf%n_beta,dp) + one)
!
   end function get_exact_s2_uhf
!
!
   subroutine print_summary_uhf(wf, print_mo_info)
!!
!!    Print Summary
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none 
!
      class(uhf), intent(inout) :: wf
!
      logical, intent(in) :: print_mo_info
!
      call output%printf('m', '- Summary of '// &
                         &trim(convert_to_uppercase(wf%name_))// ' wavefunction &
                         &energetics (a.u.):', fs='(/t3,a)')
      call wf%print_energy()
!
      if (print_mo_info) call wf%save_orbital_info()
!
      call output%printf('m', '- '// &
                         &trim(convert_to_uppercase(wf%name_))// ' wavefunction spin expectation values:', fs='(/t3,a)')     
!
      call wf%print_spin()
!
   end subroutine print_summary_uhf
!
!
   subroutine print_spin_uhf(wf)
!!
!!    Print spin
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp) :: exact_S2
      real(dp) :: contamination
!
      exact_S2      = wf%get_exact_s2()
      contamination = wf%get_spin_contamination()
!
      call output%printf('m', 'Sz:                 (f12.8)', &
                         reals=[half*real(wf%n_alpha - wf%n_beta,dp)], fs='(/t6,a)')
!
      call output%printf('m', 'Sz(Sz + 1):         (f12.8)', &
                         reals=[exact_S2], fs='(t6,a)')
!
      call output%printf('m', 'S^2:                (f12.8)', &
                         reals=[exact_S2 + contamination], fs='(t6,a)')
!
      call output%printf('m', 'Spin contamination: (f12.8)', &
                         reals=[contamination], fs='(t6,a)')
!
   end subroutine print_spin_uhf
!
!
end module uhf_class
