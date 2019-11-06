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
!  Unrestricted Hartree-Fock wavefunction
!
   type, extends(hf) :: uhf
!
      integer :: n_alpha
      integer :: n_beta
!
      real(dp), dimension(:,:), allocatable :: ao_density_a
      real(dp), dimension(:,:), allocatable :: ao_density_b
!
      real(dp), dimension(:,:), allocatable :: ao_fock_a
      real(dp), dimension(:,:), allocatable :: ao_fock_b
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_a
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_b
!
      real(dp), dimension(:), allocatable :: orbital_energies_a
      real(dp), dimension(:), allocatable :: orbital_energies_b
!
      logical :: fractional_uniform_valence = .false.
!
   contains
!
!     Preparation routines
!
      procedure :: prepare                              => prepare_uhf
      procedure :: determine_n_alpha_and_n_beta         => determine_n_alpha_and_n_beta_uhf
      procedure :: read_settings                        => read_settings_uhf
      procedure :: read_uhf_settings                    => read_uhf_settings_uhf
!
!     AO Fock and energy related routines
!
      procedure :: initialize_fock                      => initialize_fock_uhf
      procedure :: construct_ao_spin_fock               => construct_ao_spin_fock_uhf
      procedure :: calculate_uhf_energy                 => calculate_uhf_energy_uhf
      procedure :: update_fock_and_energy_no_cumulative => update_fock_and_energy_no_cumulative_uhf
      procedure :: update_fock_and_energy_cumulative    => update_fock_and_energy_cumulative_uhf
      procedure :: update_fock_and_energy               => update_fock_and_energy_uhf
      procedure :: update_fock_mm                       => update_fock_mm_uhf
      procedure :: update_fock_pcm                      => update_fock_pcm_uhf
      procedure :: set_ao_fock                          => set_ao_fock_uhf
      procedure :: get_ao_fock                          => get_ao_fock_uhf
!
!     AO Density related routines
!
      procedure :: initialize_density                   => initialize_density_uhf
      procedure :: set_initial_ao_density_guess         => set_initial_ao_density_guess_uhf
      procedure :: save_ao_density                      => save_ao_density_uhf
      procedure :: update_ao_density                    => update_ao_density_uhf
      procedure :: form_ao_density                      => form_ao_density_uhf
      procedure :: construct_ao_spin_density            => construct_ao_spin_density_uhf
      procedure :: set_ao_density_to_core_guess         => set_ao_density_to_core_guess_uhf
      procedure :: get_homo_degeneracy                  => get_homo_degeneracy_uhf
      procedure :: get_ao_density_sq                    => get_ao_density_sq_uhf
!
      procedure :: construct_mo_fock                    => construct_mo_fock_uhf
!
!     MO orbital related routines
!
      procedure :: initialize_orbitals                  => initialize_orbitals_uhf
      procedure :: roothan_hall_update_orbitals         => roothan_hall_update_orbitals_uhf
      procedure :: print_orbital_energies               => print_orbital_energies_uhf
      procedure :: print_energy                         => print_energy_uhf
      procedure :: print_orbitals                       => print_orbitals_uhf
      procedure :: save_orbital_coefficients            => save_orbital_coefficients_uhf
      procedure :: read_orbital_coefficients            => read_orbital_coefficients_uhf
      procedure :: save_orbital_energies                => save_orbital_energies_uhf
      procedure :: read_orbital_energies                => read_orbital_energies_uhf
!
!     Gradients and Hessians (todo)
!
      procedure :: get_packed_roothan_hall_gradient     => get_packed_roothan_hall_gradient_uhf
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density_a              => initialize_ao_density_a_uhf
      procedure :: initialize_ao_density_b              => initialize_ao_density_b_uhf
!
      procedure :: initialize_ao_fock_a                 => initialize_ao_fock_a_uhf
      procedure :: initialize_ao_fock_b                 => initialize_ao_fock_b_uhf
!
      procedure :: initialize_orbital_coefficients_a    => initialize_orbital_coefficients_a_uhf
      procedure :: initialize_orbital_coefficients_b    => initialize_orbital_coefficients_b_uhf
!
      procedure :: initialize_orbital_energies_a        => initialize_orbital_energies_a_uhf
      procedure :: initialize_orbital_energies_b        => initialize_orbital_energies_b_uhf
!
      procedure :: destruct_ao_density                  => destruct_ao_density_uhf
      procedure :: destruct_ao_density_a                => destruct_ao_density_a_uhf
      procedure :: destruct_ao_density_b                => destruct_ao_density_b_uhf
!
      procedure :: destruct_fock                        => destruct_fock_uhf
      procedure :: destruct_ao_fock_a                   => destruct_ao_fock_a_uhf
      procedure :: destruct_ao_fock_b                   => destruct_ao_fock_b_uhf
!
      procedure :: destruct_orbital_coefficients        => destruct_orbital_coefficients_uhf
      procedure :: destruct_orbital_coefficients_a      => destruct_orbital_coefficients_a_uhf
      procedure :: destruct_orbital_coefficients_b      => destruct_orbital_coefficients_b_uhf
!
      procedure :: destruct_orbital_energies            => destruct_orbital_energies_uhf
      procedure :: destruct_orbital_energies_a          => destruct_orbital_energies_a_uhf
      procedure :: destruct_orbital_energies_b          => destruct_orbital_energies_b_uhf
!
      procedure :: set_n_mo                              => set_n_mo_uhf
!
   end type uhf
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
      call wf%read_settings()
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
      wf%name_ = 'uhf'
!
      wf%n_ao = wf%system%get_n_aos()
!
      call wf%set_n_mo()
!
      if (wf%fractional_uniform_valence) then
!
         call output%printf('Requested fractional uniform valence. Valence electrons will be  &
                            &distributed evenly in the highest molecular orbitals (if plural).', &
                            pl='minimal', ffs='(/t3,a)')
!
      endif
!
      wf%orbital_coefficients_file = sequential_file('orbital_coefficients')
      wf%orbital_energies_file = sequential_file('orbital_energies')
      wf%restart_file = sequential_file('scf_restart_file')
!
      call wf%initialize_sp_eri_schwarz() 
      call wf%initialize_sp_eri_schwarz_list()
!
      call wf%construct_sp_eri_schwarz()
!
      if(wf%system%mm_calculation)  call wf%initialize_mm_matrices()
      if(wf%system%pcm_calculation) call wf%initialize_pcm_matrices()
!
      wf%frozen_core = .false.
      wf%frozen_hf_mos = .false.
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
      implicit none
!
      class(uhf) :: wf
!
      character(len=*) :: guess
!
      real(dp), dimension(:,:), allocatable :: h_wx
!
      if (trim(guess) == 'sad' .or. trim(guess) == 'SAD') then
!
         call wf%set_ao_density_to_sad()
!
         wf%ao_density_a = (real(wf%n_alpha, kind=dp)/real(wf%n_alpha + wf%n_beta, kind=dp))*wf%ao_density
         wf%ao_density_b = (real(wf%n_beta, kind=dp)/real(wf%n_alpha + wf%n_beta, kind=dp))*wf%ao_density
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
   end subroutine set_initial_ao_density_guess_uhf
!
!
   subroutine construct_mo_fock_uhf(wf)
!!
!!    Construct MO Fock
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
!!    Give notice to user that it does not yet exist.
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      call output%printf('Requested MO transformation of Fock matrix, but this is not yet implemented for (a0).', &
                        pl='minimal', ffs='(/t3,a)', chars=[trim(wf%name_)] )
!
   end subroutine construct_mo_fock_uhf
!
!
   subroutine print_orbital_energies_uhf(wf)
!!
!!    Print orbital energies
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints the current orbital energies to output.
!!
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Use new output%print_vector instead of deprecated print_vector.
!!    Removed indent from input.
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      call output%print_vector('normal', '- Alpha orbital energies', wf%n_ao, wf%orbital_energies_a, &
                              fs='(f16.12)', columns=4)
!
      call output%print_vector('normal', '- Beta orbital energies',  wf%n_ao, wf%orbital_energies_b, &
                              fs='(f16.12)', columns=4)
!
   end subroutine print_orbital_energies_uhf
!
!
   subroutine save_orbital_coefficients_uhf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
!
      call wf%orbital_coefficients_file%open_('write', 'rewind')
!
      call wf%orbital_coefficients_file%write_(wf%orbital_coefficients_a, wf%n_ao*wf%n_mo)
      call wf%orbital_coefficients_file%write_(wf%orbital_coefficients_b, wf%n_ao*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
!
   end subroutine save_orbital_coefficients_uhf
!
!
   subroutine read_orbital_coefficients_uhf(wf)
!!
!!    Read orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      call wf%is_restart_safe('ground state')
!
      call wf%orbital_coefficients_file%open_('read', 'rewind')
!
      call wf%orbital_coefficients_file%read_(wf%orbital_coefficients_a, wf%n_ao*wf%n_mo)
      call wf%orbital_coefficients_file%read_(wf%orbital_coefficients_b, wf%n_ao*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
!
   end subroutine read_orbital_coefficients_uhf
!
!
   subroutine save_orbital_energies_uhf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      call wf%orbital_energies_file%open_('write', 'rewind')
!
      call wf%orbital_energies_file%write_(wf%orbital_energies_a, wf%n_mo)
      call wf%orbital_energies_file%write_(wf%orbital_energies_b, wf%n_mo)
!
      call wf%orbital_energies_file%close_
!
   end subroutine save_orbital_energies_uhf
!
!
   subroutine read_orbital_energies_uhf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      call wf%is_restart_safe('ground state')
!
      call wf%orbital_energies_file%open_('read', 'rewind')
!
      call wf%orbital_energies_file%read_(wf%orbital_energies_a, wf%n_mo)
      call wf%orbital_energies_file%read_(wf%orbital_energies_b, wf%n_mo)
!
      call wf%orbital_energies_file%close_
!
   end subroutine read_orbital_energies_uhf
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
      real(dp), dimension(wf%n_ao*(wf%n_ao - 1)/2, wf%n_densities), intent(inout) :: G
!
      real(dp), dimension(:,:), allocatable :: G_sq
      real(dp), dimension(:), allocatable :: G_pck
      real(dp), dimension(:,:), allocatable :: Po, Pv
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
      call mem%alloc(G_pck, wf%n_ao*(wf%n_ao - 1)/2)
      call mem%alloc(G_sq, wf%n_ao, wf%n_ao)
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
      call mem%dealloc(G_pck, wf%n_ao*(wf%n_ao - 1)/2)
      call mem%dealloc(G_sq, wf%n_ao, wf%n_ao)
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
      wf%fractional_uniform_valence = input%requested_keyword_in_section('fractional uniform valence', 'hf')
!
   end subroutine read_uhf_settings_uhf
!
!
   subroutine initialize_orbitals_uhf(wf)
!!
!!    Initialize orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the arrays associated with the orbital
!!    coefficients. In spin-unrestricted wavefunctions, this
!!    will include alpha and beta coefficients, though these
!!    are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%initialize_orbital_coefficients()
!
      call wf%initialize_orbital_coefficients_a()
      call wf%initialize_orbital_coefficients_b()
!
      call wf%initialize_orbital_energies()
      call wf%initialize_orbital_energies_a()
      call wf%initialize_orbital_energies_b()
!
      wf%orbital_coefficients   = zero
      wf%orbital_coefficients_a = zero
      wf%orbital_coefficients_b = zero
!
      wf%orbital_energies   = zero
      wf%orbital_energies_a = zero
      wf%orbital_energies_b = zero
!
   end subroutine initialize_orbitals_uhf
!
!
   subroutine initialize_density_uhf(wf)
!!
!!    Initialize density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO density (or densities).
!!    In spin-unrestricted wavefunctions, this alpha and beta densities,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%initialize_ao_density()
!
      call wf%initialize_ao_density_a()
      call wf%initialize_ao_density_b()
!
      wf%ao_density   = zero
      wf%ao_density_a = zero
      wf%ao_density_b = zero
!
   end subroutine initialize_density_uhf
!
!
   subroutine initialize_fock_uhf(wf)
!!
!!    Initialize Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%initialize_ao_fock()
!
      call wf%initialize_ao_fock_a()
      call wf%initialize_ao_fock_b()
!
      wf%ao_fock   = zero
      wf%ao_fock_a = zero
      wf%ao_fock_b = zero
!
   end subroutine initialize_fock_uhf
!
!
   subroutine destruct_fock_uhf(wf)
!!
!!    Destruct Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%destruct_ao_fock()
!
      call wf%destruct_ao_fock_a()
      call wf%destruct_ao_fock_b()
!
   end subroutine destruct_fock_uhf
!
!
   subroutine update_fock_and_energy_no_cumulative_uhf(wf, h_wx)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified for QM/MM by Tommaso Giovannini, July 2019 
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!      
      real(dp), dimension(:, :), allocatable            :: h_wx_eff
!
      call mem%alloc(h_wx_eff, wf%n_ao,wf%n_ao)
!      
      h_wx_eff = h_wx
! 
      if(wf%system%mm_calculation.and.wf%system%mm%forcefield.eq.'non-polarizable') &
         call wf%update_h_wx_mm(h_wx_eff)
!         
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha', h_wx_eff)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta', h_wx_eff)
!
      call mem%dealloc(h_wx_eff, wf%n_ao, wf%n_ao) 
!      
      if(wf%system%mm_calculation.and.wf%system%mm%forcefield.ne.'non-polarizable') &
         call wf%update_fock_mm()
!         
      if(wf%system%pcm_calculation) call wf%update_fock_pcm()
!      
      call wf%calculate_uhf_energy(h_wx)
!
   end subroutine update_fock_and_energy_no_cumulative_uhf
!
!
   subroutine update_fock_and_energy_uhf(wf, h_wx, prev_ao_density)
!!
!!    Update Fock and energy wrapper
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Wrapper for cumulative or no_cumulative subroutines
!!    depending on the path
!!
      implicit none
!
      class(uhf) :: wf
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
   end subroutine update_fock_and_energy_uhf
!
!
   subroutine set_ao_fock_uhf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
      real(dp), dimension(:), allocatable :: F_sigma
!
      call mem%alloc(F_sigma, wf%n_ao*(wf%n_ao + 1)/2)
!
!     Alpha Fock
!
      call dcopy(wf%n_ao*(wf%n_ao + 1)/2, F, 1, F_sigma, 1)
      call squareup(F_sigma, wf%ao_fock_a, wf%n_ao)
!
!     Beta Fock
!
      call dcopy(wf%n_ao*(wf%n_ao + 1)/2, F(1, 2), 1, F_sigma, 1)
      call squareup(F_sigma, wf%ao_fock_b, wf%n_ao)
!
      call mem%dealloc(F_sigma, wf%n_ao*(wf%n_ao + 1)/2)
!
   end subroutine set_ao_fock_uhf
!
!
   subroutine get_ao_fock_uhf(wf, F)
!!
!!    Get AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns the AO Fock
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
      real(dp), dimension(:), allocatable :: F_sigma
!
      call mem%alloc(F_sigma, wf%n_ao*(wf%n_ao + 1)/2)
!
!     Alpha Fock
!
      call packin(F_sigma, wf%ao_fock_a, wf%n_ao)
      call dcopy(wf%n_ao*(wf%n_ao + 1)/2, F_sigma, 1, F, 1)
!
!     Beta Fock
!
      call packin(F_sigma, wf%ao_fock_b, wf%n_ao)
      call dcopy(wf%n_ao*(wf%n_ao + 1)/2, F_sigma, 1, F(1, 2), 1)
!
      call mem%dealloc(F_sigma, wf%n_ao*(wf%n_ao + 1)/2)
!
   end subroutine get_ao_fock_uhf
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
      wf%orbital_coefficients_a = zero
      wf%orbital_energies_a     = zero
!
      wf%orbital_coefficients_b = zero
      wf%orbital_energies_b     = zero
!
      call wf%do_roothan_hall(wf%ao_fock_a, wf%orbital_coefficients_a, wf%orbital_energies_a)
      call wf%do_roothan_hall(wf%ao_fock_b, wf%orbital_coefficients_b, wf%orbital_energies_b)
!
   end subroutine roothan_hall_update_orbitals_uhf
!
!
   subroutine print_energy_uhf(wf)
!!
!!    Print wavefunction summary
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
         call output%printf('HOMO-LUMO gap (alpha):     (f19.12)', pl='minimal', fs='(/t6,a)', reals=[homo_lumo_gap_a])
!
      endif 
!
      if (wf%n_beta > 0 .and. wf%n_beta < wf%n_mo) then 
!
         homo_lumo_gap_b = wf%orbital_energies_b(wf%n_beta + 1) - wf%orbital_energies_b(wf%n_beta)
         call output%printf('HOMO-LUMO gap (beta):      (f19.12)', pl='minimal', fs='(t6,a)',  reals=[homo_lumo_gap_b])
!
      endif
!
      nuclear_repulsion = wf%system%get_nuclear_repulsion()
!
      if(wf%system%mm_calculation.and.wf%system%mm%forcefield.eq.'non-polarizable') then
!
         nuclear_repulsion = nuclear_repulsion + wf%system%get_nuclear_repulsion_mm()
!
      endif
!
      call output%printf('Nuclear repulsion energy:  (f19.12)', pl='minimal', fs='(t6,a)',  reals=[nuclear_repulsion])
      call output%printf('Electronic energy:         (f19.12)', pl='minimal', fs='(t6,a)',  reals=[wf%energy - nuclear_repulsion])
      call output%printf('Total energy:              (f19.12)', pl='minimal', fs='(t6,a)',  reals=[wf%energy])
!
      if(wf%system%mm_calculation) call wf%print_energy_mm()
      if(wf%system%pcm_calculation) call wf%print_energy_pcm()
!      
!
   end subroutine print_energy_uhf
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
   subroutine save_ao_density_uhf(wf)
!!
!!    Save AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Save the AO density (or densities, if unrestricted) based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(uhf) :: wf
!
      type(sequential_file) :: ao_density_file
      type(sequential_file) :: ao_density_file_a
      type(sequential_file) :: ao_density_file_b
!
      ao_density_file   = sequential_file('ao_density')
      ao_density_file_a = sequential_file('ao_density_a')
      ao_density_file_b = sequential_file('ao_density_b')
!
      call ao_density_file%open_('write', 'rewind')
      call ao_density_file%write_(wf%ao_density, wf%n_ao*wf%n_ao)
      call ao_density_file%close_
!
      call ao_density_file_a%open_('write', 'rewind')
      call ao_density_file_a%write_(wf%ao_density_a, wf%n_ao*wf%n_ao)
      call ao_density_file_a%close_
!
      call ao_density_file_b%open_('write', 'rewind')
      call ao_density_file_b%write_(wf%ao_density_b, wf%n_ao*wf%n_ao)
      call ao_density_file_b%close_
!
   end subroutine save_ao_density_uhf
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
      wf%ao_fock = h_wx
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies)
!
      wf%orbital_coefficients_a = wf%orbital_coefficients
      wf%orbital_coefficients_b = wf%orbital_coefficients
!
      wf%orbital_energies_a = wf%orbital_energies
      wf%orbital_energies_b = wf%orbital_energies
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
   subroutine get_ao_density_sq_uhf(wf, D)
!!
!!    Get AO density squared
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Returns the unpacked AO density matrix D
!!    (or density matrices in descendants, see overwriting routines)
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(inout) :: D
!
      call dcopy(wf%n_ao**2, wf%ao_density_a, 1, D, 1)
      call dcopy(wf%n_ao**2, wf%ao_density_b, 1, D(1, 2), 1)
!
   end subroutine get_ao_density_sq_uhf
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
            wf%ao_density_a = zero
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
            wf%ao_density_a = zero
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
            wf%ao_density_b = zero
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
            wf%ao_density_b = zero
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
      wf%ao_density = wf%ao_density_a + wf%ao_density_b
!
   end subroutine form_ao_density_uhf
!
!
   subroutine update_fock_and_energy_cumulative_uhf(wf, prev_ao_density, h_wx)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in) :: prev_ao_density
!
      logical :: cumulative
!
!     Here, the previous AO density is sent as [D_a D_b],
!     where each is full square
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density_a, 1)
      call daxpy(wf%n_ao**2, -one, prev_ao_density(1, 2), 1, wf%ao_density_b, 1)
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
      call daxpy(wf%n_ao**2, -one, prev_ao_density(1, 2), 1, wf%ao_density, 1)
!
      cumulative = .true.
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha', &
                      h_wx, cumulative)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta', &
                      h_wx, cumulative)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density_a, 1)
      call daxpy(wf%n_ao**2, one, prev_ao_density(1, 2), 1, wf%ao_density_b, 1)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)
      call daxpy(wf%n_ao**2, one, prev_ao_density(1, 2), 1, wf%ao_density, 1)
!
      call wf%calculate_uhf_energy(h_wx)
!
   end subroutine update_fock_and_energy_cumulative_uhf
!
!
   subroutine construct_ao_spin_fock_uhf(wf, D, D_sigma, sigma, &
                     h_wx, cumulative)
!!
!!    Construct AO spin Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    The routine computes the alpha or beta Fock matrix, depending
!!    on the value of the spin 'sigma' (='alpha' or 'beta'):
!!
!!       F_αβ^alpha = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^alpha
!!       F_αβ^beta  = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^beta
!!
!!    Here the superscript refers to the spin function, while the subscripts
!!    are AO indices. In contrast to the restricted routine, this one does
!!    not calculate the energy - a separate call is required to get the
!!    unrestricted Hartree-Fock energy.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_sigma
!
      character(len=*), intent(in) :: sigma
!
!
      logical, intent(in), optional :: cumulative
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      integer :: thread = 0, n_threads = 1
      logical :: local_cumulative
!
      real(dp), dimension(:,:), allocatable :: F, sp_density_schwarz
!
      real(dp) :: coulomb_thr, exchange_thr, precision_thr    ! Actual thresholds
!
      integer :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: scaled_D_sigma ! = 2 * D_sigma
!
!     Set thresholds to ignore Coulomb and exchange terms,
!     as well as the desired Libint integral precision,
!     and determine whether the construction should be
!     cumulative or not
!
      coulomb_thr   = wf%coulomb_threshold
      exchange_thr  = wf%exchange_threshold
      precision_thr = wf%libint_epsilon
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
!     Compute number of significant ERI shell pairs (the Fock construction
!     only loops over these shell pairs) and the maximum element
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
      max_eri_schwarz = get_abs_max(wf%sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
      call wf%construct_sp_density_schwarz(sp_density_schwarz, D)
      max_D_schwarz = get_abs_max(sp_density_schwarz, wf%system%n_s**2)
!
!$      n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      F = zero
!
      call wf%construct_coulomb_ao_G(F, D, n_threads, max_D_schwarz, max_eri_schwarz,         &
                                                sp_density_schwarz,                           &
                                                n_sig_sp, coulomb_thr, precision_thr,         &
                                                wf%system%shell_limits)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
      scaled_D_sigma = two*D_sigma
!
      call wf%construct_sp_density_schwarz(sp_density_schwarz, scaled_D_sigma)
      max_D_schwarz = get_abs_max(sp_density_schwarz, wf%system%n_s**2)
!
      call wf%construct_exchange_ao_G(F, scaled_D_sigma, n_threads, max_D_schwarz, max_eri_schwarz, &
                                                   sp_density_schwarz,   &
                                                   n_sig_sp, exchange_thr, precision_thr,                 &
                                                   wf%system%shell_limits)
!
      call mem%dealloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
      call mem%dealloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
!
!     Add the accumulated Fock matrix F into the correct Fock matrix
!     (i.e., either the alpha or beta Fock matrix )
!
      if (trim(sigma) == 'alpha') then
!
         if (.not. local_cumulative) wf%ao_fock_a = zero
         do thread = 1, n_threads
!
            call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock_a, 1)
!
         enddo
!
         call symmetric_sum(wf%ao_fock_a, wf%n_ao)
         wf%ao_fock_a = wf%ao_fock_a*half
!
         if (.not. local_cumulative) wf%ao_fock_a = wf%ao_fock_a + h_wx
!
      elseif (trim(sigma) == 'beta') then
!
         if (.not. local_cumulative) wf%ao_fock_b = zero
         do thread = 1, n_threads
!
            call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock_b, 1)
!
         enddo
!
         call symmetric_sum(wf%ao_fock_b, wf%n_ao)
         wf%ao_fock_b = wf%ao_fock_b*half
!
         if (.not. local_cumulative) wf%ao_fock_b = wf%ao_fock_b + h_wx
!
      else
!
         call output%error_msg('Did not recognize spin variable in construct_ao_fock:' // trim(sigma))
!
      endif
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads)
!
   end subroutine construct_ao_spin_fock_uhf
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
   subroutine initialize_ao_density_a_uhf(wf)
!!
!!    Initialize AO density alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%ao_density_a)) call mem%alloc(wf%ao_density_a, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_a_uhf
!
!
   subroutine destruct_ao_density_uhf(wf)
!!
!!    Destruct AO density
!!    Written by Linda Goletto, Oct 2019
!!    Based on the work of Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO density matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%hf%destruct_ao_density()
!
      call wf%destruct_ao_density_a()
      call wf%destruct_ao_density_b()
!
   end subroutine destruct_ao_density_uhf
!
!
   subroutine destruct_ao_density_a_uhf(wf)
!!
!!    Destruct AO density alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%ao_density_a)) call mem%dealloc(wf%ao_density_a, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_a_uhf
!
!
   subroutine initialize_ao_density_b_uhf(wf)
!!
!!    Initialize AO density beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%ao_density_b)) call mem%alloc(wf%ao_density_b, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_b_uhf
!
!
   subroutine destruct_ao_density_b_uhf(wf)
!!
!!    Destruct AO density beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%ao_density_b)) call mem%dealloc(wf%ao_density_b, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_b_uhf
!
!
   subroutine initialize_ao_fock_a_uhf(wf)
!!
!!    Initialize AO Fock alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%ao_fock_a)) call mem%alloc(wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_fock_a_uhf
!
!
   subroutine destruct_ao_fock_a_uhf(wf)
!!
!!    Destruct AO Fock alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%ao_fock_a)) call mem%dealloc(wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_a_uhf
!
!
   subroutine initialize_ao_fock_b_uhf(wf)
!!
!!    Initialize AO Fock beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%ao_fock_b)) call mem%alloc(wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_fock_b_uhf
!
!
   subroutine destruct_ao_fock_b_uhf(wf)
!!
!!    Destruct AO Fock beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%ao_fock_b)) call mem%dealloc(wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_b_uhf
!
!
   subroutine destruct_orbital_coefficients_uhf(wf)
!!
!!    Destruct orbital coefficients
!!    Written by Linda Goletto, Oct 2019
!!    Based on the work of Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the orbital coefficients matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%hf%destruct_orbital_coefficients()
!
      call wf%destruct_orbital_coefficients_a()
      call wf%destruct_orbital_coefficients_b()
!
   end subroutine destruct_orbital_coefficients_uhf
!
!
   subroutine initialize_orbital_coefficients_a_uhf(wf)
!!
!!    Initialize orbital coefficients alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%orbital_coefficients_a)) call mem%alloc(wf%orbital_coefficients_a, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_a_uhf
!
!
   subroutine destruct_orbital_coefficients_a_uhf(wf)
!!
!!    Destruct orbital coefficients alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%orbital_coefficients_a)) call mem%dealloc(wf%orbital_coefficients_a, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_a_uhf
!
!
   subroutine initialize_orbital_coefficients_b_uhf(wf)
!!
!!    Initialize orbital coefficients beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%orbital_coefficients_b)) call mem%alloc(wf%orbital_coefficients_b, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_b_uhf
!
!
   subroutine destruct_orbital_coefficients_b_uhf(wf)
!!
!!    Destruct orbital coefficients beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%orbital_coefficients_b)) call mem%dealloc(wf%orbital_coefficients_b, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_b_uhf
!
!
   subroutine initialize_orbital_energies_a_uhf(wf)
!!
!!    Initialize orbital energies alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%orbital_energies_a)) call mem%alloc(wf%orbital_energies_a, wf%n_mo)
!
   end subroutine initialize_orbital_energies_a_uhf
!
!
   subroutine destruct_orbital_energies_uhf(wf)
!!
!!    Destruct orbital energies
!!    Written by Linda Goletto, Oct 2019
!!    Based on the work of Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the orbital energies array (or arrays).
!!    In spin-unrestricted wavefunctions, this alpha and beta arrays,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
      call wf%hf%destruct_orbital_energies()
!
      call wf%destruct_orbital_energies_a()
      call wf%destruct_orbital_energies_b()
!
   end subroutine destruct_orbital_energies_uhf
!
!
   subroutine destruct_orbital_energies_a_uhf(wf)
!!
!!    Initialize orbital energies alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%orbital_energies_a)) call mem%dealloc(wf%orbital_energies_a, wf%n_mo)
!
   end subroutine destruct_orbital_energies_a_uhf
!
!
   subroutine initialize_orbital_energies_b_uhf(wf)
!!
!!    Initialize orbital energies beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (.not. allocated(wf%orbital_energies_b)) call mem%alloc(wf%orbital_energies_b, wf%n_mo)
!
   end subroutine initialize_orbital_energies_b_uhf
!
!
   subroutine destruct_orbital_energies_b_uhf(wf)
!!
!!    Destruct orbital energies beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      if (allocated(wf%orbital_energies_b)) call mem%dealloc(wf%orbital_energies_b, wf%n_mo)
!
   end subroutine destruct_orbital_energies_b_uhf
!
!
   subroutine update_fock_mm_uhf(wf)
!!
!!    Update alpha and beta Fock matrices with polarizable QM/MM terms
!!    For now: QM/FQ model (see mm_class and output file)
!!    Written by Tommaso Giovannini, July 2019 for QM/MM
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(:),    allocatable                  :: potential_points
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
                     wf%system%mm%pol_emb_rhs,         &
                     wf%system%mm%n_variables,  &
                     zero,                      &
                     wf%system%mm%pol_emb_lhs,            &
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
         wf%ao_fock_a = wf%ao_fock_a + half*wf%pol_emb_fock
         wf%ao_fock_b = wf%ao_fock_b + half*wf%pol_emb_fock
!
!
         call output%print_matrix('debug', 'Total QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM Density - Spin alpha', &
                                  wf%ao_density_a, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM Density - Spin beta', &
                                  wf%ao_density_b, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'FQ Fock - Spin alpha + beta', &
                                  wf%pol_emb_fock, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM/FQ Fock - Spin alpha', wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM/FQ Fock - Spin beta', wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
!
         call mem%dealloc(potential_points, wf%system%mm%n_atoms)
!         
      else
!      
         call output%error_msg('The only available polarizable force field is fq')
!         
      endif
!
!
   end subroutine update_fock_mm_uhf
!
!
   subroutine update_fock_pcm_uhf(wf)
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
      class(uhf) :: wf
!
      character(kind=c_char, len=*), parameter :: mep_lbl = 'NucMEP'
      character(kind=c_char, len=*), parameter :: asc_lbl = 'NucASC'
      integer :: i
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
      call pcmsolver_set_surface_function(wf%system%pcm%pcm_context, int(wf%system%pcm%n_tesserae,kind=c_int), &
                                         -wf%system%pcm%pcm_rhs, pcmsolver_fstring_to_carray(mep_lbl))
!                                          
      call pcmsolver_compute_asc(wf%system%pcm%pcm_context, &
                                 pcmsolver_fstring_to_carray(mep_lbl), &
                                 pcmsolver_fstring_to_carray(asc_lbl), &
                                 irrep=0_c_int)
!                                 
      call pcmsolver_get_surface_function(wf%system%pcm%pcm_context, int(wf%system%pcm%n_tesserae,kind=c_int), &
                                          wf%system%pcm%charges, pcmsolver_fstring_to_carray(asc_lbl))
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
      wf%ao_fock_a = wf%ao_fock_a + half*wf%pcm_fock
      wf%ao_fock_b = wf%ao_fock_b + half*wf%pcm_fock
!
      call output%print_matrix('debug', 'Total QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM Density - Spin alpha', &
                               wf%ao_density_a, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM Density - Spin beta', &
                               wf%ao_density_b, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'PCM Fock - Spin alpha + beta', &
                               wf%pcm_fock, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM/PCM Fock - Spin alpha', wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM/PCM Fock - Spin beta', wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
!
   end subroutine update_fock_pcm_uhf
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
      call output%printf('- Cholesky decomposition of AO overlap to get linearly independent orbitals:', &
                         pl='n', fs='(/t3,a)', ll=100)
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap()
!
      wf%n_densities = 2
!
      call wf%determine_n_alpha_and_n_beta()
!
      call output%printf('Number of alpha electrons:        (i8)', ints=[wf%n_alpha], fs='(/t6,a)',pl='minimal')
      call output%printf('Number of beta electrons:         (i8)', ints=[wf%n_beta], fs='(t6,a)',pl='minimal')
!
      call output%printf('Number of virtual alpha orbitals: (i8)', ints=[wf%n_mo - wf%n_alpha], fs='(t6,a)',pl='minimal')
      call output%printf('Number of virtual beta orbitals:  (i8)', ints=[wf%n_mo - wf%n_beta], fs='(t6,a)',pl='minimal')
!
      call output%printf('Number of molecular orbitals:     (i8)', ints=[wf%n_mo], fs='(t6,a)',pl='minimal')
      call output%printf('Number of atomic orbitals:        (i8)', ints=[wf%n_ao], fs='(t6,a)',pl='minimal')
!
      if (wf%n_mo .lt. wf%n_ao) &
         call output%printf('Removed (i0) AOs due to linear dep.', ints=[wf%n_ao - wf%n_mo], fs='(/t6,a)', pl='minimal')
!
   end subroutine set_n_mo_uhf
!
!
   subroutine print_orbitals_uhf(wf, orbital_list)
!!
!!    Print orbitals
!!    Written by Eirik F. Kjønstad and Tor S. Haugland, Oct 2019
!!
!!    Prints the orbitals with atom & orbital information given.
!!
!!    orbital_list: A list of integers determining which MO to print
!!                     [1,3,20] -> print MO 1, 3 and 20
!!
      implicit none
!
      class(uhf),            intent(in)           :: wf
      integer, dimension(:), intent(in), optional :: orbital_list
!
      integer, dimension(:), allocatable :: orbital_list_alpha, orbital_list_beta
!
      integer :: n_orbitals_alpha, n_orbitals_beta
      integer :: mo, first_mo, last_mo
!
!     Set defaults
!
      if (.not. present(orbital_list)) then
!
!        Default: orbital_list_alpha = [n_alpha - 4, ..., n_alpha + 5]
!
         first_mo = max(1, wf%n_alpha - 4)
         last_mo = min(wf%n_mo, wf%n_alpha + 5)
!
         n_orbitals_alpha = last_mo - first_mo + 1
!
         call mem%alloc(orbital_list_alpha, n_orbitals_alpha)
!
         do mo = 1, n_orbitals_alpha
!
            orbital_list_alpha(mo) = first_mo + mo - 1
!
         enddo
!
!        Default: orbital_list_beta = [n_beta - 4, ..., n_beta + 5]
!
         first_mo = max(1, wf%n_beta - 4)
         last_mo = min(wf%n_mo, wf%n_beta + 5)
!
         n_orbitals_beta = last_mo - first_mo + 1
!
         call mem%alloc(orbital_list_beta, n_orbitals_beta)
!
         do mo = 1, n_orbitals_beta
!
            orbital_list_beta(mo) = first_mo + mo - 1
!
         enddo
!
      else
!
!        Input: orbital_list
!
         n_orbitals_alpha = size(orbital_list)
         n_orbitals_beta = size(orbital_list)
!
         call mem%alloc(orbital_list_alpha, n_orbitals_alpha)
         call mem%alloc(orbital_list_beta,  n_orbitals_beta)
!
         orbital_list_alpha = orbital_list
         orbital_list_beta = orbital_list
!
      endif
!
      call output%printf('- Printing alpha molecular orbitals ((i0) from total (i0))', pl='normal', fs='(/t3,a)', &
                        ints=[n_orbitals_alpha, wf%n_mo])
!
      call wf%print_orbitals_from_coefficients(orbital_list_alpha, wf%orbital_coefficients_a)
!
      call output%printf('- Printing beta molecular orbitals ((i0) from total (i0))', pl='normal', fs='(/t3,a)', &
                        ints=[n_orbitals_beta, wf%n_mo])
!
      call wf%print_orbitals_from_coefficients(orbital_list_beta, wf%orbital_coefficients_b)
!
      call mem%dealloc(orbital_list_alpha, n_orbitals_alpha)
      call mem%dealloc(orbital_list_beta,  n_orbitals_beta)
!
   end subroutine print_orbitals_uhf
!
end module uhf_class
