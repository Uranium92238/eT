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
   use reordering
   use array_utilities
   use array_analysis
   use interval_class
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
      procedure :: prepare                           => prepare_uhf
!
      procedure :: determine_n_alpha_and_n_beta      => determine_n_alpha_and_n_beta_uhf
      procedure :: read_settings                     => read_settings_uhf
      procedure :: read_uhf_settings                 => read_uhf_settings_uhf
      procedure :: print_wavefunction_summary        => print_wavefunction_summary_uhf
!
!     AO Fock and energy related routines 
!
      procedure :: initialize_fock                   => initialize_fock_uhf
      procedure :: destruct_fock                     => destruct_fock_uhf 
      procedure :: construct_ao_spin_fock            => construct_ao_spin_fock_uhf
      procedure :: calculate_uhf_energy              => calculate_uhf_energy_uhf
      procedure :: update_fock_and_energy            => update_fock_and_energy_uhf
      procedure :: update_fock_and_energy_cumulative => update_fock_and_energy_cumulative_uhf
      procedure :: set_ao_fock                       => set_ao_fock_uhf
      procedure :: get_ao_fock                       => get_ao_fock_uhf
!
!     AO Density related routines 
!
      procedure :: initialize_density                => initialize_density_uhf
      procedure :: set_initial_ao_density_guess      => set_initial_ao_density_guess_uhf
      procedure :: save_ao_density                   => save_ao_density_uhf
      procedure :: update_ao_density                 => update_ao_density_uhf
      procedure :: form_ao_density                   => form_ao_density_uhf
      procedure :: construct_ao_spin_density         => construct_ao_spin_density_uhf
      procedure :: set_ao_density_to_core_guess      => set_ao_density_to_core_guess_uhf
      procedure :: get_homo_degeneracy               => get_homo_degeneracy_uhf
      procedure :: get_ao_density_sq                 => get_ao_density_sq_uhf
!
      procedure :: construct_mo_fock                 => construct_mo_fock_uhf
!
!     MO orbital related routines 
!
      procedure :: initialize_orbitals               => initialize_orbitals_uhf
      procedure :: roothan_hall_update_orbitals      => roothan_hall_update_orbitals_uhf
      procedure :: print_orbital_energies            => print_orbital_energies_uhf
      procedure :: save_orbital_coefficients         => save_orbital_coefficients_uhf
      procedure :: read_orbital_coefficients         => read_orbital_coefficients_uhf
      procedure :: save_orbital_energies             => save_orbital_energies_uhf
      procedure :: read_orbital_energies             => read_orbital_energies_uhf
!
!     Gradients and Hessians (todo)
!
      procedure :: get_packed_roothan_hall_gradient  => get_packed_roothan_hall_gradient_uhf
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density_a           => initialize_ao_density_a_uhf
      procedure :: initialize_ao_density_b           => initialize_ao_density_b_uhf
!
      procedure :: initialize_ao_fock_a              => initialize_ao_fock_a_uhf
      procedure :: initialize_ao_fock_b              => initialize_ao_fock_b_uhf
!
      procedure :: initialize_orbital_coefficients_a => initialize_orbital_coefficients_a_uhf
      procedure :: initialize_orbital_coefficients_b => initialize_orbital_coefficients_b_uhf
!
      procedure :: initialize_orbital_energies_a     => initialize_orbital_energies_a_uhf
      procedure :: initialize_orbital_energies_b     => initialize_orbital_energies_b_uhf
!
      procedure :: destruct_ao_density_a             => destruct_ao_density_a_uhf
      procedure :: destruct_ao_density_b             => destruct_ao_density_b_uhf
!
      procedure :: destruct_ao_fock_a                => destruct_ao_fock_a_uhf
      procedure :: destruct_ao_fock_b                => destruct_ao_fock_b_uhf
!
      procedure :: destruct_orbital_coefficients_a   => destruct_orbital_coefficients_a_uhf
      procedure :: destruct_orbital_coefficients_b   => destruct_orbital_coefficients_b_uhf
!
      procedure :: destruct_orbital_energies_a       => destruct_orbital_energies_a_uhf
      procedure :: destruct_orbital_energies_b       => destruct_orbital_energies_b_uhf
!
   end type uhf
!
!
contains 
!
!
   subroutine prepare_uhf(wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      wf%name_ = 'UHF'
!
      call wf%system%prepare()
!
      wf%n_ao = wf%system%get_n_aos()
!
      call initialize_coulomb_c()
      call initialize_kinetic_c()
      call initialize_nuclear_c()
      call initialize_overlap_c()
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap() 
!
      wf%n_o         = (wf%system%get_n_electrons())/2
      wf%n_v         = wf%n_mo - wf%n_o
      wf%n_densities = 2
!
      call wf%determine_n_alpha_and_n_beta()
      call wf%read_settings()
!
      if (wf%fractional_uniform_valence) then 
!
         write(output%unit, '(/t3,a)') 'Requested fractional uniform valence. Valence electrons will be'
         write(output%unit, '(t3,a)')  'distributed evenly in the highest molecular orbitals (if plural).'
!
      endif
!
      call wf%orbital_coefficients_file%init('orbital_coefficients', 'sequential', 'unformatted')
      call wf%orbital_energies_file%init('orbital_energies', 'sequential', 'unformatted')
      call wf%restart_file%init('hf_restart_file', 'sequential', 'unformatted')
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
      write(output%unit, '(/t3,a,a,a)') 'Requested MO transformation of Fock matrix, but this ', &
                                          'is not yet implemented for ', wf%name_
!
   end subroutine construct_mo_fock_uhf
!
!
   subroutine print_orbital_energies_uhf(wf, indentation)
!!
!!    Print orbital energies  
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints the current orbital energies to output
!!    in a hopefully readable way.
!!
      implicit none 
!
      class(uhf), intent(in) :: wf 
!
      character(len=*), optional :: indentation
!
      character(len=40) :: indent 
!
      indent = '6'
      if (present(indentation)) indent = trim(indentation)
!
      write(output%unit, '(/t' // trim(indent) // ',a)') 'Alpha orbital energies:'
!
      call print_vector(wf%orbital_energies_a, wf%n_ao, indent)
!
      write(output%unit, '(/t' // trim(indent) // ',a)') 'Beta orbital energies:'
!
      call print_vector(wf%orbital_energies_b, wf%n_ao, indent)
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
      call disk%open_file(wf%orbital_coefficients_file, 'write', 'rewind')
!
      write(wf%orbital_coefficients_file%unit) wf%orbital_coefficients_a 
      write(wf%orbital_coefficients_file%unit) wf%orbital_coefficients_b
!
      call disk%close_file(wf%orbital_coefficients_file)
!
   end subroutine save_orbital_coefficients_uhf
!
!
   subroutine read_orbital_coefficients_uhf(wf)
!!
!!    Save orbital coefficients 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(uhf), intent(inout) :: wf 
!
      call wf%is_restart_safe()
!
      call disk%open_file(wf%orbital_coefficients_file, 'read', 'rewind')
!
      read(wf%orbital_coefficients_file%unit) wf%orbital_coefficients_a 
      read(wf%orbital_coefficients_file%unit) wf%orbital_coefficients_b
!
      call disk%close_file(wf%orbital_coefficients_file)
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
      call disk%open_file(wf%orbital_energies_file, 'write', 'rewind')
!
      write(wf%orbital_energies_file%unit) wf%orbital_energies_a 
      write(wf%orbital_energies_file%unit) wf%orbital_energies_b
!
      call disk%close_file(wf%orbital_energies_file)
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
      call wf%is_restart_safe()
!
      call disk%open_file(wf%orbital_energies_file, 'read', 'rewind')
!
      read(wf%orbital_energies_file%unit) wf%orbital_energies_a 
      read(wf%orbital_energies_file%unit) wf%orbital_energies_b
!
      call disk%close_file(wf%orbital_energies_file)
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
      integer :: n_records, i 
!
      character(len=100) :: line, value 
!
      if (requested_section('hf')) then ! User has requested something 
!
         call move_to_section('hf', n_records)
!
         do i = 1, n_records
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:27) == 'fractional uniform valence:') then
!
               value = line(28:100)
               value = remove_preceding_blanks(value)
!
               if (trim(value) == 'true') then 
!
                  wf%fractional_uniform_valence = .true.
!
               endif
!
               cycle
!
            endif
!
         enddo
!
      endif 
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
   subroutine update_fock_and_energy_uhf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!!
!!    Update Fock and energy 
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
      integer, intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer, dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha',    &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta',     &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)     
!
      call wf%calculate_uhf_energy(h_wx)
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
!!    Set AO Fock 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns the AO Fock
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp), dimension(:,:), intent(inout) :: F ! Packed
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
      call dcopy(wf%n_ao*(wf%n_ao + 1)/2, F_sigma, 1, F(wf%n_ao*(wf%n_ao + 1)/2 + 1, 1), 1)   
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
   subroutine print_wavefunction_summary_uhf(wf)
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
      class(uhf), intent(in) :: wf 
!
      real(dp) :: homo_lumo_gap_a
      real(dp) :: homo_lumo_gap_b
!
      write(output%unit, '(/t3,a,a,a)') '- Summary of ', trim(wf%name_), ' wavefunction energetics (a.u.):'
!
      homo_lumo_gap_a = wf%orbital_energies_a(wf%n_alpha + 1) - wf%orbital_energies_a(wf%n_alpha)
      homo_lumo_gap_b = wf%orbital_energies_b(wf%n_beta + 1) - wf%orbital_energies_b(wf%n_beta)
!
      write(output%unit, '(/t6,a26,f19.12)') 'HOMO-LUMO gap (alpha):    ', homo_lumo_gap_a
      write(output%unit, '(t6,a26,f19.12)')  'HOMO-LUMO gap (beta):     ', homo_lumo_gap_b
      write(output%unit, '(t6,a26,f19.12)')  'Nuclear repulsion energy: ', wf%system%get_nuclear_repulsion()
      write(output%unit, '(t6,a26,f19.12)')  'Electronic energy:        ', wf%energy - wf%system%get_nuclear_repulsion()
      write(output%unit, '(t6,a26,f19.12)')  'Total energy:             ', wf%energy
!
      call wf%print_orbital_energies('3')
!
   end subroutine print_wavefunction_summary_uhf
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
      type(file) :: ao_density 
!
      type(file) :: ao_density_a
      type(file) :: ao_density_b
!
      call ao_density%init('ao_density', 'sequential', 'formatted')
      call ao_density_a%init('ao_density_a', 'sequential', 'formatted')
      call ao_density_b%init('ao_density_b', 'sequential', 'formatted')
!
      call disk%open_file(ao_density, 'readwrite', 'rewind')
      write(ao_density%unit, *) wf%ao_density
      call disk%close_file(ao_density)
!
      call disk%open_file(ao_density_a, 'readwrite', 'rewind')
      write(ao_density_a%unit, *) wf%ao_density_a
      call disk%close_file(ao_density_a)
!
      call disk%open_file(ao_density_b, 'readwrite', 'rewind')
      write(ao_density_b%unit, *) wf%ao_density_b
      call disk%close_file(ao_density_b)
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
      real(dp), dimension(:,:), intent(inout) :: D
!
      call dcopy(wf%n_ao**2, wf%ao_density_a, 1, D, 1)
      call dcopy(wf%n_ao**2, wf%ao_density_b, 1, D(wf%n_ao**2 + 1, 1), 1)
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
      real(dp), parameter :: threshold = 1.0D-6 
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
      if (abs(energies(homo) - energies(homo + 1)) .le. threshold .or. &
          abs(energies(homo) - energies(homo - 1)) .le. threshold) then ! HOMO is degenerate 
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
      endif
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
   subroutine update_fock_and_energy_cumulative_uhf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, prev_ao_density, h_wx)
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
      integer, intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in) :: prev_ao_density
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer, dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
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
                     sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx, cumulative) 
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta', &
                     sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx, cumulative) 
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
                     sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx, cumulative)
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
      integer, intent(in) :: n_s
!
      logical, intent(in), optional :: cumulative 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer, dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
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
      call wf%get_n_sig_eri_sp(n_sig_sp, sp_eri_schwarz)
      max_eri_schwarz = get_abs_max(sp_eri_schwarz, n_s*(n_s + 1)/2)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors 
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(sp_density_schwarz, n_s, n_s)
      call wf%construct_sp_density_schwarz(sp_density_schwarz, D)
      max_D_schwarz = get_abs_max(sp_density_schwarz, n_s**2)
!
!$      n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      F = zero 
!
      call wf%ao_fock_coulomb_construction_loop(F, D, n_threads, max_D_schwarz, max_eri_schwarz,         & 
                                                sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list, &
                                                n_s, n_sig_sp, coulomb_thr, precision_thr,               &
                                                wf%system%shell_limits)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors 
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
      scaled_D_sigma = two*D_sigma 
!
      call wf%construct_sp_density_schwarz(sp_density_schwarz, scaled_D_sigma)
      max_D_schwarz = get_abs_max(sp_density_schwarz, n_s**2)      
!
      call wf%ao_fock_exchange_construction_loop(F, scaled_D_sigma, n_threads, max_D_schwarz, max_eri_schwarz, & 
                                                   sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list,    &
                                                   n_s, n_sig_sp, exchange_thr, precision_thr,                 &
                                                   wf%system%shell_limits)      
!
      call mem%dealloc(sp_density_schwarz, n_s, n_s)
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
end module uhf_class
