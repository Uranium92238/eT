!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module fci_class
!
!!
!! Full Configuration Interaction (fci) class module
!! Written by Enrico Ronca, 2020
!!
!
   use parameters
   use global_out,           only: output
   use global_in,            only: input
   use memory_manager_class, only: mem
   use wavefunction_class,   only: wavefunction
   use timings_class,        only: timings
   use stream_file_class,    only : stream_file
!
   implicit none
!
   type, extends(wavefunction) :: fci
!
      integer :: n_alpha_strings, n_beta_strings, n_determinants
      integer :: n_states
      integer :: multiplicity
      integer :: n_alpha, n_beta
      integer :: n_v_alpha, n_v_beta
      integer :: n_alpha_excitations, n_beta_excitations
!
      type(stream_file), dimension(:), allocatable :: fci_files
!
      integer, dimension(:), allocatable :: alpha_strings, beta_strings
!
!     The excitation maps store information about the string generated by
!     acting with a_p^dagger a_q on some string (e.g., the address of the generated strings)
      integer, dimension(:,:,:), allocatable :: excitation_maps_alpha, excitation_maps_beta
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs
      real(dp), dimension(:,:),     allocatable :: h_pq
!
      real(dp), dimension(:,:,:,:), allocatable :: effective_2e_hamiltonian
!
      real(dp), dimension(:), allocatable :: energies
!
      real(dp), dimension(:,:,:), allocatable :: ci_coefficients
!
      real(dp), dimension(:,:), allocatable :: density
!
   contains
!
      procedure, public :: save_fci_state &
                        => save_fci_state_fci
!
      procedure, public :: read_fci_state &
                        => read_fci_state_fci
!
      procedure, public :: initialize_fci_state_files &
                        => initialize_fci_state_files_fci
!
      procedure, public :: hamiltonian_transformation_no_spin_symmetry &
                        => hamiltonian_transformation_no_spin_symmetry_fci
!
      procedure, public :: calculate_spin_squared &
                        => calculate_spin_squared_fci
!
      procedure, public :: calculate_spin_multiplicities &
                        => calculate_spin_multiplicities_fci
!
      procedure, public :: print_dominant_coefficients_amplitudes &
                        => print_dominant_coefficients_amplitudes_fci
!
      procedure, public :: cleanup &
                        => cleanup_fci
!
      procedure, public :: print_banner &
                        => print_banner_fci
!
      procedure, public :: print_fci_summary
!
      procedure, public :: initialize &
                        => initialize_fci
!
      procedure, public :: get_fci_start_vector  &
                        => get_fci_start_vector_fci
!
      procedure, public :: construct_h_diagonal
!
      procedure, private :: consistency_checks
!
      procedure, private :: initialize_strings
      procedure, private :: destruct_strings
!
      procedure, private :: initialize_excitation_maps
      procedure, private :: destruct_excitation_maps
!
      procedure, private :: initialize_hamiltonian_integrals
      procedure, private :: destruct_hamiltonian_integrals
!
      procedure, private :: initialize_energies
      procedure, private :: destruct_energies
!
      procedure, private :: initialize_ci_coefficients
      procedure, private :: destruct_ci_coefficients
!
      procedure :: initialize_gs_density &
                => initialize_gs_density_fci
      procedure :: destruct_gs_density &
                => destruct_gs_density_fci
!
      procedure, private :: general_fci_preparations
      procedure, private :: set_variables_from_template_wf
      procedure, private :: read_settings
!
      procedure, private :: generate_strings
      procedure, private :: construct_excitation_maps
      procedure, private :: get_occupied_and_virtuals_in_string
      procedure, private, nopass :: calculate_occupied_indices_in_strings
!
      procedure, private :: construct_destruction_strings_and_signs
      procedure, private :: construct_creation_strings_and_signs
!
      procedure, private :: construct_hamiltonian_integrals
!
      procedure, private :: construct_h_pq
      procedure, private :: construct_g_pqrs
      procedure, private :: construct_effective_2e_hamiltonian
!
      procedure, private :: construct_D
      procedure, private :: construct_sigma
!
      procedure, private :: print_dominant_coefficients
      procedure, private :: get_determinant_string
      procedure, private :: sm_sp_expectation_value
!
!     Properties
!
      procedure :: construct_gs_density &
                => construct_gs_density_fci
!
      procedure, private :: construct_density
      procedure, private :: add_alpha_density
      procedure, private :: add_beta_density
!
      procedure :: get_electronic_dipole &
                => get_electronic_dipole_fci
      procedure :: get_electronic_quadrupole &
                => get_electronic_quadrupole_fci
      procedure :: calculate_expectation_value &
                => calculate_expectation_value_fci
!
   end type fci
!
   interface
!
      include "initialize_destruct_fci_interface.F90"
      include "strings_fci_interface.F90"
      include "preconditioning_fci_interface.F90"
      include "integrals_fci_interface.F90"
      include "contract_fci_interface.F90"
      include "spin_operators_fci_interface.F90"
      include "file_handling_fci_interface.F90"
      include "properties_fci_interface.F90"
!
   end interface
!
contains
!
!
   subroutine initialize_fci(wf, template_wf)
!!
!!    Initialize
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'fci'
!
      call wf%general_fci_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      call wf%consistency_checks()
!
      call wf%initialize_strings()
!
      call wf%generate_strings()
!
      call wf%initialize_excitation_maps()
!
      call wf%construct_excitation_maps(wf%n_alpha,               &
                                        wf%n_v_alpha,             &
                                        wf%n_alpha_strings,       &
                                        wf%alpha_strings,         &
                                        wf%excitation_maps_alpha)
!
      call wf%construct_excitation_maps(wf%n_beta,                &
                                        wf%n_v_beta,              &
                                        wf%n_beta_strings,        &
                                        wf%beta_strings,          &
                                        wf%excitation_maps_beta)
!
      call wf%initialize_hamiltonian_integrals()
      call wf%construct_hamiltonian_integrals()
!
      call wf%initialize_energies()
      call wf%initialize_ci_coefficients()
!
      call wf%initialize_fci_state_files()
!
   end subroutine initialize_fci
!
!
   subroutine general_fci_preparations(wf)
!!
!!    General FCI preparations
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci) :: wf
!
      call wf%read_settings()
!
   end subroutine general_fci_preparations
!
!
   subroutine set_variables_from_template_wf(wf, template_wf)
!!
!!    Set variables from template wavefunction
!!    Written by Enrico Ronca, 2020
!!
!!    Transfers the wavefunction variables needed to start a
!!    FCI calculation from either a HF wavefunction
!!
!
      use math_utilities, only : binomial
!
      implicit none
!
      class(fci) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      call wf%prepare_ao_tool(template = template_wf%ao)
      call wf%prepare_embedding()
!
      wf%n_alpha = (wf%multiplicity - 1 + wf%ao%get_n_electrons())/2
      wf%n_beta  = wf%ao%get_n_electrons() - wf%n_alpha
!
      wf%n_o = max(wf%n_alpha,wf%n_beta)
      wf%n_v = template_wf%n_mo - wf%n_o
!
      wf%ao%n = template_wf%ao%n
      wf%n_mo = template_wf%n_mo
!
      wf%n_v_alpha = wf%n_mo - wf%n_alpha
      wf%n_v_beta = wf%n_mo - wf%n_beta
!
!     Number of single excitations:
!     a_p^dagger a_q with q = p (n_o) and with q /= p (n_o * n_v)
!
      wf%n_alpha_excitations = wf%n_alpha + wf%n_alpha*wf%n_v_alpha
      wf%n_beta_excitations  = wf%n_beta  + wf%n_beta *wf%n_v_beta
!
      wf%hf_energy = template_wf%hf_energy
!
      wf%n_alpha_strings = binomial(wf%n_mo, wf%n_alpha)
      wf%n_beta_strings  = binomial(wf%n_mo, wf%n_beta)
!
      wf%n_determinants = wf%n_alpha_strings * wf%n_beta_strings
!
!     Set orbital coefficients and energies
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call dcopy(wf%n_mo, template_wf%orbital_energies, 1, wf%orbital_energies, 1)
!
      call dcopy(wf%ao%n*wf%n_mo, template_wf%orbital_coefficients, 1, wf%orbital_coefficients, 1)
!
      if (template_wf%exists_frozen_fock_terms) then
         call output%error_msg("Frozen core and frozen HF have not been implemented for (a0) yet.", &
                               chars=[wf%name_])
      end if
!
   end subroutine set_variables_from_template_wf
!
!
   subroutine print_banner_fci(wf)
!!
!!    Print banner
!!    Written by Enrico Ronca, 2020
!!
!
      use string_utilities, only : convert_to_uppercase
!
      implicit none
!
      class(fci), intent(in) :: wf
!
      character(len=200) :: name_
!
      name_ = trim(convert_to_uppercase(wf%name_)) // ' wavefunction'
!
      call output%printf('m', ':: (a0)', chars=[name_], fs='(//t3,a)')
      call output%print_separator('m', len_trim(name_) + 6, '=')
!
      call output%printf('m', ' - Number of orbitals:', fs='(/t3,a)')
!
      call output%printf('m', 'Alpha electrons:          (i0)', ints=[wf%n_alpha], fs='(/t6,a)')
      call output%printf('m', 'Virtual alpha orbitals:   (i0)', ints=[wf%n_v_alpha], fs='(t6,a)')
      call output%printf('m', 'Beta electrons:           (i0)', ints=[wf%n_beta], fs='(t6,a)')
      call output%printf('m', 'Virtual beta orbitals:    (i0)', ints=[wf%n_v_beta], fs='(t6,a)')
      call output%printf('m', 'Molecular orbitals:       (i0)', ints=[wf%n_mo], fs='(t6,a)')
!
      call output%printf('m', ' - Number of determinants:', fs='(/t3,a)')
!
      call output%printf('m', 'Alpha Determinants:       (i0)', &
                         ints=[wf%n_alpha_strings], fs='(t6,a)')
!
      call output%printf('m', 'Beta Determinants:        (i0)', &
                         ints=[wf%n_beta_strings], fs='(t6,a)')
!
   end subroutine print_banner_fci
!
!
   subroutine cleanup_fci(wf)
!!
!!    Cleanup
!!    Written by Enrico Ronca , 2020
!!
      implicit none
!
      class(fci) :: wf
!
      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
      call wf%destruct_strings()
      call wf%destruct_excitation_maps()
      call wf%destruct_hamiltonian_integrals()
      call wf%destruct_ci_coefficients()
      call wf%destruct_energies()
      call wf%destruct_gs_density()
!
      deallocate(wf%ao)
!
      call output%printf('v', '- Cleaning up (a0) wave function', &
                         chars=[trim(wf%name_)], fs='(/t3,a)')
!
   end subroutine cleanup_fci
!
!
   subroutine read_settings(wf)
!!
!!    Read settings
!!    Written by Enrico Ronca, Aug 2020
!!
!!    Reads the fci-section of the input file
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call input%get_required_keyword('multiplicity', 'system', wf%multiplicity)
      call input%get_required_keyword('states', 'solver fci', wf%n_states)
!
   end subroutine read_settings
!
!
   subroutine get_fci_start_vector_fci(wf, trial, c, restart)
!!
!!    Get FCI start vectors
!!    Written by Enrico Ronca, 2022
!!
      use array_utilities, only: get_l2_norm
      use array_initialization, only:  zero_array
!
      implicit none
!
      class(fci) :: wf
!
      integer, intent(in) :: trial
!
      real(dp) :: c_norm
!
      real(dp), dimension(wf%n_determinants), intent(out) :: c
!
      logical, intent(in) :: restart
!
      character(len=200) :: guess_wf
!
!     Default guess wavefunction
      guess_wf = 'single determinant'
!
      call input%get_keyword('start guess', 'solver fci', guess_wf)
!
      if (restart) then
!
         if (wf%fci_files(trial)%exists()) then
!
            call output%printf('m', 'Requested restart. Reading in solution from file.', &
                               fs='(/t3,a)')
!
            call wf%read_fci_state(c, trial)
            return
!
         endif
!
         call output%warning_msg('asked for restart of FCI energies, &
                                 &but no vectors found on disk')
!
      endif
!
!     Sets whatever remains using the start vector tool
!
      if (guess_wf == 'single determinant') then
!
         call zero_array(c, wf%n_determinants)
         c(trial) = one
!
      else if (guess_wf == 'random') then
!
         call random_number(c)
!
         c_norm = get_l2_norm(c, wf%n_determinants)
         call dscal(wf%n_determinants, one/c_norm, c, 1)
!
      else
!
         call output%error_msg('Unknown guess for wave function.')
!
      end if
!
   end subroutine get_fci_start_vector_fci
!
!
   subroutine print_dominant_coefficients_amplitudes_fci(wf, c)
!!
!!    Print dominant coefficients amplitudes
!!    Written by Enrico Ronca, 2020
!!
!!    Prints the dominant amplitudes in the FCI vector.
!!
      implicit none
!
      class(fci), intent(in) :: wf
!
      real(dp), dimension(wf%n_determinants), intent(in) :: c
!
      call wf%print_dominant_coefficients(c)
!
   end subroutine print_dominant_coefficients_amplitudes_fci
!
!
   subroutine print_dominant_coefficients(wf, c)
!!
!!    Print dominant coefficients
!!    Written by Enrico Ronca, 2020
!!
      use math_utilities,  only: convert_integer_to_binary_string
      use array_utilities, only: get_n_highest
      use string_utilities, only: remove_delimiter_from_string
      use index_invert, only: invert_compound_index
!
      implicit none
!
      class(fci), intent(in) :: wf
      real(dp), dimension(wf%n_determinants), intent(in) :: c
!
      integer :: n_elements, j
!
      integer,  dimension(:), allocatable :: dominant_indices
      real(dp), dimension(:), allocatable :: dominant_values
!
      integer,  dimension(:), allocatable :: alpha_occupation_vector
      integer,  dimension(:), allocatable :: beta_occupation_vector
      integer,  dimension(:), allocatable :: occupation_vector
!
      integer :: index_alpha_string, index_beta_string
      integer :: alpha_string, beta_string, n_separators
!
      character(len=64) :: determinant_string
!
      n_separators = wf%n_mo + 20
!
!     Sort according to largest contributions
!
      n_elements = min(10, wf%n_determinants)
!
      call mem%alloc(dominant_indices, n_elements)
      call mem%alloc(dominant_values, n_elements)
      call mem%alloc(alpha_occupation_vector, wf%n_mo)
      call mem%alloc(beta_occupation_vector, wf%n_mo)
      call mem%alloc(occupation_vector, wf%n_mo)
!
      call get_n_highest(n_elements, wf%n_determinants, &
                         abs(c), dominant_values, dominant_indices)
!
!     Print largest contributions
!
      call output%printf('m', 'Largest CI amplitudes:', fs='(/t6,a)')
      call output%print_separator('m', n_separators, '-', fs='(t6,a)')
!
      do j = 1, n_elements
!
         call invert_compound_index(dominant_indices(J), index_alpha_string, &
                                    index_beta_string, wf%n_alpha_strings, wf%n_beta_strings)
!
         alpha_string = wf%alpha_strings(index_alpha_string)
         beta_string = wf%beta_strings(index_beta_string)
!
         determinant_string = wf%get_determinant_string(alpha_string, beta_string)
!
         call output%printf('m', '(a0) (f19.12)', chars=[trim(determinant_string)], &
                            reals=[c(dominant_indices(j))], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('m', n_separators, '-', fs='(t6,a)')
!
      call mem%dealloc(dominant_indices, n_elements)
      call mem%dealloc(dominant_values, n_elements)
      call mem%dealloc(alpha_occupation_vector, wf%n_mo)
      call mem%dealloc(beta_occupation_vector, wf%n_mo)
      call mem%dealloc(occupation_vector, wf%n_mo)
!
   end subroutine print_dominant_coefficients
!
!
   function get_determinant_string(wf, alpha_index, beta_index) result(string)
!!
!!    Get determinant string
!!    Written by Alexander C. Paul
!!
!!    Makes character string of length n_mo, where the last (most right)
!!    character represents the occupation of the first (lowest energy) orbital.
!!    If the chacater is:
!!    0: unoccupied orbital
!!    a: occupied by an alpha electron
!!    b: occupied by a beta electron
!!    2: doubly occupied
!!
      use math_utilities,  only: convert_integer_to_binary_string
!
      implicit none
!
      class(fci), intent(in) :: wf
!
      integer, intent(in) :: alpha_index, beta_index
!
      character(len=64) :: string, temp
!
      integer,  dimension(:), allocatable :: alpha_occupation_vector
      integer,  dimension(:), allocatable :: beta_occupation_vector
      integer :: length, i
!
      call mem%alloc(alpha_occupation_vector, wf%n_mo)
      call mem%alloc(beta_occupation_vector, wf%n_mo)
!
      call convert_integer_to_binary_string(alpha_occupation_vector, alpha_index, wf%n_mo)
      call convert_integer_to_binary_string(beta_occupation_vector, beta_index, wf%n_mo)
!
      length = 0
      do i = 1, wf%n_mo
!
         if (alpha_occupation_vector(i) + beta_occupation_vector(i) == 2) then
!
            write(temp, '(a,i1)') string(1:length), 2
!
         else if (alpha_occupation_vector(i) == 1) then
!
            write(temp, '(a,a1)') string(1:length), 'a'
!
         else if (beta_occupation_vector(i) == 1) then
!
            write(temp, '(a,a1)') string(1:length), 'b'
!
         else
!
            write(temp, '(a,i1)') string(1:length), 0
!
         end if
!
         length = length + 1
         string = temp
!
      end do
!
      call mem%dealloc(alpha_occupation_vector, wf%n_mo)
      call mem%dealloc(beta_occupation_vector, wf%n_mo)
!
   end function get_determinant_string
!
!
   subroutine print_fci_summary(wf)
!!
!!    Print summary
!!    Written by Enrico Ronca, 2022
!!
!!    Prints summary of FCI. Lists the dominant amplitudes and
!!    the energies.
!!
      implicit none
!
      class(fci), intent(in) :: wf
!
      integer :: state
!
      real(dp), dimension(:), allocatable :: multiplicities
!
      call mem%alloc(multiplicities, wf%n_states)
      call wf%calculate_spin_multiplicities(multiplicities)
!
!     Print FCI vectors
!
      call output%printf('n', '- CI amplitudes:', fs='(/t3,a)')
!
      do state = 1, wf%n_states
!
         call output%printf('n', 'Electronic state nr. (i0)', ints=[state], fs='(/t6,a)')
!
         call output%printf('n', 'Energy (Hartree):  (f19.12)', &
                            reals=[wf%energies(state)], fs='(/t6,a)')
         call output%printf('n', 'Spin Multiplicity: (f19.12)', &
                            reals=[multiplicities(state)], fs='(t6,a)')
!
         call wf%print_dominant_coefficients_amplitudes(wf%ci_coefficients(:,:,state))
!
      enddo
!
!     Print FCI energies
!
      call output%printf('m', '- Electronic energies:', fs='(/t3,a)')
!
      call output%printf('m', 'Energy', fs='(/t39,a)')
      call output%print_separator('m', 42, '-', fs='(t27,a)')
      call output%printf('m', ' State                (Hartree)             (eV)', fs='(t6,a)')
      call output%print_separator('m', 63, '-', fs='(t6,a)')
!
      do state = 1, wf%n_states
!
         call output%printf('m', '(i4)             (f19.12)   (f19.12)',         &
                            ints=[state], reals=[wf%energies(state),   &
                            wf%energies(state)*Hartree_to_eV], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('m', 63, '-', fs='(t6,a)')
      call output%printf('m', 'eV/Hartree (CODATA 2014): (f11.8)', &
                         reals=[Hartree_to_eV], fs='(t6,a)')
!
      call mem%dealloc(multiplicities, wf%n_states)
!
   end subroutine print_fci_summary
!
!
   subroutine consistency_checks(wf)
!!
!!    Consistency check for FCI
!!    Written by Enrico Ronca, 2022
!!
!!    Check whether the numeber of orbitals is smaller than the precision of the integers
!!    to build the strings and if there are electrons in the system.
!!
      implicit none
!
      class(fci), intent(in) :: wf
!
      if (wf%n_mo > storage_size(wf%n_mo)) &
         call output%error_msg('The number of orbitals is larger than the precision of the integer')
!
      if (wf%n_v_alpha == 0 .or. wf%n_v_beta == 0) &
         call output%error_msg('Currently only possible to run systems with no virtual &
                                 &alpha or beta orbitals.')
!
   end subroutine consistency_checks
!
!
end module fci_class
