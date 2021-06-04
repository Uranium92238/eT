!
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
module ao_tool_class
!
!!
!!    Atomic orbital tool class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2020
!!
!!    Class that stores information of the atomic orbital (AO) basis
!!    and can be used to calculate AO integrals. In particular, it contains
!!    a set of centers (atomic, non-atomic), each of which consists of a set
!!    of shells.
!!
!
   use kinds
   use parameters
   use named_range_class
   use iso_c_binding,         only: c_int
!
   use global_in,             only: input
   use global_out,            only: output
   use memory_manager_class,  only: mem
!
   use atomic_center_class,   only: atomic_center
   use point_charges_class,   only: point_charges
!
   implicit none
!
!  Type definition
!
   type :: ao_tool
!
!     Number of AOs, number of shells, and the maximum shell size
!
      integer, public :: n, n_sh, max_sh_size
!
!     Vector of shell intervals (first, last, length)
!
      type(range_), dimension(:), allocatable, public :: shells
!
!     Various useful index mappings
!
      integer, dimension(:), allocatable, public :: shell_to_center
      integer, dimension(:), allocatable, public :: ao_to_center
      integer, dimension(:), allocatable, public :: ao_to_shell
!
!     One-electron integral (OEI) arrays
!
      real(dp), dimension(:,:), allocatable, public :: h ! Hamiltonian one-electron integral matrix
      real(dp), dimension(:,:), allocatable, public :: s ! Overlap integral matrix
      real(dp), dimension(:,:), allocatable, public :: v ! Coulomb external charges integral matrix
!
!     Cauchy-Schwarz electron repulsion integral (CS-ERI) screenings
!
      integer, public :: n_sig_eri_shp ! Number of AB for which (AB|AB) exceeds ERI cutoff
!
      real(dp), dimension(:,:), allocatable, public :: cs_eri_max         ! Contains (AB|AB)^1/2
      integer, dimension(:,:), allocatable, public  :: cs_eri_max_indices ! Sorting indices
!
      integer, private :: n_centers
      type(atomic_center), dimension(:), allocatable, private :: centers ! AO basis centers
!
      integer, private  :: n_orthonormal_ao ! Orthonormal AOs (OAOs)
!     Threshold for removing linear dependency between AOs n_orthonormal_ao
      real(dp), private :: lindep_threshold
!
      real(dp), dimension(:,:), allocatable, public :: P, L ! AO-to-OAO transformation matrices
!
      integer, private :: n_center_subsets
      type(named_range), dimension(:), allocatable, private :: center_subsets
!
      real(dp), private :: libint_epsilon ! Default Libint ERI precision
!
      real(dp), private :: eri_cutoff
      real(dp), private :: oei_cutoff
!
      character(len=200) :: basis_type_ ! standard, spherical, or Gaussian
!
      integer :: charge
!
   contains
!
      procedure, public :: initialize &
                        => initialize_ao_tool
!
      procedure, public :: initialize_ao_tool_from_template
!
      generic :: initialize_external_charges &
              => initialize_external_charges_ao_tool, &
                 initialize_external_charges_from_point_charges_ao_tool
!
      procedure, nopass, public :: initialize_external_charges_ao_tool
      procedure, nopass, public :: initialize_external_charges_from_point_charges_ao_tool
!
      procedure, nopass, public :: initialize_external_unit_charges &
                                => initialize_external_unit_charges_ao_tool
!
      procedure, public :: get_oei &
                        => get_oei_ao_tool ! OEI = one-electron integral
!
      procedure, public :: get_oei_1der &
                        => get_oei_1der_ao_tool ! 1der = first derivative
!
      procedure, public :: get_eri &
                        => get_eri_ao_tool ! ERI = electron repulsion integral
!
      procedure, public :: get_eri_1der &
                        => get_eri_1der_ao_tool
!
      procedure, nopass, public :: get_n_oei_components &
                                => get_n_oei_components_ao_tool ! e.g. 3 for dipole (x,y,z)
!
      procedure, public :: initialize_oei &
                        => initialize_oei_ao_tool
!
      procedure, public :: construct_stored_oei &
                        => construct_stored_oei_ao_tool
!
      procedure, public :: set_atomic_centers &
                        => set_atomic_centers_ao_tool
!
      procedure, public :: get_sad_guess &
                        => get_sad_guess_ao_tool ! SAD = superposition of atomic densities
!
      procedure, public :: get_aos_in_subset &
                        => get_aos_in_subset_ao_tool ! AO indices for a specific center subset
!
      procedure, public :: get_center &
                        => get_center_ao_tool
!
      procedure, public :: get_n_centers &
                        => get_n_centers_ao_tool
!
      procedure, public :: is_center_subset &
                        => is_center_subset_ao_tool
!
      procedure, public :: get_n_centers_in_subset &
                        => get_n_centers_in_subset_ao_tool
!
      procedure, public :: set_libint_epsilon &       ! epsilon = precision of electron repulsion
                        => set_libint_epsilon_ao_tool !           integrals
!
      procedure, public :: set_eri_cutoff &
                        => set_eri_cutoff_ao_tool
!
      procedure, public :: set_oei_cutoff &
                        => set_oei_cutoff_ao_tool
!
      procedure, public :: get_libint_epsilon &
                        => get_libint_epsilon_ao_tool
!
      procedure, public :: get_eri_cutoff &
                        => get_eri_cutoff_ao_tool
!
      procedure, public :: get_oei_cutoff &
                        => get_oei_cutoff_ao_tool
!
      procedure, public :: get_reduced_ao_metric &
                        => get_reduced_ao_metric_ao_tool
!
      procedure, public :: orthonormal_ao_pivot_basis_transformation &
                        => orthonormal_ao_pivot_basis_transformation_ao_tool ! Y = P^T X P, where P are pivots
!
      procedure, public :: orthonormal_ao_pivot_transformation &
                        => orthonormal_ao_pivot_transformation_ao_tool ! Y = P X
!
      procedure, public :: get_n_orthonormal_ao &
                        => get_n_orthonormal_ao_ao_tool
!
      procedure, public :: print_ao_vectors &
                        => print_ao_vectors_ao_tool
!
      procedure, public :: get_frozen_cores_and_n_frozen_orbitals &
                        => get_frozen_cores_and_n_frozen_orbitals_ao_tool
!
      procedure, public :: print_centers &
                        => print_centers_ao_tool
!
      procedure, public :: print_z_matrix &
                        => print_z_matrix_ao_tool
!
      procedure, public :: evaluate_aos_at_point &
                        => evaluate_aos_at_point_ao_tool
!
      procedure, public :: get_center_coordinates &
                        => get_center_coordinates_ao_tool
!
      procedure, public :: get_center_symbols &
                        => get_center_symbols_ao_tool
!
      procedure, public :: get_point_charges &
                        => get_point_charges_ao_tool
!
      procedure, public :: export_centers_to_libint &
                        => export_centers_to_libint_ao_tool
!
      procedure, public :: initialize_libint_integral_engines &
                        => initialize_libint_integral_engines_ao_tool
!
      procedure, public :: get_n_electrons &
                        => get_n_electrons_ao_tool
!
      procedure, public :: get_subset_point_charges &
                        => get_subset_point_charges_ao_tool
!
      procedure, public :: has_ghost_atoms &
                        => has_ghost_atoms_ao_tool
!
!     Routines to write a molden file
!
      procedure, public :: get_molden_ao_indices &
                        => get_molden_ao_indices_ao_tool
!
      procedure, public :: print_molden_geometry &
                        => print_molden_geometry_ao_tool
!
      procedure, public :: print_basis_set_molden &
                        => print_basis_set_molden_ao_tool
!
      procedure, public :: get_SAD_center_indices &
                        => get_SAD_center_indices_ao_tool
!
      procedure, private :: construct_cs_eri_max_screenings ! Cauchy-Schwarz (CS),
                                                            ! electron repulsion integrals (ERIs)
!
      procedure, private :: construct_oei
      procedure, private :: construct_oei_screened
      procedure, private :: construct_oei_1der
!
      procedure, private :: determine_linearly_independent_aos
!
      procedure, private :: set_up_centers
!
      procedure, private :: initialize_centers
      procedure, private :: initialize_centers_from_template
      procedure, private :: initialize_centers_from_ao_tool
      procedure, private :: initialize_centers_from_input_file
      procedure, private :: initialize_center_shells
!
      procedure, private :: set_atomic_center_positions
!
      procedure, private :: initialize_total_charge
!
      procedure, private :: calculate_n_shells
      procedure, private :: calculate_n_aos
      procedure, private :: calculate_max_shell_size
!
      procedure, private :: construct_shells
      procedure, private :: construct_shell_to_center
      procedure, private :: construct_ao_to_shell
      procedure, private :: construct_ao_to_center
!
      procedure, private :: calculate_geometry_dependent_variables ! h, s, ERI screening lists, ...
      procedure, private :: copy_geometry_dependent_variables
!
      procedure, private :: get_subset_index
      procedure, private :: shp_on_same_atom
!
      final :: destructor
!
   end type ao_tool
!
   include "../libint/atom_init_cdef.F90"
   include "../libint/oei_cdef.F90"
   include "../libint/eri_cdef.F90"
   include "../libint/libint_initialization_cdef.F90"
!
!  Interface for constructor(s)
!
   interface ao_tool
!
      procedure :: new_ao_tool
!
   end interface ao_tool
!
!
contains
!
!
   function new_ao_tool() result(ao)
!!
!!    New AO tool
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      type(ao_tool) :: ao
!
      ao%lindep_threshold = 1.0d-6
      ao%eri_cutoff       = 1.0d-12
      ao%oei_cutoff       = 1.0d-17
!
      ao%libint_epsilon = ao%eri_cutoff**2
!
      call input%get_keyword('integral precision',  &
                                        'solver scf',          &
                                        ao%libint_epsilon)
!
      call output%printf('v', 'Libint electron repulsion integral precision: (e11.4)', &
                               reals=[ao%libint_epsilon], fs='(/t3,a)')
!
!     Enforce Cartesian or spherical Gaussians, or use default for the given basis?
!
      ao%basis_type_ = 'default'
!
      if (input%is_keyword_present('cartesian gaussians', 'system')) then
!
         ao%basis_type_ = 'cartesian'
!
      elseif (input%is_keyword_present('pure gaussians', 'system')) then
!
         ao%basis_type_ = 'spherical'
!
      endif
!
   end function new_ao_tool
!
!
   subroutine initialize_ao_tool(ao, centers, charge)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Makes the AO tool ready to use.
!!    Centers are by default read from input. Uses the optional 'centers' instead if given.
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers
      integer, intent(in), optional :: charge
!
      logical :: print_z_matrix
!
      call ao%set_up_centers(centers)
      call ao%initialize_total_charge(charge)
!
      call ao%print_centers('angstrom')
      call ao%print_centers('bohr')
!
      print_z_matrix = input%is_keyword_present('z-matrix', 'print')
      if (print_z_matrix) call ao%print_z_matrix()
!
      call ao%calculate_n_aos()
      call ao%calculate_n_shells()
!
      allocate(ao%shells(ao%n_sh))
!
      call mem%alloc(ao%shell_to_center, ao%n_sh)
      call mem%alloc(ao%ao_to_center, ao%n)
      call mem%alloc(ao%ao_to_shell, ao%n)
!
      call ao%construct_shells()
      call ao%construct_shell_to_center()
      call ao%construct_ao_to_center()
      call ao%construct_ao_to_shell()
!
      call ao%calculate_max_shell_size()
!
      call mem%alloc(ao%cs_eri_max,         ao%n_sh*(ao%n_sh + 1)/2, 2)
      call mem%alloc(ao%cs_eri_max_indices, ao%n_sh*(ao%n_sh + 1)/2, 3)
!
      call ao%construct_stored_oei('overlap')
!
      call ao%determine_linearly_independent_aos()
      call ao%construct_cs_eri_max_screenings()
!
   end subroutine initialize_ao_tool
!
!
   subroutine initialize_ao_tool_from_template(ao, template)
!!
!!    Initialize ao_tool from template
!!    Written by Alexander C. Paul, Feb. 2021
!!
!!    Makes the AO tool ready to use, copying from an ao_tool template
!!
      use array_utilities, only: copy_integer
!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
      class(ao_tool), intent(in)    :: template

      call ao%set_up_centers(ao_template=template)
!
      ao%n      = template%n
      ao%n_sh   = template%n_sh
      ao%charge = template%charge
!
      ao%eri_cutoff = template%eri_cutoff
      ao%oei_cutoff = template%oei_cutoff
!
      allocate(ao%shells(ao%n_sh))
      ao%shells = template%shells
!
      call mem%alloc(ao%shell_to_center, ao%n_sh)
      call copy_integer(template%shell_to_center, ao%shell_to_center, ao%n_sh)
!
      call mem%alloc(ao%ao_to_center, ao%n)
      call copy_integer(template%ao_to_center, ao%ao_to_center, ao%n)
!
      call mem%alloc(ao%ao_to_shell, ao%n)
      call copy_integer(template%ao_to_shell, ao%ao_to_shell, ao%n)
!
      ao%max_sh_size = template%max_sh_size
!
      call ao%initialize_oei('hamiltonian')
      call ao%initialize_oei('overlap')
!
      call mem%alloc(ao%cs_eri_max,         ao%n_sh*(ao%n_sh + 1)/2, 2)
      call mem%alloc(ao%cs_eri_max_indices, ao%n_sh*(ao%n_sh + 1)/2, 3)
!
      call ao%copy_geometry_dependent_variables(template)
!
   end subroutine initialize_ao_tool_from_template
!
!
   subroutine set_up_centers(ao, centers, ao_template)
!!
!!    Set up centers
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Makes the centers fully ready for use.
!!
!!    Reads centers from input, unless optional arguments are passed:
!!
!!       - centers:     copy from an existing set of centers
!!       - ao_template: copy from centers in another ao_tool
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers
!
      type(ao_tool), optional, intent(in) :: ao_template
!
      call ao%initialize_centers(centers, ao_template)
      call ao%export_centers_to_libint()
      call ao%initialize_center_shells()
!
   end subroutine set_up_centers
!
!
   subroutine initialize_centers(ao, centers, ao_template)
!!
!!    Initialize centers
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Allocates and sets the centers as well as the center subsets.
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers
!
      type(ao_tool), optional, intent(in) :: ao_template
!
      if (present(centers)) then
!
         call ao%initialize_centers_from_template(centers)
!
      else if (present(ao_template)) then
!
         call ao%initialize_centers_from_ao_tool(ao_template)
!
      else
!
         call ao%initialize_centers_from_input_file()
!
      endif
!
   end subroutine initialize_centers
!
!
   subroutine initialize_centers_from_input_file(ao)
!!
!!    Initialize centers from input file
!!    Written by Eirik F. Kjønstad, 2021
!!
      use atomic_center_reader_class, only: atomic_center_reader
!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      class(atomic_center_reader), allocatable :: center_reader
!
      center_reader = atomic_center_reader(ao%basis_type_)
!
      call center_reader%read_centers()
!
      ao%n_centers = center_reader%get_n_centers()
      ao%n_center_subsets = center_reader%get_n_center_subsets()
!
      allocate(ao%centers(ao%n_centers))
      allocate(ao%center_subsets(ao%n_center_subsets))
!
      call center_reader%get_centers(ao%centers)
      call center_reader%get_center_subsets(ao%center_subsets)
!
      deallocate(center_reader)
!
   end subroutine initialize_centers_from_input_file
!
!
   subroutine initialize_centers_from_template(ao, centers)
!!
!!    Initialize centers from template
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      class(atomic_center), dimension(:), intent(in) :: centers
!
      integer :: I
!
      ao%n_centers = size(centers)
!
      allocate(ao%centers(ao%n_centers))
!
      do I = 1, ao%n_centers
!
         ao%centers(I) = atomic_center(I,                         &
                                       -1,                        &
                                       centers(I)%symbol,         &
                                       centers(I)%coordinates,    &
                                       centers(I)%basis,          &
                                       ao%basis_type_,            &
                                       centers(I)%is_ghost())
!
      enddo
!
      ao%n_center_subsets = 1
      allocate(ao%center_subsets(ao%n_center_subsets))
!
      ao%center_subsets(1) = named_range('unclassified',  &
                                         1,               &
                                         ao%n_centers)
!
   end subroutine initialize_centers_from_template
!
!
   subroutine initialize_centers_from_ao_tool(ao, ao_template)
!!
!!    Initialize centers from ao tool
!!    Written by Alexander C. Paul, Feb 2021
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
      class(ao_tool), intent(in)    :: ao_template
!
      integer :: I
!
      ao%n_centers = ao_template%n_centers
!
      allocate(ao%centers(ao%n_centers))
!
      do I = 1, ao%n_centers
!
         ao%centers(I) = atomic_center(ao_template%centers(i)%libint_number, &
                                       ao_template%centers(i)%input_number,  &
                                       ao_template%centers(i)%symbol,        &
                                       ao_template%centers(i)%coordinates,   &
                                       ao_template%centers(i)%basis,         &
                                       ao_template%basis_type_,              &
                                       ao_template%centers(i)%is_ghost())
!
      enddo
!
      ao%n_center_subsets = ao_template%n_center_subsets
      allocate(ao%center_subsets(ao%n_center_subsets))
!
      ao%center_subsets = ao_template%center_subsets
!
   end subroutine initialize_centers_from_ao_tool
!
!
   subroutine export_centers_to_libint_ao_tool(ao)
!!
!!    Export centers to Libint
!!    Written by Rolf H. Myhre, Mar. 2020
!!
      use libint_initialization, only: export_geometry_and_basis_to_libint
!
      implicit none
!
      class(ao_tool) :: ao
!
      call export_geometry_and_basis_to_libint(ao%centers)
!
      call ao%initialize_libint_integral_engines()
!
   end subroutine export_centers_to_libint_ao_tool
!
!
   subroutine initialize_center_shells(ao)
!!
!!    Initialize shells
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: I
!
      do I = 1, ao%n_centers
!
         call ao%centers(I)%initialize_shells()
!
      enddo
!
   end subroutine initialize_center_shells
!
!
   subroutine initialize_total_charge(ao, charge)
!!
!!    Initialize total_charge
!!    Written by Sarai D. Folkestad, 2021
!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
      integer, optional, intent(in) :: charge
!
      if (present(charge)) then
!
         ao%charge = charge
         return
!
      endif
!
      ao%charge = 0
      call input%get_keyword('charge', 'system', ao%charge)
!
   end subroutine initialize_total_charge
!
!
   subroutine calculate_n_shells(ao)
!!
!!    Calculate number of shells
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: I
!
      ao%n_sh = 0
!
      do I = 1, ao%n_centers
!
         ao%n_sh = ao%n_sh + ao%centers(I)%n_shells
!
      enddo
!
   end subroutine calculate_n_shells
!
!
   subroutine calculate_n_aos(ao)
!!
!!    Caculate number of AOs
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: I
!
      ao%n = 0
!
      do I = 1, ao%n_centers
!
         ao%n = ao%n + ao%centers(I)%n_ao
!
      enddo
!
   end subroutine calculate_n_aos
!
!
   subroutine initialize_libint_integral_engines_ao_tool(ao)
!!
!!    Initialize Libint integral engines
!!    Written by Eirik F. Kjønstad, June 2019
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      call initialize_eri_c(eri_precision = ao%libint_epsilon)
      call initialize_kinetic_c()
      call initialize_nuclear_c()
      call initialize_overlap_c()
      call initialize_dipole_c()
      call initialize_quadrupole_c()
!
   end subroutine initialize_libint_integral_engines_ao_tool
!
!
   subroutine construct_ao_to_center(ao)
!!
!!    Construct AO to center
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: I, A, B, w
!
      A = 0
!
      do I = 1, ao%n_centers
         do B = 1, ao%centers(I)%n_shells
!
            A = A + 1
!
            do w = ao%shells(A)%first, ao%shells(A)%get_last()
!
               ao%ao_to_center(w) = I
!
            enddo
         enddo
      enddo
!
   end subroutine construct_ao_to_center
!
!
   subroutine construct_ao_to_shell(ao)
!!
!!    Construct AO to shell
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: A, w
!
      A = 0
!
      do A = 1, ao%n_sh
         do w = ao%shells(A)%first, ao%shells(A)%get_last()
!
            ao%ao_to_shell(w) = A
!
         enddo
      enddo
!
   end subroutine construct_ao_to_shell
!
!
   subroutine construct_shell_to_center(ao)
!!
!!    Construct shell to center
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: I, A, B
!
      A = 0
!
      do I = 1, ao%n_centers
         do B = 1, ao%centers(I)%n_shells
!
            A = A + 1
!
            ao%shell_to_center(A) = I
!
         enddo
      enddo
!
   end subroutine construct_shell_to_center
!
!
   subroutine construct_shells(ao)
!!
!!    Construct shells
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Constructs 'shells', which is a vector containing the AO intervals for each shell.
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: I, A, B
!
      A = 0
!
      do I = 1, ao%n_centers
         do B = 1, ao%centers(I)%n_shells
!
            A = A + 1
!
            ao%shells(A) = range_(ao%centers(I)%shells(B))
!
         enddo
      enddo
!
   end subroutine construct_shells
!
!
   subroutine calculate_max_shell_size(ao)
!!
!!    Calculate maximum shell size
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: A
!
      ao%max_sh_size = 0
!
      do A = 1, ao%n_sh
!
         if (ao%shells(A)%length .gt. ao%max_sh_size) then
!
            ao%max_sh_size = ao%shells(A)%length
!
         endif
!
      enddo
!
   end subroutine calculate_max_shell_size
!
!
   function get_n_oei_components_ao_tool(oei_type) result(n_components)
!!
!!    Get n OEI components
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns the number of one-electron integral components for the
!!    integral of type 'oei_type'.
!!
      implicit none
!
      character(len=*), intent(in) :: oei_type
!
      integer :: n_components
!
      if (trim(oei_type) == 'hamiltonian') then
!
         n_components = 1
!
      elseif (trim(oei_type) == 'electrostatic potential') then
!
         n_components = 1
!
      elseif (trim(oei_type) == 'electrostatic potential unit') then
!
         n_components = 1
!
      elseif (trim(oei_type) == 'overlap') then
!
         n_components = 1
!
      elseif (trim(oei_type) == 'dipole') then
!
         n_components = 3
!
      elseif (trim(oei_type) == 'quadrupole') then
!
         n_components = 6
!
      else
!
         n_components = 0
!
         call output%error_msg('Could not recognize one-electron integral (a0)', &
                               chars=[trim(oei_type)])
!
      endif
!
   end function get_n_oei_components_ao_tool
!
!
   subroutine get_oei_ao_tool(ao, oei_type, oei, screening)
!!
!!    Get OEI (one-electron integral)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Constructs one-electron integrals in the AO basis and places the result in 'oei'
!!
!!       oei:      integral array (n_ao x n_ao x n_components)
!!       oei_type: which integral to calculate
!!
!!    'oei_type' can be:
!!
!!       - 'hamiltonian'  One-electron Hamiltonian (h)  n_components = 1
!!       - 'overlap'      AO overlap (S)                n_components = 1
!!       - 'dipole'       Dipole moment (mu)            n_components = 3 (mu_x, mu_y, mu_z)
!!       - 'quadrupole'   Quadrupole moment (q)         n_components = 6 (q_xx, q_xy, q_xz,
!!                                                                        q_yy, q_yz, q_zz)
      use timings_class, only: timings
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      character(len=*), intent(in) :: oei_type
!
      real(dp), dimension(*), intent(out) :: oei
!
      logical, intent(in), optional :: screening
!
      integer :: n_components
!
      type(timings), allocatable :: timer
!
      logical :: screening_local
!
      timer = timings('One-electron integrals (' // trim(oei_type) // ')', 'm')
      call timer%turn_on()
!
      screening_local = .false.
      if (present(screening)) screening_local = screening
!
      n_components = ao%get_n_oei_components(oei_type)
!
      if (screening_local) then
!
         call ao%construct_oei_screened(trim(oei_type), oei, n_components)
!
      else
!
         call ao%construct_oei(trim(oei_type), oei, n_components)
!
      endif
!
      call timer%turn_off()
!
   end subroutine get_oei_ao_tool
!
!
   subroutine construct_oei(ao, oei_type, x, n_components)
!!
!!    Construct OEI (one-electron integral)
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Computes the one-electron AO integral matrix x using
!!    the C++ construction routine get_oei_c. This C++ routine constructs
!!    the x_ABk contributions to x for the shells A and B, where k
!!    denotes the components (e.g., x, y, z, for dipole integrals)
!!
!!    Modified by EFK and SDF, 2020: generalized to any number of components.
!
      use iso_c_binding
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(ao%n, ao%n, n_components), intent(out) :: x
!
      character(len=*) :: oei_type
!
      integer :: w, y, k, w_f, y_f
!
      integer(c_int) :: A, B
!
      real(dp), dimension(:,:,:), pointer                           :: x_ABk_p
      real(dp), dimension((ao%max_sh_size)**2*n_components), target :: x_ABk
!
      character(len=100, kind=c_char) :: oei_type_c
!
      oei_type_c = trim(oei_type) // c_null_char
!
!$omp parallel do private(A, B, x_ABk, x_ABk_p, w, y, w_f, y_f, k) schedule(dynamic)
      do A = 1, int(ao%n_sh, kind=c_int)
         do B = 1, A
!
            call get_oei_c(oei_type_c, x_ABk, A, B)
!
            x_ABk_p(1 : ao%shells(A)%length, &
                    1 : ao%shells(B)%length, &
                    1 : n_components)        &
            => x_ABk(1 : ao%shells(A)%length*&
                         ao%shells(B)%length*&
                         n_components)
!
            do k = 1, n_components
               do w = 1, ao%shells(A)%length
                  do y = 1, ao%shells(B)%length
!
                     w_f = ao%shells(A)%first - 1 + w
                     y_f = ao%shells(B)%first - 1 + y
!
                     x(w_f, y_f, k) = x_ABk_p(w, y, k)
                     x(y_f, w_f, k) = x_ABk_p(w, y, k)
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_oei
!
!
   subroutine get_oei_1der_ao_tool(ao,        &
                                   oei_type,  &
                                   oei)
!!
!!    Get OEI 1der (one-electron integral)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Constructs first derivative of one-electron integrals in the AO basis and places
!!    the rersult in 'oei'.
!!
!!       oei:      integral array (n_ao x n_ao x 3 * n_centers)
!!       oei_type: integral type
!!
!!    'oei_type' can be:
!!
!!       - 'kinetic'      Kinetic contribution to one-electron Hamiltonian (h)
!!       - 'nuclear'      Nuclear attraction contribution to one-electron Hamiltonian (h)
!!       - 'overlap'      AO overlap (S)
!!
      use array_utilities, only: zero_array
      use timings_class, only: timings
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      character(len=*), intent(in) :: oei_type
!
      real(dp), dimension(*), intent(out) :: oei
!
      integer(c_int) :: A, B, n_ao
!
      type(timings), allocatable :: timer
!
      timer = timings('First-derivative, one-electron integrals (' // trim(oei_type) // ')', 'v')
      call timer%turn_on()
!
      call output%printf('v',                                                             &
                         'Constructing first-derivative one-electron integrals ((a0))',   &
                         chars=[oei_type],                                                &
                         fs='(/t3,a)')
!
!     Call appropriate construction routine
!
      if (trim(oei_type) == 'kinetic' .or. &
          trim(oei_type) == 'overlap') then
!
!        Operator is independent of nuclear centers => call general routine
!
         call ao%construct_oei_1der(trim(oei_type), oei)
!
      elseif (trim(oei_type) == 'nuclear') then
!
!        Operator depends on nuclear centers
!
         call zero_array(oei, 3 * ao%n_centers * ao%n**2)
!
         n_ao = int(ao%n, c_int)
!
!$omp parallel do private(A, B)
         do A = 1, int(ao%n_sh, c_int)
            do B = 1, A
!
               call add_nuclear_h_1der_c(oei, A, B, n_ao)
!
            enddo
         enddo
!$omp end parallel do
!
      else
!
!        Could not recognize integral type, or not yet implemented; give error
!
         call output%error_msg('Could not recognize integral type (a0) in get_oei_1der!', &
                               chars=[trim(oei_type)])
!
      endif
!
      call timer%turn_off()
!
   end subroutine get_oei_1der_ao_tool
!
!
   subroutine construct_oei_1der(ao, oei_type, x)
!!
!!    Construct OEI 1 der (one-electron integral)
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Computes the first derivative of one-electron AO integral matrix x.
!!
!!    Uses get_oei_c routine, which gives x_ABqk contributions to x for the shells A and B,
!!    where k refers to the shell centers (k = 1, 2 refers to A, B) and q to their xyz-coord.
!!
      use array_utilities, only: zero_array
      use iso_c_binding
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp), dimension(ao%n, ao%n, 3, ao%n_centers), intent(out) :: x
!
      character(len=*) :: oei_type
!
      integer :: w, y, q, w_f, y_f
!
      integer(c_int) :: A, B
!
      real(dp), dimension(:,:,:,:), pointer              :: x_ABqk_p
      real(dp), dimension(6*(ao%max_sh_size)**2), target :: x_ABqk
!
      character(len=100, kind=c_char) :: oei_type_c
!
      integer :: A_atom, B_atom
!
      oei_type_c = trim(oei_type) // c_null_char
!
      call zero_array(x, 3*(ao%n_centers)*(ao%n)**2)
!
!$omp parallel do private(A, B, A_atom, B_atom, x_ABqk, x_ABqk_p, w, y, w_f, y_f, q) schedule(dynamic)
      do A = 1, int(ao%n_sh, kind=c_int)
         do B = 1, A
!
            call get_oei_1der_c(oei_type_c, x_ABqk, A, B)
!
            x_ABqk_p(1 : ao%shells(A)%length, 1 : ao%shells(B)%length, 1 : 3, 1 : 2) &
                                       => x_ABqk(1 : 6 * ao%shells(A)%length * ao%shells(B)%length)
!
            A_atom = ao%shell_to_center(A)
            B_atom = ao%shell_to_center(B)
!
            do q = 1, 3
               do w = 1, ao%shells(A)%length
                  do y = 1, ao%shells(B)%length
!
                     w_f = ao%shells(A)%first - 1 + w
                     y_f = ao%shells(B)%first - 1 + y
!
                     x(w_f, y_f, q, A_atom) = x(w_f, y_f, q, A_atom) + x_ABqk_p(w, y, q, 1)
                     x(y_f, w_f, q, A_atom) = x(y_f, w_f, q, A_atom) + x_ABqk_p(w, y, q, 1)
!
                     x(w_f, y_f, q, B_atom) = x(w_f, y_f, q, B_atom) + x_ABqk_p(w, y, q, 2)
                     x(y_f, w_f, q, B_atom) = x(y_f, w_f, q, B_atom) + x_ABqk_p(w, y, q, 2)
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_oei_1der
!
!
   subroutine get_eri_ao_tool(ao, g, A, B, C, D, precision_, skip)
!!
!!    Get ERI (electron repulsion integral)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Wrapper for constructing electron repulsion integrals g for the
!!    shell quartet (A, B, C, D).
!!
!!    The integrals are placed in g.
!!
!!    Two optional arguments:
!!
!!       precision_:  (intent in) Double precision real corresponding to the Libint precision
!!                   'epsilon' to use when calculating the integral. Does not guarantee a precision
!!                   to the given value and should therefore be selected conservatively.
!!
!!       skip:       (intent out) If present, this integer will be 1 if Libint decided not to
!!                   calculate the integral; it will be zero otherwise. If it is present, g will
!!                   not be zeroed out if Libint decides not to calculate g. Thus, only pass 'skip'
!!                   to the routine if you wish to avoid zeroing out elements that are negligible.
!!
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer :: A, B, C, D
!
      real(dp), dimension(ao%shells(A)%length * ao%shells(B)%length * &
                          ao%shells(C)%length * ao%shells(D)%length), intent(out) :: g
!
      real(dp), optional, intent(in) :: precision_
      integer, optional, intent(out) :: skip
!
      integer(c_int) :: A_, B_, C_, D_
      integer(c_int) :: n_A, n_B, n_C, n_D
!
      integer(c_int) :: skip_local
!
      real(dp) :: precision_local
!
!     Convert to C integers
!
      A_ = int(A, c_int)
      B_ = int(B, c_int)
      C_ = int(C, c_int)
      D_ = int(D, c_int)
!
      n_A = int(ao%shells(A)%length, c_int)
      n_B = int(ao%shells(B)%length, c_int)
      n_C = int(ao%shells(C)%length, c_int)
      n_D = int(ao%shells(D)%length, c_int)
!
!     Set precision to non-default value if requested
!
      if (present(precision_)) then
!
         precision_local = precision_
!
      else
!
         precision_local = ao%libint_epsilon
!
      endif
!
!     Get the integrals g from Libint
!
      call get_eri_c(g,                    &
                     A_, B_, C_, D_,       &
                     precision_local,      &
                     skip_local,           &
                     n_A, n_B, n_C, n_D)
!
!     Return skip if requested; otherwise, zero g if skip_local = 1
!
      if (present(skip)) then
!
         skip = skip_local
!
      else
!
         if (skip_local .eq. 1) call zero_array(g, int(n_A * n_B * n_C * n_D))
!
      endif
!
   end subroutine get_eri_ao_tool
!
!
   subroutine get_eri_1der_ao_tool(ao, g, A, B, C, D)
!!
!!    Get ERI 1der (electron repulsion integral)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Wrapper for constructing first derivative electron repulsion integrals g for the
!!    shell quartet (A, B, C, D). The result is stored in g, in the order of ABCDqk, where
!!    q = 1,2,3 are Cartesian coordinates (x,y,z) and k = 1,2,3,4 are the atomic centers
!!    corresponding to the atoms where A,B,C,D are located.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer :: A, B, C, D
!
      real(dp), dimension(ao%shells(A)%length * ao%shells(B)%length * &
                          ao%shells(C)%length * ao%shells(D)%length * 12), intent(out) :: g
!
      integer(c_int) :: A_, B_, C_, D_
!
!     Convert to C integers
!
      A_ = int(A, c_int)
      B_ = int(B, c_int)
      C_ = int(C, c_int)
      D_ = int(D, c_int)
!
!     Get the integrals g from Libint
!
      call get_eri_1der_c(g, A_, B_, C_, D_)
!
   end subroutine get_eri_1der_ao_tool
!
!
   subroutine initialize_external_charges_from_point_charges_ao_tool(pc)
!!
!!    Initialize external charges from point_charges
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Initializes the external charge engines that are needed for
!!    embedding (QM/MM)
!!
      use libint_initialization, only: initialize_coulomb_external_charges_c
!
      use iso_c_binding
!
      implicit none
!
      type(point_charges), intent(in) :: pc
!
      call initialize_coulomb_external_charges_c(pc%q,      &
                                  pc%r*angstrom_to_bohr,    &
                                  int(pc%n_charges, kind=c_int))
!
   end subroutine initialize_external_charges_from_point_charges_ao_tool
!
!
   subroutine initialize_external_charges_ao_tool(n_charges, positions, charges)
!!
!!    Initialize external charges
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Initializes the external charge engines that are needed for
!!    embedding (QM/MM)
!!
      use libint_initialization, only: initialize_coulomb_external_charges_c

!
      use iso_c_binding
!
      implicit none
!
      integer,                            intent(in) :: n_charges
      real(dp), dimension(3, n_charges),  intent(in) :: positions
      real(dp), dimension(n_charges),     intent(in) :: charges
!
      call initialize_coulomb_external_charges_c(charges,      &
                                  positions*angstrom_to_bohr,  &
                                  int(n_charges, kind=c_int))
!
   end subroutine initialize_external_charges_ao_tool
!
!
   subroutine initialize_external_unit_charges_ao_tool(n_charges, positions)
!!
!!    Initialize external unit charges
!!    Written by Sarai D. Folkestad, 2020
!!
      use libint_initialization, only: initialize_coulomb_external_unit_charges_c
!
      use iso_c_binding
!
      implicit none
!
      integer,                            intent(in) :: n_charges
      real(dp), dimension(3, n_charges),  intent(in) :: positions
!
      call initialize_coulomb_external_unit_charges_c(positions*angstrom_to_bohr,  &
                                                      int(n_charges, kind=c_int))
!
   end subroutine initialize_external_unit_charges_ao_tool
!
!
   subroutine initialize_oei_ao_tool(ao, oei_type)
!!
!!    Initialize OEI
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Routine to allocate one-electron integrals that are
!!    to be kept in memory in the AO tool.
!!
!!    Valid oei_type values:
!!
!!       - 'hamiltonian'  One-electron Hamiltonian (h)
!!       - 'overlap'      AO overlap (S)
!!       - 'electrostatic potential'      Coulombic interaction of electrons
!!                                 with the external charges (v)
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      character(len=*), intent(in) :: oei_type
!
      if (trim(oei_type) == 'hamiltonian') then
!
         if (.not. allocated(ao%h)) call mem%alloc(ao%h, ao%n, ao%n)
!
      elseif (trim(oei_type) == 'overlap') then
!
         if (.not. allocated(ao%s)) call mem%alloc(ao%s, ao%n, ao%n)
!
      elseif (trim(oei_type) == 'electrostatic potential') then
!
         if (.not. allocated(ao%v)) call mem%alloc(ao%v, ao%n, ao%n)
!
      else
!
         call output%error_msg('Did not recognize one-electron integral (a0) &
                                 &to initialize. Some integrals (e.g. dipole) can &
                                 &only be constructed on-the-fly - see the get_oei routine.')
!
      endif
!
   end subroutine initialize_oei_ao_tool
!
!
   subroutine construct_stored_oei_ao_tool(ao, oei_type, screening)
!!
!!    Construct stored OEI
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    This routine only constructs integrals that are internally stored.
!!    For computing integrals on-the-fly, use the get_oei routine.
!!
!!    'oei_type' can be:
!!
!!       - 'hamiltonian'           One-electron Hamiltonian (h)
!!       - 'overlap'               AO overlap (S)
!!       - 'electrostatic potential'      Coulombic interaction of electrons
!!                                 with the external charges (v)
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      character(len=*), intent(in) :: oei_type
!
      logical, intent(in), optional :: screening
!
      if (trim(oei_type) == 'hamiltonian') then
!
         call ao%initialize_oei('hamiltonian')
         call ao%get_oei('hamiltonian', ao%h, screening)
!
      elseif (trim(oei_type) == 'overlap') then
!
         call ao%initialize_oei('overlap')
         call ao%get_oei('overlap', ao%s)
!
      elseif (trim(oei_type) == 'electrostatic potential') then
!
         call ao%initialize_oei('electrostatic potential')
         call ao%get_oei('electrostatic potential', ao%v, screening)
!
      else
!
         call output%error_msg('Did not recognize one-electron integral (a0) &
                                 &to construct. Some integrals (e.g. dipole) can &
                                 &only be constructed on-the-fly - see the get_oei routine.')
!
      endif
!
   end subroutine construct_stored_oei_ao_tool
!
!
   subroutine print_centers_ao_tool(ao, units)
!!
!!    Print centers
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints xyz of the molecular geometry in the unit specified ('angstrom' or 'bohr').
!!
!!    Modified by Tor S. Haugland, Oct 2019
!!    Printf and print_separator. Modified to take unit as input
!!
!!    Modifified by Alexander C. Paul, Dec 2019
!!    Print of atom numbers and separate lines for basis sets
!!
!!    Modified by Eirik F. Kjønstad, Sep 2020
!!    Renamed some variables, moved and adapted routine to AO tool
!!
!!    Modified by Tor S. Haugland, May 2021
!!    Print ghost atoms
!!
      implicit none
!
      class(ao_tool), intent(in)   :: ao
      character(len=*), intent(in) :: units
!
      real(dp) :: conversion_factor
!
      integer :: I, line_length
!
      logical :: print_basis, print_ghost, print_separator
!
      real(dp), dimension(3) :: position_
!
!     Print header for coordinate print
!
      line_length = 78
!
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(/t6,a)')
!
      call output%printf('m', 'Geometry ((a0))', chars=[get_units_label(units)], fs='(t38,a)')
!
      call output%print_separator(pl='minimal', symbol='=', n=line_length, fs='(t6,a)')
!
      call output%printf('m', 'Atom           X                  Y                  &
                         &Z         # in input', ll=100, fs='(t9,a)')
!
      call output%print_separator(pl='minimal', symbol='=', n=line_length, fs='(t6,a)')
!
!     Print actual coordinates
!
      conversion_factor = get_conversion_factor(from='angstrom', to=units)
!
      do I = 1, ao%n_centers
!
         if (I == 1) then
            print_basis = .true.
            print_ghost = ao%centers(I)%is_ghost()
            print_separator = .false.
         else
            print_basis = ao%centers(I-1)%basis .ne. ao%centers(I)%basis
            print_ghost = (.not. ao%centers(I-1)%is_ghost()) .and. ao%centers(I)%is_ghost()
            print_separator = ao%centers(I-1)%is_ghost() .and. (.not. ao%centers(I)%is_ghost())
         endif
!
         if (print_ghost) then
            call output%printf('m', &
               '=============================== Ghost atoms ==================================',&
               fs='(t6,a)')
         endif
!
         if (print_separator) then
            call output%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
         endif
!
         if (print_basis) then
            call output%printf('m', &
                               'Basis: ' // trim(ao%centers(I)%get_basis_set_name()), &
                               fs='(t9,a)')
         end if
!
         position_ = ao%centers(I)%coordinates * conversion_factor
!
         call output%printf('m', '(i4) (a2) (f18.12) (f18.12) (f18.12)  (i7)',   &
                            chars=[ao%centers(I)%symbol],                        &
                            ints=[I, ao%centers(I)%input_number],                &
                            reals=position_, ll=100, fs='(t6,a)')
!
      enddo
!
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
!
   end subroutine print_centers_ao_tool
!
!
   subroutine set_atomic_centers_ao_tool(ao, R, units)
!!
!!    Set atomic centers
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Updates the position of the atomic centers and updates internal variables
!!    that depend on these positions.
!!
!!       R:     New position for the centers (3 x n_atoms)
!!       units: 'angstrom' or 'bohr', unit of R
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      real(dp), dimension(3, ao%n_centers), intent(in) :: R
!
      character(len=*), intent(in) :: units
!
      call ao%set_atomic_center_positions(R, units)
!
      call ao%export_centers_to_libint()
!
      call ao%calculate_geometry_dependent_variables()
!
   end subroutine set_atomic_centers_ao_tool
!
!
   subroutine calculate_geometry_dependent_variables(ao)
!!
!!    Calculates geometry-dependent variables
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Calculates any variable in the object that will depends on the position of
!!    the centers. Should be called on initialization of the object and whenever the
!!    center positions are altered.
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      call ao%construct_stored_oei('overlap')
      call ao%construct_stored_oei('hamiltonian', screening = .true.)
!
      call ao%determine_linearly_independent_aos()
!
      call ao%construct_cs_eri_max_screenings()
!
   end subroutine calculate_geometry_dependent_variables
!
!
   subroutine set_atomic_center_positions(ao, R, units)
!!
!!    Set atomic center positions
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      real(dp), dimension(3, ao%n_centers), intent(in) :: R
!
      character(len=*), intent(in) :: units
!
      integer :: I
!
      real(dp) :: conversion
!
      conversion = get_conversion_factor(from=units, to='angstrom')
!
      do I = 1, ao%n_centers
!
         ao%centers(I)%coordinates(:) = R(:,I) * conversion
!
      enddo
!
   end subroutine set_atomic_center_positions
!
!
   subroutine determine_linearly_independent_aos(ao, s)
!!
!!    Determine linearly independent orbitals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    This routine determines P and L, and the number of linearly
!!    independent/orthonormal AOs (n_orthonormal_ao), resulting from a Cholesky decomposition
!!    of the AO overlap matrix S to within the linear dependency threshold:
!!
!!       P^T S P = L L^T.
!!
!!    If the optional argument 's' is not present, the routine assumes that the internally
!!    stored overlap matrix (ao%s) has been allocated and constructed.
!!
!!    Modified by Rolf H. Myhre and Alexander C. Paul, Oct 2019. Included sanity check.
!!    Modified by Eirik F. Kjønstad, Sep 2020. Introduced optional S.
!!
      use timings_class, only: timings
!
      use array_utilities, only: zero_array,       &
                                 zero_array_int,   &
                                 full_cholesky_decomposition
!
      implicit none
!
      class(ao_tool), intent(inout), target :: ao
!
      real(dp), dimension(ao%n, ao%n), optional, intent(in), target :: s
!
      real(dp), dimension(:,:), pointer :: s_p        ! Points to s or ao%s
!
      integer, dimension(:), allocatable    :: pivots ! Vector containing selected pivots
      real(dp), dimension(:,:), allocatable :: L      ! Cholesky factor, full dimensionality
!
      integer :: i, j
!
      type(timings), allocatable :: timer
!
      timer = timings('Cholesky decomposition of AO overlap', 'normal')
      call timer%turn_on()
!
      call output%printf('n', '- Cholesky decomposition of AO overlap to get &
                                &linearly independent AOs:', fs='(/t3,a)')
!
      if (present(s)) then
!
         s_p => s
!
      else
!
         s_p => ao%s
!
      endif
!
      if (allocated(ao%L)) call mem%dealloc(ao%L, ao%n_orthonormal_ao, ao%n_orthonormal_ao)
      if (allocated(ao%P)) call mem%dealloc(ao%P, ao%n, ao%n_orthonormal_ao)
!
!     Decompose AO overlap
!
      call mem%alloc(L, ao%n, ao%n)
      call mem%alloc(pivots, ao%n)
!
      call zero_array(L, ao%n**2)
      call zero_array_int(pivots, ao%n)
!
      call full_cholesky_decomposition(s_p,                   &
                                       L,                     &
                                       ao%n,                  &
                                       ao%n_orthonormal_ao,   &
                                       ao%lindep_threshold,   &
                                       pivots)
!
      if (    ao%n_orthonormal_ao   .gt. ao%n     &
         .or. ao%n_orthonormal_ao   .le. 0        &
         .or. any(pivots .gt. ao%n)    &
         .or. any(pivots .le. 0)     ) then
!
         call output%printf('m', 'Something went wrong when decomposing the AO overlap.', &
                            fs='(/t3,a)')
!
         call output%printf('m', 'Did you compile with wrong type of integers in setup? &
                            &For example system native BLAS with default 64-bit integers.', &
                            ffs='(/t3,a)')
!
         call output%printf('m', 'If that is the case, use setup with --int32 or install MKL.')
!
         call output%error_msg('Failed to decompose AO overlap.')
!
      end if
!
!     Set Cholesky factor L and pivot matrix P
!
      call mem%alloc(ao%L, ao%n_orthonormal_ao, ao%n_orthonormal_ao)
!
!$omp parallel do private(j, i)
      do j = 1, ao%n_orthonormal_ao
         do i = 1, ao%n_orthonormal_ao
!
            ao%L(i, j) = L(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(L, ao%n, ao%n)
!
      call mem%alloc(ao%P, ao%n, ao%n_orthonormal_ao)
!
      call zero_array(ao%P, ao%n * ao%n_orthonormal_ao)
!
!$omp parallel do private(j)
      do j = 1, ao%n_orthonormal_ao
!
         ao%P(pivots(j), j) = one
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(pivots, ao%n)
!
      call output%printf('m', 'Linear dependence threshold:             (e8.2)', &
                         reals=[ao%lindep_threshold], fs='(/t6,a)')
!
      call output%printf('m', 'Number of atomic orbitals:               (i0)', &
                         ints=[ao%n], fs='(t6,a)')
!
      call output%printf('m', 'Number of orthonormal atomic orbitals:   (i0)', &
                         ints=[ao%n_orthonormal_ao], fs='(t6,a)')
!
      if (ao%n_orthonormal_ao .lt. ao%n) &
         call output%printf('m', 'Removed (i0) AOs due to linear dep.', &
                            ints=[ao%n - ao%n_orthonormal_ao], fs='(/t6,a)')
!
      call timer%turn_off()
!
   end subroutine determine_linearly_independent_aos
!
!
   subroutine get_sad_guess_ao_tool(ao, D)
!!
!!    Get SAD guess
!!    Written by Eirik F. Kjønstad, Aug-Sep 2018 and Sep 2020
!!
!!    Constructs the block-diagonal SAD density D from atomic UHF densities stored on file.
!!
!!    The UHF calculations have been performed with valence electrons evenly spread out in
!!    the degenerate HOMO orbitals to ensure a spherically symmetric AO density.
!!
!!    Moved and adapted to AO tool, removed loop allocations, Eirik F. Kjønstad, 2020
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ao_tool) :: ao
!
      real(dp), dimension(ao%n, ao%n), intent(out) :: D
!
      real(dp), dimension(:), allocatable, target   :: D_I     ! Holds Ith atomic density
      real(dp), dimension(:,:), pointer             :: D_I_p   ! Pointer to Ith atomic density
!
      integer :: I, w, x, w_full, x_full
!
      type(range_), allocatable :: aos_I
!
      call zero_array(D, ao%n**2)
!
      allocate(aos_I)
!
      call mem%alloc(D_I, ao%n**2)
!
      do I = 1, ao%n_centers
!
         if (ao%centers(I)%is_ghost()) cycle
!
!        Read the atomic density at the Ith center and set pointer
!        to the relevant portion of the array
!
         call ao%centers(I)%read_atomic_uhf_density(D_I)
!
         aos_I = ao%centers(I)%get_ao_range()
!
         D_I_p(1 : aos_I%length, 1 : aos_I%length) => D_I(1 : aos_I%length**2)
!
!        Copy the atomic density into the full space matrix D
!
         do w = 1, aos_I%length
            do x = 1, aos_I%length
!
               w_full = w + aos_I%first - 1
               x_full = x + aos_I%first - 1
!
               D(w_full, x_full) = D_I_p(w, x)
!
            enddo
         enddo
!
      enddo
!
      call mem%dealloc(D_I, ao%n**2)
!
      deallocate(aos_I)
!
   end subroutine get_sad_guess_ao_tool
!
!
   subroutine get_aos_in_subset_ao_tool(ao, subset, first, last)
!!
!!    Get AOs in subset
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Gets AOs in center subset specified by 'subset' and optionally returns the first
!!    and last AO indices.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      character(len=*), intent(in) :: subset
!
      integer, intent(out), optional :: first
      integer, intent(out), optional :: last
!
      integer :: set, I
!
      class(range_), allocatable :: aos
!
!     Test if there's nothing to do; if so, tell the developer not to be silly
!
      if (.not. present(first) .and. .not. present(last)) &
         call output%error_msg('You must specify either first or last in get_aos_in_subset.')
!
!     Get the center subset associated with 'subset'
!
      do set = 1, ao%n_center_subsets
!
         if (trim(subset) .eq. trim(ao%center_subsets(set)%get_name())) then
!
!           Return the requested AO indices
!
            if (present(first)) then
!
               I   = ao%center_subsets(set)%first
               aos = ao%centers(I)%get_ao_range()
!
               first = aos%first
!
            endif
!
            if (present(last)) then
!
               I   = ao%center_subsets(set)%get_last()
               aos = ao%centers(I)%get_ao_range()
!
               last = aos%get_last()
!
            endif
!
            return
!
         endif
!
      enddo
!
!     Haven't returned? Subset not found!
!
      call output%error_msg('Could not recognize subset (a0)!', chars=[subset])
!
   end subroutine get_aos_in_subset_ao_tool
!
!
   function is_center_subset_ao_tool(ao, subset) result(is_center_subset)
!!
!!    Is center subset?
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Checks whether 'subset' corresponds to one of center subsets.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      character(len=*), intent(in) :: subset
!
      logical :: is_center_subset
!
      integer :: set
!
      is_center_subset = .false.
!
      do set = 1, ao%n_center_subsets
!
         if (trim(subset) .eq. trim(ao%center_subsets(set)%get_name())) then
!
            is_center_subset = .true.
            return
!
         endif
!
      enddo
!
   end function is_center_subset_ao_tool
!
!
   function get_n_centers_in_subset_ao_tool(ao, subset) result(n_centers)
!!
!!    Get number of centers in subset
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns the number of centers in subset with name 'subset'.
!!    Returns 0 if the set does not exist.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      character(len=*), intent(in) :: subset
!
      integer :: n_centers, I
!
      n_centers = 0
!
      do I = 1, ao%n_center_subsets
!
         if (trim(ao%center_subsets(I)%get_name()) == trim(subset)) then
!
            n_centers = ao%center_subsets(I)%length
            return
!
         endif
!
      enddo
!
   end function get_n_centers_in_subset_ao_tool
!
!
   subroutine set_libint_epsilon_ao_tool(ao, epsilon_)
!!
!!    Set Libint epsilon
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Sets the default Libint epsilon value and updates Libint. Note: it is not a precise
!!    measure of ERI accuracy and should always be chosen conservatively.
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      real(dp), intent(in) :: epsilon_
!
      ao%libint_epsilon = epsilon_
      call set_eri_precision_c(ao%libint_epsilon)
!
   end subroutine set_libint_epsilon_ao_tool
!
!
   subroutine set_eri_cutoff_ao_tool(ao, eri_cutoff)
!!
!!    Set ERI cutoff
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      real(dp), intent(in) :: eri_cutoff
!
      ao%eri_cutoff = eri_cutoff
!
   end subroutine set_eri_cutoff_ao_tool
!
!
   subroutine set_oei_cutoff_ao_tool(ao, oei_cutoff)
!!
!!    Set one-electron integral (oei) cutoff
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      real(dp), intent(in) :: oei_cutoff
!
      ao%oei_cutoff = oei_cutoff
!
   end subroutine set_oei_cutoff_ao_tool
!
!
   pure function get_libint_epsilon_ao_tool(ao) result(epsilon_)
!!
!!    Get Libint epsilon
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp) :: epsilon_
!
      epsilon_ = ao%libint_epsilon
!
   end function get_libint_epsilon_ao_tool
!
!
   pure function get_eri_cutoff_ao_tool(ao) result(eri_cutoff)
!!
!!    Get ERI cutoff
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp) :: eri_cutoff
!
      eri_cutoff = ao%eri_cutoff
!
   end function get_eri_cutoff_ao_tool
!
!
   pure function get_oei_cutoff_ao_tool(ao) result(oei_cutoff)
!!
!!    Get one-electron integral (oei) cutoff
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp) :: oei_cutoff
!
      oei_cutoff = ao%oei_cutoff
!
   end function get_oei_cutoff_ao_tool
!
!
   subroutine get_center_ao_tool(ao, I, center)
!!
!!    Get center
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Sets 'center' equal to the Ith atomic center according
!!    to the AO tool.
!!
!!    Note: the returned 'center' object does not have initialized shells.
!!          The initialization of shells is performed during the initialization
!!          of the AO tool based on center objects, e.g. 'center'.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: I
!
      type(atomic_center), intent(out) :: center
!
      if (I .gt. ao%n_centers) &
         call output%error_msg('Center index (i0) greater than total number of centers (i0)', &
                               ints=[I, ao%n_centers])
!
      center = ao%centers(I)
!
   end subroutine get_center_ao_tool
!
!
   pure function get_n_centers_ao_tool(ao) result(n_centers)
!!
!!    Get number of centers
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns ao%n_centers.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer :: n_centers
!
      n_centers = ao%n_centers
!
   end function get_n_centers_ao_tool
!
!
   subroutine get_reduced_ao_metric_ao_tool(ao, S)
!!
!!    Get reduced AO metric
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2020
!!
!!    Constructs the AO overlap matrix in the
!!    linearly independent ao basis
!!
!!       S = L L^T
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
      real(dp), dimension(ao%n_orthonormal_ao, ao%n_orthonormal_ao), intent(out)  :: S
!
      call dgemm('N','T',              &
                  ao%n_orthonormal_ao, &
                  ao%n_orthonormal_ao, &
                  ao%n_orthonormal_ao, &
                  one,                 &
                  ao%L,                &
                  ao%n_orthonormal_ao, &
                  ao%L,                &
                  ao%n_orthonormal_ao, &
                  zero,                &
                  S,                   &
                  ao%n_orthonormal_ao)
!
   end subroutine get_reduced_ao_metric_ao_tool
!
!
   subroutine orthonormal_ao_pivot_basis_transformation_ao_tool(ao, X, Y)
!!
!!    OAO pivot basis transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Constructs
!!
!!       Y = P^T X P,
!!
!!    where P is the pivot matrix from Cholesky decomposition of the AO overlap matrix.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp), dimension(ao%n, ao%n), intent(in) :: X
!
      real(dp), dimension(ao%n_orthonormal_ao, ao%n_orthonormal_ao), intent(out) :: Y
!
      real(dp), dimension(:,:), allocatable :: XP
!
      call mem%alloc(XP, ao%n, ao%n_orthonormal_ao)
!
      call dgemm('N', 'N',  &
                  ao%n,     &
                  ao%n_orthonormal_ao, &
                  ao%n,     &
                  one,      &
                  X,        &
                  ao%n,     &
                  ao%P,     &
                  ao%n,     &
                  zero,     &
                  XP,       &
                  ao%n)
!
      call dgemm('T', 'N',  &
                  ao%n_orthonormal_ao, &
                  ao%n_orthonormal_ao, &
                  ao%n,     &
                  one,      &
                  ao%P,     &
                  ao%n,     &
                  XP,       &
                  ao%n,     &
                  zero,     &
                  Y,        &
                  ao%n_orthonormal_ao)
!
      call mem%dealloc(XP, ao%n, ao%n_orthonormal_ao)
!
   end subroutine orthonormal_ao_pivot_basis_transformation_ao_tool
!
!
   subroutine orthonormal_ao_pivot_transformation_ao_tool(ao, X, Y)
!!
!!    OAO pivot transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Constructs
!!
!!       Y = P X,
!!
!!    where P is the pivot matrix from Cholesky decomposition of the AO overlap matrix.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp), dimension(ao%n_orthonormal_ao, ao%n_orthonormal_ao), intent(in) :: X
!
      real(dp), dimension(ao%n, ao%n_orthonormal_ao), intent(out) :: Y
!
      call dgemm('N', 'N',  &
                  ao%n,     &
                  ao%n_orthonormal_ao, &
                  ao%n_orthonormal_ao, &
                  one,      &
                  ao%P,     &
                  ao%n,     &
                  X,        &
                  ao%n_orthonormal_ao, &
                  zero,     &
                  Y,        &
                  ao%n)
!
   end subroutine orthonormal_ao_pivot_transformation_ao_tool
!
!
   function get_n_orthonormal_ao_ao_tool(ao) result(n_orthonormal_ao)
!!
!!    Get number of OAOs
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer :: n_orthonormal_ao
!
      n_orthonormal_ao = ao%n_orthonormal_ao
!
   end function get_n_orthonormal_ao_ao_tool
!
!
   subroutine print_ao_vectors_ao_tool(ao, C, out_file, m, offset)
!!
!!    Print AO vector
!!    Written by Eirik F. Kjønstad, Tor S. Haugland, and Alexander C. Paul, 2019-2020
!!
!!    This routine contains parts of a routine, written & modified by Eirik F. Kjønstad,
!!    Tor S. Haugland, and Alexander C. Paul, that was originally part of the HF wavefunction.
!!    Moved to AO tool by Eirik F. Kjønstad, 2020.
!!
!!    Prints the array 'C' to the output file 'out'
!!
!!       C:        (n_ao x m) real double precision vector to print
!!       out_file: output file where C is to be printed
!!       m:        Number of vectors; i.e., second dimension of C
!!       offset:   Number/label for the vectors in C. In particular, k + offset
!!                 is the label given to kth column of C.
!!
!!    The C-amplitudes along each AO is printed with detailed information
!!    of the AO in question (center, angular momentum and projection, l and m_l,
!!    and so on).
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: m
!
      integer, intent(in) :: offset
!
      real(dp), dimension(ao%n, m), intent(in) :: C
!
      type(output_file), intent(in) :: out_file
!
      integer :: w, A, w_red, p, I
!
      logical :: adv
!
      character(len=2)  :: symbol
      character(len=8)  :: ang_mom
      character(len=50) :: n_format_string
      character(len=7)  :: format_string = '(f10.6)'
!
      call out_file%printf('n', '  AO    Center  l m_l', fs='(/t3,a)', adv=.false.)
!
      do p = 1, m
!
         adv = (p == m)
!
         call out_file%printf('n', '(i4)', adv=adv, ints=[p + offset], fs='(8x,a)')
!
      enddo
!
      call out_file%print_separator(pl='normal', n=87, symbol='-')
!
      do I = 1, ao%n_centers
!
         symbol  = trim(ao%centers(I)%symbol)
!
         do A = 1, ao%centers(I)%n_shells
!
            do w = ao%centers(I)%shells(A)%first, ao%centers(I)%shells(A)%get_last()
!
               w_red = w - ao%centers(I)%shells(A)%first + 1
!
               ang_mom = ao%centers(I)%shells(A)%get_angular_momentum_label(w_red, &
                                                            ao%centers(I)%cartesian)
!
!              Setup the string for printf to print the right number of reals
!
               n_format_string = repeat('  ' // format_string, m)
!
               call out_file%printf('n', '(i4) (i4) (b4)  '// ang_mom // n_format_string,   &
                                     reals=[C(w, 1 : m)],                                   &
                                     ints=[w, I], chars=[symbol], ll=87)
!
            enddo
!
         enddo
!
      enddo
!
      call out_file%print_separator(pl='normal', n=87, symbol='-')
!
   end subroutine print_ao_vectors_ao_tool
!
!
   subroutine construct_cs_eri_max_screenings(ao)
!!
!!    Construct Cauchy-Schwarz electron repulsion integral (CS-ERI) max screenings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of g_wxwx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively. This array is sorted from largest to smallest and the
!!    sorting indices stored for convenient access.
!!
!!    In addition, the routine uses the array to determine the number of
!!    shell pairs that are significant (with respect to the ERI cutoff value).
!!    This is based on the Cauchy-Schwarz inequality
!!
!!       abs ( g_wxyz ) <= g_wxwx^1/2 * g_yzyz^1/2 <= g_wxwx^1/2 * (max g)
!!
      use timings_class, only: timings
      use array_utilities, only: quicksort_with_index_descending, get_abs_max
!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
!
      integer :: A, B, AB
!
      real(dp) :: maximum
!
      real(dp), dimension(ao%max_sh_size**4) :: g
!
      type(timings), allocatable :: timer, timer_sort
!
      timer = timings('Construct Cauchy-Schwarz ERI screenings (total time)', 'v')
      call timer%turn_on()
!
!     Compute the maximum element in each shell pair, and store indices
!
!$omp parallel do private(A, B, AB, g, maximum) schedule(dynamic)
      do A = 1, ao%n_sh
         do B = 1, A
!
            AB = (max(A,B)*(max(A,B)-3)/2) + A + B
!
            call ao%get_eri(g, A, B, A, B)
!
            maximum = get_abs_max(g, ( ao%shells(A)%length * ao%shells(B)%length )**2 )
!
            ao%cs_eri_max(AB, 2) = sqrt(maximum)
!
            ao%cs_eri_max_indices(AB, 1) = A
            ao%cs_eri_max_indices(AB, 2) = B
!
         enddo
      enddo
!$omp end parallel do
!
!     Sort the Cauchy-Schwarz maximum elements from largest to smallest, and store sorting indices
!
      call dcopy(ao%n_sh*(ao%n_sh + 1)/2,    &
                 ao%cs_eri_max(:, 2),        &
                 1,                          &
                 ao%cs_eri_max(:, 1),        &
                 1)
!
      timer_sort = timings('Construct Cauchy-Schwarz ERI screenings (time to sort)','v')
      call timer_sort%turn_on()
!
      call quicksort_with_index_descending(ao%cs_eri_max(:, 1),         & ! On exit, sorted array
                                           ao%cs_eri_max_indices(:, 3), & ! On exit, sorting indices
                                           ao%n_sh*(ao%n_sh + 1)/2)
!
      call timer_sort%turn_off()
!
!     Count the number of significant ERI shell-pairs
!
      ao%n_sig_eri_shp = 0
!
      do AB = 1, ao%n_sh*(ao%n_sh + 1)/2
!
!        Check if the upper bound of (AB / CD) integrals is below ERI cutoff; if so,
!        we do not consider AB to be a significant shell pair
!
         if (ao%cs_eri_max(AB, 1) * ao%cs_eri_max(1, 1) .lt. ao%eri_cutoff) then
!
            exit
!
         else
!
            ao%n_sig_eri_shp = ao%n_sig_eri_shp + 1
!
         endif
!
      enddo
!
      call timer%turn_off()
!
   end subroutine construct_cs_eri_max_screenings
!
!
   subroutine get_frozen_cores_and_n_frozen_orbitals_ao_tool(ao,           &
                                                             frozen,       &
                                                             n_frozen_ao,  &
                                                             first,        &
                                                             last)
!!
!!    Get frozen cores and number of frozen orbitals
!!    Written by Sarai D. Folkestad, Sep 2018
!!
!!    Determines which centers to freeze (given by the logical array 'frozen') and
!!    the total number of frozen AOs. It considers only centers 'first' to 'last'.
!!
!!    Extened by Alexander C. Paul, Jan 2021
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: first, last
!
      logical, dimension(last - first + 1), intent(out) :: frozen
!
      integer, intent(out) :: n_frozen_ao
!
      integer :: I
!
      frozen = .false.
!
      n_frozen_ao = 0
!
      do I = first, last
!
         if (ao%centers(I)%number_ .ge. 5 .and. ao%centers(I)%number_ .le. 12) then
!
            n_frozen_ao = n_frozen_ao + 1
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 13 .and. ao%centers(I)%number_ .le. 30) then
!
            n_frozen_ao = n_frozen_ao + 5
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 31 .and. ao%centers(I)%number_ .le. 38) then
!
            n_frozen_ao = n_frozen_ao + 9
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 39 .and. ao%centers(I)%number_ .le. 48) then
!
            n_frozen_ao = n_frozen_ao + 14
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 49 .and. ao%centers(I)%number_ .le. 70) then
!
            n_frozen_ao = n_frozen_ao + 18
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 71 .and. ao%centers(I)%number_ .le. 80) then
!
            n_frozen_ao = n_frozen_ao + 23
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 81 .and. ao%centers(I)%number_ .le. 103) then
!
            n_frozen_ao = n_frozen_ao + 34
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .ge. 104 .and. ao%centers(I)%number_ .le. 105) then
!
            n_frozen_ao = n_frozen_ao + 50
            frozen(I)   = .true.
!
         elseif (ao%centers(I)%number_ .gt. 105) then
!
            call output%error_msg('No support for frozen core for Z > 105.')
!
         endif
!
      enddo
!
   end subroutine get_frozen_cores_and_n_frozen_orbitals_ao_tool
!
!
   subroutine evaluate_aos_at_point_ao_tool(ao, x, y, z, aos_at_point)
!!
!!    Evaluate AOs at point
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
!!    Evaluates all the aos at a point on a grid.
!!
!!    'x', 'y', 'z'  : coordinates of the grid point
!!
!!    'aos_at_point' : Value of the aos at the grid point (vector of dim n_ao)
!!
      use math_utilities, only: double_factorial, binomial, factorial
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp), intent(in) :: x, y, z
      real(dp), dimension(ao%n), intent(out) :: aos_at_point
!
      integer :: w, n_primitives, i, j, shell, center, t, u, two_v, two_v_m, two_v_max, t_max
      integer :: floored_term, l, ml, abs_ml, count_ml
!
      real(dp) :: x_rel, y_rel, z_rel, r_squared, radial_part, angular_part, N_S_lm, C_lm_tuv
      real(dp) :: coefficient_i, coefficient_j, exponent_i, exponent_j
      real(dp) :: normalization_constant_i, overlap_primitives
!
      real(dp), parameter :: sqrt_three      = sqrt(three)
      real(dp), parameter :: sqrt_three_half = sqrt(three*half)
      real(dp), parameter :: sqrt_fifteen    = sqrt(15.0d0)
      real(dp), parameter :: sqrt_five_half  = sqrt(five*half)
!
      aos_at_point = zero
!
      w = 1 ! AO index
!
!     Loop over atoms
!
      do center = 1, ao%n_centers
!
!        Find the position of the point relative to the nucleus, in atomic units
!
         x_rel = (x - ao%centers(center)%coordinates(1)) * angstrom_to_bohr
         y_rel = (y - ao%centers(center)%coordinates(2)) * angstrom_to_bohr
         z_rel = (z - ao%centers(center)%coordinates(3)) * angstrom_to_bohr
!
!        Calculate distance between the point and the nucleus
!
         r_squared = x_rel**2 + y_rel**2 + z_rel**2
!
         do shell = 1, ao%centers(center)%n_shells
!
!           Determine angular momentum
!
            l = ao%centers(center)%shells(shell)%l
!
!           Determine radial part
!
            n_primitives = ao%centers(center)%shells(shell)%get_n_primitives()
!
            overlap_primitives = zero
!
            do i = 1, n_primitives
               do j = 1, n_primitives
!
                  exponent_i    = ao%centers(center)%shells(shell)%get_exponent_i(i)
                  exponent_j    = ao%centers(center)%shells(shell)%get_exponent_i(j)
                  coefficient_i = ao%centers(center)%shells(shell)%get_coefficient_i(i)
                  coefficient_j = ao%centers(center)%shells(shell)%get_coefficient_i(j)
!
                  overlap_primitives = overlap_primitives                  &
                                       + coefficient_i*coefficient_j       &
                                       *sqrt(sqrt(exponent_i*exponent_j)  &
                                       /(exponent_i + exponent_j))**(3 + 2*l)
!
               enddo
            enddo
!
            overlap_primitives = overlap_primitives*(two**(three*half + real(l, dp)))
!
            radial_part = zero
!
            do i = 1, n_primitives
!
!              Evaluate linear combination of normalized primitives at point
!
               exponent_i    = ao%centers(center)%shells(shell)%get_exponent_i(i)
               coefficient_i = ao%centers(center)%shells(shell)%get_coefficient_i(i)
!
!              Normalization constant for primitive containing a Racah's normalized angular part
!              (from eqn. (94) without Racah's normalization constant and (95)
!              in Giesea, T. J. HSERILib: Gaussian integral evaluation)
!
               normalization_constant_i = (four*exponent_i)**(real(l, dp)*half + three*quarter)
!
               radial_part = radial_part &
                             + normalization_constant_i*coefficient_i*exp(-exponent_i*r_squared)
!
            enddo
!
            radial_part = radial_part/(sqrt(overlap_primitives*((two*pi)**(three*half)&
                                       *(real(double_factorial(2*l-1), dp)))))
!
!           Construct Racah's normalized orbitals
!
            if (l == 0) then
!
!              Cartesian s-shell
!
               aos_at_point(w) = radial_part
!
            elseif (l == 1) then
!
!              Cartesian p-shell
!
!              p_x, p_y, p_z (first, second, third in CCA standard)
!
               aos_at_point(w)   = radial_part*x_rel
               aos_at_point(w+1) = radial_part*y_rel
               aos_at_point(w+2) = radial_part*z_rel
!
            elseif (l == 2) then
!
!              ml = -2, angular part: sqrt(3)xy
!
               aos_at_point(w)   = sqrt_three*x_rel*y_rel*radial_part
!
!              ml = -1, angular part: sqrt(3)yz
!
               aos_at_point(w+1) = sqrt_three*y_rel*z_rel*radial_part
!
!              ml = 0, angular part: 1.5z^2 - 0.5r^2
!
               aos_at_point(w+2) = half*(three*z_rel**2 - r_squared)*radial_part
!
!              ml = 1, angular part: sqrt(3)xz
!
               aos_at_point(w+3) = sqrt_three*x_rel*z_rel*radial_part
!
!              ml = 2, angular part: 1/2 sqrt(3)(x^2 - y^2)
!
               aos_at_point(w+4) = half*sqrt_three*(x_rel**2 - y_rel**2)*radial_part
!
            elseif (l == 3) then
!
!              ml = -3, angular part: 1/2 sqrt(5/2) (3x^2 - y^2)y
!
               aos_at_point(w) = half*sqrt_five_half*(three*x_rel**2 - y_rel**2)*y_rel*radial_part
!
!              ml = -2, angular part: sqrt(15) xyz
!
               aos_at_point(w+1) = sqrt_fifteen*x_rel*y_rel*z_rel*radial_part
!
!              ml = -1, angular part: 1/2 sqrt(3/2) (5*z**2 - r**2)y
!              ml =  1, angular part: 1/2 sqrt(3/2) (5*z**2 - r**2)x
!
               aos_at_point(w+2) = half*sqrt_three_half*(five*z_rel**2 - r_squared)*radial_part
               aos_at_point(w+4) = aos_at_point(w+2)*x_rel
               aos_at_point(w+2) = aos_at_point(w+2)*y_rel
!
!              ml = 0, angular part: 1/2(5z**2 - 3r^2)z
!
               aos_at_point(w+3) = half*(five*z_rel**2 - three*r_squared)*z_rel*radial_part
!
!              ml = 2, angular_part: 1/2 sqrt(15)(x^2-y^2)z
!
               aos_at_point(w+5) = half*sqrt_fifteen*(x_rel**2 - y_rel**2)*z_rel*radial_part
!
!              ml = 3, angular part: 1/2 sqrt(5/2)(x^2-3y^2)x
!
               aos_at_point(w+6) = half*sqrt_five_half*(x_rel**2 - three*y_rel**2)&
                                    *x_rel*radial_part
!
            elseif (l > 3) then
!
!              Solid harmonic shell
!
               count_ml = 0
!
               do ml = -l, l, 1
!
                  abs_ml = abs(ml)
!
!                 Solid harmonic from Molecular electronic structure theory eqn. (6.4.47)
!
                  t_max = (l - abs_ml)/2 ! floors (l - abs_ml)/2, since abs_ml <= the non-negative l
!
                  if (ml .ge. 0) then
                     two_v_m = 0
                  else
                     two_v_m = 1
                  endif
!
                  floored_term = (abs_ml - two_v_m)/2 ! floors (abs_ml - two_v_m)/2,
!                                                       since 2*v_m <= the non-negative abs_ml
                  two_v_max = 2*floored_term + two_v_m
!
!                 Calculate value of solid harmonic at point
!
                  angular_part = zero
!
                  do t = 0, t_max
                     do u = 0, t
                        do two_v = two_v_m, two_v_max, 2
!
!                          C_lm_tuv from Molecular electronic structure theory eqn.
!                          (6.4.48), where (two_v - two_v_m)/2 is an int
!
                           C_lm_tuv = (-one)**(t + (two_v - two_v_m)/2)*quarter**t &
                                      *real(binomial(l, t)*binomial(l-t, abs_ml+t)*&
                                          binomial(t, u)*binomial(abs_ml, two_v), dp)
!
                           angular_part = angular_part                     &
                                 + C_lm_tuv*x_rel**(2*t+abs_ml-2*u-two_v)  &
                                 *y_rel**(2*u+two_v)                       &
                                 *z_rel**(l-2*t-abs_ml)

!
                        enddo
                     enddo
                  enddo
!
!                 Multiply by N_S_lm from Molecular electronic structure theory eqn. (6.4.49)
!
                  N_S_lm = one/real(2**abs_ml*factorial(l), dp)*sqrt(real(2*factorial(l + abs_ml)&
                                       *factorial(l - abs_ml), dp))
!
                  if (ml == 0) N_S_lm = N_S_lm/sqrt(two)
!
                  angular_part = N_S_lm*angular_part
!
                  aos_at_point(w + count_ml) = radial_part*angular_part
!
                  count_ml = count_ml + 1
!
               enddo
!
            endif
!
            w = w + 2*l + 1
!
         enddo ! end loop over shells
!
      enddo ! end loop over atoms
!
   end subroutine evaluate_aos_at_point_ao_tool
!
!
   function get_center_coordinates_ao_tool(ao) result(R_qk)
!!
!!    Get center coordinates
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns the center coordinates R_qk, where q = 1, 2, 3 (corresponding to x, y, and z)
!!    and k = 1, 2, 3, ..., n_centers.
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      real(dp), dimension(3, ao%n_centers) :: R_qk
!
      integer :: I
!
      do I = 1, ao%n_centers
!
         R_qk(:, I) = ao%centers(I)%coordinates
!
      enddo
!
   end function get_center_coordinates_ao_tool
!
!
   subroutine get_center_symbols_ao_tool(ao, symbols)
!!
!!    Get center symbols
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Returns the center symbols
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      character(len=2), dimension(ao%n_centers), intent(out) :: symbols
!
      integer :: I
!
      do I = 1, ao%n_centers
!
         symbols(I) = ao%centers(I)%symbol
!
      enddo
!
   end subroutine get_center_symbols_ao_tool
!
!
   subroutine get_point_charges_ao_tool(ao, pc)
!!
!!    Get point charges
!!    Written by Sarai D.Folkestad, 2020
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      type(point_charges), intent(out) :: pc
!
      integer :: I
!
      pc = point_charges(ao%n_centers)
      call pc%initialize()
!
      pc%r = ao%get_center_coordinates()
!
      do I = 1, ao%n_centers
!
         pc%q(I) = real(ao%centers(I)%nuclear_charge,kind=dp)
!
      enddo
!
   end subroutine get_point_charges_ao_tool
!
!
   subroutine get_subset_point_charges_ao_tool(ao, pc, name_, include_higher_priority)
!!
!!    Get subset point charges
!!    Written by Sarai D.Folkestad, 2020
!!
!!    Returns the point charges of a subset
!!    or the union of subsets
!!
!!    name_                   : Name of the subset
!!
!!    include_higher_priority : Determines if subsets of higher priority
!!                            should also be included.
!!
!!                            If include_higher_priority = .true.,
!!                            The named subset is the last (lowest priority)
!!                            subset that is included
!!
!!                            If include_higher_priority = .false.,
!!                            The named subset is the only subset that is included
!!
!!
      implicit none
!
      class(ao_tool), intent(in)       :: ao
      character(len=*), intent(in)     :: name_
      logical, intent(in)              :: include_higher_priority
!
      type(point_charges), intent(out) :: pc
!
      integer                          :: I, first, last, subset_index
!
      subset_index = ao%get_subset_index(name_)
!
      if (subset_index == 0) call output%error_msg('did not find subset of name '// trim(name_))
!
      if (include_higher_priority) then
         first = 1
      else
         first = ao%center_subsets(subset_index)%first
      endif
!
      last = ao%center_subsets(subset_index)%get_last()
!
      pc = point_charges(last - first + 1)
      call pc%initialize()
!
      do I = first, last
!
         pc%r(:,I) = ao%centers(I)%coordinates(:)
         pc%q(I)   = real(ao%centers(I)%nuclear_charge, kind = dp)
!
      enddo
!
   end subroutine get_subset_point_charges_ao_tool
!
!
   pure function get_subset_index(ao, name_) result(subset_index)
!!
!!    Get subset index
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Returns the index of the subset with name = name_
!!
!!    Returns 0, if no such set exists
!!
      implicit none
!
      class(ao_tool),   intent(in) :: ao
      character(len=*), intent(in) :: name_
      integer                      :: subset_index
      integer                      :: I
!
      subset_index = 0
!
      do I = 1, ao%n_center_subsets
!
         if (ao%center_subsets(I)%get_name() == trim(name_)) subset_index = I
!
      enddo
!
   end function get_subset_index
!
!
   function get_n_electrons_ao_tool(ao) result(n)
!!
!!    Get number of electrons
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Gets the number of electrons based on the total charge specified and the
!!    charges at the atomic centers
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer :: n
!
      integer :: I
!
      n = 0
!
      do I = 1, ao%n_centers
!
         n = n + ao%centers(I)%nuclear_charge
!
      enddo
!
      n = n - ao%charge
!
   end function get_n_electrons_ao_tool
!
!
   subroutine copy_geometry_dependent_variables(ao, template)
!!
!!    Copy geometry dependent variables
!!    Written by Alexander C. Paul, Feb 2021
!!
      use array_utilities, only: copy_integer
!
      implicit none
!
      class(ao_tool), intent(inout) :: ao
      type(ao_tool),  intent(in)    :: template
!
      call dcopy(ao%n**2, template%h, 1, ao%h, 1)
      call dcopy(ao%n**2, template%s, 1, ao%s, 1)
!
      ao%n_orthonormal_ao = template%n_orthonormal_ao
!
      call mem%alloc(ao%L, ao%n_orthonormal_ao, ao%n_orthonormal_ao)
      call mem%alloc(ao%P, ao%n, ao%n_orthonormal_ao)
!
      call dcopy(ao%n_orthonormal_ao**2, template%L, 1, ao%L, 1)
      call dcopy(ao%n*ao%n_orthonormal_ao, template%P, 1, ao%P, 1)
!
      ao%n_sig_eri_shp = template%n_sig_eri_shp
!
      call dcopy(2*ao%n_sh*(ao%n_sh + 1)/2, template%cs_eri_max, 1, &
                 ao%cs_eri_max, 1)
!
      call copy_integer(template%cs_eri_max_indices, ao%cs_eri_max_indices, &
                        3*ao%n_sh*(ao%n_sh + 1)/2)
!
   end subroutine copy_geometry_dependent_variables
!
!
   subroutine destructor(ao)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      type(ao_tool) :: ao
!
      integer :: I
!
      if (allocated(ao%h)) call mem%dealloc(ao%h, ao%n, ao%n)
      if (allocated(ao%s)) call mem%dealloc(ao%s, ao%n, ao%n)
      if (allocated(ao%v)) call mem%dealloc(ao%v, ao%n, ao%n)
!
      if (allocated(ao%centers)) then
!
         do I = 1, ao%n_centers
!
            call ao%centers(I)%cleanup()
!
         enddo
!
         deallocate(ao%centers)
!
      endif
!
      if (allocated(ao%shell_to_center))  call mem%dealloc(ao%shell_to_center, ao%n_sh)
      if (allocated(ao%ao_to_center))     call mem%dealloc(ao%ao_to_center, ao%n)
      if (allocated(ao%ao_to_shell))      call mem%dealloc(ao%ao_to_shell, ao%n)
!
      if (allocated(ao%L)) call mem%dealloc(ao%L, ao%n_orthonormal_ao, ao%n_orthonormal_ao)
      if (allocated(ao%P)) call mem%dealloc(ao%P, ao%n, ao%n_orthonormal_ao)
!
      if (allocated(ao%cs_eri_max)) &
         call mem%dealloc(ao%cs_eri_max, ao%n_sh*(ao%n_sh + 1)/2, 2)
!
      if (allocated(ao%cs_eri_max_indices)) &
         call mem%dealloc(ao%cs_eri_max_indices, ao%n_sh*(ao%n_sh + 1)/2, 3)
!
   end subroutine destructor
!
!
   subroutine construct_oei_screened(ao, oei_type, x, n_components)
!!
!!    Construct OEI screened (one-electron integral)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2021
!
      use iso_c_binding
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(ao%n, ao%n, n_components), intent(out) :: x
!
      character(len=*) :: oei_type
!
      integer :: w, y, k, w_f, y_f
!
      integer :: A, B, n_shp, AB
      integer(c_int) :: A_c, B_c
!
      real(dp) :: max_s
!
      real(dp), dimension(:,:,:), pointer                           :: x_ABk_p
      real(dp), dimension((ao%max_sh_size)**2*n_components), target :: x_ABk
!
      character(len=100, kind=c_char) :: oei_type_c
!
      integer, dimension(:,:), allocatable :: shp_list
!
      if (trim(oei_type) == 'overlap') &
         call output%error_msg('cannot use overlap screening on overlap')
!
      call zero_array(x, ao%n**2*n_components)
!
!     Count the number of significant shell pairs
!
      call mem%alloc(shp_list, ao%n_sh**2, 2)
!
      n_shp = 0
!
      do A = 1, ao%n_sh
         do B = 1, A
!
            if (ao%shp_on_same_atom(A,B)) then
!
               n_shp = n_shp + 1
               shp_list(n_shp, 1) = A
               shp_list(n_shp, 2) = B
!
            else
!
               max_s = maxval(abs(ao%s(ao%shells(A)%first:ao%shells(A)%get_last(),&
                                       ao%shells(B)%first:ao%shells(B)%get_last())))
!
               if (max_s .gt. ao%oei_cutoff) then
!
                  n_shp = n_shp + 1
                  shp_list(n_shp, 1) = A
                  shp_list(n_shp, 2) = B
!
               endif
!
            endif
!
         enddo
      enddo
!
      oei_type_c = trim(oei_type) // c_null_char
!
!$omp parallel do private(AB, A_c, B_c, x_ABk, x_ABk_p, w, y, w_f, y_f, k) schedule(dynamic)
      do AB = 1, n_shp
!
            A_c = int(shp_list(AB, 1), kind=c_int)
            B_c = int(shp_list(AB, 2), kind=c_int)
!
           call get_oei_c(oei_type_c, x_ABk, A_c, B_c)
!
           x_ABk_p(1 : ao%shells(A_c)%length, &
                   1 : ao%shells(B_c)%length, &
                   1 : n_components)          &
           => x_ABk(1 : ao%shells(A_c)%length * &
                        ao%shells(B_c)%length * &
                        n_components)
!
           do k = 1, n_components
              do w = 1, ao%shells(A_c)%length
                 do y = 1, ao%shells(B_c)%length
!
                    w_f = ao%shells(A_c)%first - 1 + w
                    y_f = ao%shells(B_c)%first - 1 + y
!
                    x(w_f, y_f, k) = x_ABk_p(w, y, k)
                    x(y_f, w_f, k) = x_ABk_p(w, y, k)
!
                 enddo
              enddo
           enddo
         enddo
!$omp end parallel do
!
      call mem%dealloc(shp_list, ao%n_sh**2, 2)
!
   end subroutine construct_oei_screened
!
!
   function shp_on_same_atom(ao, A, B) result(on_same_atom)
!!
!!    Shell pair on same atom
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Checks if shells A and B are on the same atom
!!
      implicit none

      class(ao_tool), intent(in) :: ao
      integer, intent(in) :: A, B
      logical :: on_same_atom
!
      on_same_atom = (ao%shell_to_center(A) == ao%shell_to_center(B))
!
   end function shp_on_same_atom
!
!
   subroutine print_z_matrix_ao_tool(ao)
!!
!!    Print Z-matrix
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Prints the Z-matrix in angstrom.
!!
!
      use z_matrix_tool_class, only: z_matrix_tool
!
      implicit none
!
      class(ao_tool), intent(in)   :: ao
!
      type(z_matrix_tool), allocatable :: z_matrix
!
      real(dp), dimension(:,:), allocatable :: R
      character(len=2), dimension(ao%n_centers) :: symbols
!
      z_matrix = z_matrix_tool(ao%n_centers)
      call z_matrix%initialize()
!
      call mem%alloc(R, 3, ao%n_centers)
!
      R =  ao%get_center_coordinates()
      call ao%get_center_symbols(symbols)
!
      call z_matrix%construct(R, symbols)
      call z_matrix%print_(output)
!
      call mem%dealloc(R, 3, ao%n_centers)
!
      call z_matrix%cleanup_z_matrix()
!
   end subroutine print_z_matrix_ao_tool
!
!
   pure function has_ghost_atoms_ao_tool(ao) result(has_ghosts)
!!
!!    Has ghost atoms?
!!    Written by Tor S. Haugland, May 2021
!!
!!    Are ghost atoms present in system?
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      logical :: has_ghosts
      integer :: I
!
      has_ghosts = .false.
      do I = 1, ao%get_n_centers()
         has_ghosts = has_ghosts .or. ao%centers(I)%is_ghost()
      enddo
!
   end function has_ghost_atoms_ao_tool
!
!
   function get_molden_ao_indices_ao_tool(ao) result(map)
!!
!!    Get Molden AO indices
!!    Written by Alexander C. Paul, May 2021
!!
!!    Molden expects e.g. the d-orbitals to be ordered
!!    as xx, yy, zz, xy, xz, yz, similarly for f and g orbitals
!!    (c.f. angular_momentum.F90)
!!    This routine creates a map so that map(10) returns the index
!!    of the 10th AO as expected by molden
!!
      implicit none
!
      class(ao_tool), intent(in) :: ao
!
      integer, dimension(ao%n) :: map
!
      type(range_), allocatable :: aos
      integer :: c
!
      do c = 1, ao%n_centers
!
         aos = ao%centers(c)%get_ao_range()
         map(aos%first:aos%get_last()) = &
            ao%centers(c)%get_ao_molden_order(aos%first, aos%get_last())
!
      end do
!
   end function get_molden_ao_indices_ao_tool
!
!
   subroutine print_molden_geometry_ao_tool(ao, file_)
!!
!!    Print molden geometry
!!    Written by Alexander C. Paul, May 2021
!!
!!    Geometry in bohr in the format:
!!    Element  number  "Atomic number"  x  y  z
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
      type(output_file), intent(inout) :: file_
!
      integer :: counter, c
!
      counter = 0
      do c = 1, ao%n_centers
         counter = counter + 1
         call file_%printf('m', '(b2)   (i5)   (i3)   (f15.10)  (f15.10)  (f15.10)', &
                            chars=[ao%centers(c)%symbol], fs='(t1,a)',               &
                            reals=[ao%centers(c)%coordinates*angstrom_to_bohr],      &
                            ints=[counter, ao%centers(c)%number_])
      end do
!
   end subroutine print_molden_geometry_ao_tool
!
!
   subroutine print_basis_set_molden_ao_tool(ao, file_)
!!
!!    Print basis set molden
!!    Written by Alexander C. Paul, May 2021
!!
!!    Prints basis set for a molden file:
!!    Atom number
!!    shell-label n_primitives
!!    exponent coefficient
!!    "empty line"
!!
!!    The atom number is the number in eT (as printed in print_molden_geometry)
!!    The shell label is the angular momentum converted to a letter (s,p,d,f,g)
!!
!!    If a spherical basis set is used we need to add [5D7F] and [9G].
!!    NB: Mixed spherical and cartesian basis sets are not supported
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(ao_tool), intent(in) :: ao
      type(output_file), intent(inout) :: file_
!
      integer :: counter, c, s, i, n
      logical :: cartesian
!
      counter = 0
      do c = 1, ao%n_centers
!
         counter = counter + 1
!
         call file_%printf('m', '(i0)', fs='(t1,a)', ints=[counter])
!
         do s = 1, ao%centers(c)%n_shells
!
            n = ao%centers(c)%shells(s)%get_n_primitives()
            call file_%printf('m', '(a0) (i0)', ints=[n], fs='(t1,a)', &
                 chars=[ao%centers(c)%shells(s)%get_angular_momentum()])
!
            do i = 1, n
               call file_%printf('m', '(f16.7)   (f15.10)', fs='(t1,a)', &
                                  reals=[ao%centers(c)%shells(s)%get_exponent_i(i), &
                                         ao%centers(c)%shells(s)%get_coefficient_i(i)])
            end do
         end do
!
         call file_%printf('m', '')
!
      end do
!
!     Either only spherical or only cartesian
!
      cartesian = ao%centers(1)%cartesian
!
      do c = 1, ao%n_centers
         if (ao%centers(1)%cartesian .neqv. ao%centers(c)%cartesian) then
            call output%error_msg('Mixed cartesian and spherical basis sets &
                                  &are not supported in molden-file')
         end if
      end do
!
      if (.not. cartesian) then
         call file_%printf('m', '[5D7F]', fs='(t1,a)')
         call file_%printf('m', '[9G]', fs='(t1,a)')
      end if
!
   end subroutine print_basis_set_molden_ao_tool
!
!
   subroutine get_SAD_center_indices_ao_tool(ao, center_indices)
!!
!!    Get SAD center indices
!!    Written by Tor S. Haugland and Sarai D. Folkestad, Jun 2021
!!
!!    Collects the indices of unique atoms, excluding ghost atoms.
!!
!
      use string_utilities, only: index_of_unique_strings
!
      implicit none
!
      class(ao_tool),                           intent(in) :: ao
      integer, dimension(ao%n_centers),         intent(out) :: center_indices
!
      integer :: n_centers, I
!
      character(len=50), dimension(ao%n_centers) :: atom_and_basis
      integer, dimension(ao%n_centers) :: unique_center_indices
!
      do I = 1, ao%n_centers
!
         atom_and_basis(I) = ao%centers(I)%get_identifier_string()
!
      enddo
!
      call index_of_unique_strings(unique_center_indices, ao%n_centers, atom_and_basis)
!
      center_indices = 0
      n_centers = 0
!
      do I = 1, ao%n_centers
!
         if (all(unique_center_indices /= I)) cycle
         if (ao%centers(I)%is_ghost()) cycle ! No density to generate for ghost atom
!
         n_centers = n_centers + 1
         center_indices(n_centers) = I
!
      enddo
!
   end subroutine get_SAD_center_indices_ao_tool
!
!
end module ao_tool_class
