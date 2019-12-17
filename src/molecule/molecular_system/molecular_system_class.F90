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
module molecular_system_class
!
!!
!!    Molecule class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use parameters
   use libint_initialization
   use active_atoms_class, only : active_atoms
!
   use global_in, only : input
   use global_out, only : output
   use sequential_file_class, only : sequential_file
   use direct_file_class, only : direct_file
   use output_file_class, only : output_file
   use memory_manager_class, only : mem
   use string_utilities, only : convert_to_lowercase
   use interval_class, only : interval
   use atomic_class, only : atomic
   use mm_class, only : mm
   use pcm_class, only : pcm
!
   implicit none
!
   include "../../libint/atom_init_cdef.F90"
!
   type :: molecular_system
!
      character(len=200) :: name_
      character(len=100), dimension(:), allocatable :: basis_sets
!
      character(len=100) :: coordinate_units = 'angstrom'
!
      integer :: n_atoms
      integer :: n_basis_sets 
      integer :: charge
      integer :: multiplicity
      integer :: n_electrons 
      integer :: n_s
!
      integer :: cartesian_gaussians_int
!
      type(atomic), dimension(:), allocatable :: atoms
!
      type(interval), dimension(:), allocatable :: shell_limits 
!
      integer :: n_active_atom_spaces                                    ! Number of active spaces
      type(active_atoms), dimension(:),allocatable :: active_atom_spaces ! Array of active spaces
!
      integer, dimension(:), allocatable :: shell2atom 
!
      integer :: max_shell_size
!     
!     CGTO information                                           
!                                                              
      integer :: n_cart_basis = 0                              
      integer :: n_pure_basis = 0                              
      integer :: n_primitives_cart = 0                         
      logical :: cartesian_basis = .false. 
      logical :: mm_calculation = .false.
      logical :: pcm_calculation = .false.
!
      type(mm) :: mm
      type(pcm) :: pcm
!
!     AO Cholesky vectors
!
      type(direct_file) :: ao_cholesky_file
!
      integer :: n_J = 0
!
   contains
!
      procedure, private :: prepare                         => prepare_molecular_system
      procedure :: cleanup                                  => cleanup_molecular_system
!
      procedure, private :: write_libint_files              => write_libint_files_molecular_system
      procedure, private :: delete_libint_files             => delete_libint_files_molecular_system
!
      procedure, private :: read_settings                   => read_settings_molecular_system
      procedure, private :: read_system                     => read_system_molecular_system
      procedure, private :: read_geometry                   => read_geometry_molecular_system
      procedure, private :: read_active_atoms               => read_active_atoms_molecular_system
!
      procedure :: print_system                             => print_system_molecular_system
      procedure :: print_geometry                           => print_geometry_molecular_system
      procedure :: print_active_atoms                       => print_active_atoms_molecular_system
!
      procedure :: get_geometry                             => get_geometry_molecular_system
      procedure :: set_geometry                             => set_geometry_molecular_system
!
      procedure :: get_total_nuclear_repulsion              => get_total_nuclear_repulsion_molecular_system
!
      procedure :: get_nuclear_repulsion                    => get_nuclear_repulsion_molecular_system
      procedure :: get_nuclear_repulsion_1der_numerical     => get_nuclear_repulsion_1der_numerical_molecular_system
      procedure :: get_nuclear_repulsion_1der               => get_nuclear_repulsion_1der_molecular_system
!
      procedure :: get_n_electrons                          => get_n_electrons_molecular_system
      procedure :: get_nuclear_dipole                       => get_nuclear_dipole_molecular_system
      procedure :: get_nuclear_quadrupole                   => get_nuclear_quadrupole_molecular_system
      procedure :: get_nuclear_dipole_active                => get_nuclear_dipole_active_molecular_system
      procedure :: get_nuclear_quadrupole_active            => get_nuclear_quadrupole_active_molecular_system
!
      procedure :: get_n_aos                                => get_n_aos_molecular_system
      procedure :: get_n_shells                             => get_n_shells_molecular_system
      procedure :: get_shell_limits                         => get_shell_limits_molecular_system
      procedure :: basis2shell                              => basis2shell_molecular_system
      procedure :: get_max_shell_size                       => get_max_shell_size_molecular_system
!
      procedure :: shell_to_atom                            => shell_to_atom_molecular_system
!
      procedure :: check_if_basis_present_and_pure             => check_if_basis_present_and_pure_molecular_system
      procedure :: initialize_libint_atoms_and_bases           => initialize_libint_atoms_and_bases_molecular_system
      procedure, nopass :: initialize_libint_integral_engines  => initialize_libint_integral_engines_molecular_system
!
      procedure :: initialize_basis_sets                    => initialize_basis_sets_molecular_system
      procedure :: initialize_atoms                         => initialize_atoms_molecular_system
      procedure :: initialize_shell_limits                  => initialize_shell_limits_molecular_system
      procedure :: initialize_shell2atom                    => initialize_shell2atom_molecular_system
!
      procedure :: destruct_basis_sets                      => destruct_basis_sets_molecular_system
      procedure :: destruct_atoms                           => destruct_atoms_molecular_system
      procedure :: destruct_shell_limits                    => destruct_shell_limits_molecular_system
      procedure :: destruct_shell2atom                      => destruct_shell2atom_molecular_system
!
      procedure :: translate_from_input_order_to_eT_order   => translate_from_input_order_to_eT_order_molecular_system
!
      procedure :: construct_ao_v_wx                        => construct_ao_v_wx_molecular_system      
      procedure :: construct_ao_h_wx                        => construct_ao_h_wx_molecular_system     
      procedure :: construct_ao_h_wx_kinetic_1der           => construct_ao_h_wx_kinetic_1der_molecular_system     
      procedure :: construct_ao_h_wx_nuclear_1der           => construct_ao_h_wx_nuclear_1der_molecular_system
! 
      procedure :: construct_ao_g_wxyz_1der                 => construct_ao_g_wxyz_1der_molecular_system
      procedure :: construct_ao_g_wxyz                      => construct_ao_g_wxyz_molecular_system  
      procedure, nopass :: construct_ao_g_wxyz_epsilon      => construct_ao_g_wxyz_epsilon_molecular_system
!      
      procedure :: construct_ao_s_wx                        => construct_ao_s_wx_molecular_system
      procedure :: construct_ao_s_wx_1der                   => construct_ao_s_wx_1der_molecular_system
!   
      procedure :: construct_ao_mu_wx                       => construct_ao_mu_wx_molecular_system     
      procedure :: construct_ao_q_wx                        => construct_ao_q_wx_molecular_system   
!
      procedure :: read_n_active_atoms_for_method           => read_n_active_atoms_for_method_molecular_system 
!
      procedure :: first_and_last_ao_active_space           => first_and_last_ao_active_space_molecular_system 
      procedure :: first_ao_active_space                    => first_ao_active_space_molecular_system
      procedure :: last_ao_active_space                     => last_ao_active_space_molecular_system
!
      procedure :: set_basis_info                           => set_basis_info_molecular_system
      procedure :: set_shell_basis_info                     => set_shell_basis_info_molecular_system
      procedure :: check_convert_pure_to_cartesian_basis    => check_convert_pure_to_cartesian_basis_molecular_system
      procedure :: normalize_raw_primitives                 => normalize_raw_primitives_molecular_system
      procedure :: get_nuclear_repulsion_mm                 => get_nuclear_repulsion_mm_molecular_system
!
      procedure :: basis2atom                               => basis2atom_molecular_system
!
      procedure :: evaluate_aos_at_point     => evaluate_aos_at_point_molecular_system
!
      procedure :: rename_core_valence_dunning_sets   &
                                       => rename_core_valence_dunning_sets_molecular_system
!
   end type molecular_system
!
!
   interface
!
      include "ao_integrals_interface.F90"
!
   end interface 
!
!
   interface molecular_system
!
      procedure :: new_molecular_system
      procedure :: new_molecular_system_from_parameters
!
   end interface
!
!
contains
!
!
   function new_molecular_system() result(molecule)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(molecular_system) :: molecule
!
      molecule%charge = 0
      molecule%multiplicity = 1
      molecule%cartesian_gaussians_int = 0      
      molecule%n_active_atom_spaces = 0
!
      call molecule%read_settings()
!
      call molecule%prepare()
!
   end function new_molecular_system
!
!
   function new_molecular_system_from_parameters(atoms, name_, charge, multiplicity, mm_calculation, &
                                                 pcm_calculation) result(molecule)
!!
!!    Initialize From Parameters
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      type(molecular_system) :: molecule
!
      type(atomic), dimension(:), intent(in) :: atoms
      character(len=100),         intent(in) :: name_
      integer,                    intent(in) :: charge
      integer,                    intent(in) :: multiplicity
      logical,                    intent(in) :: mm_calculation
      logical,                    intent(in) :: pcm_calculation
!
      allocate(molecule%atoms(size(atoms)))
!
      molecule%atoms           = atoms
      molecule%n_atoms         = size(atoms)
      molecule%name_           = name_
      molecule%charge          = charge
      molecule%multiplicity    = multiplicity
      molecule%mm_calculation  = mm_calculation
      molecule%pcm_calculation = pcm_calculation
!
      molecule%cartesian_gaussians_int = 0
      molecule%n_active_atom_spaces = 0
!
      call molecule%prepare()
!
   end function new_molecular_system_from_parameters
!
!
   subroutine read_settings_molecular_system(molecule)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      call molecule%read_system()
      call molecule%read_geometry()
      call molecule%read_active_atoms()
!
      molecule%mm_calculation = input%requested_section('molecular mechanics')
      molecule%pcm_calculation = input%requested_section('pcm')
!
   end subroutine read_settings_molecular_system
!
!
   subroutine prepare_molecular_system(molecule)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use shell_class, only: shell
!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: s 
!
      integer :: n_s, i, j
      integer(i6) :: k
!
      integer(i6), dimension(:), allocatable :: n_shells_on_atoms
      integer(i6), dimension(:), allocatable :: n_basis_in_shells
      integer(i6), dimension(:), allocatable :: first_ao_in_shells
      integer(i6), dimension(:), allocatable :: shell_numbers
!
      real(dp), dimension(:,:), allocatable :: qm_coordinates
      real(dp), dimension(:), allocatable :: qm_charges
!
!     First have a look to the basis set infos
!
      call molecule%rename_core_valence_dunning_sets()
      call molecule%check_if_basis_present_and_pure()
!
!     Initialize libint with atoms and basis sets,
!     then initialize the integral engines 
!
      call molecule%initialize_libint_atoms_and_bases()
      call molecule%initialize_libint_integral_engines()
!
!     Initialize atoms and shells for eT
!
      allocate(n_shells_on_atoms(molecule%n_atoms))
      n_shells_on_atoms = 0
!
      call get_n_shells_on_atoms_c(n_shells_on_atoms)
!
      do k = 1, int(molecule%n_atoms,4) ! Loop over atoms
!
         call molecule%atoms(k)%set_number()
!
!        Allocate and initialize the corresponding shells
!
         molecule%atoms(k)%n_shells = n_shells_on_atoms(k)
         call molecule%atoms(k)%initialize_shells()
!
!        Then determine the number of basis functions in each shell
!        and save number of AOs per atom
!
         allocate(n_basis_in_shells(n_shells_on_atoms(k)))
         allocate(shell_numbers(n_shells_on_atoms(k)))
         allocate(first_ao_in_shells(n_shells_on_atoms(k)))
!
         call get_n_basis_in_shells_c(k, n_basis_in_shells)
         call get_shell_numbers_c(k, shell_numbers)
         call get_first_ao_in_shells_c(k, first_ao_in_shells)
!
         molecule%atoms(k)%n_ao = 0
!
         do j = 1, n_shells_on_atoms(k) 
!
            molecule%atoms(k)%shells(j) = shell(first=int(first_ao_in_shells(j)), &
                                                length=int(n_basis_in_shells(j)), &
                                                number_=int(shell_numbers(j)))
!
            molecule%atoms(k)%n_ao = molecule%atoms(k)%n_ao + &
                                       molecule%atoms(k)%shells(j)%length
!
         enddo
!
         deallocate(n_basis_in_shells)
         deallocate(shell_numbers)
         deallocate(first_ao_in_shells)
!
      enddo
!
      deallocate(n_shells_on_atoms)
!
      call molecule%initialize_shell_limits()
!
      n_s = molecule%get_n_shells()
      do s = 1, n_s 
!  
         molecule%shell_limits(s) = molecule%get_shell_limits(s)
!
      enddo
!
      molecule%n_electrons = molecule%get_n_electrons()
      molecule%n_s         = molecule%get_n_shells()
!
!     Some sanity checks and stops
!
      do i = 1, molecule%n_atoms
!
         if (molecule%atoms(i)%n_ao == 0) then
!   
            call output%error_msg('Is basis defined for all atoms in input?')
!
         endif
!
      enddo
!
      call molecule%set_basis_info                          
      call molecule%check_convert_pure_to_cartesian_basis    
      call molecule%normalize_raw_primitives                 
!
      call molecule%print_system()
!
      call molecule%get_max_shell_size(molecule%max_shell_size)
!
      if (molecule%mm_calculation) then
!
         call molecule%mm%prepare()
!         
      endif
!
!     pcm
!
      if (molecule%pcm_calculation) then
!      
         call mem%alloc(qm_coordinates, 3, molecule%n_atoms)
         call mem%alloc(qm_charges, molecule%n_atoms)
!      
         do i = 1, molecule%n_atoms 
!      
            qm_coordinates(1,i) = molecule%atoms(i)%x * angstrom_to_bohr
            qm_coordinates(2,i) = molecule%atoms(i)%y * angstrom_to_bohr
            qm_coordinates(3,i) = molecule%atoms(i)%z * angstrom_to_bohr
            qm_charges(i)       = molecule%atoms(i)%number_
!      
         enddo
!      
         call molecule%pcm%prepare(molecule%n_atoms, qm_coordinates, qm_charges)
!      
         call mem%dealloc(qm_coordinates, 3, molecule%n_atoms)
         call mem%dealloc(qm_charges, molecule%n_atoms)
!         
      endif
!
      call initialize_shell2atom_c()
!
      call molecule%initialize_shell2atom()
!
      do i = 1, molecule%n_s 
!
         molecule%shell2atom(i) = molecule%shell_to_atom(i)
!
      enddo
!
!
   end subroutine prepare_molecular_system
!
!
   subroutine cleanup_molecular_system(molecule)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: i
!
      do i = 1, molecule%n_atoms 
!
         call molecule%atoms(i)%cleanup()
!
      enddo
!
      call molecule%destruct_atoms()
      call molecule%destruct_basis_sets()
      call molecule%destruct_shell_limits()
      call molecule%destruct_shell2atom()
      call molecule%delete_libint_files()
      if (molecule%mm_calculation) call molecule%mm%cleanup()
      if (molecule%pcm_calculation) call molecule%pcm%cleanup()
!
   end subroutine cleanup_molecular_system
!
!
   subroutine initialize_libint_atoms_and_bases_molecular_system(molecule)
!!
!!    Initialize Libint atoms and bases 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad 2018-2019
!!
      implicit none 
!
      class(molecular_system) :: molecule 
!           
      character(len=100) :: temp_name
!
      integer :: i
!
      call molecule%write_libint_files()
!
      call initialize_atoms(molecule%name_)
!
      call reset_basis_c()
!
      do i = 1, molecule%n_basis_sets
!
         write(temp_name, '(a, a1, i4.4)') trim(molecule%name_), '_', i
!
         call initialize_basis(molecule%basis_sets(i), temp_name, molecule%cartesian_gaussians_int)
!
      enddo
!
   end subroutine initialize_libint_atoms_and_bases_molecular_system
!
!
   subroutine read_system_molecular_system(molecule)
!!
!!    Read system
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the multiplicity charge and name of the molecule.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      call input%get_required_keyword_in_section('name', 'system', molecule%name_)
!
      call input%get_keyword_in_section('charge', 'system', molecule%charge)
      call input%get_keyword_in_section('multiplicity', 'system', molecule%multiplicity)
!
   end subroutine read_system_molecular_system
!
!
   subroutine set_geometry_molecular_system(molecule, R_qk)
!!
!!    Set geometry 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Sets the molecular geometry. R_qk is the qth coordinate (q = 1(x), 2(y), 3(z))
!!    of the kth atom (k = 1, 2, 3, ..., n_atoms).
!!
      implicit none 
!
      class(molecular_system), intent(inout) :: molecule 
!
      real(dp), dimension(3, molecule%n_atoms), intent(in) :: R_qk 
!
      integer :: k
!
!     Update geometry on eT side
!
      do k = 1, molecule%n_atoms
!
         molecule%atoms(k)%x = (R_qk(1,k))/angstrom_to_bohr
         molecule%atoms(k)%y = (R_qk(2,k))/angstrom_to_bohr
         molecule%atoms(k)%z = (R_qk(3,k))/angstrom_to_bohr
!
      enddo
!
!     Write Libint files, then update the atoms and bases 
!     on the Libint side.
!
      call molecule%initialize_libint_atoms_and_bases()
!
!     Finally, reinitialize the Libint integral engines 
!
      call molecule%initialize_libint_integral_engines()
!
   end subroutine set_geometry_molecular_system
!
!
   function get_geometry_molecular_system(molecule) result(R_qk)
!!
!!    Get geometry 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Gets the molecular geometry. R_qk is the qth coordinate (q = 1(x), 2(y), 3(z))
!!    of the kth atom (k = 1, 2, 3, ..., n_atoms).
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp), dimension(3, molecule%n_atoms) :: R_qk 
!
      integer :: k
!
      do k = 1, molecule%n_atoms
!
         R_qk(1,k) = (molecule%atoms(k)%x)*angstrom_to_bohr
         R_qk(2,k) = (molecule%atoms(k)%y)*angstrom_to_bohr
         R_qk(3,k) = (molecule%atoms(k)%z)*angstrom_to_bohr
!
      enddo
!
   end function get_geometry_molecular_system
!
!
   subroutine initialize_libint_integral_engines_molecular_system()
!!
!!    Initialize Libint integral engines 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      call initialize_coulomb_c()
      call initialize_kinetic_c()
      call initialize_nuclear_c()
      call initialize_overlap_c()
      call initialize_dipole_c()
      call initialize_quadrupole_c()
!
   end subroutine initialize_libint_integral_engines_molecular_system
!
!
   subroutine read_geometry_molecular_system(molecule)
!!
!!    Read geometry
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    Modified by Åsmund H. Tveten, Oct 2019. Ensuring that
!!    coordinates are read in units of Ångström.
!!
!!    Read atoms and their coordinates.
!!    The coordinates are assumed to be in units of Ångström by
!!    default; if units of Bohr are specified, they will be
!!    converted to Ångström.
!!
!!    In eT, atoms are ordered after basis set, with the first basis
!!    set specified on input first.
!!
!!    If there are active atoms, the active atoms are first in eT,
!!    irrespective of basis set.
!!
!!    This routine handles reading of the geometry and from eT.inp 
!!    and orders the atoms according to basis set.
!!    If active atoms are specified, they will be read and 
!!    reordered in read_active_atoms_molecular_system
!!
      implicit none
!
      class(molecular_system) :: molecule
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable :: positions
!
      character(len=2), dimension(:), allocatable :: symbols
      character(len=100), dimension(:), allocatable :: basis_sets
!
      logical, dimension(:), allocatable :: selected_atom
!
      integer :: atom, current_atom
!
      character(len=100) :: current_basis
!
      logical :: units_angstrom
!
      real(dp) :: angstrom_conversion
!
      units_angstrom = .true.
      angstrom_conversion = one
!
      molecule%n_atoms = input%get_n_atoms()
!
      call mem%alloc(positions, molecule%n_atoms, 3)
!
      allocate(symbols(molecule%n_atoms))
      allocate(basis_sets(molecule%n_atoms))
!
      call input%get_geometry(molecule%n_atoms, symbols, positions, basis_sets, units_angstrom)
!
      call molecule%initialize_atoms()
!
!     Set coordinate conversion factor if geometry is given in units of Bohr
!
      if(.not. units_angstrom) then
!
         angstrom_conversion = bohr_to_angstrom
         molecule%coordinate_units = 'bohr'
!
      endif
!
!     1. Place the first atom in atoms
!
      current_atom = 1
!
      molecule%atoms(current_atom)%symbol       = symbols(current_atom)
      molecule%atoms(current_atom)%basis        = basis_sets(current_atom)
      molecule%atoms(current_atom)%x            = positions(current_atom,1)*angstrom_conversion
      molecule%atoms(current_atom)%y            = positions(current_atom,2)*angstrom_conversion
      molecule%atoms(current_atom)%z            = positions(current_atom,3)*angstrom_conversion
!
      molecule%atoms(current_atom)%input_number = current_atom
!
      current_basis = molecule%atoms(current_atom)%basis
!
      allocate(selected_atom(molecule%n_atoms))
      selected_atom = .false.
!
      selected_atom(current_atom) = .true.
!
!     2. Place rest of atoms such that they come in order of basis set
!
      do while (current_atom < molecule%n_atoms)
!
         do atom = 1, molecule%n_atoms
!
            if (trim(basis_sets(atom)) == trim(current_basis) .and. .not. selected_atom(atom)) then
!
               current_atom = current_atom + 1
!
               molecule%atoms(current_atom)%symbol       = symbols(atom)
               molecule%atoms(current_atom)%basis        = basis_sets(atom)
               molecule%atoms(current_atom)%x            = positions(atom,1)*angstrom_conversion
               molecule%atoms(current_atom)%y            = positions(atom,2)*angstrom_conversion
               molecule%atoms(current_atom)%z            = positions(atom,3)*angstrom_conversion
!
               molecule%atoms(current_atom)%input_number = atom
!
               selected_atom(atom) = .true.      
!
            endif
!
         enddo
!
         do atom = 1, molecule%n_atoms
!
            if (.not. selected_atom(atom)) then
!
               current_atom = current_atom + 1
!
               molecule%atoms(current_atom)%symbol       = symbols(atom)
               molecule%atoms(current_atom)%basis        = basis_sets(atom)
               molecule%atoms(current_atom)%x            = positions(atom,1)*angstrom_conversion
               molecule%atoms(current_atom)%y            = positions(atom,2)*angstrom_conversion
               molecule%atoms(current_atom)%z            = positions(atom,3)*angstrom_conversion
!
               molecule%atoms(current_atom)%input_number = atom
!
               current_basis = molecule%atoms(current_atom)%basis
!
               selected_atom(atom) = .true.
!
               exit             
!
            endif
!
         enddo         
!
      enddo
!
      call mem%dealloc(positions, molecule%n_atoms, 3)
!
      deallocate(symbols)
      deallocate(basis_sets)
      deallocate(selected_atom)
!
   end subroutine read_geometry_molecular_system
!
!
   subroutine rename_core_valence_dunning_sets_molecular_system(molecule)
!!
!!    Rename core valence Dunning sets 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Renames the basis sets (if any) of the core-valence type:
!!
!!       aug-cc-pCVXZ -> _aug-cc-pCVXZ
!!
!!    Note to developers: This is necessary because the Libint files must be named 
!!    with an "_" prefix to avoid it looking for an "augmentation" file. 
!!
      implicit none 
!
      class(molecular_system) :: molecule 
!
      integer :: I, k
!
      integer, parameter :: n_renamings = 4
!
      character(len=20), dimension(n_renamings), parameter :: original_names = ['aug-cc-pcvdz', &
                                                                                'aug-cc-pcvtz', &
                                                                                'aug-cc-pcvqz', &
                                                                                'aug-cc-pcv5z']
!
      character(len=20), dimension(n_renamings), parameter :: new_names = ['_aug-cc-pcvdz', &
                                                                           '_aug-cc-pcvtz', &
                                                                           '_aug-cc-pcvqz', &
                                                                           '_aug-cc-pcv5z']
!
      do I = 1, molecule%n_atoms 
!
         do k = 1, n_renamings
!
            if (trim(molecule%atoms(I)%basis) == trim(original_names(k))) then 
!
               molecule%atoms(I)%basis = trim(new_names(k))
!
            endif
!
         enddo
!
      enddo
!
   end subroutine rename_core_valence_dunning_sets_molecular_system
!
!
   subroutine write_libint_files_molecular_system(molecule)
!!
!!    Write LibInt files
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Writes to disk an xyz file, which is used by the libint package
!!    to calculate integrals.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: atom = 0
!
      character(len=100) temp_name
      character(len=100) current_basis
!
      type(output_file) :: mol_file
      type(output_file) :: basis_file
!
      integer :: basis_set_counter, atom_offset, current_basis_nbr, i
!
      integer, dimension(:), allocatable :: n_atoms_in_basis
!
!     Write atom file
!
      mol_file = output_file(trim(molecule%name_) // '.xyz')
      call mol_file%open_('rewind')
!
      call mol_file%printf('m', '(i5)', ints=[molecule%n_atoms], fs='(a/)')
!
      do atom = 1, molecule%n_atoms
!
         call mol_file%printf('m', '(a2)   (f21.16)   (f21.16)   (f21.16)', &
                              chars = [molecule%atoms(atom)%symbol],        &
                              reals = [molecule%atoms(atom)%x,              &
                                       molecule%atoms(atom)%y,              &
                                       molecule%atoms(atom)%z],             &
                              fs ='(a)', ll=80)
!
      enddo
!
      call mol_file%close_
!
!     Count number of basis sets
!
      current_basis = molecule%atoms(1)%basis
      molecule%n_basis_sets = 1
!
      do atom = 2, molecule%n_atoms
!
         if (molecule%atoms(atom)%basis .ne. current_basis) then
!
            molecule%n_basis_sets = molecule%n_basis_sets + 1
            current_basis = molecule%atoms(atom)%basis
!
         endif
!
      enddo
!
      call molecule%initialize_basis_sets()
!
      call mem%alloc(n_atoms_in_basis, molecule%n_basis_sets)
!
      n_atoms_in_basis = 1
      basis_set_counter = 1
      molecule%basis_sets(basis_set_counter) = molecule%atoms(1)%basis
!
!     Count number of atoms in each basis
!
      do atom = 2, molecule%n_atoms
!
         if (molecule%atoms(atom)%basis .ne. molecule%basis_sets(basis_set_counter)) then
!
            basis_set_counter = basis_set_counter + 1
            molecule%basis_sets(basis_set_counter) = molecule%atoms(atom)%basis
!
         else
!
            n_atoms_in_basis(basis_set_counter) = n_atoms_in_basis(basis_set_counter) + 1
!
         endif
!
      enddo
!
      atom_offset = 0
!
      do current_basis_nbr = 1, molecule%n_basis_sets
!
         write(temp_name, '(a, a1, i4.4, a4)')trim(molecule%name_), '_', current_basis_nbr,  '.xyz'
!
         basis_file = output_file(trim(temp_name))
         call basis_file%open_('rewind')
!
         call basis_file%printf('m', '(i5)', ints=[n_atoms_in_basis(current_basis_nbr)], &
                                fs='(a/)')
!
         do i = 1, n_atoms_in_basis(current_basis_nbr)
!
            atom = i + atom_offset
!
            call basis_file%printf('m', '(a2)   (f21.16)   (f21.16)   (f21.16)', &
                                   chars = [molecule%atoms(atom)%symbol],        &
                                   reals = [molecule%atoms(atom)%x,              &
                                            molecule%atoms(atom)%y,              &
                                            molecule%atoms(atom)%z],             &
                                   fs ='(a)', ll=80)
!
         enddo
!
         call basis_file%close_
!
         atom_offset = atom_offset + n_atoms_in_basis(current_basis_nbr)
!
      enddo
!
      call mem%dealloc(n_atoms_in_basis, molecule%n_basis_sets)
!
   end subroutine write_libint_files_molecular_system
!
!
   subroutine delete_libint_files_molecular_system(molecule)
!!
!!    Delete LibInt Files
!!    Written by Tor S. Haugland, 2019
!!
!!    Deletes .xyz files used by LibInt to generate integrals.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      type(sequential_file) :: xyz_file
      integer :: current_basis_nbr
      character(len=100) temp_name
!
      xyz_file = sequential_file(trim(molecule%name_) // '.xyz', 'formatted')
!
      call xyz_file%delete_()
!
      do current_basis_nbr = 1, molecule%n_basis_sets
!
         write(temp_name, '(a, a1, i4.4, a4)') trim(molecule%name_), '_', current_basis_nbr, '.xyz'
!
         xyz_file = sequential_file(trim(temp_name), 'formatted')
!
         call xyz_file%delete_()
!
      enddo
!
   end subroutine delete_libint_files_molecular_system
!
!
   subroutine read_active_atoms_molecular_system(molecule)
!!
!!    Read_active_atoms
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    In eT, atoms are ordered after basis set, with the first basis
!!    set specified on input first.
!!
!!    If there are active atoms, the active atoms are first in eT,
!!    irrespective of basis set.
!!
!!    This routine handles reading of the active atoms 
!!    and their reordering to the begining of the atoms
!!    array.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: method, offset, atom, lower_level_method, space
!
      integer :: n_total, central_atom, i, j
!
      integer, dimension(:), allocatable :: n_active
      real(dp), dimension(:), allocatable :: radius
!
      integer, dimension(:), allocatable :: active_atoms, atoms_tmp, counter
!
      character(len=4), dimension(:), allocatable :: methods
!
      character(len=200), dimension(:), allocatable :: basis
!
      character(len=200) :: selection_type, inactive_basis
!
      real(dp) :: x, y, z
!
      type(atomic), dimension(:), allocatable :: atoms_copy
!
      logical :: found, inactive_basis_found
!
      if (.not. input%requested_section('active atoms')) return
!
!     Possible methods for active atoms ( ordered after level !! )
!
      methods = (/   'ccsd', &
                     'cc2 ', &
                     'ccs ', &
                     'hf  '/)  
!
      call mem%alloc(n_active, size(methods))
!
      n_active = 0
      n_total  = 0
!
!     Based on 'selection type' different approaches are used to construct the 
!     active atoms array
!
!     Determine selection type
! 
      call input%get_keyword_in_section('selection type', 'active atoms', selection_type)
!
      if (trim(selection_type) == 'range' .or. trim(selection_type) == 'list') then
!
!        Set the number of actives for each method
!
         do method = 1, size(methods)
!
            n_active(method) = input%get_n_elements_for_keyword_in_section(trim(methods(method)), 'active atoms')
            n_total = n_total + n_active(method)
!
         enddo
!
         call mem%alloc(active_atoms, n_total)
!
         offset = 0
!
         do method = 1, size(methods)
!
            if (n_active(method) == 0) cycle
!
            call mem%alloc(atoms_tmp, n_active(method))
!
            call input%get_array_for_keyword_in_section(trim(methods(method)), 'active atoms', n_active(method), atoms_tmp)
!
!           Translate atom list to current eT ordering (based on basis set)
!
            call molecule%translate_from_input_order_to_eT_order(n_active(method), atoms_tmp, &
                                                      active_atoms(offset + 1: offset + n_active(method)))
!
            call mem%dealloc(atoms_tmp, n_active(method))
!
            offset = offset + n_active(method)
!
         enddo
!
!        Sanity check for duplicates
!
         do i = 1, n_total
            do j = i + 1, n_total
!
             if (active_atoms(i) == active_atoms(j)) &
                call output%error_msg('duplicate in active atoms, check input for overlapping lists or ranges.')
!
            enddo
         enddo
!
      elseif  (trim(selection_type) == 'central atom') then
!
!        Set central atom (using current eT ordering)
!
         call input%get_keyword_in_section('central atom', 'active atoms', central_atom)
!
         do atom = 1, molecule%n_atoms
!
            if (molecule%atoms(atom)%input_number == central_atom) then
!
               central_atom = atom
               exit
!
            endif
!
         enddo
!
         call mem%alloc(radius, size(methods))
!
!        Set radius
!
         radius = zero
!
         do method = 1, size(methods)
!
            call input%get_keyword_in_section(trim(methods(method)), 'active atoms', radius(method))
!
         enddo
!
!        Sanity check on radii ( r_cc3 < r_ccsd < r_cc2 < r_ccs < r_hf)
!
         do method = 1, size(methods) - 1
!
            if (radius(method) == zero) cycle
!
            do lower_level_method = method + 1, size(methods)
!
               if (radius(lower_level_method) == 0) cycle
!
               if (radius(method) .gt. radius(lower_level_method)) &
                  call output%error_msg('Active atoms by radius, but ' // trim(methods(method)) //&
                     ' radius is greater than '// trim(methods(lower_level_method)) //' radius.')
!
            enddo
!
         enddo
!
!        Get number of active atoms
!
         do atom = 1, molecule%n_atoms 
!
            x = (molecule%atoms(central_atom)%x - molecule%atoms(atom)%x)
            y = (molecule%atoms(central_atom)%y - molecule%atoms(atom)%y)
            z = (molecule%atoms(central_atom)%z - molecule%atoms(atom)%z)
!
            do method = 1, size(methods)
!
               if (radius(method) == zero) cycle
!
               if (sqrt(x**2 + y**2 + z**2) .lt. radius(method)) then 
!
                  n_active(method) = n_active(method) + 1
                  exit 
!
               endif
!
            enddo
!
         enddo
!
         do method = 1, size(methods)
!
            n_total = n_total + n_active(method)
!
         enddo
!
         call mem%alloc(active_atoms, n_total)
!
         call mem%alloc(counter, size(methods))
         counter = 0
!
         do atom = 1, molecule%n_atoms 
!
            x = (molecule%atoms(central_atom)%x - molecule%atoms(atom)%x)
            y = (molecule%atoms(central_atom)%y - molecule%atoms(atom)%y)
            z = (molecule%atoms(central_atom)%z - molecule%atoms(atom)%z)
!
            offset = 0
!
            do method = 1, size(methods)
!
               if (radius(method) == zero) cycle
!
               if (sqrt(x**2 + y**2 + z**2) .lt. radius(method)) then 
!
                  counter(method) = counter(method) + 1
                  active_atoms(offset + counter(method)) = atom

                  exit 
!
               endif
!
               offset = offset + n_active(method)
!
            enddo
!
         enddo
!
         call mem%dealloc(radius, size(methods))
         call mem%dealloc(counter, size(methods))
!
      endif
!
      allocate(basis(size(methods)))
!
      do method = 1, size(methods)
!
         basis(method) = repeat(' ', 200)
         call input%get_keyword_in_section(trim(methods(method))// ' basis', 'active atoms', basis(method))
!
      enddo
!
      molecule%n_active_atom_spaces = 0
!
      do method = 1, size(methods)
!
         if (n_active(method) > 0) molecule%n_active_atom_spaces = molecule%n_active_atom_spaces + 1
!
      enddo
!
      if (molecule%n_active_atom_spaces == 0) call output%error_msg('Requested active atoms, but no active spaces found.')
!
      allocate(molecule%active_atom_spaces(molecule%n_active_atom_spaces))
!
!     Reorder atoms
!
      allocate(atoms_copy(molecule%n_atoms))
      atoms_copy = molecule%atoms
!
      atom  = 0
      space = 0
      offset = 0
!
      do method = 1, size(methods)
         do i = 1, n_active(method)
!
            atom = atom + 1
            molecule%atoms(atom) = atoms_copy(active_atoms(atom))
!
            if (adjustl(basis(method)) /= '') molecule%atoms(atom)%basis = trim(adjustl(basis(method)))
!
         enddo        
!
         if (n_active(method) > 0) then
! 
            space = space + 1
!
            molecule%active_atom_spaces(space)%level = trim(methods(method))
!
            molecule%active_atom_spaces(space)%first_atom = offset + 1
            molecule%active_atom_spaces(space)%last_atom = offset + n_active(method)
!
            offset = offset + n_active(method)
!
         endif
!
      enddo
!
      inactive_basis_found = input%requested_keyword_in_section('inactive basis', 'active atoms')
!
      if (inactive_basis_found) then
!
         call input%get_keyword_in_section('inactive basis', 'active atoms', inactive_basis)
!
      endif 
!
      atom = n_total + 1
!
      do i = 1, molecule%n_atoms
!
         found = .false.
!
         do j = 1, n_total
!
            if (i == active_atoms(j)) then
!
               found = .true.
               exit
!
            endif
!
         enddo
!
         if (.not. found) then
!
            molecule%atoms(atom) = atoms_copy(i)
!
            if (inactive_basis_found) molecule%atoms(atom)%basis = trim(adjustl(inactive_basis))
!
            atom = atom + 1
!
         endif 
!
      enddo
!
      deallocate(atoms_copy)
!
      call mem%dealloc(n_active, size(methods))
      call mem%dealloc(active_atoms, n_total)
      deallocate(basis)
!
   end subroutine read_active_atoms_molecular_system
!
!
   function get_nuclear_repulsion_1der_numerical_molecular_system(molecule, dx) result(h_nuc_qk)
!!
!!    Get nuclear repulsion 1der numerical 
!!    Written by Eirik F. Kjønstad, June 2019
!!
!!    Computes derivative numerically.
!!
      implicit none 
!
      class(molecular_system), intent(inout) :: molecule 
!
      real(dp), intent(in) :: dx
!
      real(dp), dimension(3, molecule%n_atoms) :: h_nuc_qk 
!
      real(dp) :: energy_x, energy_xdx
!
      real(dp), dimension(3, molecule%n_atoms) :: R_qk, R_qk_displaced
!
      integer :: k, q
!
      h_nuc_qk = zero
!
      energy_x = molecule%get_nuclear_repulsion()
! 
      R_qk = molecule%get_geometry()
!
      do k = 1, molecule%n_atoms
         do q = 1, 3
!
            R_qk_displaced = R_qk 
            R_qk_displaced(q,k) = R_qk_displaced(q,k) + dx 
!
            call molecule%set_geometry(R_qk_displaced)
!
            energy_xdx = molecule%get_nuclear_repulsion()
!
            h_nuc_qk(q,k) = (energy_xdx-energy_x)/dx
!
         enddo
      enddo
!
      call molecule%set_geometry(R_qk)
!
   end function get_nuclear_repulsion_1der_numerical_molecular_system
!
!
   function get_nuclear_repulsion_1der_molecular_system(molecule) result(h_nuc_qk)
!!
!!    Get nuclear repulsion 1der 
!!    Written by Eirik F. Kjønstad, June 2019
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp), dimension(3, molecule%n_atoms) :: h_nuc_qk 
!
      real(dp) :: x_ij, y_ij, z_ij, r_ij_3
      integer :: i, j 
!
      h_nuc_qk = zero
!
      do i = 1, molecule%n_atoms 
         do j = 1, i - 1
!
            x_ij = (molecule%atoms(i)%x - molecule%atoms(j)%x)*angstrom_to_bohr
            y_ij = (molecule%atoms(i)%y - molecule%atoms(j)%y)*angstrom_to_bohr 
            z_ij = (molecule%atoms(i)%z - molecule%atoms(j)%z)*angstrom_to_bohr
!
            r_ij_3 = sqrt(x_ij**2 + y_ij**2 + z_ij**2)**3
!
            h_nuc_qk(1, i) = h_nuc_qk(1, i) - x_ij * molecule%atoms(i)%number_*molecule%atoms(j)%number_/r_ij_3
            h_nuc_qk(2, i) = h_nuc_qk(2, i) - y_ij * molecule%atoms(i)%number_*molecule%atoms(j)%number_/r_ij_3
            h_nuc_qk(3, i) = h_nuc_qk(3, i) - z_ij * molecule%atoms(i)%number_*molecule%atoms(j)%number_/r_ij_3
!
            h_nuc_qk(1, j) = h_nuc_qk(1, j) + x_ij * molecule%atoms(i)%number_*molecule%atoms(j)%number_/r_ij_3
            h_nuc_qk(2, j) = h_nuc_qk(2, j) + y_ij * molecule%atoms(i)%number_*molecule%atoms(j)%number_/r_ij_3
            h_nuc_qk(3, j) = h_nuc_qk(3, j) + z_ij * molecule%atoms(i)%number_*molecule%atoms(j)%number_/r_ij_3
!
         enddo
      enddo
!
   end function get_nuclear_repulsion_1der_molecular_system
!
!
   function get_total_nuclear_repulsion_molecular_system(molecule) result(nuclear_repulsion)
!!
!!    Get total nuclear repulsion 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Returns the total nuclear repulsion energy, accounting for the MM contribution 
!!    if QM/MM calculation.
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp) :: nuclear_repulsion
!
      if (molecule%mm_calculation .and. molecule%mm%forcefield .eq. 'non-polarizable') then 
!
         nuclear_repulsion = molecule%get_nuclear_repulsion() + molecule%get_nuclear_repulsion_mm()
!
      else
!
         nuclear_repulsion = molecule%get_nuclear_repulsion()
!
      endif 
!
   end function get_total_nuclear_repulsion_molecular_system
!
!
   function get_nuclear_repulsion_molecular_system(molecule)
!!
!!    Get nuclear repulsion
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates, and returns, the nuclear repulsion term for the molecule,
!!    in units of Hartree. Makes use of the Ångstrøm to Bohr conversion
!!    factor defined in the parameters module.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      real(dp) :: get_nuclear_repulsion_molecular_system
!
      integer :: i = 0, j = 0
!
      real(dp) :: x_ij, y_ij, z_ij, r_ij
!
      get_nuclear_repulsion_molecular_system = zero
!
      do i = 1, molecule%n_atoms
         do j = i + 1, molecule%n_atoms
!
            x_ij = molecule%atoms(i)%x - molecule%atoms(j)%x
            y_ij = molecule%atoms(i)%y - molecule%atoms(j)%y
            z_ij = molecule%atoms(i)%z - molecule%atoms(j)%z
!
            r_ij = sqrt(x_ij**2 + y_ij**2 + z_ij**2)
!
            r_ij = angstrom_to_bohr*r_ij
!
            if (abs(r_ij) .lt. 1.0D-7) then
!
               call output%error_msg('two atoms are placed on top of each other.')
!
            endif
!
            get_nuclear_repulsion_molecular_system = get_nuclear_repulsion_molecular_system &
                  + ((molecule%atoms(i)%number_)*(molecule%atoms(j)%number_))/r_ij
!
         enddo
      enddo
!
  end function get_nuclear_repulsion_molecular_system
!
!
   function get_n_electrons_molecular_system(molecule)
!!
!!    Get number of electrons
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates and returns the number of electrons.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: get_n_electrons_molecular_system
!
      integer :: i = 0
!
      get_n_electrons_molecular_system = 0
!
      do i = 1, molecule%n_atoms
!
         get_n_electrons_molecular_system = get_n_electrons_molecular_system + molecule%atoms(i)%number_
!
      enddo
!
      get_n_electrons_molecular_system = get_n_electrons_molecular_system - molecule%charge
!
   end function get_n_electrons_molecular_system
!
!
   integer function get_n_shells_molecular_system(molecule)
!!
!!    Get number of shells
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: I = 0
!
      get_n_shells_molecular_system = 0
!
      do I = 1, molecule%n_atoms
!
         get_n_shells_molecular_system = get_n_shells_molecular_system + molecule%atoms(I)%n_shells
!
      enddo
!
   end function get_n_shells_molecular_system
!
!
   integer function get_n_aos_molecular_system(molecule)
!!
!!    Get number of AOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: I = 0, J = 0
!
      get_n_aos_molecular_system = 0
!
      do I = 1, molecule%n_atoms
!
         do J = 1, molecule%atoms(I)%n_shells
!
            get_n_aos_molecular_system = get_n_aos_molecular_system &
                           + molecule%atoms(I)%shells(J)%length
!
         enddo
!
      enddo
!
   end function get_n_aos_molecular_system
!
!
   subroutine get_max_shell_size_molecular_system(molecule, max_shell_size)
!!
!!    Get maximum shell size 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Loops through the shells and determines the size of the largest shell. 
!!    This can be useful to preallocate a matrix able to hold shell pair 
!!    or shell quadruple vector without having to allocate inside the 
!!    loop over such shell pairs or quadruples or what have you.
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      integer, intent(inout) :: max_shell_size
!
      integer :: s
!
      max_shell_size = 0
      do s = 1, molecule%n_s 
!
         if (max_shell_size .lt. molecule%shell_limits(s)%length) then
! 
            max_shell_size = molecule%shell_limits(s)%length
!
         endif
!
      enddo
!
   end subroutine get_max_shell_size_molecular_system
!
!
   function get_shell_limits_molecular_system(molecule, A) result(the_interval)
!!
!!    Get shell limits
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer, intent(in) :: A
!
      integer :: I, J
!
      type(interval) :: the_interval
!
      do I = 1, molecule%n_atoms
         do J = 1, molecule%atoms(I)%n_shells
!
            if (A .eq. molecule%atoms(I)%shells(J)%number_) then
!
               the_interval = interval(molecule%atoms(I)%shells(J)%first, &
                                       molecule%atoms(I)%shells(J)%last)
!
            endif
!
         enddo
      enddo
!
       if (the_interval%length == 0) call output%error_msg('in get_shell_limits.')
!
   end function get_shell_limits_molecular_system
!
!
   function basis2shell_molecular_system(molecule, basis_function)
!!
!!    Basis2shell
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer, intent(in) :: basis_function
!
      integer :: basis2shell_molecular_system
!
      integer :: I, J
!
      basis2shell_molecular_system = 0
!
      do I = 1, molecule%n_atoms
         do J = 1, molecule%atoms(I)%n_shells
!
            if (molecule%atoms(I)%shells(J)%last  .ge. basis_function .and. &
                molecule%atoms(I)%shells(J)%first .le. basis_function) then
!
               basis2shell_molecular_system = &
               molecule%atoms(I)%shells(J)%number_
!
            endif
!
         enddo
      enddo
!
      if (basis2shell_molecular_system == 0) call output%error_msg('in basis2shell.')
!
   end function basis2shell_molecular_system
!
!
   function basis2atom_molecular_system(molecule, basis_function)
!!
!!    Basis2atom
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer, intent(in) :: basis_function
!
      integer :: basis2atom_molecular_system
!
      integer :: I, J
!
      basis2atom_molecular_system = 0
!
      do I = 1, molecule%n_atoms
         do J = 1, molecule%atoms(I)%n_shells
!
            if (molecule%atoms(I)%shells(J)%last  .ge. basis_function .and. &
                molecule%atoms(I)%shells(J)%first .le. basis_function) then
!
               basis2atom_molecular_system = I
!
            endif
!
         enddo
      enddo
!
      if (basis2atom_molecular_system == 0) call output%error_msg('in basis2atom.')
!
   end function basis2atom_molecular_system
!
!
   function shell_to_atom_molecular_system(molecule, shell)
!!
!!
!!
      implicit none
!
      class(molecular_system) :: molecule
      integer, intent(in) :: shell
!
      integer :: I, accumulated_shells
!
      integer :: shell_to_atom_molecular_system
!
      shell_to_atom_molecular_system = 0
!
      if (shell .le. 0) then
!
         call output%error_msg('shell number has illegal value 0.')
!
      elseif (shell .gt. molecule%get_n_shells()) then
!
         call output%error_msg('shell number exceeds total number of shells.')
!
      endif
!
      accumulated_shells = 0
!
      do I = 1, molecule%n_atoms
!
         accumulated_shells = accumulated_shells + molecule%atoms(I)%n_shells
!
         if (shell .le. accumulated_shells) then
!
            shell_to_atom_molecular_system = I
            return
!
         endif    
!
      enddo
!     
      if (shell_to_atom_molecular_system == 0) call output%error_msg('in shell_to_atom.')
!
   end function shell_to_atom_molecular_system
!
!
   subroutine initialize_atoms_molecular_system(molecule)
!!
!!    Initialize atoms
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (.not. allocated(molecule%atoms)) allocate(molecule%atoms(molecule%n_atoms))
!
   end subroutine initialize_atoms_molecular_system
!
!
   subroutine destruct_atoms_molecular_system(molecule)
!!
!!    destruct atoms
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (allocated(molecule%atoms)) deallocate(molecule%atoms)
!
   end subroutine destruct_atoms_molecular_system
!
!
   subroutine initialize_basis_sets_molecular_system(molecule)
!!
!!    Initialize basis sets
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (.not. allocated(molecule%basis_sets)) allocate(molecule%basis_sets(molecule%n_basis_sets))
!
   end subroutine initialize_basis_sets_molecular_system
!
!
   subroutine destruct_basis_sets_molecular_system(molecule)
!!
!!    Destruct basis sets
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (allocated(molecule%basis_sets)) deallocate(molecule%basis_sets)
!
   end subroutine destruct_basis_sets_molecular_system
!
!
   subroutine initialize_shell_limits_molecular_system(molecule)
!!
!!    Initialize shell limits 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: n_s 
!
      n_s = molecule%get_n_shells()
!
      if (.not. allocated(molecule%shell_limits)) allocate(molecule%shell_limits(n_s))
!
   end subroutine initialize_shell_limits_molecular_system
!
!
   subroutine initialize_shell2atom_molecular_system(molecule)
!!
!!    Initialize shell to atom 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (.not. allocated(molecule%shell2atom)) call mem%alloc(molecule%shell2atom, molecule%n_s)
!
   end subroutine initialize_shell2atom_molecular_system
!
!
   subroutine destruct_shell2atom_molecular_system(molecule)
!!
!!    Destruct shell to atom 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (allocated(molecule%shell2atom)) call mem%dealloc(molecule%shell2atom, molecule%n_s)
!
   end subroutine destruct_shell2atom_molecular_system
!
!
   subroutine destruct_shell_limits_molecular_system(molecule)
!!
!!    Destruct basis sets
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      if (allocated(molecule%shell_limits)) deallocate(molecule%shell_limits)
!
   end subroutine destruct_shell_limits_molecular_system
!
!
   subroutine print_geometry_molecular_system(molecule, units)
!!
!!    Print geometry 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints xyz of the molecular geometry in the unit specified.
!!
!!    units:   'angstrom' or 'bohr'
!!
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Printf and print_separator. Modified to take unit as input
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule  
      character(len=*),        intent(in) :: units
!
      real(dp)           :: scaling
      character(len=100) :: local_units
!
      integer            :: I, line_length, start_integer
      type(atomic)       :: atom
!
!     Choose unit
!
      if (trim(units) == 'angstrom') then
!
         scaling = one
         local_units = 'Å'
!
      elseif (trim(units) == 'bohr') then
!
         scaling = angstrom_to_bohr
         local_units = 'a.u.'
!
      else
!
         scaling = zero
         call output%error_msg("Could not find unit: "// units // ". Use angstrom/bohr")
!
      endif
!
!     Print header
!
      line_length = 73
!
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(/t5,a)')
      call output%printf('m', 'Geometry((a0))', fs='(t32,a)', chars=[local_units])
      call output%print_separator(pl='minimal', symbol='=', n=line_length, fs='(t5,a)')
      call output%printf('m', 'Atom          X                Y                &
                         &Z                Basis', ll=100, fs='(t6,a)')
      call output%print_separator(pl='minimal', symbol='=', n=line_length, fs='(t5,a)')
!
!     Print geometry for every atom
!
      do I = 1, molecule%n_atoms
!
         atom = molecule%atoms(I)
!
         start_integer = 1
         if (atom%basis(1:1) == '_') start_integer = 2
!
         call output%printf('m', adjustl(atom%symbol) // &
                            ' (f17.12)(f17.12)(f17.12)  (a14)', &
                            chars=[atom%basis(start_integer:)], &
                            reals=[atom%x * scaling, atom%y * scaling, &
                            atom%z * scaling], ll=100, fs='(t7,a)')
!
      enddo 
!
      call output%print_separator(pl='minimal', symbol='=', n=line_length, fs='(t5,a)')
!
   end subroutine print_geometry_molecular_system
!
!
   subroutine print_system_molecular_system(molecule)
!!
!!    Print geometry 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!    Modified by Tommaso Giovannini, March 2019
!!    Modified by Åsmund H. Tveten, Oct 2019. Switched to printf function.
!!
      implicit none 
!
      class(molecular_system) :: molecule  
!
      if (molecule%mm_calculation) then
!
         call output%printf('m', ':: Molecular system specifications (QM)', fs='(//t3,a)')
         call output%print_separator('m', 42, '=')
!
      else
!
         call output%printf('m', ':: Molecular system specifications', fs='(//t3,a)')
         call output%print_separator('m', 37, '=')
!
      endif
!
      call output%printf('m', 'Name:             (a0)', &
                         chars=[trim(molecule%name_)], fs='(/t6,a)')
      call output%printf('m', 'Charge:         (i3)', fs='(t6,a)',  ints=[molecule%charge]) 
      call output%printf('m', 'Multiplicity:   (i3)', &
                         ints=[molecule%multiplicity], fs='(t6,a)')
      call output%printf('m', 'Coordinate units: (a0)', &
                         chars=[trim(molecule%coordinate_units)], fs='(t6,a)')
!
      if (molecule%cartesian_gaussians_int.eq.1) &
         call output%printf('m', 'Using Cartesian gaussians.', fs='(/t6,a)')
!
      call output%printf('m', 'Pure basis functions:      (i5)', &
                         ints=[molecule%n_pure_basis], fs='(/t6,a)')
      call output%printf('m', 'Cartesian basis functions: (i5)', &
                         ints=[molecule%n_cart_basis], fs='(t6,a)')
      call output%printf('m', 'Primitive basis functions: (i5)', &
                         ints=[molecule%n_primitives_cart], fs='(t6,a)')
!
      call output%printf('m', 'Nuclear repulsion energy (a.u.):   (f25.12)', &
                         reals=[molecule%get_nuclear_repulsion()], fs='(/t6,a)')
      call output%printf('m', 'Bohr/angstrom value (CODATA 2010): (f25.12)', &
                         reals=[bohr_to_angstrom], fs='(t6,a)')
!
      call molecule%print_active_atoms()
!
      call molecule%print_geometry('angstrom')
      call molecule%print_geometry('bohr')
!
   end subroutine print_system_molecular_system
!
!
   subroutine translate_from_input_order_to_eT_order_molecular_system(molecule, n_elements,  &
                                                         array_input_ordering, array_eT_ordering)
!!
!!    Translate from input order to eT order
!!    Written by Sarai D. Folkestad, Mar 2019
!!
!!    In eT, atoms are ordered after basis set, with the first basis
!!    set specified on input first.
!!
!!    If there are active atoms, the active atoms are first in eT,
!!    irrespective of basis set.
!!
!!    This routine translates an array of input indices,
!!    which corresponds to the order of the atoms in the
!!    input file to an array of indices corresponding to atom 
!!    order in eT.
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: n_elements
!
      integer, dimension(n_elements), intent(in) :: array_input_ordering
      integer, dimension(n_elements), intent(out) :: array_eT_ordering
!  
!     Local variables
!
      integer :: i, j
!
      do i = 1, n_elements
         do j = 1, molecule%n_atoms
!
            if (array_input_ordering(i) == molecule%atoms(j)%input_number) then
!
               array_eT_ordering(i) = j
! 
            endif
!
         enddo
      enddo
!
   end subroutine translate_from_input_order_to_eT_order_molecular_system
!
!
   subroutine get_nuclear_dipole_molecular_system(molecule, mu_k)
!!
!!    Get nuclear dipole moment 
!!    Written by Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp), dimension(3), intent(out) :: mu_k 
!
      integer :: j 
!
      mu_k = zero 
!
      do j = 1, molecule%n_atoms 
!
         mu_k(1) = mu_k(1) + (molecule%atoms(j)%x)*(molecule%atoms(j)%number_)*(angstrom_to_bohr)
         mu_k(2) = mu_k(2) + (molecule%atoms(j)%y)*(molecule%atoms(j)%number_)*(angstrom_to_bohr)
         mu_k(3) = mu_k(3) + (molecule%atoms(j)%z)*(molecule%atoms(j)%number_)*(angstrom_to_bohr)
!
      enddo
!
   end subroutine get_nuclear_dipole_molecular_system
!
!
   subroutine get_nuclear_quadrupole_molecular_system(molecule, Q_k)
!!
!!    Get nuclear dipole moment 
!!    Written by Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp), dimension(6), intent(out) :: Q_k 
!
      integer :: j 
!
      Q_k = zero 
!
      do j = 1, molecule%n_atoms 
!
         Q_k(1) = Q_k(1) + (molecule%atoms(j)%x**2)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(2) = Q_k(2) + (molecule%atoms(j)%x*molecule%atoms(j)%y)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(3) = Q_k(3) + (molecule%atoms(j)%x*molecule%atoms(j)%z)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(4) = Q_k(4) + (molecule%atoms(j)%y**2)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(5) = Q_k(5) + (molecule%atoms(j)%y*molecule%atoms(j)%z)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(6) = Q_k(6) + (molecule%atoms(j)%z**2)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
!
      enddo
!
   end subroutine get_nuclear_quadrupole_molecular_system
!
!
   function read_n_active_atoms_for_method_molecular_system(molecule, method) result(number_)
!!
!!    Read number of active for method
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      character(len=*) :: method
!
      integer :: number_
!
!     Local variables
!
      character(len=200) :: selection_type
!
      integer :: central_atom, i
!
      real(dp) :: radius, x, y, z
!
      number_ = 0
!
      call input%get_keyword_in_section('selection type', 'active atoms', selection_type)
!
      if (trim(selection_type) == 'range' .or. trim(selection_type) == 'list') then
!
         number_ = input%get_n_elements_for_keyword_in_section(trim(method), 'active atoms')
!
      elseif (trim(selection_type) == 'central atom') then
!
         call input%get_keyword_in_section('central atom', 'active atoms', central_atom)
         call input%get_keyword_in_section(trim(method), 'active atoms', radius)
!
         do i = 1, molecule%n_atoms
!
            if (molecule%atoms(i)%input_number == central_atom) then
!
               central_atom = i
               exit
!
            endif
!
         enddo
!
!        Set active atoms
!
         number_ = 0
!
         do i = 1, molecule%n_atoms 
!
            x = (molecule%atoms(central_atom)%x - molecule%atoms(i)%x)
            y = (molecule%atoms(central_atom)%y - molecule%atoms(i)%y)
            z = (molecule%atoms(central_atom)%z - molecule%atoms(i)%z)
!
            if (sqrt(x**2 + y**2 + z**2) .lt. radius) number_ = number_ + 1
!
         enddo
!
      else
!
         call output%error_msg('Could not recognize active atom selection type.')
!
      endif
!
   end function read_n_active_atoms_for_method_molecular_system
!
!
   subroutine first_and_last_ao_active_space_molecular_system(molecule, level, first, last)
!!
!!    First and last ao in active space
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      character(len=*), intent(in) :: level
!
      integer, intent(out) :: first, last
!
      call molecule%first_ao_active_space(level, first)
      call molecule%last_ao_active_space(level, last)
!
   end subroutine first_and_last_ao_active_space_molecular_system
!
!
   subroutine first_ao_active_space_molecular_system(molecule, level, first)
!!
!!    First ao in active space
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      character(len=*), intent(in) :: level
!
      integer, intent(out) :: first
!
      integer :: i, first_atom 
!  
      first_atom = 0
!
      do i = 1, molecule%n_active_atom_spaces
!
         if (trim(molecule%active_atom_spaces(i)%level) == trim(level)) then
!
            first_atom = molecule%active_atom_spaces(i)%first_atom
            exit
!
         endif
!
      enddo
!
      if (first_atom == 0) &
         call output%warning_msg('Could not find the requested active space in molecular system.')
!
      first = 1
!
      do i = 1, first_atom - 1
!
         first = first + molecule%atoms(i)%n_ao
!
      enddo
!
   end subroutine first_ao_active_space_molecular_system
!
!
   subroutine last_ao_active_space_molecular_system(molecule, level, last)
!!
!!    Last ao in active space
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      character(len=*), intent(in) :: level
!
      integer, intent(out) :: last
!
      integer :: i, last_atom 
!  
      last_atom = 0
!
      do i = 1, molecule%n_active_atom_spaces
!
         if (trim(molecule%active_atom_spaces(i)%level) == trim(level)) then
!
            last_atom = molecule%active_atom_spaces(i)%last_atom
            exit
!
         endif
!
      enddo
!
      if (last_atom == 0) &
         call output%warning_msg('Could not find the requested active space in molecular system.')
!
      last = 0
!
      do i = 1, last_atom
!
         last = last + molecule%atoms(i)%n_ao
!
      enddo
!
   end subroutine last_ao_active_space_molecular_system
!
!
  subroutine set_basis_info_molecular_system(molecule)
!!
!!    Set basis info
!!    Written by Sarai D. Folkestad, Dec 2018
!!
!!    Opens and reads the basis set file
!!    for each atom in the molecule.
!!
!!    The information for each shell 
!!    (n_primitives, exponents, and coefficients)
!!    is then set by the set_shell_basis_info routine
!!
      implicit none
!
      class(molecular_system) :: molecule
      integer :: atom_index, shell
      character(len=100) :: basis
      character(len=100) :: libint_path
!
      type(sequential_file), allocatable :: basis_set_file
!
!     Loop over atoms 
!
      do atom_index = 1, molecule%n_atoms
!
         shell = 0
!
!        Deal with non-augmented part first
!
         if (molecule%atoms(atom_index)%basis(1:3) .ne. 'aug') then
            write(basis,'(a)') molecule%atoms(atom_index)%basis
         else
            write(basis,'(a)') molecule%atoms(atom_index)%basis(5:100)
         endif
!
         call convert_to_lowercase(basis)
         call get_environment_variable("LIBINT_DATA_PATH",libint_path)
         basis_set_file = sequential_file(trim(libint_path) // '/' // trim(basis) // '.g94', 'formatted')
!
         call basis_set_file%open_('read', 'rewind')  
!
         call molecule%set_shell_basis_info(basis_set_file, atom_index, shell)
!
         call basis_set_file%close_
!
!        Get the augmented part
         if (molecule%atoms(atom_index)%basis(1:3) .eq. 'aug') then
!
            basis = 'augmentation-'//trim(basis)
            basis_set_file = sequential_file(trim(libint_path) // '/' // trim(basis) // '.g94', 'formatted')
!
            call basis_set_file%open_('read', 'rewind')  
!
            call molecule%set_shell_basis_info(basis_set_file, atom_index, shell)
!
            call basis_set_file%close_
!
         endif
!
      enddo 
!
   end subroutine set_basis_info_molecular_system
!
!
   subroutine set_shell_basis_info_molecular_system(molecule, basis_set_file, atom_index, shell)
!!
!!    Set shell basis info
!!    Written by Sarai D. Folkestad, Dec 2018
!!    Modified by Marco Scavino, 2019
!!
!!    Sets the number of primitives, and the 
!!    coefficient and exponents of the primitives in 
!!    the basis_detail variable for each shell
!!
!!    Modified by Marco Scavino, 2019
!!    Extended to read the SP primitives from
!!    gaussian format
!!
      implicit none
!
      class(molecular_system), intent(inout) :: molecule
!
      integer, intent(inout)            :: shell
      integer, intent(in)               :: atom_index
!
      type(sequential_file), intent(in)                 :: basis_set_file
!
      character(len=200)   :: line
!
      integer         :: n_primitive, primitive
!
      character(len=2)     :: ang_mom
!
      logical              :: elm_found
!
      real(dp)             :: coefficient, coefficient_2, exponent_
!
!     Find position of element in file
!     Look for:  ****
!     Symbol
!
      elm_found = .false.
!
      call basis_set_file%read_(line,'(a200)')
!
      do while (.not. elm_found)
!
         if (trim(line) == '****') then
!
            call basis_set_file%read_(line,'(a200)')
!
            if ((line(1:2)) == molecule%atoms(atom_index)%symbol) then
!
               elm_found = .true.
!
            endif
!
         endif
!
         if (.not. elm_found) then
!
            call basis_set_file%read_(line,'(a200)')
!
         endif
!
      enddo
!
!     Read angular momentum and number of primitives
!
      call basis_set_file%read_(line,'(a200)')
!
      do while (trim(line) .ne. '****') ! Loop over AO
!
         shell = shell + 1
!
         read(line, *) ang_mom, n_primitive
!
!        Sanity check -> does shell exceed number of shells on atom
!
         if (shell .gt. molecule%atoms(atom_index)%n_shells) &
            call output%error_msg('Mismatch in number of shells in set_shell_basis_info')
!
!        Sanity check -> does shell have the correct angular momentum

         if (angular_momentum_from_symbol(ang_mom) .ne. molecule%atoms(atom_index)%shells(shell)%l) &
            call output%error_msg('Mismatch in angular momentum in set_shell_basis_info')
!
!        Set number of primitives and initialize exponents and coefficient array
!
         call molecule%atoms(atom_index)%shells(shell)%set_n_primitives(n_primitive)
         call molecule%atoms(atom_index)%shells(shell)%initialize_exponents()
         call molecule%atoms(atom_index)%shells(shell)%initialize_coefficients()

!
!        Loop over primitives and set coefficient and exponent
!
         if(ang_mom == "SP") then
!
!           In case of "SP" shell, split S and P coefficients
!
            call molecule%atoms(atom_index)%shells(shell+1)%set_n_primitives(n_primitive)
            call molecule%atoms(atom_index)%shells(shell+1)%initialize_exponents()
            call molecule%atoms(atom_index)%shells(shell+1)%initialize_coefficients()
!
            do primitive = 1, n_primitive
!
!              Read coefficient and exponent
!
               call basis_set_file%read_(line,'(a200)')
               read(line, *) exponent_, coefficient, coefficient_2
!
               call molecule%atoms(atom_index)%shells(shell)%set_exponent_i(primitive, exponent_)
               call molecule%atoms(atom_index)%shells(shell+1)%set_exponent_i(primitive, exponent_)
               call molecule%atoms(atom_index)%shells(shell)%set_coefficient_i(primitive, coefficient)
               call molecule%atoms(atom_index)%shells(shell+1)%set_coefficient_i(primitive, coefficient_2)
!
            enddo
!
            shell = shell + 1
!
         else
!
            do primitive = 1, n_primitive
!
!              Read coefficient and exponent
!
               call basis_set_file%read_(line,'(a200)')
               read(line, *) exponent_, coefficient
!
               call molecule%atoms(atom_index)%shells(shell)%set_exponent_i(primitive, exponent_)
               call molecule%atoms(atom_index)%shells(shell)%set_coefficient_i(primitive, coefficient)

!
            enddo
!
         end if
!
         call basis_set_file%read_(line,'(a200)')
!
      enddo
!
   end subroutine set_shell_basis_info_molecular_system
!
!
   function angular_momentum_from_symbol(angular_momentum_symbol)
!!
!!    Angular momentum from symbol
!!    Written by Sarai D. Folkestad, Dec 2018
!!
      implicit none
!
      character(len=1), intent(in) :: angular_momentum_symbol
!
      integer :: angular_momentum_from_symbol
!
      angular_momentum_from_symbol = 0
!
      if (angular_momentum_symbol == 'S') then
!
         angular_momentum_from_symbol = 0
!
      elseif (angular_momentum_symbol == 'P') then
!
         angular_momentum_from_symbol = 1
!
      elseif (angular_momentum_symbol == 'D') then
!
         angular_momentum_from_symbol = 2
!
      elseif (angular_momentum_symbol == 'F') then
!
         angular_momentum_from_symbol = 3
!
      elseif (angular_momentum_symbol == 'G') then
!
         angular_momentum_from_symbol = 4
!
      elseif (angular_momentum_symbol == 'H') then
!
         angular_momentum_from_symbol = 5
!
      elseif (angular_momentum_symbol == 'I') then
!
         angular_momentum_from_symbol = 6
!
      else
!
         call output%error_msg('no support for visualization of orbitals with l > 6')
!
      endif
!
   end function angular_momentum_from_symbol
!
!
   subroutine check_convert_pure_to_cartesian_basis_molecular_system(molecule)
!!   
!!    Check convert pure to cartesian
!!    Written by Tommaso Giovannini, 2019
!!
!!    Check if it is necessary to convert from pure to cartesian
!!
      implicit none
!
      class(molecular_system) :: molecule
      integer :: i,j, angmom, n_primitives, n_functions
!
      do i = 1, molecule%n_atoms
!  
         do j = 1, molecule%atoms(i)%n_shells
!  
            angmom                  = molecule%atoms(i)%shells(j)%l
            n_primitives            = molecule%atoms(i)%shells(j)%get_n_primitives()
            n_functions             = molecule%atoms(i)%shells(j)%length
            molecule%n_pure_basis   = molecule%n_pure_basis + n_functions 
!  
             if(angmom .ge. 1) then
!  
                molecule%atoms(i)%shells(j)%size_cart = n_functions
!  
                if(angmom .ge. 2) then
!  
                   if(angmom .eq. 2 .and. n_functions .eq. 6) continue
!  
                   if(angmom .eq. 2 .and. n_functions .eq. 5) then
!  
                      molecule%atoms(i)%shells(j)%size_cart = n_functions + 1
!  
                   endIf
!  
                   if(angmom .eq. 3 .and. n_functions .eq.10) continue
!  
                   if(angmom .eq. 3 .and. n_functions .eq. 7) then
!  
                      molecule%atoms(i)%shells(j)%size_cart = n_functions + 3
!  
                   endIf
!  
                   if(angmom .eq. 4 .and. n_functions .eq. 15) continue
!  
                   if(angmom .eq. 4 .and. n_functions .eq. 9) then
!  
                      molecule%atoms(i)%shells(j)%size_cart = n_functions + 6
!  
                   endIf
!  
                   if(angmom .gt. 4) then
!  
                      call output%error_msg('Cartesian G functions not yet implemented')
!  
                   endIf
!  
                endif
!  
            else 
!  
               molecule%atoms(i)%shells(j)%size_cart = n_functions 
!  
            endIf
!  
            molecule%n_cart_basis      = molecule%n_cart_basis + &
                                          molecule%atoms(i)%shells(j)%size_cart
!
            molecule%n_primitives_cart = molecule%n_primitives_cart + &
                                          molecule%atoms(i)%shells(j)%size_cart*n_primitives
!  
         enddo
!  
      enddo
!  
      if(molecule%n_cart_basis.eq.molecule%n_pure_basis) molecule%cartesian_basis = .True.
!  
!
   end subroutine check_convert_pure_to_cartesian_basis_molecular_system
!
!
   subroutine normalize_raw_primitives_molecular_system(molecule)
!!   
!!    Normalize raw primitives coefficients
!!    Written by Tommaso Giovannini
!!
      use math_utilities, only: double_factorial
!!
      implicit none
!
      class(molecular_system) :: molecule
      integer :: i, j, k, l, angmom, n_primitives
      real(kind=dp) :: alpha_1, coeff_1, pi, overlap_kk, overlap_kl, alpha_2, coeff_2
      real(kind=dp) :: sum_
!
         pi = 4.0d0*atan(1.0d0)
!  
         do i = 1, molecule%n_atoms
!  
            do j = 1, molecule%atoms(i)%n_shells
!  
                angmom = molecule%atoms(i)%shells(j)%l
                n_primitives = molecule%atoms(i)%shells(j)%get_n_primitives()
!  
                do k = 1, n_primitives
!  
                  alpha_1 = molecule%atoms(i)%shells(j)%get_exponent_i(k)
                  coeff_1 = molecule%atoms(i)%shells(j)%get_coefficient_i(k)
!
                  overlap_kk = ((pi/(2*alpha_1))**1.5d0)/(4.0d0*(alpha_1))**angmom
                  overlap_kk = double_factorial(2*angmom-1) * overlap_kk
!
                  coeff_1 = coeff_1 / sqrt(overlap_kk)
                  call molecule%atoms(i)%shells(j)%set_coefficient_i(k, coeff_1)
!  
                enddo
!  
                sum_ = 0.0d0
!  
                do k = 1, n_primitives
!  
                   alpha_1 = molecule%atoms(i)%shells(j)%get_exponent_i(k)
                   coeff_1 = molecule%atoms(i)%shells(j)%get_coefficient_i(k)
!
                   overlap_kk = ((pi/(2*alpha_1))**1.5d0)/(4.0d0*(alpha_1))**angmom
                   overlap_kk = double_factorial(2*angmom-1) * overlap_kk
                   sum_ = sum_ + coeff_1*coeff_1*overlap_kk
!  
                   do l = 1,k-1
!  
                     alpha_2 = molecule%atoms(i)%shells(j)%get_exponent_i(l)
                     coeff_2 = molecule%atoms(i)%shells(j)%get_coefficient_i(l)
!
                     overlap_kl = ((pi/(alpha_1+alpha_2))**1.5d0)/(2.0d0*(alpha_1+alpha_2))**angmom
                     overlap_kl = double_factorial(2*angmom-1) * overlap_kl
                     sum_ = sum_ + 2.0d0*coeff_1*coeff_2*overlap_kl
!  
                   enddo
!  
                enddo
!  
                do k =1, n_primitives
!  
                  coeff_1 = molecule%atoms(i)%shells(j)%get_coefficient_i(k)
                  coeff_1 = coeff_1/sqrt(sum_)
                  call molecule%atoms(i)%shells(j)%set_coefficient_i(k, coeff_1)
!  
                enddo
!  
            enddo
!  
         enddo
!  
!
   end subroutine normalize_raw_primitives_molecular_system
!
!
   function get_nuclear_repulsion_mm_molecular_system(molecule)
!!
!!    Get nuclear repulsion if non-polarizable QM/MM
!!    Written by Tommaso Giovannini, 2019
!!
!!    Calculates, and returns, the nuclear repulsion term for electrostatic
!!    embedding, in units of Hartree. Makes use of the Ångstrøm to Bohr conversion
!!    factor defined in the parameters module.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      real(dp) :: get_nuclear_repulsion_mm_molecular_system
!
      integer :: i = 0, j = 0
!
      real(dp) :: x_ij, y_ij, z_ij, r_ij
!
      get_nuclear_repulsion_mm_molecular_system = zero
!
      do i = 1, molecule%mm%n_atoms
         do j = i + 1, molecule%mm%n_atoms
!
            x_ij = molecule%mm%coordinates(1,i) - molecule%mm%coordinates(1,j)
            y_ij = molecule%mm%coordinates(2,i) - molecule%mm%coordinates(2,j)
            z_ij = molecule%mm%coordinates(3,i) - molecule%mm%coordinates(3,j)
!
            r_ij = sqrt(x_ij**2 + y_ij**2 + z_ij**2)
!
            r_ij = angstrom_to_bohr*r_ij
!
            if (abs(r_ij) .lt. 1.0D-7) then
!
               call output%error_msg('two atoms are placed on top of each other.')
!
            endif
!
            get_nuclear_repulsion_mm_molecular_system = get_nuclear_repulsion_mm_molecular_system &
                  + ((molecule%mm%charge(i))*(molecule%mm%charge(j)))/r_ij
!
         enddo
      enddo
!
!
      do i = 1, molecule%n_atoms
         do j = 1, molecule%mm%n_atoms
!
            x_ij = molecule%atoms(i)%x - molecule%mm%coordinates(1,j)
            y_ij = molecule%atoms(i)%y - molecule%mm%coordinates(2,j)
            z_ij = molecule%atoms(i)%z - molecule%mm%coordinates(3,j)
!
            r_ij = sqrt(x_ij**2 + y_ij**2 + z_ij**2)
!
            r_ij = angstrom_to_bohr*r_ij
!
            if (abs(r_ij) .lt. 1.0D-7) then
!
               call output%error_msg('two atoms are placed on top of each other.')
!
            endif
!
            get_nuclear_repulsion_mm_molecular_system = get_nuclear_repulsion_mm_molecular_system &
                  + ((molecule%atoms(i)%number_)*(molecule%mm%charge(j)))/r_ij
!
         enddo
      enddo
!
  end function get_nuclear_repulsion_mm_molecular_system
!
!
   subroutine check_if_basis_present_and_pure_molecular_system(molecule)
!!
!!    Check Basis Set Info
!!    Written by Tommaso Giovannini, Oct 2019
!!
!!    The subroutine looks for the basis set inside LIBINT_DATA_PATH
!!    if the basis file.g94 is not found than the program stops
!!
!!    Then it checks whether the basis set is defined in cartesian 
!!    or in pure.
!!    
!!    The default is pure.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: atom_index
      character(len=100) :: basis
!      
      character(len=300) :: libint_path
      character(len=300) :: file_name
!      
      logical :: exists_in_libint_path
!
!     Loop over atoms 
!
      do atom_index = 1, molecule%n_atoms
!
         call get_environment_variable("LIBINT_DATA_PATH",libint_path)
!         
         if (molecule%atoms(atom_index)%basis(1:3) .ne. 'aug') then
!        
            write(basis,'(a)') molecule%atoms(atom_index)%basis
            call convert_to_lowercase(basis)
!            
         else
!         
            write(basis,'(a)') molecule%atoms(atom_index)%basis(5:100)
            call convert_to_lowercase(basis)
            basis = 'augmentation-'//trim(basis)
!
         endif
!
         file_name = trim(libint_path) // '/' // trim(basis) // '.g94'
!
         inquire(file=file_name, exist=exists_in_libint_path)
!         
         if(exists_in_libint_path) then
!          
           continue
!
         else
!           
            call output%error_msg(trim(basis)//' basis set has not been found'//&
                                & new_line('a') //&
                                & '  Maybe you forgot -basis --basis_dir?')
!            
         endif
!
      enddo
!
      if(trim(basis).eq.'3-21g'.or.         &
         trim(basis).eq.'6-31g'.or.         &
         trim(basis).eq.'6-31+g'.or.        &
         trim(basis).eq.'6-31++g'.or.       &
         trim(basis).eq.'6-31g*'.or.        &
         trim(basis).eq.'6-31g**'.or.       &
         trim(basis).eq.'6-31+g*'.or.       &
         trim(basis).eq.'6-31+g**'.or.      &
         trim(basis).eq.'6-31++g**'.or.     &
         trim(basis).eq.'6-31g(d,p)'.or.    &
         trim(basis).eq.'6-31g(2df,p)'.or.  &
         trim(basis).eq.'6-31g(3df,3pd)') molecule%cartesian_gaussians_int = 1
!
      if (input%requested_keyword_in_section('cartesian gaussians', 'system')) then 
!
         molecule%cartesian_gaussians_int = 1
! 
      else if(input%requested_keyword_in_section('pure gaussians','system')) then
!
         molecule%cartesian_gaussians_int = 0
!
      endif 

   end subroutine check_if_basis_present_and_pure_molecular_system
!
!
   subroutine get_nuclear_dipole_active_molecular_system(molecule, mu_k)
!!
!!    Get nuclear dipole moment 
!!    Written by Eirik F. Kjønstad, Apr 2019
!!
!!    Modified by Tommaso Giovannini, Linda Goletto 
!!    and Anders Hutcheson, Oct 2019 
!!
!!    Used when an active space exists, only loops over the active atoms 
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp), dimension(3), intent(out) :: mu_k 
!
      integer :: j 
      integer :: n_active_atoms, n_active_spaces
!
      n_active_spaces  = molecule%n_active_atom_spaces
      n_active_atoms = molecule%active_atom_spaces(n_active_spaces)%last_atom
!
      mu_k = zero 
!
      do j = 1, n_active_atoms
!
         mu_k(1) = mu_k(1) + (molecule%atoms(j)%x)*(molecule%atoms(j)%number_)*(angstrom_to_bohr)
         mu_k(2) = mu_k(2) + (molecule%atoms(j)%y)*(molecule%atoms(j)%number_)*(angstrom_to_bohr)
         mu_k(3) = mu_k(3) + (molecule%atoms(j)%z)*(molecule%atoms(j)%number_)*(angstrom_to_bohr)
!
      enddo
!
   end subroutine get_nuclear_dipole_active_molecular_system
!
!
   subroutine get_nuclear_quadrupole_active_molecular_system(molecule, Q_k)
!!
!!    Get nuclear dipole moment 
!!    Written by Eirik F. Kjønstad, Apr 2019
!!
!!    Modified by Tommaso Giovannini, Linda Goletto 
!!    and Anders Hutcheson, Oct 2019 
!!
!!    Used when an active space exists, only loops over the active atoms
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      real(dp), dimension(6), intent(out) :: Q_k 
!
      integer :: j 
      integer :: n_active_atoms, n_active_spaces
!
      n_active_spaces  = molecule%n_active_atom_spaces
      n_active_atoms = molecule%active_atom_spaces(n_active_spaces)%last_atom
!      
      Q_k = zero 
!
      do j = 1, n_active_atoms
!
         Q_k(1) = Q_k(1) + (molecule%atoms(j)%x**2)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(2) = Q_k(2) + (molecule%atoms(j)%x*molecule%atoms(j)%y)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(3) = Q_k(3) + (molecule%atoms(j)%x*molecule%atoms(j)%z)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(4) = Q_k(4) + (molecule%atoms(j)%y**2)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(5) = Q_k(5) + (molecule%atoms(j)%y*molecule%atoms(j)%z)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
         Q_k(6) = Q_k(6) + (molecule%atoms(j)%z**2)*(molecule%atoms(j)%number_)*(angstrom_to_bohr**2)
!
      enddo
!
   end subroutine get_nuclear_quadrupole_active_molecular_system
!
!
   subroutine evaluate_aos_at_point_molecular_system(molecule, x, y, z, aos_at_point, n_ao)
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
!!    'n_ao'         : number of aos
!!
!
      use math_utilities, only: double_factorial, binomial, factorial
!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      integer :: n_ao
!
      real(dp), intent(in) :: x, y, z
      real(dp), dimension(n_ao), intent(out) :: aos_at_point
!
      integer :: ao, n_primitives, i, j, shell, atom, t, u, two_v, two_v_m, two_v_max, t_max
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
      ao = 1
!
!     Loop over atoms
!
      do atom = 1, molecule%n_atoms
!
!        Find the position of the point relative to the nucleus, in atomic units
!
         x_rel = (x - molecule%atoms(atom)%x)*angstrom_to_bohr
         y_rel = (y - molecule%atoms(atom)%y)*angstrom_to_bohr
         z_rel = (z - molecule%atoms(atom)%z)*angstrom_to_bohr
!
!        Calculate distance between the point and the nucleus
!
         r_squared = x_rel**2 + y_rel**2 + z_rel**2
!
         do shell = 1, molecule%atoms(atom)%n_shells
!
!           Determine angular momentum
!
            l = molecule%atoms(atom)%shells(shell)%l 
!
!           Determine radial part
!
            n_primitives = molecule%atoms(atom)%shells(shell)%get_n_primitives()
!
            overlap_primitives = zero
!
            do i = 1, n_primitives
               do j = 1, n_primitives
!
                  exponent_i = molecule%atoms(atom)%shells(shell)%get_exponent_i(i)
                  exponent_j = molecule%atoms(atom)%shells(shell)%get_exponent_i(j)
                  coefficient_i = molecule%atoms(atom)%shells(shell)%get_coefficient_i(i)
                  coefficient_j = molecule%atoms(atom)%shells(shell)%get_coefficient_i(j)
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
               exponent_i = molecule%atoms(atom)%shells(shell)%get_exponent_i(i)
               coefficient_i = molecule%atoms(atom)%shells(shell)%get_coefficient_i(i)
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
               aos_at_point(ao) = radial_part
!
            elseif (l == 1) then
!
!              Cartesian p-shell
!
!              p_x, p_y, p_z (first, second, third in CCA standard)
!
               aos_at_point(ao)   = radial_part*x_rel
               aos_at_point(ao+1) = radial_part*y_rel
               aos_at_point(ao+2) = radial_part*z_rel
!
            elseif (l == 2) then
!
!              ml = -2, angular part: sqrt(3)xy
!
               aos_at_point(ao)   = sqrt_three*x_rel*y_rel*radial_part
!
!              ml = -1, angular part: sqrt(3)yz
!
               aos_at_point(ao+1) = sqrt_three*y_rel*z_rel*radial_part
!
!              ml = 0, angular part: 1.5z^2 - 0.5r^2
!
               aos_at_point(ao+2) = half*(three*z_rel**2 - r_squared)*radial_part
!
!              ml = 1, angular part: sqrt(3)xz
!
               aos_at_point(ao+3) = sqrt_three*x_rel*z_rel*radial_part
!
!              ml = 2, angular part: 1/2 sqrt(3)(x^2 - y^2)
!  
               aos_at_point(ao+4) = half*sqrt_three*(x_rel**2 - y_rel**2)*radial_part
!
            elseif (l == 3) then
!
!              ml = -3, angular part: 1/2 sqrt(5/2) (3x^2 - y^2)y
!
               aos_at_point(ao) = half*sqrt_five_half*(three*x_rel**2 - y_rel**2)*y_rel*radial_part
!
!              ml = -2, angular part: sqrt(15) xyz
!
               aos_at_point(ao+1) = sqrt_fifteen*x_rel*y_rel*z_rel*radial_part
!
!              ml = -1, angular part: 1/2 sqrt(3/2) (5*z**2 - r**2)y
!              ml =  1, angular part: 1/2 sqrt(3/2) (5*z**2 - r**2)x
!
               aos_at_point(ao+2) = half*sqrt_three_half*(five*z_rel**2 - r_squared)*radial_part
               aos_at_point(ao+4) = aos_at_point(ao+2)*x_rel
               aos_at_point(ao+2) = aos_at_point(ao+2)*y_rel
!
!              ml = 0, angular part: 1/2(5z**2 - 3r^2)z
!
               aos_at_point(ao+3) = half*(five*z_rel**2 - three*r_squared)*z_rel*radial_part
!
!              ml = 2, angular_part: 1/2 sqrt(15)(x^2-y^2)z
!
               aos_at_point(ao+5) = half*sqrt_fifteen*(x_rel**2 - y_rel**2)*z_rel*radial_part
!
!              ml = 3, angular part: 1/2 sqrt(5/2)(x^2-3y^2)x
!
               aos_at_point(ao+6) = half*sqrt_five_half*(x_rel**2 - three*y_rel**2)&
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
                  if (ml == 0) N_S_lm = N_S_lm/sqrt(two)
!
                  angular_part = N_S_lm*angular_part
!
                  aos_at_point(ao + count_ml) = radial_part*angular_part
!
                  count_ml = count_ml + 1
!
               enddo
!
            endif
!
            ao = ao + 2*l + 1
!
         enddo ! end loop over shells
!
      enddo ! end loop over atoms
!          
   end subroutine evaluate_aos_at_point_molecular_system
!
!
   subroutine print_active_atoms_molecular_system(molecule)
!!
!!    Print active atoms
!!    Written by Sarai D. Folkestad, Dec 2019
!!
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer :: i, atom, first, last, n_total
!
      if (molecule%n_active_atom_spaces == 0) return
!
      call output%printf('m', 'Active atoms:', fs='(/t6, a)')
!
      call output%print_separator('m', 38,'=', fs='(/t6, a)')
      call output%printf('m',' Atom   Symbol       Basis     Method', fs='(t6, a)')
      call output%print_separator('m', 38,'=', fs='(t6, a)')
!
      n_total = 0
!
      do i = 1, molecule%n_active_atom_spaces
!
         first = molecule%active_atom_spaces(i)%first_atom
         last  = molecule%active_atom_spaces(i)%last_atom
         n_total = n_total + last - first + 1
!
         do atom = first, last
!
            call output%printf('m','(i4)      '// molecule%atoms(atom)%symbol //'      '      &
                              // trim(molecule%atoms(atom)%basis)                         &
                              // '       '// trim(molecule%active_atom_spaces(i)%level), &
                              ints=[molecule%atoms(atom)%input_number],                   &
                              fs='(t6,a)')
         enddo
      enddo
!
      call output%print_separator('m', 38,'=', fs='(t6, a)')
      call output%printf('m','Total number of active atoms: (i0)', ints=[n_total], fs='(t6, a)')
      call output%printf('m','OBS: Atoms will be reordered, active atoms first', fs='(t6, a)')
!
   end subroutine print_active_atoms_molecular_system
!
!
end module molecular_system_class
