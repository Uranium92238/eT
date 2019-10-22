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
      integer :: n_atoms
      integer :: n_basis_sets 
      integer :: charge
      integer :: multiplicity
      integer :: n_electrons 
      integer :: n_s
!
      integer(i6) :: cartesian_gaussians_int
!
      type(atomic), dimension(:), allocatable :: atoms
!
      type(interval), dimension(:), allocatable :: shell_limits 
!
      integer :: n_active_atoms_spaces
      type(active_atoms), dimension(:),allocatable :: active_atoms_spaces
!
      integer, dimension(:), allocatable :: shell2atom 
!
      logical :: active_atoms = .false.
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
!
      type(mm) :: mm
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
!
      procedure :: get_n_aos                                => get_n_aos_molecular_system
      procedure :: get_n_shells                             => get_n_shells_molecular_system
      procedure :: get_shell_limits                         => get_shell_limits_molecular_system
      procedure :: basis2shell                              => basis2shell_molecular_system
      procedure :: get_max_shell_size                       => get_max_shell_size_molecular_system
!
      procedure :: shell_to_atom                            => shell_to_atom_molecular_system
!
      procedure :: check_basis_set                             => check_basis_set_molecular_system
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
!
      procedure :: set_basis_info                           => set_basis_info_molecular_system
      procedure :: set_shell_basis_info                     => set_shell_basis_info_molecular_system
      procedure :: check_convert_pure_to_cartesian_basis    => check_convert_pure_to_cartesian_basis_molecular_system
      procedure :: normalize_raw_primitives                 => normalize_raw_primitives_molecular_system
      procedure :: get_nuclear_repulsion_mm                 => get_nuclear_repulsion_mm_molecular_system
!
      procedure :: basis2atom                               => basis2atom_molecular_system
!
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
!
      call molecule%read_settings()
!
      call molecule%prepare()
!
   end function new_molecular_system
!
!
   function new_molecular_system_from_parameters(atoms, name, charge, multiplicity, mm_calculation) result(molecule)
!!
!!    Initialize From Parameters
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      type(molecular_system) :: molecule
!
      type(atomic), dimension(:), intent(in) :: atoms
      character(len=100),         intent(in) :: name
      integer,                    intent(in) :: charge
      integer,                    intent(in) :: multiplicity
      logical,                    intent(in) :: mm_calculation
!
      allocate(molecule%atoms(size(atoms)))
!
      molecule%atoms          = atoms
      molecule%n_atoms        = size(atoms)
      molecule%name_          = name
      molecule%charge         = charge
      molecule%multiplicity   = multiplicity
      molecule%mm_calculation = mm_calculation
!
      molecule%cartesian_gaussians_int = 0
      molecule%n_active_atoms_spaces = 0
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
!     First have a look to the basis set
!
      call molecule%check_basis_set()
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
!!      if (input%requested_keyword_in_section('cartesian gaussians', 'system')) then 
!!!
!!         molecule%cartesian_gaussians_int = 1
!!!
!!      endif 
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
!!
!!    Read atoms and their coordinates, assumed to be in units of Ångstrøm.
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
      molecule%n_atoms = input%get_n_atoms()
!
      call mem%alloc(positions, molecule%n_atoms, 3)
!
      allocate(symbols(molecule%n_atoms))
      allocate(basis_sets(molecule%n_atoms))
!
      call input%get_geometry(molecule%n_atoms, symbols, positions, basis_sets)
!
      call molecule%initialize_atoms()
!
!     1. Place the first atom in atoms
!
      current_atom = 1
!
      molecule%atoms(current_atom)%symbol       = symbols(current_atom)
      molecule%atoms(current_atom)%basis        = basis_sets(current_atom)
      molecule%atoms(current_atom)%x            = positions(current_atom,1)
      molecule%atoms(current_atom)%y            = positions(current_atom,2)
      molecule%atoms(current_atom)%z            = positions(current_atom,3)
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
               molecule%atoms(current_atom)%x            = positions(atom,1)
               molecule%atoms(current_atom)%y            = positions(atom,2)
               molecule%atoms(current_atom)%z            = positions(atom,3)
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
               molecule%atoms(current_atom)%x            = positions(atom,1)
               molecule%atoms(current_atom)%y            = positions(atom,2)
               molecule%atoms(current_atom)%z            = positions(atom,3)
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
      call mol_file%printf('(i5)',ints=[molecule%n_atoms], fs='(a/)')
!
      do atom = 1, molecule%n_atoms
!
         call mol_file%printf('(a2)   (f21.16)   (f21.16)   (f21.16)', &
                              chars = [molecule%atoms(atom)%symbol],   &
                              reals = [molecule%atoms(atom)%x,         &
                                       molecule%atoms(atom)%y,         &
                                       molecule%atoms(atom)%z],        &
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
         call basis_file%printf('(i5)',ints=[n_atoms_in_basis(current_basis_nbr)], fs='(a/)')
!
         do i = 1, n_atoms_in_basis(current_basis_nbr)
!
            atom = i + atom_offset
!
            call basis_file%printf('(a2)   (f21.16)   (f21.16)   (f21.16)', &
                                   chars = [molecule%atoms(atom)%symbol],   &
                                   reals = [molecule%atoms(atom)%x,         &
                                            molecule%atoms(atom)%y,         &
                                            molecule%atoms(atom)%z],        &
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
      molecule%n_active_atoms_spaces = 0
!
      do method = 1, size(methods)
!
         if (n_active(method) > 0) molecule%n_active_atoms_spaces = molecule%n_active_atoms_spaces + 1
!
      enddo
!
      if (molecule%n_active_atoms_spaces == 0) call output%error_msg('Requested active atoms, but no active spaces found.')
!
      allocate(molecule%active_atoms_spaces(molecule%n_active_atoms_spaces))
!
      call output%printf('Active atoms:', fs='(/t6, a)', pl='minimal')
!
      call output%printf('--------------------------------------', fs='(t6, a)', pl='minimal')
      call output%printf(' Atom   Symbol       Basis     Method', fs='(t6, a)', pl='minimal')
      call output%printf('--------------------------------------', fs='(t6, a)', pl='minimal')
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
            call output%printf('(i4)      '// molecule%atoms(atom)%symbol //'      '                     &
                               // trim(molecule%atoms(atom)%basis)// '       '// trim(methods(method)),  &
                     ints=[molecule%atoms(atom)%input_number], fs='(t6,a)', pl='minimal')
!
         enddo        
!
         if (n_active(method) > 0) then
! 
            space = space + 1
!
            molecule%active_atoms_spaces(space)%level = trim(methods(method))
!
            molecule%active_atoms_spaces(space)%first_atom = offset + 1
            molecule%active_atoms_spaces(space)%last_atom = offset + n_active(method)
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
      call output%printf('--------------------------------------', fs='(t6, a)', pl='minimal')
      call output%printf('Total number of active atoms: ', ints=[n_total], fs='(t6, a)', pl='minimal')
      call output%printf('OBS: Atoms will be reordered, active atoms first', fs='(t6, a)', pl='minimal')
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
   subroutine print_geometry_molecular_system(molecule)
!!
!!    Print geometry 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(molecular_system) :: molecule  
!
      integer :: I 
!
      character(len=12) :: frmt0
      character(len=46) :: frmt1
!
      frmt0="(t5,68('='))"
      frmt1="(t5,'Atom',9x,'X',16x,'Y',16x,'Z',13x,'Basis')"
      write(output%unit,'(/)')
      write(output%unit,frmt0 )
      write(output%unit,'(t32,a)') 'Geometry (Å)'
      write(output%unit,frmt0 )
      write(output%unit,frmt1 )
      write(output%unit,frmt0 )
!
      do I = 1, molecule%n_atoms 
!
         write(output%unit, '(t6, a2, f17.12, f17.12, f17.12, 3x, a11)')  molecule%atoms(I)%symbol, &
                                                                          molecule%atoms(I)%x,      &
                                                                          molecule%atoms(I)%y,      &
                                                                          molecule%atoms(I)%z,      &
                                                                          adjustr(trim(molecule%atoms(I)%basis))
!  
         flush(output%unit)
!
      enddo 
!
      write(output%unit,frmt0 )
      write(output%unit,'(t30,a)') 'Geometry (a.u.)'
      write(output%unit,frmt0 )
      write(output%unit,frmt1 )
      write(output%unit,frmt0 )
      do I = 1, molecule%n_atoms 
!
         write(output%unit, '(t6, a2, f17.12, f17.12, f17.12, 3x, a11)')  molecule%atoms(I)%symbol, &
                                                                          angstrom_to_bohr*molecule%atoms(I)%x,      &
                                                                          angstrom_to_bohr*molecule%atoms(I)%y,      &
                                                                          angstrom_to_bohr*molecule%atoms(I)%z,      &
                                                                          adjustr(trim(molecule%atoms(I)%basis))
!  
         flush(output%unit)
!
      enddo 
!      
      write(output%unit,frmt0 )
!
   end subroutine print_geometry_molecular_system
!
!
   subroutine print_system_molecular_system(molecule)
!!
!!    Print geometry 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!    Modified by Tommaso Giovannini, March 2019
!!
      implicit none 
!
      class(molecular_system) :: molecule  
!
      write(output%unit, '(/t3,a)')       '- Molecular system specifications:'
!
      write(output%unit, '(/t6,a14,a)')      'Name:         ', trim(molecule%name_)
      write(output%unit, '(t6,a14,i1)')      'Charge:       ', molecule%charge 
      write(output%unit, '(t6,a14,i1)')      'Multiplicity: ', molecule%multiplicity 
!
      if (molecule%cartesian_gaussians_int.eq.1) call output%printf('Using Cartesian gaussians.', pl='m', fs='(/t6,a)')
!
      write(output%unit, '(/t6,a27,i5)')     'Pure basis functions:      ', molecule%n_pure_basis
      write(output%unit, '(t6,a27,i5)')      'Cartesian basis functions: ', molecule%n_cart_basis
      write(output%unit, '(t6,a27,i5)')      'Primitive basis functions: ', molecule%n_primitives_cart
!
      write(output%unit, '(/t6,a35,f25.12)') 'Nuclear repulsion energy (a.u.):   ', molecule%get_nuclear_repulsion()
      write(output%unit, '(t6,a35,f25.12)')  'Bohr/angstrom value (CODATA 2010): ', bohr_to_angstrom
!
      flush(output%unit)
!
      call molecule%print_geometry()
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
      integer :: i, first_atom, last_atom
!
      first_atom = 0
      last_atom = 0 
!
      do i = 1, molecule%n_active_atoms_spaces
!
         if (trim(molecule%active_atoms_spaces(i)%level) == trim(level)) then
!
            first_atom = molecule%active_atoms_spaces(i)%first_atom
            last_atom = molecule%active_atoms_spaces(i)%last_atom
            exit
!
         endif
!
      enddo
!
      if (first_atom == 0 .or. last_atom == 0) &
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
      last = first - 1
!
      do i = first_atom, last_atom
!
         last = last + molecule%atoms(i)%n_ao
!
      enddo
!
   end subroutine first_and_last_ao_active_space_molecular_system
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
         call molecule%atoms(atom_index)%shells(shell)%basis_details%set_n_primitives(n_primitive)
         call molecule%atoms(atom_index)%shells(shell)%basis_details%initialize_exponents()
         call molecule%atoms(atom_index)%shells(shell)%basis_details%initialize_coefficients()
!
!        Loop over primitives and set coefficient and exponent
!
         if(ang_mom == "SP") then
!
!           In case of "SP" shell, split S and P coefficients
!
            call molecule%atoms(atom_index)%shells(shell+1)%basis_details%set_n_primitives(n_primitive)
            call molecule%atoms(atom_index)%shells(shell+1)%basis_details%initialize_exponents()
            call molecule%atoms(atom_index)%shells(shell+1)%basis_details%initialize_coefficients()
!
            do primitive = 1, n_primitive
!
!              Read coefficient and exponent
!
               call basis_set_file%read_(line,'(a200)')
               read(line, *) exponent_, coefficient, coefficient_2
!
               call molecule%atoms(atom_index)%shells(shell)%basis_details%set_exponent_i(primitive, exponent_)
               call molecule%atoms(atom_index)%shells(shell+1)%basis_details%set_exponent_i(primitive, exponent_)
               call molecule%atoms(atom_index)%shells(shell)%basis_details%set_coefficient_i(primitive, coefficient)
               call molecule%atoms(atom_index)%shells(shell+1)%basis_details%set_coefficient_i(primitive, coefficient_2)
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
               call molecule%atoms(atom_index)%shells(shell)%basis_details%set_exponent_i(primitive, exponent_)
               call molecule%atoms(atom_index)%shells(shell)%basis_details%set_coefficient_i(primitive, coefficient)
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
      integer :: i,j, angmom1, n_prim1, n_func
!
      do i = 1, molecule%n_atoms
!  
         do j = 1, molecule%atoms(i)%n_shells
!  
             angmom1 = molecule%atoms(i)%shells(j)%l
             n_prim1 = int(molecule%atoms(i)%shells(j)%basis_details%n_primitives,kind(n_prim1))
             n_func  = molecule%atoms(i)%shells(j)%length
             molecule%n_pure_basis = molecule%n_pure_basis + n_func 
!  
             if(angmom1.ge.1) then
!  
                molecule%atoms(i)%shells(j)%size_cart = n_func
!  
                if(angmom1.ge.2) then
!  
                   if(angmom1.eq.2.and.n_func.eq.6) continue
!  
                   if(angmom1.eq.2.and.n_func.eq.5) then
!  
                      molecule%atoms(i)%shells(j)%size_cart = n_func + 1
!  
                   endIf
!  
                   if(angmom1.eq.3.and.n_func.eq.10) continue
!  
                   if(angmom1.eq.3.and.n_func.eq.7) then
!  
                      molecule%atoms(i)%shells(j)%size_cart = n_func + 3
!  
                   endIf
!  
                   if(angmom1.eq.4.and.n_func.eq.15) continue
!  
                   if(angmom1.eq.4.and.n_func.eq.9) then
!  
                      molecule%atoms(i)%shells(j)%size_cart = n_func + 6
!  
                   endIf
!  
                   if(angmom1.gt.4) then
!  
                      write(output%unit,'(a)') 'Cartesian G functions NYI'
                      Stop
!  
                   endIf
!  
                endif
!  
            else 
!  
               molecule%atoms(i)%shells(j)%size_cart = n_func 
!  
            endIf
!  
            molecule%n_cart_basis      = molecule%n_cart_basis + molecule%atoms(i)%shells(j)%size_cart
            molecule%n_primitives_cart = molecule%n_primitives_cart + molecule%atoms(i)%shells(j)%size_cart*n_prim1
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
      integer :: i,j,k,l, angmom, n_prim
      real(kind=dp) :: alpha1, coeff1, pi, overlap_kk, overlap_kl,alpha2,coeff2
      real(kind=dp) :: sum_
!
!
         pi = 4.0d0*atan(1.0d0)
!  
         do i = 1, molecule%n_atoms
!  
            do j = 1, molecule%atoms(i)%n_shells
!  
                angmom = molecule%atoms(i)%shells(j)%l
                n_prim = int(molecule%atoms(i)%shells(j)%basis_details%n_primitives,kind(n_prim))
!  
                do k = 1, n_prim
!  
                   alpha1 = molecule%atoms(i)%shells(j)%basis_details%exponents(k)
                   coeff1 = molecule%atoms(i)%shells(j)%basis_details%coefficients(k)
                   overlap_kk = ((pi/(2*alpha1))**1.5d0)/(4.0d0*(alpha1))**angmom
                   overlap_kk = double_factorial(2*angmom-1) * overlap_kk
                   coeff1 = coeff1 / sqrt(overlap_kk)
                   molecule%atoms(i)%shells(j)%basis_details%coefficients(k) = coeff1
!  
                enddo
!  
                sum_ = 0.0d0
!  
                do k = 1, n_prim
!  
                   alpha1 = molecule%atoms(i)%shells(j)%basis_details%exponents(k)
                   coeff1 = molecule%atoms(i)%shells(j)%basis_details%coefficients(k)
                   overlap_kk = ((pi/(2*alpha1))**1.5d0)/(4.0d0*(alpha1))**angmom
                   overlap_kk = double_factorial(2*angmom-1) * overlap_kk
                   sum_ = sum_ + coeff1*coeff1*overlap_kk
!  
                   do l = 1,k-1
!  
                      alpha2 = molecule%atoms(i)%shells(j)%basis_details%exponents(l)
                      coeff2 = molecule%atoms(i)%shells(j)%basis_details%coefficients(l)
                      overlap_kl = ((pi/(alpha1+alpha2))**1.5d0)/(2.0d0*(alpha1+alpha2))**angmom
                      overlap_kl = double_factorial(2*angmom-1) * overlap_kl
                      sum_ = sum_ + 2.0d0*coeff1*coeff2*overlap_kl
!  
                   enddo
!  
                enddo
!  
                do k =1, n_prim
!  
                   coeff1 = molecule%atoms(i)%shells(j)%basis_details%coefficients(k)
                   coeff1 = coeff1/sqrt(sum_)
                   molecule%atoms(i)%shells(j)%basis_details%coefficients(k) = coeff1
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
   subroutine check_basis_set_molecular_system(molecule)
!!
!!    Check Basis Path
!!    Written by Tommaso Giovannini, Oct 2019
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
                                & '  Maybe you forgot -ldd path_to_basis?')
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

   end subroutine check_basis_set_molecular_system
!
!
end module molecular_system_class
