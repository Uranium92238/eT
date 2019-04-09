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
   use atomic_class
   use io_utilities
   use interval_class
   use libint_initialization
   use ao_integral_tool_class
   use active_atoms_info_class
!
   implicit none
!
   include "../../libint/atom_init_cdef.F90"
!
   type :: molecular_system
!
      character(len=200) :: name
      character(len=100), dimension(:), allocatable :: basis_sets
!
      integer :: n_atoms
      integer :: n_basis_sets 
      integer :: charge
      integer :: multiplicity
      integer :: n_electrons 
      integer :: n_s
!
      type(atomic), dimension(:), allocatable :: atoms
!
      type(ao_integral_tool) :: ao_integrals
!
      type(interval), dimension(:), allocatable :: shell_limits 
!
      logical :: active_atoms = .false.
!
      integer :: n_active_atoms = 0
!
   contains
!
      procedure :: prepare                                  => prepare_molecular_system
      procedure :: cleanup                                  => cleanup_molecular_system
!
      procedure, private :: write_libint_files              => write_libint_files_molecular_system
!
      procedure, private :: read_settings                   => read_settings_molecular_system
      procedure, private :: read_system                     => read_system_molecular_system
      procedure, private :: read_geometry                   => read_geometry_molecular_system
      procedure, private :: read_active_atoms               => read_active_atoms_molecular_system
!
      procedure :: print_system                             => print_system_molecular_system
      procedure :: print_geometry                           => print_geometry_molecular_system
!
      procedure :: get_nuclear_repulsion                    => get_nuclear_repulsion_molecular_system
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
      procedure :: initialize_basis_sets                    => initialize_basis_sets_molecular_system
      procedure :: initialize_atoms                         => initialize_atoms_molecular_system
      procedure :: initialize_shell_limits                  => initialize_shell_limits_molecular_system
!
      procedure :: destruct_basis_sets                      => destruct_basis_sets_molecular_system
      procedure :: destruct_atoms                           => destruct_atoms_molecular_system
      procedure :: destruct_shell_limits                    => destruct_shell_limits_molecular_system
!
      procedure :: translate_from_input_order_to_eT_order   => translate_from_input_order_to_eT_order_molecular_system
!
      procedure :: construct_ao_h_wx                        => construct_ao_h_wx_molecular_system      
!
   end type molecular_system
!
contains
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
   end subroutine read_settings_molecular_system
!
!
   subroutine prepare_molecular_system(molecule)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      character(len=100) :: temp_name
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
!     Read eT.inp and write files for Libint
!
      molecule%charge = 0
      molecule%multiplicity = 1
!
      call molecule%read_settings()
      call molecule%write_libint_files()
!
!     Initialize libint with atoms and basis sets
!
      call initialize_atoms(molecule%name)
!
      do i = 1, molecule%n_basis_sets
!
         write(temp_name, '(a, a1, i4.4)') trim(molecule%name), '_', i
!
         call initialize_basis(molecule%basis_sets(i), temp_name)
!
      enddo
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
         call get_n_basis_in_shells_c(k, n_basis_in_shells)
!
         molecule%atoms(k)%n_ao = 0
!
         do j = 1, n_shells_on_atoms(k)

            molecule%atoms(k)%shells(j)%size = n_basis_in_shells(j)
!
            molecule%atoms(k)%n_ao = molecule%atoms(k)%n_ao + n_basis_in_shells(j)

         enddo
!
         deallocate(n_basis_in_shells)
!
!        Get shell numbers
!
         allocate(shell_numbers(n_shells_on_atoms(k)))
         call get_shell_numbers_c(k, shell_numbers)
!
         do j = 1, n_shells_on_atoms(k)
!
            molecule%atoms(k)%shells(j)%number_ = shell_numbers(j)
!
         enddo
!
         deallocate(shell_numbers)
!
!        And the first AO index in each shell
!
         allocate(first_ao_in_shells(n_shells_on_atoms(k)))
         call get_first_ao_in_shells_c(k, first_ao_in_shells)
!
         do j = 1, n_shells_on_atoms(k)

            molecule%atoms(k)%shells(j)%first = first_ao_in_shells(j)
!
         enddo
!
         deallocate(first_ao_in_shells)
!
!        Then determine the angular momentum of shells & the last AO index
!
         do j = 1, n_shells_on_atoms(k)
!
            call molecule%atoms(k)%shells(j)%determine_angular_momentum()
            call molecule%atoms(k)%shells(j)%determine_last_ao_index()
!
         enddo
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
      if (molecule%charge .ne. 0) then
!
         call output%error_msg('SAD not yet implemented for charged species!')
!
      endif
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
      call molecule%print_system()
!
   end subroutine prepare_molecular_system
!
!
   subroutine cleanup_molecular_system(molecule)
!!
!!    Initialize
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
!
   end subroutine cleanup_molecular_system
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
!
      call input%get_required_keyword_in_section('name', 'system', molecule%name)
!
      call input%get_keyword_in_section('charge', 'system', molecule%charge)
      call input%get_keyword_in_section('multiplicity', 'system', molecule%multiplicity)
!
   end subroutine read_system_molecular_system
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
      molecule%atoms(current_atom)%input_nbr = current_atom
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
               molecule%atoms(current_atom)%input_nbr = atom
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
               molecule%atoms(current_atom)%input_nbr = atom
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
      type(file) :: mol_file, basis_file
!
      integer :: basis_set_counter, atom_offset, current_basis_nbr, i
!
      integer, dimension(:), allocatable :: n_atoms_in_basis
!
!     Write atom file
!
      call mol_file%init(trim(molecule%name) // '.xyz', 'sequential', 'formatted')
      call disk%open_file(mol_file, 'write', 'rewind')
!
      write(mol_file%unit, '(i5/)') molecule%n_atoms
!
      do atom = 1, molecule%n_atoms
!
         write(mol_file%unit, '(a2, 3x, f21.16, 3x, f21.16, 3x, f21.16)')  &
                                    molecule%atoms(atom)%symbol,           &
                                    molecule%atoms(atom)%x,                &
                                    molecule%atoms(atom)%y,                &
                                    molecule%atoms(atom)%z
!
      enddo
!
      call disk%close_file(mol_file)
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
         write(temp_name, '(a, a1, i4.4, a4)')trim(molecule%name), '_', current_basis_nbr,  '.xyz'
!
         call basis_file%init(trim(temp_name), 'sequential', 'formatted')
         call disk%open_file(basis_file, 'write', 'rewind')
!
         write(basis_file%unit, '(i5/)') n_atoms_in_basis(current_basis_nbr)
!
         do i = 1, n_atoms_in_basis(current_basis_nbr)
!
            atom = i + atom_offset
!
            write(basis_file%unit, '(a2, 3x, f21.16, 3x, f21.16, 3x, f21.16)') &
                                       molecule%atoms(atom)%symbol,           &
                                       molecule%atoms(atom)%x,                &
                                       molecule%atoms(atom)%y,                &
                                       molecule%atoms(atom)%z
!
         enddo
!
         call disk%close_file(basis_file)
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
      character(len=200) :: active_basis, selection_type
!
      integer :: i, j, active_atom_counter, atom_counter !ioerror = 0, first, last
      integer :: central_atom
!
      integer, dimension(:), allocatable :: active_atoms, active_atoms_copy
!
      type(atomic), dimension(:), allocatable :: atoms_copy
!
      logical :: found
!
      real(dp) :: hf_radius, x, y, z
!
!     Do nothing if active atoms not requsted
!
      if (.not. input%requested_section('active atoms')) return
!
!     Find selection type
! 
      call input%get_keyword_in_section('selection type', 'active atoms', selection_type)
!
!     For the given selection type get the active atoms
!
      if (trim(selection_type) == 'range' .or. trim(selection_type) == 'list') then
!
!        Get the nuumber of elements in the list/range
!
         molecule%n_active_atoms = input%get_n_elements_for_keyword_in_section('hf', 'active atoms')
!
         call mem%alloc(active_atoms, molecule%n_active_atoms)
!
!        Get the active atoms
!
         call input%get_array_for_keyword_in_section('hf', 'active atoms', &
                                                      molecule%n_active_atoms, active_atoms)
!
!        Translate atom list to current eT ordering (based on basis set)
!
         call mem%alloc(active_atoms_copy, molecule%n_active_atoms)
         active_atoms_copy = active_atoms
!
         call molecule%translate_from_input_order_to_eT_order(molecule%n_active_atoms, active_atoms_copy, active_atoms)
!
         call mem%dealloc(active_atoms_copy, molecule%n_active_atoms)
!
      elseif (selection_type == 'central atom') then
!
         call input%get_keyword_in_section('central atom', 'active atoms', central_atom)
         call input%get_keyword_in_section('hf', 'active atoms', hf_radius)
!
!        Central atom in current eT ordering:
!
         do i = 1, molecule%n_atoms
!
            if (molecule%atoms(i)%input_nbr == central_atom) then
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
         molecule%n_active_atoms = 0
!
         do i = 1, molecule%n_atoms 
!
            x = (molecule%atoms(central_atom)%x - molecule%atoms(i)%x)
            y = (molecule%atoms(central_atom)%y - molecule%atoms(i)%y)
            z = (molecule%atoms(central_atom)%z - molecule%atoms(i)%z)
!
            if (sqrt(x**2 + y**2 + z**2) .le. hf_radius) &
                     molecule%n_active_atoms = molecule%n_active_atoms + 1
!
         enddo
!
         call mem%alloc(active_atoms, molecule%n_active_atoms)
!
         active_atom_counter = 0
!
         do i = 1, molecule%n_atoms 
!
           x = (molecule%atoms(central_atom)%x - molecule%atoms(i)%x)
           y = (molecule%atoms(central_atom)%y - molecule%atoms(i)%y)
           z = (molecule%atoms(central_atom)%z - molecule%atoms(i)%z)
!
           if (sqrt(x**2 + y**2 + z**2) .le. hf_radius) then 
!
              active_atom_counter = active_atom_counter + 1
!
              active_atoms(active_atom_counter) = i
!
           endif
!
         enddo
!
      else
!
         call output%error_msg('did not recognize selection type for active atoms.')
!
      endif
!
!     Find and set active basis
!
      if (input%requested_keyword_in_section('active basis', 'active atoms')) then
!
         call input%get_keyword_in_section('active basis', 'active atoms', active_basis)
!
         do i = 1, molecule%n_active_atoms
!
            molecule%atoms(active_atoms(i))%basis = trim(active_basis)
!
         enddo
!
      endif
!
!     Print active atoms
!
      
      write(output%unit, '(/t6, a)')'Active atoms:'
!
      write(output%unit, '(t6, a18)')'------------------'
      write(output%unit, '(t6, a18)')' Atom      Symbol '
      write(output%unit, '(t6, a18)')'------------------'
!
      do i = 1, molecule%n_active_atoms
!
         write(output%unit, '(t6, i5, 11x, a2)') molecule%atoms(active_atoms(i))%input_nbr, molecule%atoms(active_atoms(i))%symbol
!
      enddo
!
      write(output%unit, '(t6, a18)')'------------------'
      write(output%unit, '(t6, a30, i4)')'Total number of active atoms: ', molecule%n_active_atoms
      write(output%unit, '(t6, a/)')'OBS: Atoms will be reordered, active atoms first.'
      flush(output%unit)
!
!     Reorder atoms
!
      allocate(atoms_copy(molecule%n_atoms))
      atoms_copy = molecule%atoms
!
      do i = 1, molecule%n_active_atoms
!
         molecule%atoms(i) = atoms_copy(active_atoms(i))
!
      enddo
!
      atom_counter = molecule%n_active_atoms + 1
!
      do i = 1, molecule%n_atoms
!
         found = .false.
!
         do j = 1, molecule%n_active_atoms
!
            if(i == active_atoms(j)) then
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
            molecule%atoms(atom_counter) = atoms_copy(i)
            atom_counter = atom_counter + 1
!
         endif 
!
      enddo
!
      deallocate(atoms_copy)
!
   end subroutine read_active_atoms_molecular_system
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
                           + molecule%atoms(I)%shells(J)%size
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
      integer :: s1 
!
      max_shell_size = 0
      do s1 = 1, molecule%n_s 
!
         if (max_shell_size .lt. molecule%shell_limits(s1)%size) max_shell_size = molecule%shell_limits(s1)%size
!
      enddo
!
   end subroutine get_max_shell_size_molecular_system
!
!
   function get_shell_limits_molecular_system(molecule, A)
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
      type(interval) :: get_shell_limits_molecular_system
!
      get_shell_limits_molecular_system%first = 0
      get_shell_limits_molecular_system%last  = 0
      get_shell_limits_molecular_system%size  = 0
!
      do I = 1, molecule%n_atoms
!
         do J = 1, molecule%atoms(I)%n_shells
!
            if (A .eq. molecule%atoms(I)%shells(J)%number_) then
!
               get_shell_limits_molecular_system%first = molecule%atoms(I)%shells(J)%first
               get_shell_limits_molecular_system%last  = molecule%atoms(I)%shells(J)%last
               get_shell_limits_molecular_system%size  = molecule%atoms(I)%shells(J)%size
!
            endif
!
         enddo
!
      enddo
!
       if (get_shell_limits_molecular_system%size == 0) call output%error_msg('in get_shell_limits.')
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
!!    destruct basis sets
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
!!    Initialize basis sets
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
      write(output%unit, *)
!
      do I = 1, molecule%n_atoms 
!
         write(output%unit, '(t6, a2, f17.12, f17.12, f17.12, 3x, a11)')  molecule%atoms(I)%symbol, &
                                                                          molecule%atoms(I)%x,      &
                                                                          molecule%atoms(I)%y,      &
                                                                          molecule%atoms(I)%z,      &
                                                                          molecule%atoms(I)%basis 
!  
         flush(output%unit)
!
      enddo 
!
   end subroutine print_geometry_molecular_system
!
!
   subroutine print_system_molecular_system(molecule)
!!
!!    Print geometry 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(molecular_system) :: molecule  
!
      write(output%unit, '(/t3,a)')       '- Molecular system specifications:'
!
      write(output%unit, '(/t6,a14,a)')      'Name:         ', trim(molecule%name)
      write(output%unit, '(t6,a14,i1)')      'Charge:       ', molecule%charge 
      write(output%unit, '(t6,a14,i1)')      'Multiplicity: ', molecule%multiplicity 
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
            if (array_input_ordering(i) == molecule%atoms(j)%input_nbr) then
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
end module molecular_system_class
