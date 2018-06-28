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
   use atom_init
   use libint_initialization
!
   implicit none
!
   type :: molecular_system
!
      character(len=100) :: name
!
      integer(i15) :: n_atoms
      integer(i15) :: charge
!
      type(atomic), dimension(:,:), allocatable :: atoms
!
   contains
!
      procedure :: initialize => initialize_molecular_system
      procedure :: write      => write_molecular_system
!
      procedure, private :: read_info     => read_info_molecular_system
      procedure, private :: read_geometry => read_geometry_molecular_system
!
      procedure :: get_nuclear_repulsion => get_nuclear_repulsion_molecular_system
      procedure :: get_n_electrons       => get_n_electrons_molecular_system
!
      procedure :: get_n_aos => get_n_aos_molecular_system
      procedure :: get_n_shells => get_n_shells_molecular_system
      procedure :: get_shell_limits => get_shell_limits_molecular_system
      procedure :: basis2shell => basis2shell_molecular_system
!
   end type molecular_system
!
contains
!
!
   subroutine initialize_molecular_system(molecule)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(kind=4) :: i = 0, j = 0
!
      integer(kind=4), dimension(:,:), allocatable :: n_shells_on_atoms
      integer(kind=4), dimension(:,:), allocatable :: n_basis_in_shells
      integer(kind=4), dimension(:,:), allocatable :: first_ao_in_shells
      integer(kind=4), dimension(:,:), allocatable :: shell_numbers
!
      call molecule%read_info
!
      write(output%unit, *) 'N atoms? ', molecule%n_atoms
      flush(output%unit)
!
      allocate(molecule%atoms(molecule%n_atoms, 1))
!
      allocate(n_shells_on_atoms(molecule%n_atoms,1))
      n_shells_on_atoms = 0
!
      call molecule%read_geometry
!
      do i = 1, molecule%n_atoms
!
         call molecule%atoms(i, 1)%set_number()
!
      enddo
!
      call molecule%write      ! Write an xyz file for the read geometry
!
         write(output%unit, *) 'hei???'
         flush(output%unit)
!
      call initialize_basis()
      call get_n_shells_on_atoms(n_shells_on_atoms)
         write(output%unit, *) 'here maybe???'
         flush(output%unit)
!
      do i = 1, molecule%n_atoms ! Loop over atoms
!
!        Set atomic number
!
       !  call molecule%atoms(i, 1)%set_number()
!
!        Allocate and initialize the corresponding shells
!
         molecule%atoms(i,1)%n_shells = n_shells_on_atoms(i,1)
         write(output%unit, *) 'n_shells', molecule%atoms(i,1)%n_shells
         flush(output%unit)
         allocate(molecule%atoms(i,1)%shells(molecule%atoms(i,1)%n_shells, 1))
!
!        Then determine the number of basis functions in each shell
!
         allocate(n_basis_in_shells(n_shells_on_atoms(i,1), 1))
         call get_n_basis_in_shells(i, n_basis_in_shells)

         do j = 1, n_shells_on_atoms(i,1)

            molecule%atoms(i,1)%shells(j,1)%size = n_basis_in_shells(j,1)

         enddo

         deallocate(n_basis_in_shells)
!
!        Get shell numbers
!
         allocate(shell_numbers(n_shells_on_atoms(i,1), 1))
         call get_shell_numbers(i, shell_numbers)
!
         do j = 1, n_shells_on_atoms(i,1)
!
            molecule%atoms(i,1)%shells(j,1)%number = shell_numbers(j, 1)
            write(output%unit, *) 'The ', j, 'th shell on atom ', i, ' has shell nr. ', molecule%atoms(i,1)%shells(j,1)%number
!
         enddo
!
         deallocate(shell_numbers)
!
!        And the first AO index in each shell
!
         allocate(first_ao_in_shells(n_shells_on_atoms(i,1), 1))
         call get_first_ao_in_shells(i, first_ao_in_shells)
!
         do j = 1, n_shells_on_atoms(i,1)

            molecule%atoms(i,1)%shells(j,1)%first = first_ao_in_shells(j,1)
!
         enddo
!
         deallocate(first_ao_in_shells)
!
!        Then determine the angular momentum of shells & the last AO index
!
         do j = 1, n_shells_on_atoms(i,1)
!
            call molecule%atoms(i,1)%shells(j,1)%determine_angular_momentum()
            call molecule%atoms(i,1)%shells(j,1)%determine_last_ao_index()
!
         enddo
!
      enddo
!
   end subroutine initialize_molecular_system
!
!
   subroutine read_info_molecular_system(molecule)
!!
!!    Read information
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the number of atoms, charge and name of the molecule.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      character(len=100) :: line
      character(len=100) :: current_basis
!
      integer(i15) :: i = 0
!
      call input%init('eT.inp', 'sequential', 'formatted')
      call disk%open_file(input, 'read')
      rewind(input%unit)
!
      molecule%n_atoms = 0
!
      read(input%unit,'(a)') line
      line = remove_preceding_blanks(line)
!
      do while (trim(line) .ne. 'end geometry')

         if (line(1:6) == 'basis:' .or.  &
             line(1:6) == 'Basis:' .or.  &
             line(1:6) == 'BASIS:' ) then

            read(input%unit,'(a)') line
            line = remove_preceding_blanks(line)

            do while (trim(line) .ne. 'end geometry'  .and.  &
                       line(1:6) .ne. 'basis:'        .and.  &
                       line(1:6) .ne. 'Basis:'        .and.  &
                       line(1:6) .ne. 'BASIS:')

               molecule%n_atoms = molecule%n_atoms + 1

               read(input%unit,'(a)') line
               line = remove_preceding_blanks(line)

            enddo

            backspace(input%unit)

            elseif (line(1:5) == 'name:' .or.  &
                    line(1:5) == 'Name:' .or.  &
                    line(1:5) == 'NAME:' ) then

               molecule%name = trim(line(6:100))
               molecule%name = remove_preceding_blanks(molecule%name)

            elseif (line(1:7) == 'charge:' .or.  &
                    line(1:7) == 'Charge:' .or.  &
                    line(1:7) == 'CHARGE:' ) then

               read(line(8:100),*) molecule%charge

         endif

         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)

      enddo

      call disk%close_file(input)
!
   end subroutine read_info_molecular_system
!
!
   subroutine read_geometry_molecular_system(molecule)
!!
!!    Read
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Read atoms and their coordinates, assumed to be in units of Ångstrøm.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      character(len=100) :: line
      character(len=100) :: current_basis
!
      integer(i15) :: i = 0, current_atom = 0
!
      integer(i15) :: cursor
      character(len=100) :: coordinate
!
      call disk%open_file(input, 'read')
      rewind(input%unit)
!
      read(input%unit,'(a)') line
      line = remove_preceding_blanks(line)
!
      current_atom = 0
!
      do while (trim(line) .ne. 'end geometry')

         if (line(1:6) == 'basis:' .or.  &
             line(1:6) == 'Basis:' .or.  &
             line(1:6) == 'BASIS:' ) then
!
            current_basis = trim(line(7:100))
            current_basis = remove_preceding_blanks(current_basis)
!
            read(input%unit,'(a)') line
            line = remove_preceding_blanks(line)


           do while (trim(line) .ne. 'end geometry'  .and.  &
                      line(1:6) .ne. 'basis:'        .and.  &
                      line(1:6) .ne. 'Basis:'        .and.  &
                      line(1:6) .ne. 'BASIS:')
!
               current_atom = current_atom + 1
!
               molecule%atoms(current_atom, 1)%basis = current_basis
               molecule%atoms(current_atom, 1)%symbol = trim(line(1:2))
!
               line = line(3:100)
               line = remove_preceding_blanks(line)
!
               cursor = 1
!
               do
                  if (line(cursor:cursor) .eq. ' ') then
                     exit
                  else
                     cursor = cursor + 1
                     cycle
                  endif
               enddo
!
               coordinate = line(1:cursor)
               read(coordinate, '(f30.25)') molecule%atoms(current_atom,1)%x
!
               cursor = cursor + 1
!
               line = line(cursor:100)
               line = remove_preceding_blanks(line)
!
               cursor = 1 ! Initial value
!
               do
                  if (line(cursor:cursor) .eq. ' ') then
                     exit
                  else
                     cursor = cursor + 1
                     cycle
                  endif
               enddo
!
               coordinate = line(1:cursor)
               read(coordinate, '(f30.25)') molecule%atoms(current_atom,1)%y
!
               cursor = cursor + 1
!
               coordinate = line(cursor:100)
               coordinate = remove_preceding_blanks(coordinate)
!
               read(coordinate, '(f30.25)') molecule%atoms(current_atom,1)%z
!
              read(input%unit,'(a)') line
              line = remove_preceding_blanks(line)

           enddo
!
           backspace(input%unit)
!
         endif

         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
!
      enddo

      call disk%close_file(input)
!
   end subroutine read_geometry_molecular_system
!
!
   subroutine write_molecular_system(molecule)
!!
!!    Write
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Writes to disk an xyz file, which is used by the libint package
!!    to calculate integrals.
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(i15) :: atom = 0
!
      type(file) :: mol_file
!
      call mol_file%init(trim(molecule%name) // '.xyz', 'sequential', 'formatted')
      call disk%open_file(mol_file, 'write', 'rewind')
!
      write(mol_file%unit, '(i3/)') molecule%n_atoms
!
      do atom = 1, molecule%n_atoms
!
         write(mol_file%unit, '(a3, 3x, f30.25, 3x, f30.25, 3x, f30.25)')  &
                                 molecule%atoms(atom, 1)%symbol,           &
                                 molecule%atoms(atom, 1)%x,                &
                                 molecule%atoms(atom, 1)%y,                &
                                 molecule%atoms(atom, 1)%z
!
      enddo
!
      call disk%close_file(mol_file)
!
   end subroutine write_molecular_system
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
      integer(i15) :: i = 0, j = 0
!
      real(dp) :: x_ij, y_ij, z_ij, r_ij
!
      get_nuclear_repulsion_molecular_system = zero
!
      do i = 1, molecule%n_atoms
         do j = i + 1, molecule%n_atoms
!
            x_ij = molecule%atoms(i, 1)%x - molecule%atoms(j, 1)%x
            y_ij = molecule%atoms(i, 1)%y - molecule%atoms(j, 1)%y
            z_ij = molecule%atoms(i, 1)%z - molecule%atoms(j, 1)%z
!
            r_ij = sqrt(x_ij**2 + y_ij**2 + z_ij**2)
!
            r_ij = angstrom_to_bohr*r_ij
!
            if (abs(r_ij) .lt. 1.0D-7) then
!
               write(output%unit,'(/t3,a)') 'Error: Two atoms are placed on top of each other'
               stop
!
            endif
!
            get_nuclear_repulsion_molecular_system = get_nuclear_repulsion_molecular_system &
                  + ((molecule%atoms(i, 1)%number)*(molecule%atoms(j, 1)%number))/r_ij
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
      integer(i15) :: get_n_electrons_molecular_system
!
      integer(i15) :: i = 0
!
      get_n_electrons_molecular_system = 0
!
      do i = 1, molecule%n_atoms
!
         get_n_electrons_molecular_system = get_n_electrons_molecular_system + molecule%atoms(i,1)%number
!
      enddo
!
      get_n_electrons_molecular_system = get_n_electrons_molecular_system - molecule%charge
!
  end function get_n_electrons_molecular_system
!
!
   integer(i15) function get_n_shells_molecular_system(molecule)
!!
!!    Get number of shells
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(i15) :: I = 0
!
      get_n_shells_molecular_system = 0
!
      do I = 1, molecule%n_atoms
!
         get_n_shells_molecular_system = get_n_shells_molecular_system + molecule%atoms(I,1)%n_shells
!
      enddo
!
   end function get_n_shells_molecular_system
!
!
   integer(i15) function get_n_aos_molecular_system(molecule)
!!
!!    Get number of AOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(i15) :: I = 0, J = 0
!
      get_n_aos_molecular_system = 0
!
      do I = 1, molecule%n_atoms
!
         do J = 1, molecule%atoms(I,1)%n_shells
!
            get_n_aos_molecular_system = get_n_aos_molecular_system &
                           + molecule%atoms(I,1)%shells(J,1)%size
!
         enddo
!
      enddo
!
   end function get_n_aos_molecular_system
!
!
   type(interval) function get_shell_limits_molecular_system(molecule, A)
!!
!!    Get shell limits
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(kind=8) :: A
!
      integer(kind=8) :: I, J
!
      do I = 1, molecule%n_atoms
!
         do J = 1, molecule%atoms(I,1)%n_shells
!
            if (A .eq. molecule%atoms(I,1)%shells(J,1)%number) then
!
               get_shell_limits_molecular_system%first = molecule%atoms(I,1)%shells(J,1)%first
               get_shell_limits_molecular_system%last  = molecule%atoms(I,1)%shells(J,1)%last
               get_shell_limits_molecular_system%size  = molecule%atoms(I,1)%shells(J,1)%size
!
            endif
!
         enddo
!
      enddo
!
   end function get_shell_limits_molecular_system
!
!
   integer(i15) function basis2shell_molecular_system(molecule, basis_function)
!!
!!    Basis2shell
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(i15), intent(in) :: basis_function
!
      integer(i15) :: I, J
!
      do I = 1, molecule%n_atoms
         do J = 1, molecule%atoms(I,1)%n_shells
!
            if (molecule%atoms(I,1)%shells(J,1)%last  .ge. basis_function .and. &
                molecule%atoms(I,1)%shells(J,1)%first .le. basis_function) then
!
               basis2shell_molecular_system = molecule%atoms(I,1)%shells(J,1)%number
!
            endif
!
         enddo
      enddo
!
   end function basis2shell_molecular_system
!
!
end module molecular_system_class
