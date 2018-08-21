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
   use ao_integral_manager_class
   use active_atoms_info_class
!
   implicit none
!
   type :: molecular_system
!
      character(len=100) :: name
      character(len=40), dimension(:), allocatable :: basis_sets
!
      integer(i15) :: n_atoms
      integer(i15) :: n_basis_sets 
      integer(i15) :: charge
!
      type(atomic), dimension(:), allocatable :: atoms
!
      type(ao_integral_manager) :: ao_integrals
!
      type(interval), dimension(:), allocatable :: shell_limits 
!
      logical :: active_atoms = .false.
!
      integer(i15) :: n_active_atoms = 0
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
      procedure :: SAD => SAD_molecular_system
!
      procedure :: shell_to_atom => shell_to_atom_molecular_system
!
      procedure :: reorder_atoms => reorder_atoms_molecular_system
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
      character(len=100) :: temp_name
!
      integer(kind=4) :: i = 0, j = 0, n_atoms_libint
!
      integer(i15) :: s, n_s 
!
      integer(kind=4), dimension(:,:), allocatable :: n_shells_on_atoms
      integer(kind=4), dimension(:,:), allocatable :: n_basis_in_shells
      integer(kind=4), dimension(:,:), allocatable :: first_ao_in_shells
      integer(kind=4), dimension(:,:), allocatable :: shell_numbers
!
      call molecule%read_info()
!
      allocate(molecule%atoms(molecule%n_atoms))
!
      allocate(n_shells_on_atoms(molecule%n_atoms, 1))
      n_shells_on_atoms = 0
!
      call molecule%read_geometry()
!
      if (requested_section('active atoms')) then
!
         call molecule%reorder_atoms()
!
      endif
!
      call molecule%write()
!
      do i = 1, molecule%n_atoms
!
         call molecule%atoms(i)%set_number()
!
      enddo
!
      call initialize_atoms(molecule%name)
!
      do i = 1, molecule%n_basis_sets ! Loop over atoms 
         write(temp_name, '(a, a1, i4.4)')trim(molecule%name), '_', i
         call initialize_basis(molecule%basis_sets(i), temp_name) 
      enddo
!
      call get_n_shells_on_atoms(n_shells_on_atoms)
!
      do i = 1, molecule%n_atoms ! Loop over atoms
!
!        Allocate and initialize the corresponding shells
!
         molecule%atoms(i)%n_shells = n_shells_on_atoms(i, 1)
!
         allocate(molecule%atoms(i)%shells(molecule%atoms(i)%n_shells))
!
!        Then determine the number of basis functions in each shell
!        and save number of aos per atom

!
         allocate(n_basis_in_shells(n_shells_on_atoms(i,1), 1))
         call get_n_basis_in_shells(i, n_basis_in_shells)
!
         molecule%atoms(i)%n_ao = 0
!
         do j = 1, n_shells_on_atoms(i, 1)

            molecule%atoms(i)%shells(j)%size = n_basis_in_shells(j,1)
!
            molecule%atoms(i)%n_ao = molecule%atoms(i)%n_ao + n_basis_in_shells(j,1)

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
            molecule%atoms(i)%shells(j)%number = shell_numbers(j, 1)
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

            molecule%atoms(i)%shells(j)%first = first_ao_in_shells(j,1)
!
         enddo
!
         deallocate(first_ao_in_shells)
!
!        Then determine the angular momentum of shells & the last AO index
!
         do j = 1, n_shells_on_atoms(i,1)
!
            call molecule%atoms(i)%shells(j)%determine_angular_momentum()
            call molecule%atoms(i)%shells(j)%determine_last_ao_index()
!
         enddo
!
      enddo
!
      if (molecule%charge .ne. 0) then
!
         write(output%unit) 'Error: SAD not yet implemented for charged species!'
         stop
!
      endif
!
      do i = 1, molecule%n_atoms
!
         if (molecule%atoms(i)%n_ao == 0) then
!   
            write(output%unit,*)'Error: Is basis defined for all atoms in input?'
            stop
!
         endif
!
      enddo
!
!     Allocate and set shell limits vector 
!
      n_s = molecule%get_n_shells()
      allocate(molecule%shell_limits(n_s))
      do s = 1, n_s 
!  
         molecule%shell_limits(s) = molecule%get_shell_limits(s)
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
!!    Determine if we have active space 
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
      rewind(input%unit)
!
      molecule%n_atoms = 0
!
      read(input%unit,'(a)') line
      line = remove_preceding_blanks(line)
!
      molecule%n_basis_sets = 0
!
      do while (trim(line) .ne. 'end geometry')

         if (line(1:6) == 'basis:' .or.  &
             line(1:6) == 'Basis:' .or.  &
             line(1:6) == 'BASIS:' ) then
!
            molecule%n_basis_sets = molecule%n_basis_sets + 1
!
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
!
         elseif (line(1:12) == 'active atoms') then
!
            molecule%active_atoms = .true.
!
         endif

         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)

      enddo
!
      if ((molecule%n_basis_sets .le. 0 ).or.( molecule%n_basis_sets .gt. molecule%n_atoms)) then
!
         write(output%unit, '(a)')'Error: Number of basis sets specified exceeds number of atoms or is zero.'
         stop
!
      endif
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
      character(len=100) :: current_basis, temp_name
!
      integer(i15) :: i = 0, current_atom = 0, current_basis_nbr = 0
      integer(i15), dimension(:,:), allocatable :: atoms_with_current_basis
!
      integer(i15) :: cursor
      character(len=100) :: coordinate
!
      type(file) :: basis_file, mol_file
!
    !  call mem%alloc_int(atoms_with_current_basis, molecule%n_basis_sets, 1)
!
     ! call mol_file%init(trim(molecule%name) // '.xyz', 'sequential', 'formatted')
     ! call disk%open_file(mol_file, 'write', 'rewind')
!
     ! write(mol_file%unit, '(i4/)') molecule%n_atoms
!
      rewind(input%unit)
!
      read(input%unit,'(a)') line
      line = remove_preceding_blanks(line)
!
      current_atom = 0
      current_basis_nbr = 0
!
      do while (trim(line) .ne. 'end geometry')

         if (line(1:6) == 'basis:' .or.  &
             line(1:6) == 'Basis:' .or.  &
             line(1:6) == 'BASIS:' ) then
!     
            current_basis_nbr = current_basis_nbr + 1
           ! atoms_with_current_basis(current_basis_nbr , 1) = 0
!
            current_basis = trim(line(7:100))
            current_basis = remove_preceding_blanks(current_basis)
           ! molecule%basis_sets(current_basis_nbr) = current_basis    
!
            read(input%unit,'(a)') line
            line = remove_preceding_blanks(line)


           do while (trim(line) .ne. 'end geometry'  .and.  &
                      line(1:6) .ne. 'basis:'        .and.  &
                      line(1:6) .ne. 'Basis:'        .and.  &
                      line(1:6) .ne. 'BASIS:')
!
              ! write(mol_file%unit, '(a)') line 
              ! atoms_with_current_basis(current_basis_nbr , 1) = atoms_with_current_basis(current_basis_nbr , 1) + 1
!
               current_atom = current_atom + 1
!
               molecule%atoms(current_atom)%basis = current_basis
               molecule%atoms(current_atom)%symbol = trim(line(1:2))
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
               read(coordinate, '(f21.16)') molecule%atoms(current_atom)%x
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
               read(coordinate, '(f21.16)') molecule%atoms(current_atom)%y
!
               cursor = cursor + 1
!
               coordinate = line(cursor:100)
               coordinate = remove_preceding_blanks(coordinate)
!
               read(coordinate, '(f21.16)') molecule%atoms(current_atom)%z
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
!
      !call disk%close_file(mol_file)
!
      !call disk%open_file(mol_file, 'read')
      !rewind(mol_file%unit)
     ! read(mol_file%unit, '(a)') line
     ! read(mol_file%unit, '(a)') line
!
     ! do current_basis_nbr = 1, molecule%n_basis_sets
!
     !    write(temp_name, '(a, a1, i4.4, a4)')trim(molecule%name), '_', current_basis_nbr,  '.xyz'
!
     !    call basis_file%init(trim(temp_name), 'sequential', 'formatted')
    !     call disk%open_file(basis_file, 'write', 'rewind')
!
     !    write(basis_file%unit, '(i4/)') atoms_with_current_basis(current_basis_nbr, 1)
     !    do i = 1, atoms_with_current_basis(current_basis_nbr, 1)
     !       read(mol_file%unit, '(a)') line 
     !       write(basis_file%unit, '(a)') line 
    !    enddo
!
    !     call disk%close_file(basis_file)
!
   !   enddo
!
   !   call mem%dealloc_int(atoms_with_current_basis, molecule%n_basis_sets, 1)
!
    !  call disk%close_file(mol_file)
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
      character(len=100) temp_name
      character(len=100) current_basis
!
      type(file) :: mol_file, basis_file
!
      integer(i15) :: basis_set_counter, atom_offset, current_basis_nbr, i
!
      integer(i15), dimension(:,:), allocatable :: n_atoms_in_basis
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
         if (molecule%atoms(1)%basis .ne. current_basis) then
!
            current_basis = molecule%atoms(atom)%basis
            molecule%n_basis_sets = molecule%n_basis_sets + 1
!
         endif
!
      enddo
!
      allocate(molecule%basis_sets(molecule%n_basis_sets))
      call mem%alloc_int(n_atoms_in_basis, molecule%n_basis_sets, 1)
!
      n_atoms_in_basis = 1
      basis_set_counter = 1
      molecule%basis_sets(basis_set_counter) = molecule%atoms(1)%basis
!
!     Count number of atoms in each basis
!
      do atom = 2, molecule%n_atoms
!
         if (molecule%atoms(1)%basis .ne. molecule%basis_sets(basis_set_counter)) then
!
            basis_set_counter = basis_set_counter + 1
            molecule%basis_sets(basis_set_counter) = molecule%atoms(atom)%basis
!
         else
!
            n_atoms_in_basis(basis_set_counter, 1) = n_atoms_in_basis(basis_set_counter, 1) + 1
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
         write(basis_file%unit, '(i5/)') n_atoms_in_basis(current_basis_nbr, 1)
!
         do i = 1, n_atoms_in_basis(current_basis_nbr, 1)
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
         atom_offset = atom_offset + n_atoms_in_basis(current_basis_nbr, 1)
!
      enddo
!
      call mem%dealloc_int(n_atoms_in_basis, molecule%n_basis_sets, 1)
!
   end subroutine write_molecular_system
!
!
   subroutine reorder_atoms_molecular_system(molecule)
!!
!!    Reorder atoms
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Reorder atoms in case of active atoms
!!
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      character(len=100) :: line
!
      integer(i15) :: i, j, active_atom_counter, ioerror = 0, first, last, atom_counter
      integer(i15) :: central_atom
!
      integer(i15), dimension(:,:), allocatable :: active_atoms
!
      type(atomic), dimension(:), allocatable :: atoms_copy
!
      logical :: found
!
      real(dp) :: hf_radius, x, y, z
!
      rewind(input%unit)
!
      read(input%unit,'(a)', iostat=ioerror) line
      line = remove_preceding_blanks(line)
!
      do while (trim(line) .ne. 'geometry')
!
         if (trim(line) == 'active atoms') then
!
            do while (trim(line) .ne. 'end active atoms')
!
               read(input%unit,'(a)', iostat=ioerror) line
               line = remove_preceding_blanks(line)
!
               if (line(1:8) == 'hf list:') then
!
                  line = line(9:100)
                  line = remove_preceding_blanks(line)
                  molecule%n_active_atoms = 0
!
                  do i = 1, 92
!
                     if (line(i:i) .ne. ' ') molecule%n_active_atoms = molecule%n_active_atoms + 1
!
                  enddo
!
                  call mem%alloc_int(active_atoms, molecule%n_active_atoms, 1)
                  read(line, *) active_atoms
                  exit
!
               elseif (line(1:9) == 'hf range:') then 
!
                  line = line(10:100)
                  line = remove_preceding_blanks(line)
!
                  if (line(1:1)=='[') then ! range given
!
                     do i = 2, 100
!
                        if (line(i:i) == ',') exit
!
                     enddo
!
                     read(line(2:i-1), *) first
!
                     do j = i, 100
!
                        if (line(j:j) == ']') exit
!
                     enddo
!
                     read(line(i+1:j-1), *) last
!
                     molecule%n_active_atoms = last - first + 1
!
                     call mem%alloc_int(active_atoms, molecule%n_active_atoms, 1)
!
                     do i = first, last
!
                        active_atoms(i - first + 1, 1) = i
!
                     enddo
!
                     exit
!
                  else 
!
                     call output%error_msg('active atom range not detected.')
!
                  endif
!
               elseif (line(1:13) == 'central atom:') then
!
                  read(line(14:100), *) central_atom
!
                  do while (trim(line) .ne. 'end active atoms')
!
                     read(input%unit,'(a)') line
                     line = remove_preceding_blanks(line)
!
                     if (line(1:10) == 'hf radius:') then
!
                        line = line(11:100)
                        line = remove_preceding_blanks(line)
                        read(line, *) hf_radius ! In Ångstom
!
                     endif
!
                     molecule%n_active_atoms = 0
!
                     do i = 1, molecule%n_atoms 
!
                        x = (molecule%atoms(central_atom)%x - molecule%atoms(i)%x)
                        y = (molecule%atoms(central_atom)%y - molecule%atoms(i)%y)
                        z = (molecule%atoms(central_atom)%z - molecule%atoms(i)%z)
!
                        if (sqrt(x**2 + y**2 + z**2) .le. hf_radius) molecule%n_active_atoms = molecule%n_active_atoms + 1
!
                     enddo
!
                     call mem%alloc_int(active_atoms, molecule%n_active_atoms, 1)
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
                           active_atoms(active_atom_counter, 1) = i
!
                        endif
!
                     enddo
!
                  enddo
                  exit
!
               else
!
                  call output%error_msg('active atom input not recognized.')
!
               endif
            enddo
         endif
!
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
!
      enddo
!
!     Reorder atoms
!
      allocate(atoms_copy(molecule%n_atoms))
      atoms_copy = molecule%atoms
!
      do i = 1, molecule%n_active_atoms
!
         molecule%atoms(i) = atoms_copy(active_atoms(i, 1))
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
            if (i == active_atoms(j, 1)) then
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
   end subroutine reorder_atoms_molecular_system
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
               write(output%unit,'(/t3,a)') 'Error: Two atoms are placed on top of each other'
               stop
!
            endif
!
            get_nuclear_repulsion_molecular_system = get_nuclear_repulsion_molecular_system &
                  + ((molecule%atoms(i)%number)*(molecule%atoms(j)%number))/r_ij
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
         get_n_electrons_molecular_system = get_n_electrons_molecular_system + molecule%atoms(i)%number
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
         get_n_shells_molecular_system = get_n_shells_molecular_system + molecule%atoms(I)%n_shells
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
   type(interval) function get_shell_limits_molecular_system(molecule, A)
!!
!!    Get shell limits
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(kind=8), intent(in) :: A
!
      integer(kind=8) :: I, J
!
      do I = 1, molecule%n_atoms
!
         do J = 1, molecule%atoms(I)%n_shells
!
            if (A .eq. molecule%atoms(I)%shells(J)%number) then
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
         do J = 1, molecule%atoms(I)%n_shells
!
            if (molecule%atoms(I)%shells(J)%last  .ge. basis_function .and. &
                molecule%atoms(I)%shells(J)%first .le. basis_function) then
!
               basis2shell_molecular_system = molecule%atoms(I)%shells(J)%number
!
            endif
!
         enddo
      enddo
!
   end function basis2shell_molecular_system
!
!
   subroutine SAD_molecular_system(molecule, n_ao, density_diagonal)
!!
!!    Superposition of atomic desities
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initial guess for HF-calculation
!!
      implicit none
!
      class(molecular_system) :: molecule
!
      integer(i15) :: n_ao
!
      real(dp), dimension(n_ao, 1) :: density_diagonal
!
      integer(i15) :: I, offset_diagonal
!
      real(dp) :: electrons
!
      real(dp), dimension(:,:), allocatable :: atom_density_diagonal
!
!     Loop over atoms and let them set their own density diagonal
!
      offset_diagonal = 0
!
      do I = 1, molecule%n_atoms
!
            call mem%alloc(atom_density_diagonal, molecule%atoms(I)%n_ao, 1)
!
            call molecule%atoms(I)%AD(atom_density_diagonal)
!
            density_diagonal(offset_diagonal + 1 : offset_diagonal + molecule%atoms(I)%n_ao, 1) = &
               atom_density_diagonal(:,1)
!
            call mem%dealloc(atom_density_diagonal, molecule%atoms(I)%n_ao, 1)
!
            offset_diagonal = offset_diagonal + molecule%atoms(I)%n_ao
!
!
      enddo
!
      electrons = 0
!
      do I = 1, n_ao
!
         electrons = electrons + density_diagonal(I, 1)
!
      enddo
!
      if (abs(electrons - molecule%get_n_electrons()) .gt. 1.0d-7) then
!
         write(output%unit, '(a)') 'Error: Mismatch in electron number SAD'
!
      endif
!
   end subroutine SAD_molecular_system
!
!
   integer function shell_to_atom_molecular_system(molecule, shell)
!!
!!
!!
      implicit none
!
      class(molecular_system) :: molecule
      integer(i15), intent(in) :: shell
!
      integer(i15) :: I, accumulated_shells
!
      if (shell .le. 0) then
!
         write(output%unit, '(a)') 'Error: shell number has illegal value 0.'
         stop
!
      elseif (shell .gt. molecule%get_n_shells()) then
!
         write(output%unit, '(a)') 'Error: shell number exceeds total number of shells.'
         stop
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
   end function shell_to_atom_molecular_system
!
!
end module molecular_system_class
