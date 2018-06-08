module molecule_class
!
!!
!!    Molecule class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use atom_class
   use io_utilities
!
   implicit none
!
   type :: molecule
!
      character(len=100) :: name
!
      integer(i15) :: n_atoms
      integer(i15) :: charge
!
      type(atom), dimension(:,:), allocatable :: atoms
!
   contains
!
      procedure :: initialize => initialize_molecule
      procedure :: write      => write_molecule
!
      procedure, private :: read_info => read_info_molecule
      procedure, private :: read_geometry => read_geometry_molecule
!
   end type molecule
!
contains
!
!
   subroutine initialize_molecule(mol)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecule) :: mol
!
      integer(i15) :: i = 0
!
      call mol%read_info
!
      allocate(mol%atoms(mol%n_atoms, 1))
!
      call mol%read_geometry
!
      do i = 1, mol%n_atoms
!
        call mol%atoms(i, 1)%set_number
        write(output%unit, *)mol%atoms(i, 1)%number
!
      enddo
!
   end subroutine initialize_molecule
!
!
   subroutine read_info_molecule(mol)
!!
!!    Set
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets number of atoms, charge and name of system  
!!
      implicit none
!
      class(molecule) :: mol
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
      mol%n_atoms = 0
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
      
               mol%n_atoms = mol%n_atoms + 1
      
               read(input%unit,'(a)') line
               line = remove_preceding_blanks(line)
      
            enddo
      
            backspace(input%unit)
      
            elseif (line(1:5) == 'name:' .or.  &
                    line(1:5) == 'Name:' .or.  &
                    line(1:5) == 'NAME:' ) then
      
               mol%name = trim(line(6:100))
               mol%name = remove_preceding_blanks(mol%name)
      
            elseif (line(1:7) == 'charge:' .or.  &
                    line(1:7) == 'Charge:' .or.  &
                    line(1:7) == 'CHARGE:' ) then
      
               read(line(8:100),*) mol%charge
      
         endif
       
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
      
      enddo
      
      call disk%close_file(input)
!
   end subroutine read_info_molecule
!
!
   subroutine read_geometry_molecule(mol)
!!
!!    Read
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Read atoms and their coordinates
!!
      implicit none
!
      class(molecule) :: mol
!
      character(len=100) :: line
      character(len=100) :: current_basis
!
      integer(i15) :: i = 0, current_atom = 0
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
              mol%atoms(current_atom, 1)%basis = current_basis
              mol%atoms(current_atom, 1)%symbol = trim(line(1:2))
!   
              line = line(3:100)
              line = remove_preceding_blanks(line)     
              read(line, *) mol%atoms(current_atom,1)%x, &
                                    mol%atoms(current_atom,1)%y, &
                                    mol%atoms(current_atom,1)%z
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
   end subroutine read_geometry_molecule
!
!
   subroutine write_molecule(mol)
!!
!!    Write
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Write xyz file for libint 
!!
      implicit none
!
      class(molecule) :: mol
      integer(i15) :: atom = 0
!
      type(file) :: mol_file
!
      call mol_file%init(trim(mol%name) // '.xyz', 'sequential', 'formatted')
      call disk%open_file(mol_file, 'write', 'rewind')
!
      write(mol_file%unit, '(i3/)') mol%n_atoms
!
      do atom = 1, mol%n_atoms
         write(mol_file%unit, '(a3, 3x, f12.5, 3x, f12.5, 3x, f12.5)')  &
                                 mol%atoms(atom, 1)%symbol,             &   
                                 mol%atoms(atom, 1)%x,                  &
                                 mol%atoms(atom, 1)%y,                  &
                                 mol%atoms(atom, 1)%z
      enddo
!
      call disk%close_file(mol_file)
!
   end subroutine write_molecule
!
!
   function nuclear_repulsion_molecule(mol)
!!
!!    Nuclear repulsion
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates and returns nuclear repulsion for the molecule
!!
      implicit none
!
      class(molecule) :: mol
!
      real(dp) :: nuclear_repulsion_molecule
!
      integer(i15) :: i = 0, j = 0
!
      real(dp) :: x_ij, y_ij, z_ij, r_ij
!
      nuclear_repulsion_molecule = 0
!
      do i = 1, mol%n_atoms
         do j = i + 1, mol%n_atoms
!
            x_ij = mol%atoms(i, 1)%x - mol%atoms(j, 1)%x
            y_ij = mol%atoms(i, 1)%y - mol%atoms(j, 1)%y
            z_ij = mol%atoms(i, 1)%x - mol%atoms(j, 1)%z
!
            r_ij = sqrt(x_ij**2 + y_ij**2 +z_ij**2)
!
             if (abs(r_ij) .lt. 1.d-7) then
                write(output%unit,'(/t3,a)') 'Error: Two atoms are placed on top of each other'
                stop
             endif
!
             nuclear_repulsion_molecule = nuclear_repulsion_molecule &
                          + (mol%atoms(i, 1)%number)*(mol%atoms(j, 1)%number)
!
         enddo
      enddo
!
  end function nuclear_repulsion_molecule
!
!
end module molecule_class

               