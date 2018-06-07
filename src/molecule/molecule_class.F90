module molecule_class
!
!!
!!    Molecule class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
   use atom_class
   use io_utilities
!
   implicit none
!
   type :: molecule
!
      character(len=40) :: name
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
      procedure, private :: set => set_molecule
      procedure, private :: read => read_molecule
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
      call mol%set
!
      allocate(mol%atoms(mol%n_atoms, 1))
!
      call mol%read
!
   end subroutine initialize_molecule
!
!
   subroutine write_molecule(mol)
!!
!!    Write
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Saves to file an xyz file with the current geometry
!!
      implicit none
!
      class(molecule) :: mol
!
   end subroutine write_molecule
!
!
   subroutine set_molecule(mol)
!!
!!    Set ....
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets number of atoms, charge and name of system  
!!
      implicit none
!
      class(molecule) :: mol
!
      character(len=40) :: line
      character(len=40) :: current_basis
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
      
               mol%name = trim(line(6:40))
               mol%name = remove_preceding_blanks(mol%name)
      
            elseif (line(1:7) == 'charge:' .or.  &
                    line(1:7) == 'Charge:' .or.  &
                    line(1:7) == 'CHARGE:' ) then
      
               read(line(8:40),*) mol%charge
      
         endif
       
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
      
      enddo
      
      call disk%close_file(input)
!
   end subroutine set_molecule
!
!
   subroutine read_molecule(mol)
!!
!!    Read
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Read atoms 
!!
      implicit none
!
      class(molecule) :: mol
!
      character(len=40) :: line
      character(len=40) :: current_basis
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
            current_basis = trim(line(7:40))
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
               mol%atoms(current_atom, 1)%type = trim(line(1:2))
               read(line(3:40), *) mol%atoms(current_atom, 1)%xyz
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
   end subroutine read_molecule
!
!
end module molecule_class
!call molecule%init('molecule.inp', 'sequential', 'formatted')
!      call disk%open_file(molecule, 'write', 'rewind')
               