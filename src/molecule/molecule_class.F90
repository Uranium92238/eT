module molecule_class
!
!!
!!    Molecule class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
	use kinds
   use file_class
   use atom_class
!
	implicit none
!
	type :: molecule
!
      character(len=40) :: name
!
      integer(i15) :: n_atoms
!
      type(atom), dimension(:,:), allocatable :: atoms
!
   contains
!
      procedure :: initialize => initialize_molecule
      procedure :: write      => write_molecule
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
      type(file) :: input
!
!
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
end module molecule_class
