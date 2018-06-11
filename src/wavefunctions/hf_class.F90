module hf_class
!
!!
!!    Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use atom_class
   use integrals_class
   use molecule_class
   use disk_manager_class
!
   implicit none
!
   type :: hf
!
      character(len=40) :: name = 'HF'
!
      integer(i15) :: n_ao
      integer(i15) :: n_mo
!
      real(dp), dimension(:,:), allocatable :: mo_coefficients
!
      type(molecule)  :: molecule
      type(integrals) :: integrals
!
	contains
!
      procedure :: initialize => initialize_hf
      procedure :: finalize   => finalize_hf
!
   end type hf
!
!
contains
!
!
   subroutine initialize_hf(wf)
!!
!! 	Initialize
!!  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: repulsion
!
!     Initialize molecule
!
      call wf%molecule%initialize ! Read atoms
      call wf%molecule%write      ! Write an xyz file for the read geometry
!
!     Get the number of atomic orbitals
!
      wf%n_ao = 0
      call get_n_aos(wf%n_ao)
!
      call wf%integrals%get_ao_xy(wf%n_ao)
!
      repulsion = wf%molecule%nuclear_repulsion()
      write(output%unit,*) 'The nuclear repulsion:', repulsion
!
   end subroutine initialize_hf
!
!
   subroutine finalize_hf(wf)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine finalize_hf
!
!
end module hf_class
