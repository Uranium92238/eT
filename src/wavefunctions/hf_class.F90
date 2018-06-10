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
      real(dp) :: n_ao
      real(dp) :: n_mo
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
      call wf%molecule%initialize
      call wf%molecule%write
!
      call wf%integrals%get_ao_xy
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
