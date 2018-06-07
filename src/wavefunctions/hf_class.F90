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
      type(molecule) :: molecule
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
      call wf%molecule%initialize()
      call wf%molecule%write()
!
      write(output%unit, '(/t6,a,a,a)') &
      'A wavefunction of type ', trim(wf%name), ' was initialized!'
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
      write(output%unit, '(/t6,a,a,a)') &
         'A wavefunction of type ', trim(wf%name), ' was finalized!'
!
   end subroutine finalize_hf
!
!
end module hf_class
