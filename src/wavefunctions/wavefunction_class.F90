module wavefunction_class
!
!!
!!    Wavefunction class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
   use molecular_system_class
!
   implicit none
!
   type, abstract :: wavefunction
!
      character(len=40) :: name
!
      integer(i15) :: n_ao
      integer(i15) :: n_mo
      integer(i15) :: n_o
      integer(i15) :: n_v
!
      type(molecular_system) :: system
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients
!
   contains
!
!
   end type wavefunction
!
!
contains
!
!
end module wavefunction_class
