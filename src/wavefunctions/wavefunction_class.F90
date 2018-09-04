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
      procedure :: initialize_mo_coefficients => initialize_mo_coefficients_wavefunction
      procedure :: destruct_mo_coefficients => destruct_mo_coefficients_wavefunction
!
   end type wavefunction
!
contains
!
!
   subroutine initialize_mo_coefficients_wavefunction(wf)
!!
!!    Initialize MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_coefficients)) call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_mo_coefficients_wavefunction
!
!
   subroutine destruct_mo_coefficients_wavefunction(wf)
!!
!!    Destruct MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_coefficients)) call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_mo_coefficients_wavefunction
!
end module wavefunction_class
