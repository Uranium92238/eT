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
      real(dp) :: energy
!
      integer(i15) :: n_ao
      integer(i15) :: n_mo
      integer(i15) :: n_o
      integer(i15) :: n_v
!
      type(molecular_system) :: system
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients
      real(dp), dimension(:,:), allocatable :: orbital_energies
!
   contains
!
      procedure :: initialize_orbital_coefficients => initialize_orbital_coefficients_wavefunction
      procedure :: initialize_orbital_energies     => initialize_orbital_energies_wavefunction
!

      procedure :: destruct_orbital_coefficients   => destruct_orbital_coefficients_wavefunction
      procedure :: destruct_orbital_energies       => destruct_orbital_energies_wavefunction
!
   end type wavefunction
!
contains
!
!
   subroutine initialize_orbital_coefficients_wavefunction(wf)
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
   end subroutine initialize_orbital_coefficients_wavefunction
!
!
   subroutine destruct_orbital_coefficients_wavefunction(wf)
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
   end subroutine destruct_orbital_coefficients_wavefunction
!
!
   subroutine initialize_orbital_energies_wavefunction(wf)
!!
!!    Initialize orbital energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_energies)) call mem%alloc(wf%orbital_energies, wf%n_mo, 1)
!
   end subroutine initialize_orbital_energies_wavefunction
!
!
   subroutine destruct_orbital_energies_wavefunction(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_energies)) call mem%dealloc(wf%orbital_energies, wf%n_mo, 1)
!
   end subroutine destruct_orbital_energies_wavefunction
!
!
end module wavefunction_class
