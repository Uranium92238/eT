!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
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
      character(len=40) :: name_
!
      real(dp) :: energy
!
      integer :: n_ao
      integer :: n_mo
      integer :: n_o
      integer :: n_v
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
      procedure :: print_wavefunction_summary      => print_wavefunction_summary_wavefunction
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
   subroutine print_wavefunction_summary_wavefunction(wf)
!!
!!    Print wavefunction summary 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Prints information related to the wavefunction,
!!    most of which is meaningful only for a properly 
!!    converged wavefunction. Should be overwritten in 
!!    descendants if more or less or other information 
!!    is present. 
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      write(output%unit, '(/t3,a,a,a)') '- Summary of ', trim(wf%name_), ' wavefunction:'
!
!
   end subroutine print_wavefunction_summary_wavefunction
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
