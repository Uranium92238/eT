!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module mm_molecule_class
!
!!
!!    MM molecule class module
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2018 and 2020
!!
!
   use mm_atom_class, only: mm_atom
   use parameters
!
   implicit none
!
   type :: mm_molecule
!
      integer                                  :: n_atoms
      type(mm_atom), dimension(:), allocatable :: atoms
!
   contains
!
      procedure :: initialize_mm_atoms &
                => initialize_mm_atoms_mm_molecule
!
      procedure :: set_chi &
                => set_chi_mm_molecule
!
      procedure :: set_eta &
                => set_eta_mm_molecule
!
      procedure :: set_q &
                => set_q_mm_molecule
!
      procedure :: set_r &
                => set_r_mm_molecule
!
      procedure :: get_chi_i &
                => get_chi_i_mm_molecule
!
      procedure :: get_eta_i &
                => get_eta_i_mm_molecule
!
      procedure :: get_q_i &
                => get_q_i_mm_molecule
!
      procedure :: get_r_i &
                => get_r_i_mm_molecule
!
   end type  mm_molecule
!
!
   interface  mm_molecule
!
      procedure :: new_mm_molecule
!
   end interface  mm_molecule
!
!
contains
!
!
   pure function new_mm_molecule(n_atoms) result(molecule)
!!
!!    New mm_molecule 
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none 
!
      type(mm_molecule) :: molecule 
!
      integer, intent(in) :: n_atoms
!
      molecule%n_atoms = n_atoms           
!
   end function new_mm_molecule
!
!
   subroutine initialize_mm_atoms_mm_molecule(molecule, R, S)
!!
!!    Initialize MM atoms
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    R : positions
!!    S : symbol in periodic table
!!
      implicit none
!
      class(mm_molecule)   :: molecule
!
      real(dp), dimension(3, molecule%n_atoms),       intent(in) :: R
      character(len=2), dimension(molecule%n_atoms),  intent(in) :: S
!
      integer  :: I
!
      allocate(molecule%atoms(molecule%n_atoms))
!
      do I = 1, molecule%n_atoms 
!
         molecule%atoms(I) = mm_atom(R(:,I), S(I))
!  
      enddo   
!
   end subroutine initialize_mm_atoms_mm_molecule
!
!
   pure subroutine set_q_mm_molecule(molecule, q)
!!
!!    Set q
!!    Written by Sarai D. Folkestad
!!
!!    Sets charge
!!    
      implicit none
!  
      class(mm_molecule), intent(inout) :: molecule
!
      real(dp), dimension(molecule%n_atoms), intent(in) :: q
!
      integer :: I
!
      do I = 1, molecule%n_atoms 
!
         call molecule%atoms(I)%set_q(q(I))
!  
      enddo
!
   end subroutine set_q_mm_molecule
!
!
   pure subroutine set_r_mm_molecule(molecule, r)
!!
!!    Set r
!!    Written by Sarai D. Folkestad
!!
!!    Sets position
!!    
      implicit none
!  
      class(mm_molecule), intent(inout) :: molecule
!
      real(dp), dimension(3, molecule%n_atoms), intent(in) :: r
!
      integer :: I
!
      do I = 1, molecule%n_atoms 
!
         call molecule%atoms(I)%set_r(r(:,I))
!  
      enddo
!
   end subroutine set_r_mm_molecule
!
!
   pure subroutine set_chi_mm_molecule(molecule, chi)
!!
!!    Set chi
!!    Written by Sarai D. Folkestad
!!
!!    Sets electronegativity
!!    
      implicit none
!  
      class(mm_molecule), intent(inout) :: molecule
!
      real(dp), dimension(molecule%n_atoms), intent(in) :: chi
!
      integer :: I
!
      do I = 1, molecule%n_atoms 
!
         call molecule%atoms(I)%set_chi(chi(I))
!  
      enddo
!
   end subroutine set_chi_mm_molecule
!
!
   pure subroutine set_eta_mm_molecule(molecule, eta)
!!
!!    Set eta
!!    Written by Sarai D. Folkestad
!!
!!    Sets chemical hardness
!!    
      implicit none
!  
      class(mm_molecule), intent(inout) :: molecule
!
      real(dp), dimension(molecule%n_atoms), intent(in) :: eta
!
      integer :: I
!
      do I = 1, molecule%n_atoms 
!
         call molecule%atoms(I)%set_eta(eta(I))
!  
      enddo
!
   end subroutine set_eta_mm_molecule
!
!
   pure function get_q_i_mm_molecule(molecule, i) result(q)
!!
!!    Get q i
!!    Written by Sarai D. Folkestad
!!
!!    Get charge of atom i
!!    
      implicit none
!  
      class(mm_molecule), intent(in) :: molecule
!
      integer, intent(in) :: i
!
      real(dp) :: q
!
      q = molecule%atoms(I)%get_q()
!
   end function get_q_i_mm_molecule
!
!
   pure function get_r_i_mm_molecule(molecule, i) result(r)
!!
!!    Get r i
!!    Written by Sarai D. Folkestad
!!
!!    Get position of atom i
!!    
      implicit none
!  
      class(mm_molecule), intent(in) :: molecule
!
      integer, intent(in) :: i
!
      real(dp), dimension(3) :: r
!
      r = molecule%atoms(I)%get_r()
!  
   end function get_r_i_mm_molecule
!
!
   pure function get_chi_i_mm_molecule(molecule, i) result(chi)
!!
!!    Get chi
!!    Written by Sarai D. Folkestad
!!
!!    Get electronegativity of atom i
!!    
      implicit none
!  
      class(mm_molecule), intent(in) :: molecule
!
      integer, intent(in) :: i
!
      real(dp) :: chi
!
      chi = molecule%atoms(I)%get_chi()
!
   end function get_chi_i_mm_molecule
!
!
   pure function get_eta_i_mm_molecule(molecule, i) result(eta)
!!
!!    Get eta
!!    Written by Sarai D. Folkestad
!!
!!    Get chemical hardness of atom i
!!    
      implicit none
!  
      class(mm_molecule), intent(in) :: molecule
!
      integer, intent(in) :: i
!
      real(dp) :: eta
!
      eta = molecule%atoms(I)%get_eta()
!
   end function get_eta_i_mm_molecule
!
!
end module mm_molecule_class
