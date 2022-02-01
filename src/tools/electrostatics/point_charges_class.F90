!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WiTHOUT ANY WArrANTY; without even the implied warranty of
!  MErCHANTABiLiTY or FiTNESS FOr A PArTiCULAr PUrPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. if not, see <https://www.gnu.org/licenses/>.
!
!
module point_charges_class
!
!!
!!    Point charges class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Oct 2020
!!
!
   use kinds
   use parameters
   use memory_manager_class, only: mem
   use global_out,           only: output
!
   implicit none 
!
   type :: point_charges
!
      integer                                :: n_charges
      real(dp), dimension(:,:), allocatable  :: r
      real(dp), dimension(:),   allocatable  :: q
!
   contains 
!
      procedure :: initialize &
                => initialize_point_charges
!
      final :: destructor
!
      generic :: get_coulomb_interaction => get_coulomb_interaction_1_point_charges, &
                                            get_coulomb_interaction_2_point_charges
!
      procedure :: get_coulomb_interaction_1_point_charges
      procedure :: get_coulomb_interaction_2_point_charges
!
      procedure :: get_dipole &
                => get_dipole_point_charges
!
      procedure :: get_quadrupole &
                => get_quadrupole_point_charges
!
      procedure :: get_coulomb_interaction_1der &
                => get_coulomb_interaction_1der_point_charges 
!
      procedure :: get_potential_at_external_points &
                => get_potential_at_external_points_point_charges
!
   end type point_charges
!
!
   interface point_charges 
!
      procedure :: new_point_charges 
!
   end interface point_charges
!
contains  
!
!
   pure function new_point_charges(n_charges) result(this)
!!
!!    New point charges
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      integer, intent(in) :: n_charges
      type(point_charges) :: this
!
      this%n_charges = n_charges
!
   end function new_point_charges
!
!
   subroutine initialize_point_charges(this)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(point_charges) :: this
!
      call mem%alloc(this%r, 3, this%n_charges)
      call mem%alloc(this%q, this%n_charges)
!
   end subroutine initialize_point_charges
!
!
   subroutine destructor(this)
!!
!!    Destructor
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      type(point_charges) :: this
!
      if (allocated(this%r)) call mem%dealloc(this%r, 3, this%n_charges)
      if (allocated(this%q)) call mem%dealloc(this%q, this%n_charges)
!
   end subroutine destructor
!
!
   function get_coulomb_interaction_1_point_charges(this) result(E)
!!
!!    Get coulomb interaction
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2020
!!
!!    Calculates 
!!
!!       E = sum_(i != j)  1/2 q_i q_j / |r_i - r_j|
!!
      implicit none
!
      class(point_charges), intent(in) :: this
!
      real(dp)                         :: E
!
      real(dp), dimension(3)           :: r_ij
      real(dp)                         :: abs_r_ij
      integer                          :: i, j
!
      E = zero
!
      do i = 1, this%n_charges
         do j = i + 1, this%n_charges
!
            r_ij = this%r(:,i) - this%r(:, j) 
!
            abs_r_ij = sqrt(r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2)*angstrom_to_bohr
!
            if (abs_r_ij .lt. 1.0D-7) &
               call output%error_msg('two charges are placed on top of each other.')
!
            E = E + this%q(i)*this%q(j)/abs_r_ij
! 
         enddo
      enddo
!
   end function get_coulomb_interaction_1_point_charges
!
!
   function get_coulomb_interaction_2_point_charges(this, that) result(E)
!!
!!    Get coulomb interaction
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2020
!!
!!    Calculate 
!!
!!       E = sum_(i != j)  1/2 q_i q_j / |r_i - r_j|
!!
!!    for two different sets of point chares   
!!
      implicit none
!
      class(point_charges), intent(in) :: this
      class(point_charges), intent(in) :: that
!
      real(dp)                         :: E
!
      real(dp), dimension(3)           :: r_ij
      real(dp)                         :: abs_r_ij
      integer                          :: i, j
!
      E = zero
!
      do i = 1, this%n_charges
         do j = 1, that%n_charges
!
            r_ij = this%r(:,i) - that%r(:, j) 
!
            abs_r_ij = sqrt(r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2)*angstrom_to_bohr
!
            if (abs_r_ij .lt. 1.0D-7) &
               call output%error_msg('two charges are placed on top of each other.')
!
            E = E + (this%q(i))*(that%q(j))/abs_r_ij
!
         enddo
      enddo
!
   end function get_coulomb_interaction_2_point_charges
!
!
   subroutine get_potential_at_external_points_point_charges(this, n_external_charges, r, potential)
!!
!!    Get potential at external points
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates and returns the potential on external charges from the nuclei
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(point_charges),                        intent(in)  :: this
      integer,                                     intent(in)  :: n_external_charges 
      real(dp), dimension(3, n_external_charges),  intent(in)  :: r
      real(dp), dimension(n_external_charges),     intent(out) :: potential
!
      integer                 :: i, p
      real(dp), dimension(3)  :: r_ip
      real(dp)                :: abs_r_ip
!
      call zero_array(potential, n_external_charges)
!
      do p = 1, n_external_charges
         do i = 1, this%n_charges
!
            r_ip = this%r(:,i) - r(:, p) 
!
            abs_r_ip = sqrt(r_ip(1)**2 + r_ip(2)**2 + r_ip(3)**2)*angstrom_to_bohr
!
            if (abs_r_ip .lt. 1.0D-7) then
!
               call output%error_msg('two atoms are placed on top of each other.')
!
            endif
!
            potential(p) = potential(p) + this%q(i)/abs_r_ip
!
         enddo
!
      enddo
!
   end subroutine get_potential_at_external_points_point_charges
!
!
   pure function get_dipole_point_charges(this) result(d)
!!
!!    Get dipole 
!!    Written by Eirik F. Kjønstad, Apr 2019 
!!
!!    Calculates the dipole moment for a set of point charges    
!!
      implicit none 
!
      class(point_charges), intent(in) :: this 
!
      real(dp), dimension(3) :: d 
!
      integer :: i 
!
      d = zero 
!
      do i = 1, this%n_charges 
!
         d(:) = d(:) + this%q(i) * this%r(:,i)
!
      enddo
!
      d = d * angstrom_to_bohr
!
   end function get_dipole_point_charges
!
!
   pure function get_quadrupole_point_charges(this) result(q)
!!
!!    Get quadrupole 
!!    Written by Eirik F. Kjønstad, Apr 2019 
!!
!!    Returns the quadrupole q ordered as xx, xy, xz, yy, yz, zz
!!
      implicit none 
!
      class(point_charges), intent(in) :: this 
!
      real(dp), dimension(6) :: Q 
!
      integer :: j
!
      Q = zero 
!
      do j = 1, this%n_charges  
!
         Q(1) = Q(1) + this%R(1, j) * this%R(1, j) * this%q(j) 
         Q(2) = Q(2) + this%R(1, j) * this%R(2, j) * this%q(j) 
         Q(3) = Q(3) + this%R(1, j) * this%R(3, j) * this%q(j)
         Q(4) = Q(4) + this%R(2, j) * this%R(2, j) * this%q(j) 
         Q(5) = Q(5) + this%R(2, j) * this%R(3, j) * this%q(j)
         Q(6) = Q(6) + this%R(3, j) * this%R(3, j) * this%q(j)
!
      enddo
!
      Q = Q * angstrom_to_bohr**2
!
   end function get_quadrupole_point_charges
!
!
   pure function get_coulomb_interaction_1der_point_charges(this) result(E)
!!
!!    Get coulomb interaction 1der
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2020
!!
!!    Calculates first derivative of Coulomb interaction energy 
!!
      implicit none
!
      class(point_charges), intent(in) :: this
!
      real(dp), dimension(3, this%n_charges) :: E
!
      real(dp), dimension(3) :: r_ij
!
      real(dp) :: r_ij_3, ZiZj
!
      integer :: i, j
!
      E = zero
!
      do i = 1, this%n_charges 
         do j = 1, i - 1
!
            r_ij = (this%r(:,i) - this%r(:,j)) * angstrom_to_bohr
!
            r_ij_3 = sqrt(r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2)**3
!
            ZiZj = this%q(i) * this%q(j)
!
            E(:, i) = E(:, i) - r_ij(:) * ZiZj / r_ij_3
            E(:, j) = E(:, j) + r_ij(:) * ZiZj / r_ij_3
!
         enddo
      enddo
!
   end function get_coulomb_interaction_1der_point_charges
!
!
end module point_charges_class
