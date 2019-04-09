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
   module subroutine construct_ao_h_wx_molecular_system(molecule, h, s1, s2)
!!
!!    Construct h_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      integer, intent(in) :: s1, s2 
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h 
!
   end subroutine construct_ao_h_wx_molecular_system
!
!
   module subroutine construct_ao_s_wx_molecular_system(molecule, s, s1, s2)
!!
!!    Construct s_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s1)%size), intent(inout) :: s
!
   end subroutine construct_ao_s_wx_molecular_system
!
!
   module subroutine construct_ao_mu_wx_molecular_system(molecule, mu_X, mu_Y, mu_Z, s1, s2)
!!
!!    Construct μ_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: mu_X ! x component
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: mu_Y ! y component 
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: mu_Z ! z component
!
      integer(i6) :: s1_4, s2_4
!
   end subroutine construct_ao_mu_wx_molecular_system
!
!
   module subroutine construct_ao_q_wx_molecular_system(molecule, q_xx, q_xy, q_xz, q_yy, q_yz, q_zz, s1, s2)
!!
!!    Construct q_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: q_xx 
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: q_xy 
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: q_xz 
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: q_yy 
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: q_yz 
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(inout) :: q_zz 
!
   end subroutine construct_ao_q_wx_molecular_system