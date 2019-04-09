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
submodule (molecular_system_class) ao_integrals
!
!!
!!    AO integrals submodule
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2019
!!
!!    Routines to handle AO integrals. Provides an interface to Libint to eT developers.
!!  
!!    Note: C++ is not fond of 64-bit integers, so ints are explicitly 
!!    translated to 32-bits here before calling Libint routines 
!!
!
   implicit none
!
   include "../../libint/h_wx_cdef.F90"
   include "../../libint/s_wx_cdef.F90"
   include "../../libint/mu_wx_cdef.F90"
   include "../../libint/q_wx_cdef.F90"
!
!
contains
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
      integer(i6) :: s1_4, s2_4 
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_h_wx_c(h, s1_4, s2_4)
!
   end subroutine construct_ao_h_wx_molecular_system
!
!
   module subroutine construct_ao_s_wx_molecular_system(molecule, s, s1, s2)
!!
!!    Construct s_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the s_αβ integral in the array s. s1 and s2 are the shells
!!    that w and x, respectively belong to.
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(molecule%shell_limits(s1)%size,molecule%shell_limits(s1)%size), intent(inout) :: s
!
      integer(i6) :: s1_4, s2_4 ! Integers that are passed to libint
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_s_wx_c(s, s1_4, s2_4) 
!
   end subroutine construct_ao_s_wx_molecular_system
!
!
   module subroutine construct_ao_mu_wx_molecular_system(molecule, mu_X, mu_Y, mu_Z, s1, s2)
!!
!!    Construct μ_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the μ_αβ integral in the array h. s1-s2 are the shells
!!    that alpha and beta belong to.
!!
!!    Note that the routine calculates the X, Y, and Z components 
!!    of the dipole simultaneously for the requested shells s1 and s2.
!!    (Because this is how Libint computes them.)
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
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_mu_wx_c(mu_X, mu_Y, mu_Z, s1_4, s2_4)
!
   end subroutine construct_ao_mu_wx_molecular_system
!
!
   module subroutine construct_ao_q_wx_molecular_system(molecule, q_xx, q_xy, q_xz, q_yy, q_yz, q_zz, s1, s2)
!!
!!    Construct q_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the q_αβ integral in the array h. s1-s2 are the shells
!!    that alpha and beta belong to.
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
      integer(i6) :: s1_4, s2_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_q_wx_c(q_xx, q_xy, q_xz, q_yy, q_yz, q_zz, s1_4, s2_4)
!
   end subroutine construct_ao_q_wx_molecular_system
!
!
end submodule ao_integrals
