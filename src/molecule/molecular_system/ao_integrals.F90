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
   include "../../libint/g_wxyz_cdef.F90"
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
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral in the array h. s1 and s2 are 
!!    the shells that alpha and beta belong to. 
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
   module subroutine construct_ao_h_wx_kinetic_1der_molecular_system(molecule, &
                           h_1x, h_1y, h_1z, h_2x, h_2y, h_2z, s1, s2) 
!!
!!    Construct h_αβ kinetic 1st-derivative
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the h_αβ (kinetic contribution) first derivative integrals 
!!    in the arrays h_1x, h_1y, ..., h_2z. 
!!
!!    s1 and s2 are the shells that w and x respectively belong to.
!! 
!!    h_1x contains the first derivative with respect to the x component of the atom that s1 is centered on,
!!    h_1y contains the first derivative with respect to the y component of the atom that s1 is centered on,
!!    and so on...
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h_1x
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h_1y
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h_1z
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h_2x
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h_2y
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h_2z
!
      integer(i6) :: s1_4, s2_4 
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_h_wx_kinetic_1der_c(h_1x, h_1y, h_1z, h_2x, h_2y, h_2z, s1_4, s2_4) 
!
   end subroutine construct_ao_h_wx_kinetic_1der_molecular_system
!
!
   module subroutine construct_and_add_ao_h_wx_nuclear_1der_molecular_system(h_wxqk, s1, s2) 
!!
!!    Construct and add h_αβ nuclear 1st-derivative
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    adds parts of the h_αβ (nuclear contribution) first derivative integrals 
!!    in the arrays. 
!!
!!    s1 and s2 are the shells that w and x respectively belong to.
!! 
!!    h_1x contains the first derivative with respect to the x component of the atom that s1 is centered on,
!!    h_1y contains the first derivative with respect to the y component of the atom that s1 is centered on,
!!    and so on...
!!
      implicit none
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(:,:,:,:), contiguous :: h_wxqk
!
      integer(i6) :: s1_4, s2_4 
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_and_add_ao_h_wx_nuclear_1der_c(h_wxqk, s1_4, s2_4) 
!
   end subroutine construct_and_add_ao_h_wx_nuclear_1der_molecular_system
!
!
   module subroutine construct_ao_g_wxyz_molecular_system(molecule, g, s1, s2, s3, s4)
!!
!!    Construct g_αβγδ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral in the array g. s1-s4 are 
!!    the shells that alpha, beta, gamma and delta belong to. 
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule 
!
      integer, intent(in) :: s1, s2, s3, s4
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size, &
                           molecule%shell_limits(s3)%size, molecule%shell_limits(s4)%size), intent(out) :: g
!
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
      s3_4 = int(s3,i6)
      s4_4 = int(s4,i6)
!
      call construct_ao_g_wxyz_c(g, s1_4, s2_4, s3_4, s4_4)
!
   end subroutine construct_ao_g_wxyz_molecular_system
!
!
   module subroutine construct_ao_g_wxyz_1der_molecular_system(molecule, g_wxyzqk, s1, s2, s3, s4, &
                                                                  a1, a2, a3, a4)
!!
!!    Construct g_αβγδ 1der
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral in the array g. s1-s4 are 
!!    the shells that alpha, beta, gamma and delta belong to. 
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule 
!
      integer, intent(in) :: s1, s2, s3, s4 ! Shells 
      integer, intent(in) :: a1, a2, a3, a4 ! Atoms 
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size,    &
                           molecule%shell_limits(s3)%size, molecule%shell_limits(s4)%size,   &
                           3, 4), intent(out) :: g_wxyzqk
!
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
      integer(i6) :: a1_4, a2_4, a3_4, a4_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
      s3_4 = int(s3,i6)
      s4_4 = int(s4,i6)
!
      a1_4 = int(a1,i6)
      a2_4 = int(a2,i6)
      a3_4 = int(a3,i6)
      a4_4 = int(a4,i6)
!
      call construct_ao_g_wxyz_1der_c(g_wxyzqk, s1_4, s2_4, s3_4, s4_4, &
                                                a1_4, a2_4, a3_4, a4_4)
!
   end subroutine construct_ao_g_wxyz_1der_molecular_system
!
!
   module subroutine construct_ao_g_wxyz_epsilon_molecular_system(g, s1, s2, s3, s4, eps, thread, skip, &
                                                            n1, n2, n3, n4)
!!
!!    Construct g_αβγδ epsilon 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral in the array g. s1-s4 are 
!!    the shells that alpha, beta, gamma and delta belong to. 
!!
!!    This is the most efficient routine to calculate these 
!!    integrals, due mostly to the precision control parameter 
!!    (epsilon) but also because we avoid computing information
!!    that might already be available (such as thread ID and 
!!    the size of each of the shells, n1-n4). 
!!
!!       - To get thread, use omp_get_thread_num()
!!       - To get n1-n4, see the shells array of the wavefunction's system object
!!       - Skip is an integer that is either 0 or 1 on exit, where a 0 means 
!!         Libint decided to skip calculating the integrals. In other words,
!!         be sure to zero g if skip is 1 and this is necessary for your 
!!         application! The reason for this integer is that it is sometimes
!!         not necessary to actually zero g (which can be a relevant penalty).
!!       - In order to determine epsilon, one useful approach is consider the 
!!         desired accuracy epsilon' = 1.0D-14 of Y, where e.g.
!!
!!             Y_αβ = Y_αβ + g_αβγδ X_γδ.
!!
!!         If X_γδ is on the order 1.0D-5, then to get 1.0D-14 accuracy in g_αβγδ X_γδ
!!         requires only that g_αβγδ is accurate to 1.0D-9 (the error is multiplied 
!!         by X_γδ). The Libint integral can be much faster for large epsilons,
!!         but care should be taken when dynamically changing epsilon.
!!
      implicit none
!
      integer, intent(in) :: s1, s2, s3, s4
      integer, intent(in) :: thread, n1, n2, n3, n4
!
      real(dp), dimension(:) :: g
!
      real(dp), intent(in) :: eps 
!
      integer, intent(inout) :: skip
!
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
      integer(i6) :: n1_4, n2_4, n3_4, n4_4
!
      integer(i6) :: thread_4, skip_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
      s3_4 = int(s3,i6)
      s4_4 = int(s4,i6)
!
      n1_4 = int(n1,i6)
      n2_4 = int(n2,i6)
      n3_4 = int(n3,i6)
      n4_4 = int(n4,i6)
!
      thread_4 = int(thread,i6)
      skip_4 = int(skip,i6)
!
      call construct_ao_g_wxyz_epsilon_c(g, s1_4, s2_4, s3_4, s4_4, eps, & 
                                       thread_4, skip_4, n1_4, n2_4, n3_4, n4_4)
!
      skip = int(skip_4)
!
   end subroutine construct_ao_g_wxyz_epsilon_molecular_system
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
      real(dp), dimension(molecule%shell_limits(s1)%size,molecule%shell_limits(s2)%size), intent(out) :: s
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
   module subroutine construct_ao_s_wx_1der_molecular_system(molecule, s_1x, s_1y, s_1z, s_2x, s_2y, s_2z, s1, s2) 
!!
!!    Construct s_αβ 1st-derivative
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the s_αβ first derivative integrals in the arrays s_1x, s_1y, ..., s_2z. 
!!
!!    s1 and s2 are the shells that w and x respectively belong to.
!! 
!!    s_1x contains the first derivative with respect to the x component of the atom that s1 is centered on,
!!    s_1y contains the first derivative with respect to the y component of the atom that s1 is centered on,
!!    and so on...
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      integer, intent(in) :: s1, s2
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: s_1x
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: s_1y
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: s_1z
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: s_2x
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: s_2y
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: s_2z
!
      integer(i6) :: s1_4, s2_4 
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_s_wx_1der_c(s_1x, s_1y, s_1z, s_2x, s_2y, s_2z, s1_4, s2_4) 
!
   end subroutine construct_ao_s_wx_1der_molecular_system
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
