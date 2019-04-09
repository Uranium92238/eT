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
module ao_integral_tool_class
!
!!
!!    AO integral tool class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Note: C++ is not fond of 64-bit integers, so ints are explicitly 
!!    translated to 32-bits here before calling Libint routines 
!!
!
!  Disk & memory class modules
!
   use file_class
   use memory_manager_class
!
!  Fortran interfaces to C++ routines
!
   use iso_c_binding
!
 !  include "../libint/h_wx_cdef.F90"
   include "../libint/s_wx_cdef.F90"
   include "../libint/mu_wx_cdef.F90"
   include "../libint/q_wx_cdef.F90"
   include "../libint/g_wxyz_cdef.F90"
!
!  Class definition
!
   type :: ao_integral_tool
!
!     No attributes yet
!
   contains
!
   !   procedure, nopass :: construct_ao_h_wx             => construct_ao_h_wx_ao_integral_tool  ! h_αβ
      procedure, nopass :: construct_ao_s_wx             => construct_ao_s_wx_ao_integral_tool  ! s_αβ
      procedure, nopass :: construct_ao_mu_wx            => construct_ao_mu_wx_ao_integral_tool ! μ_αβ
      procedure, nopass :: construct_ao_q_wx             => construct_ao_q_wx_ao_integral_tool  ! q_αβ
!
      generic :: construct_ao_g_wxyz => construct_ao_g_wxyz_2_ao_integral_tool, &
                                          construct_ao_g_wxyz_4_ao_integral_tool
!
      generic :: construct_ao_g_wxyz_epsilon => construct_ao_g_wxyz_epsilon_1_ao_integral_tool, &
                                                   construct_ao_g_wxyz_epsilon_2_ao_integral_tool, &
                                                   construct_ao_g_wxyz_epsilon_4_ao_integral_tool
!
      procedure, nopass :: construct_ao_g_wxyz_2_ao_integral_tool, &
                              construct_ao_g_wxyz_4_ao_integral_tool ! g_αβγδ
!
      procedure, nopass :: construct_ao_g_wxyz_epsilon_1_ao_integral_tool, &
                             construct_ao_g_wxyz_epsilon_2_ao_integral_tool, &
                             construct_ao_g_wxyz_epsilon_4_ao_integral_tool ! g_αβγδ
!
   end type ao_integral_tool
!
!
contains
!
!
!    subroutine construct_ao_h_wx_ao_integral_tool(h, s1, s2)
! !!
! !!    Construct h_αβ
! !!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
! !!
! !!    Fortran wrapper for the C++ routine that calculates and
! !!    saves parts of the h_αβ integral in the array h. s1-s2 are the shells
! !!    that alpha and beta belong to.
! !!
!       implicit none
! !
!       real(dp), dimension(:,:), intent(inout) :: h
!       integer, intent(in) :: s1, s2
!       integer(i6) :: s1_4, s2_4
! !
!       s1_4 = int(s1,i6)
!       s2_4 = int(s2,i6)
! !
!       call construct_ao_h_wx_c(h, s1_4, s2_4)
! !
!    end subroutine construct_ao_h_wx_ao_integral_tool
!
!
   subroutine construct_ao_s_wx_ao_integral_tool(s, s1, s2)
!!
!!    Construct s_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the s_αβ integral in the array h. s1-s2 are the shells
!!    that alpha and beta belong to.
!!
      implicit none
!
      real(dp), dimension(:,:), intent(inout) :: s
      integer, intent(in) :: s1, s2
      integer(i6) :: s1_4, s2_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_s_wx_c(s, s1_4, s2_4) 
!
   end subroutine construct_ao_s_wx_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_4_ao_integral_tool(g, s1, s2, s3, s4)
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
      real(dp), dimension(:,:,:,:), intent(inout) :: g
!
      integer, intent(in) :: s1, s2, s3, s4
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
      s3_4 = int(s3,i6)
      s4_4 = int(s4,i6)
!
      call construct_ao_g_wxyz_c(g, s1_4, s2_4, s3_4, s4_4)
!
   end subroutine construct_ao_g_wxyz_4_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_2_ao_integral_tool(g, s1, s2, s3, s4)
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
      real(dp), dimension(:,:), intent(inout) :: g
!
      integer, intent(in) :: s1, s2, s3, s4
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
      s3_4 = int(s3,i6)
      s4_4 = int(s4,i6)
!
      call construct_ao_g_wxyz_c(g, s1_4, s2_4, s3_4, s4_4)
!
   end subroutine construct_ao_g_wxyz_2_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_epsilon_1_ao_integral_tool(g, s1, s2, s3, s4, eps, thread, skip, &
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
      real(dp), dimension(:), intent(inout) :: g
!
      real(dp), intent(in) :: eps 
!
      integer, intent(in) :: s1, s2, s3, s4, thread, n1, n2, n3, n4
      integer :: skip 
!
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
      integer(i6) :: n1_4, n2_4, n3_4, n4_4
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
   end subroutine construct_ao_g_wxyz_epsilon_1_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_epsilon_2_ao_integral_tool(g, s1, s2, s3, s4, eps, thread, skip, &
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
      real(dp), dimension(:,:), intent(inout) :: g
!
      real(dp), intent(in) :: eps 
!
      integer, intent(in) :: s1, s2, s3, s4, thread, n1, n2, n3, n4
      integer :: skip 
!
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
      integer(i6) :: n1_4, n2_4, n3_4, n4_4
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
   end subroutine construct_ao_g_wxyz_epsilon_2_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_epsilon_4_ao_integral_tool(g, s1, s2, s3, s4, eps, thread, skip, &
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
      real(dp), dimension(:,:,:,:), intent(inout) :: g
!
      real(dp), intent(in) :: eps 
!
      integer, intent(in) :: s1, s2, s3, s4, thread, n1, n2, n3, n4
      integer :: skip 
!
      integer(i6) :: s1_4, s2_4, s3_4, s4_4
      integer(i6) :: n1_4, n2_4, n3_4, n4_4
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
   end subroutine construct_ao_g_wxyz_epsilon_4_ao_integral_tool
!
!
   subroutine construct_ao_mu_wx_ao_integral_tool(mu_X, mu_Y, mu_Z, s1, s2)
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
      real(dp), dimension(:,:), intent(inout) :: mu_X ! x component
      real(dp), dimension(:,:), intent(inout) :: mu_Y ! y component 
      real(dp), dimension(:,:), intent(inout) :: mu_Z ! z component
!
      integer, intent(in) :: s1, s2
      integer(i6) :: s1_4, s2_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_mu_wx_c(mu_X, mu_Y, mu_Z, s1_4, s2_4)
!
   end subroutine construct_ao_mu_wx_ao_integral_tool
!
!
   subroutine construct_ao_q_wx_ao_integral_tool(q_xx, q_xy, q_xz, q_yy, q_yz, q_zz, s1, s2)
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
      real(dp), dimension(:,:), intent(inout) :: q_xx 
      real(dp), dimension(:,:), intent(inout) :: q_xy 
      real(dp), dimension(:,:), intent(inout) :: q_xz 
      real(dp), dimension(:,:), intent(inout) :: q_yy 
      real(dp), dimension(:,:), intent(inout) :: q_yz 
      real(dp), dimension(:,:), intent(inout) :: q_zz 
!
      integer, intent(in) :: s1, s2
      integer(i6) :: s1_4, s2_4
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_q_wx_c(q_xx, q_xy, q_xz, q_yy, q_yz, q_zz, s1_4, s2_4)
!
   end subroutine construct_ao_q_wx_ao_integral_tool
!
!
end module ao_integral_tool_class

