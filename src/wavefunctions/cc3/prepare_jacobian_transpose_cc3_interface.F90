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
   module subroutine prepare_for_jacobian_transpose_cc3(wf)
!!
!!    Prepare for Jacobian transpose (CC3)
!!    Set up files containing integrals and intermediates for CC3 jacobian transpose
!!    Called from solver
!!
!!    Written by Rolf H. Myhre and Alexander Paul, April 2019
!!
   implicit none
!
   class(cc3) :: wf
!
   end subroutine prepare_for_jacobian_transpose_cc3
!
!
   module subroutine prepare_cc3_jacobian_transpose_integrals_cc3(wf)
!!
!!    Construct integrals needed in CC3 jacobian transpose and store on disk
!!    (ab|cd) ordered as abc,d
!!    (mi|lk) ordered as lm,ik
!!    (lb|kc) ordered as bcl,k
!!    (le|ck) ordered as lce,k
!!    (cd|mk) ordered as dcm,k
!!
!!    written by Rolf H. Myhre and Alexander Paul, April 2019
!!
      implicit none
!!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_jacobian_transpose_integrals_cc3
!
!
   module subroutine prepare_cc3_jacobian_transpose_intermediates_cc3(wf)
!!
!!    Construct X_abdi and Y_akil needed in CC3 jacobian transpose and store on disk
!!    For that: construct t^abc_ijk in single batches of ijk 
!!    and contract with the respective integrals
!!
!!    written by Rolf H. Myhre and Alexander Paul, April 2019
!!
      implicit none
!!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_jacobian_transpose_intermediates_cc3
!
!
   module subroutine construct_X_and_Y_cc3(wf, i, j, k, t_abc, u_abc, v_abc, Y_aijl,      &
                                           X_abdi, X_abdj, X_abdk, g_lbic, g_lbjc, g_lbkc)
!!
!!    Constructs the intermediates X_abdi and Y_akil used to compute the contributions to sigma_ai
!!
!!    X_abdi = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_kcjd
!!    Y_akil = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_lbkc
!!
!!    g_lbic, g_lbjc, g_lbkc can be used for g_pcqd as well: 
!!    The p(i,j,k) can be set in dgemm and q(i,j,k) is defined by the array used
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the loops
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout)   :: Y_aijl
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdk
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)              :: g_lbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)              :: g_lbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)              :: g_lbkc
!                       
   end subroutine construct_X_and_Y_cc3
!
!
   module subroutine sort_X_to_abid_and_write_cc3(wf)
!!
!!    Read in intermediate X_abdi from file, resort to X_baid and write to file again
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine sort_X_to_abid_and_write_cc3
