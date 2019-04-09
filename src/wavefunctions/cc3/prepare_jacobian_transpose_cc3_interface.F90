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
!!    Write some integrals and intermediates to disk
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
   module subroutine construct_X_and_Y_cc3(wf, i, j, k, t_abc, u_abc, X_abdi, X_abdj, X_abdk,   &
                                          Y_aikl, g_lbic, g_lbjc, g_lbkc)
!!
!!    Constructs the intermediates X_abdi and Y_akil used to compute the contributions to sigma_ai
!!
!!    X_abdi = sum_cjk (t^cab_ijk * g_jckd + t^acb_ijk * g_jdkc - 2 * t^acb_ijk * g_jckd)
!!           = sum_cjk (t^cab_ijk + t^abc_ijk - 2 * t^acb_ijk) * g_jckd
!!    Y_akil = sum_bjc (t^bac_ijk * g_jblc + t^abc_ijk * g_jclb - 2 * t^abc_ijk * g_jblc)
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
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdk
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout)   :: Y_aikl
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_lbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_lbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_lbkc
!                       
   end subroutine construct_X_and_Y_cc3
!
!
   module subroutine jacobian_transpose_cc3_write_X_cc3(wf, batch_x, X_abdx)
!!
!!    Write the contributions to the X_abdi intermediate to file in the respective batches
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_x%length), intent(in) :: X_abdx
!
   end subroutine jacobian_transpose_cc3_write_X_cc3
!
!
   module subroutine sort_X_to_baid_and_write_cc3(wf)
!!
!!    Read in intermediate X_abdi from file, resort to X_baid and write to file again
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine sort_X_to_baid_and_write_cc3