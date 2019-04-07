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
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transpose transformation (CC3)
!!    Alexander Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_A_cc3(wf, omega, c_ai, c_abij, sigma_ai, sigma_abij)
!!
!!    Terms of the transpose of the  CC3 Jacobi matrix
!!    Alex C. Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
   end subroutine jacobian_transpose_cc3_A_cc3
!
!
   module subroutine prepare_cc3_jacobian_transpose_integrals_cc3(wf)
!!
!!    Construct integrals needed in CC3 jacobian transpose and store on disk
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
!!    Construct some intermediates needed in CC3 jacobian transpose and store on disk
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
   module subroutine construct_X_and_Y_cc3(wf, i, j, k, t_abc, u_abc, X_acdi, X_acdj, X_acdk, Y_aikl,  &
                                             g_jbic, g_kbic, g_kbjc, g_ibjc, g_ibkc, g_jbkc,       &
                                             L_jbic, L_kbic, L_kbjc, L_ibjc, L_ibkc, L_jbkc)
!!
!!    Constructs the intermediates X_acdi and Y_akil used to compute the contributions to sigma_ai
!!
!!    X_acdi = sum_bjk (t^bac_ijk * g_jbkd - t^abc_ijk * L_jbkd)
!!    Y_akil = sum_bjc (t^bac_ijk * g_jblc - t^abc_ijk * L_jblc)
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
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_acdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_acdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_acdk
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout)   :: Y_aikl
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jbkc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!                       
   end subroutine construct_X_and_Y_cc3
!
!
   module subroutine jacobian_transpose_cc3_X_reader_cc3(wf, batch_x, X_acdx)
!!
!!    Read the X_acdx intermediate in the current batch
!!
!!    Based on omega_cc3_vvv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: X_acdx
!
   end subroutine jacobian_transpose_cc3_X_reader_cc3
!
!
   module subroutine jacobian_transpose_cc3_vvv_reader_cc3(wf, batch_x, g_bdcx)
!!
!!    Read the bdck, integrals in the current batch
!!
!!    Based on omega_cc3_vvv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_bdcx
!
   end subroutine jacobian_transpose_cc3_vvv_reader_cc3
!
!
   module subroutine jacobian_transpose_cc3_ov_vv_reader_cc3(wf, batch_y, batch_x, g_lycx, g_ybxd, L_ybxc)
!!
!!    Read the ljck, g_jbkc, L_jbkc integrals in the current batches
!!
!!    Based on omega_cc3_ov_vv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_lycx
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ybxd
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: L_ybxc
!
   end subroutine jacobian_transpose_cc3_ov_vv_reader_cc3
!
!
   module subroutine jacobian_transpose_cc3_write_intermediates_cc3(wf, batch_i, batch_j, batch_k, &
                                                                     X_acdi, X_acdj, X_acdk)
!!
!!    Write the contributions to the X_acdx intermediates to file in the respective batches
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_i, batch_j, batch_k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_i%length), intent(in) :: X_acdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_j%length), intent(in) :: X_acdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_k%length), intent(in) :: X_acdk
!
   end subroutine jacobian_transpose_cc3_write_intermediates_cc3