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
   module subroutine prepare_for_jacobian_transpose_ccsd_complex(wf)
!!
!!    Prepare for jacobian transpose
!!    Written by Tor S. Haugland, Oct 2019
!!
!!    Creates intermediates needed in the jacobian transpose calculation.
!!
!!    Based on prepare_for_jacobian_ccsd by E. F. Kjønstad and S. D. Folkestad
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_transpose_ccsd_complex
!
!
   module subroutine jacobian_transpose_transformation_ccsd_complex(wf, b)
!!
!!    Jacobian transpose transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    The transformation is performed as sigma^T = b^T A, where b is the vector
!!    sent to the routine. On exit, the vector b is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      complex(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
   end subroutine jacobian_transpose_transformation_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_d1_ccsd_complex(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD D1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D1 term,
!!
!!       - sum_ckdl (b_ckal F_id t_kl^cd + b_ckdi F_la t_kl^cd),
!!
!!    and adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_e1_ccsd_complex(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD E1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E1 term,
!!
!!       sum_ckdle (b_ckdi L_dale t_kl^ce + b_ckdl L_deia t_kl^ce)
!!      -sum_ckdlm (b_ckal L_ilmd t_km^cd + b_ckdl L_mlia t_km^cd)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
!!    The routine adds the third and forth terms first.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd_complex(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD F1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the F1 term,
!!
!!       sum_ckdlm (b_akdl t_lm^cd g_ikmc + b_ckal t_ml^cd g_mkid + b_ckdi t_ml^cd g_mkla)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_g1_ccsd_complex(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD G1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the G1 term,
!!
!!       - sum_ckdle (b_akdl t_kl^ce g_icde + b_cidl t_kl^ce g_kade + b_cldi t_kl^ce g_keda)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_g1_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_b2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the B2 term,
!!
!!       sum_c b_aicj F_cb - sum_k b_aibk F_jk + sum_ck b_aick L_ckjb
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_b2_ccsd_complex
!

!
   module subroutine jacobian_transpose_ccsd_c2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD C2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the C2 term,
!!
!!       - sum_ck (b_ajck g_ibck + b_akcj g_ikcb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_c2_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_d2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D2 term,
!!
!!       2 * sum_ckdl b_aick L_jbld t_kl^cd
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
   implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_d2_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_e2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD E2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E2 term,
!!
!!       - sum_ckdl (b_aibl t_kl^cd L_kcjd + b_aicl t_kl^cd L_jbkd + b_aicj t_kl^cd L_ldkb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_e2_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_f2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD F2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the F2 term,
!!
!!       - sum_ckdl (b_alck t_kl^cd L_jbid + b_ajck t_kl^cd L_ldib + b_djck t_kl^cd L_ialb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_f2_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_g2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD G2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the G2 term,
!!
!!       sum_ckdl (b_alcj t_kl^cd g_kbid + b_ajcl t_kl^cd g_kdib)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_g2_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_h2_ccsd_complex(wf, sigma_abij, b_abij)
!!
!!    Jacobian transpose CCSD H2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the H2 term,
!!
!!       sum_kl b_akbl g_ikjl + sum_cd b_cidj g_cadb
!!
!!    and adds it to the transformed vector sigma_abij.
!!
!!    In this routine, the b and sigma vectors are ordered as
!!
!!       b_abij = b_aibj
!!       sigma_abij = sigma_abij
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
   end subroutine jacobian_transpose_ccsd_h2_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_i2_ccsd_complex(wf, sigma_abij, b_abij)
!!
!!    Jacobian transpose CCSD I2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the I2 term,
!!
!!       sum_ckdl b_cidj t_kl^cd g_kalb + sum_ckdl b_akbl t_kl^cd g_icjd
!!
!!    and adds it to the transformed vector sigma_abij.
!!
!!    In this routine, the b and sigma vectors are ordered as
!!
!!       b_abij = b_aibj
!!       sigma_abij = sigma_abij
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
   end subroutine jacobian_transpose_ccsd_i2_ccsd_complex