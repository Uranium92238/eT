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
   module subroutine prepare_for_jacobian_transpose_ccsd_complex(wf)
!!
!!    Prepare for jacobian transpose
!!    Written by Tor S. Haugland, Andreas Skeidsvoll and Sarai D. Folkestad, Oct-Nov 2019
!!
!!    Creates intermediates needed in the jacobian transpose calculation.
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
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
   end subroutine jacobian_transpose_transformation_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_d1_intermediates_ccsd_complex(wf, t_aibj)
!!
!!    Save Jacobian transpose D1 intermediates
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and 
!!    Tor S. Haugland, Oct 2019
!!    
!!    Calculates the intermediate
!!
!!       X_lcki = sum_d t_dlck F_id
!!
!!    (E. F. K. and S. D. F 2017-2018)
!!      
!!    and write the intermediate to the file 'jacobian_transpose_d1_intermediate'.
!! 
!!    (T. S. H., Nov 2019)
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
!
   end subroutine save_jacobian_transpose_d1_intermediates_ccsd_complex
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
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Reads intermediate for term 1 from the file 'jacobian_transpose_d1_intermediate'. Removed
!!    re-order in term 2.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_e1_intermediates_ccsd_complex(wf, t_aibj, L_ilmd)
!!
!!    Save Jacobian transpose E1 intermediates
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and 
!!    Tor S. Haugland, Oct 2019
!!
!!    Calculates the intermediate
!!
!!       X_ilck = sum_md L_ilmd t_mk^dc = sum_md L_ilmd t_ckdm
!!
!!    (E. F. K and S. D. F., 2017-2018)
!!
!!    and saves it to the file 'jacobian_transpose_e1_intermediate_oovo'.
!!
!!    (T. S. H., Nov 2019)
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      complex(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_v), intent(in) :: L_ilmd
!
   end subroutine save_jacobian_transpose_e1_intermediates_ccsd_complex
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
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Reads intermediate for term 3 from the file 'jacobian_transpose_e1_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_f1_intermediates_ccsd_complex(wf, t_aibj, g_ikmc)
!!
!!    Save Jacobian transpose F1 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, and
!!    Andreas Skeidsvoll, Oct 2019
!!
!!    Construct and save the intermediates
!!
!!       X_ikdl = sum_mc t_lm^cd g_ikmc = t_mcdl g_ikmc
!!
!!    (E. F. K. and S. D. F. 2017-2018)
!! 
!!       X_lidk = sum_mc t_mk^dc g_mlic = t_mcdk g_limc
!!
!!    adds them reordered together
!!
!!       X_kdli = X_ikdl + X_lidk
!!
!!    and writes it to the file 'jacobian_transpose_f1_intermediate'.
!!
!!    (A. S. and S. D. F. Nov 2019)
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      complex(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_v), intent(in) :: g_ikmc
!
   end subroutine save_jacobian_transpose_f1_intermediates_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd_complex(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD F1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the F1 term,
!!
!!       sum_ckdlm (b_akdl t_lm^cd g_ikmc + b_dlak t_mk^dc g_mlic + b_ckdi t_ml^cd g_mkla)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
!!    Modified by Andreas Skeidsvoll, Oct 2019
!!
!!    Reads intermediate for term 1 and 2 from the file 'jacobian_transpose_f1_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_g1_intermediates_ccsd_complex(wf, t_aibj)
!!
!!    Save jacobian transpose g1 intermediates
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad and
!!    Tor S. Haugland, Oct 2019
!!
!!    Calculates the intermediate,
!!
!!       X_idkl = sum_ce t_ckel g_icde ordered as X_kdli
!!
!!    (E. F. K. and S. D. F. 2017-1028)
!!
!!    and saves it to the file 'jacobian_transpose_g1_intermediate'.
!!
!!    (T. H. S., Oct 2019)
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
!
   end subroutine save_jacobian_transpose_g1_intermediates_ccsd_complex
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
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Reads intermediate for term 1 from the file 'jacobian_transpose_g1_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
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
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_c2_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_d2_intermediates_ccsd_complex(wf, u_ckdl, L_dlbj)
!!
!!    Save Jacobian transpose D2 intermediates
!!    Written by Andreas Skeidsvoll, Tor S. Haugland,
!!    Sarai D. Folkestad and Eirik F. Kjønstad , Nov 2019
!!
!!    Constructs the intermediate
!!
!!       X_ckbj = sum_dl u_ckdl L_jbld
!!
!!    where u_ckdl = 2 t_ckdl - t_cldk.
!!
!!    and saves it to the file 'jacobian_transpose_d2_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_ckdl
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: L_dlbj ! Reordered L_jbld
!
   end subroutine save_jacobian_transpose_d2_intermediates_ccsd_complex
!
!
   module subroutine jacobian_transpose_ccsd_d2_ccsd_complex(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D2 term,
!!
!!       sum_ckdl b_aick L_jbld u_kl^cd
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
!!    Modified by Andreas Skeidsvoll and Tor S. Haugland, Nov 2019
!!
!!    Reads intermediate from the file 'jacobian_transpose_d2_intermediate'. Contribution from
!!    e2 was moved to d2, changing the term
!!       sum_ckdl b_aick L_jbld 2 t_kl^cd
!!    to
!!       sum_ckdl b_aick L_jbld (2 t_kl^cd - t_kl^dc)
!!
   implicit none
!
      class(ccsd) :: wf
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
!!       - sum_ckdl (b_aibl t_kl^cd L_kcjd + b_aicj t_kl^cd L_ldkb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
!!    Modified by Sarai D. Folkestad and Tor S. Haugland, Nov 2019
!!
!!    Reads term 1 and 2 intermediates from the files 'jacobian_transpose_e2_oo_intermediate'
!!    and 'jacobian_transpose_e2_vv_intermediate'.
!!
!!    Moved term to jacobian_transpose_d2,
!!       - sum_ckdl b_aicl t_kl^cd L_jbkd
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_e2_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_f2_intermediates_ccsd_complex(wf, t_ckdl, L_dlbi)
!!
!!    Save Jacobian transpose f2 intermediates
!!    Written by Eirik F. Kjønstad, Tor S. Haugland and 
!!    Sarai D. Folkestad, Nov 2019
!!
!!    Calculates the F2 intermediate,
!!
!!       X_ckbi = sum_dl t_ckdl L_ldib
!!
!!       (S. D. F and E. F. K. 2017-2018)
!!
!!    and writes it to the file 'jacobian_transpose_f2_intermediate'.
!!
!!       (T. S. H., Nov 2019)
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_dlbi ! Reordered L_ldib
!
   end subroutine save_jacobian_transpose_f2_intermediates_ccsd_complex
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
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Reads term 2 intermediate from the file 'jacobian_transpose_f2_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_f2_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_g2_intermediates_ccsd_complex(wf, t_aibj, g_kdib)
!!
!!    Save Jacobian transpose g2 intermediates
!!    Written by Tor S. Haugland, Eirik F. Kjønstad and 
!!    S. D. Folkestad, Nov 2019
!!
!!    Constructs intermediates 
!!
!!       X_clbi = sum_dk t_ckdl g_kbid (E. F. K. and S. D. F., 2017-2018)
!!       X_clib = sum_kd t_clkd g_kdib (T. S. H., Nov 2019)
!!
!!    and saves them to files 'jacobian_transpose_g2_intermediate' 
!!    and 'jacobian_transpose_g2_intermediate_2'
!!
!!       (T. S. H., Nov 2019)
!!
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      complex(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_kdib
!
   end subroutine save_jacobian_transpose_g2_intermediates_ccsd_complex
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
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Reads term 1 and 2 intermediates from the files 'jacobian_transpose_g2_intermediate'
!!    and 'jacobian_transpose_g2_intermediate_2'.
!!
      implicit none
!
      class(ccsd) :: wf
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
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
   end subroutine jacobian_transpose_ccsd_h2_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_i2_intermediates_ccsd_complex(wf, t_aibj, g_ovov)
!!
!!    Save Jacobian transpose i2 intermediates
!!    Written by Tor S. Haugland, Eirik F. Kjønstad and 
!!    Sarai D. Folkestad, Nov 2019
!!
!!    Construct intermediate
!!
!!       X_klij = sum_cd t_kl^cd g_icjd
!!
!!       (E. F. K. and S. D. F. 2017-2018)
!!
!!    and saves them to the file jacobian_transpose_i2_intermediate.
!!
!!       (T. S. H., Nov 2019)
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      complex(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_ovov
!
   end subroutine save_jacobian_transpose_i2_intermediates_ccsd_complex
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
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Reads term 2 intermediate from the file 'jacobian_transpose_i2_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      complex(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
   end subroutine jacobian_transpose_ccsd_i2_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_e2_oo_intermediate_ccsd_complex(wf, t_ckdl, L_ckdj)
!!
!!    Save Jacobian transpose e2 oo intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Constructs the intermediate 
!!
!!       X_jl = sum_kcd L_kcjd t_kl^cd
!!
!!    and saves it to the file 'jacobian_transpose_e2_intermediate'.
!!
!!    Modified by Sarai D. Folkestad, Nov 2019
!!
!!    Separated intermediate from jacobian_ccsd_e2_ccsd, it is now saved to file.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_ckdj ! L_kcjd
!
   end subroutine save_jacobian_transpose_e2_oo_intermediate_ccsd_complex
!
!
   module subroutine save_jacobian_transpose_e2_vv_intermediate_ccsd_complex(wf, t_ckdl, L_bkdl)
!!
!!    Save Jacobian transpose e2 vv intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Constructs the intermediate 
!!
!!       X_cd = sum_kdl t_kl^cd L_ldkb
!!
!!    and saves it to file 'jacobian_transpose_e2_intermediate_vv'.
!!
!!    Modified by Sarai D. Folkestad, Nov 2019
!!
!!    Separated intermediate from jacobian_ccsd_e2_ccsd, it is now saved to file.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_bkdl ! L_kbld
!
   end subroutine save_jacobian_transpose_e2_vv_intermediate_ccsd_complex
