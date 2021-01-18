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
   module subroutine construct_omega_ccsd_complex(wf, omega)
!!
!!    Construct omega (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_ccsd_complex
!
!
   module subroutine omega_ccsd_a2_ccsd_complex(wf, omega_abij, t_abij)
!!
!!    Omega A2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!      
!!    A2 = sum_(cd)g_acbd * t_cidj 
!!
!!    Structure: Batching over both a and b 
!!                t^+_ci_dj = t_cidj + t_di_cj
!!                t^-_ci_dj = t_cidj - t_di_cj
!!                g^+_ac_bd = g_acbd + g_bc_ad
!!                g^-_ac_bd = g_acbd - g_bc_ad
!!
!!                omega_A2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2_bj_ai
!!                omega_A2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2_bi_aj
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(inout) :: omega_abij
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
   end subroutine omega_ccsd_a2_ccsd_complex
!
!
   module subroutine omega_ccsd_b2_ccsd_complex(wf, omega_abij, t_abij)
!!
!!    Omega B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Omega B2 = g_aibj + sum_(kl) t_akbl*(g_kilj + sum_(cd) t_cidj * g_kcld)
!!
!!    Structure: g_kilj is constructed first and reordered as g_klij.
!!    Then the contraction over cd is performed, and the results added to g_klij.
!!    t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(inout):: omega_abij
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
   end subroutine omega_ccsd_b2_ccsd_complex
!
!
   module subroutine omega_ccsd_c2_ccsd_complex(wf, omega_aibj, t_aibj)
!!
!!    Omega C2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Omega C2 = -1/2 * sum_(ck) t_bk_cj*(g_kiac -1/2 sum_(dl)t_al_di * g_kdlc)
!!                    - sum_(ck) t_bk_ci*(g_kj_ac - sum_(dl)t_al_dj * g_kdlc)
!!                    - 1/2 * sum_ck u_jk^bc g_acki
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
   end subroutine omega_ccsd_c2_ccsd_complex
!
!
   module subroutine omega_ccsd_d2_ccsd_complex(wf, omega_aibj, t_aibj)
!!
!!    Omega D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D2 term,
!!
!!      D2: sum_ck u_jk^bc g_aikc
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!    where
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!    The first and second terms are referred to as D2.1 and D2.2.
!!
!!    The routine adds the terms in the following order: D2.2, D2.1
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
   end subroutine omega_ccsd_d2_ccsd_complex
!
!
   module subroutine omega_ccsd_e2_ccsd_complex(wf, omega_aibj, t_aibj)
!!
!!    Omega E2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E2 term,
!!
!!      E2: sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd)
!!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!    where
!!
!!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!!
!!    The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!!    The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!!
!!    Both are permuted added to the projection vector element omega2(ai,bj) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(ccsd) :: wf
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
   end subroutine omega_ccsd_e2_ccsd_complex
!
!
   module subroutine construct_u_aibj_ccsd_complex(wf)
!!
!!    Construct u_aibj
!!    Written by Tor S. Haugland, Nov 2019
!!
!!    Constructs 
!!
!!       u_aibj = 2t_aibj - t_ajbi
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine construct_u_aibj_ccsd_complex
