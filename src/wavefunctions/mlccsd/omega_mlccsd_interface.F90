!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   module subroutine construct_omega_mlccsd(wf, omega)
!!
!!    Construct omega (MLCCSD)
!!    Written by Sarai D. Folkestad, 2017-2019
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wf
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_mlccsd
!
!
   module subroutine omega_ccsd_a2_mlccsd(wf, omega_aibj)
!!
!!    Omega A2 term
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!      
!!    A2 = g_aibj 
!!
!!    INDEX RESTRICTIONS
!!
!!    a, b, i, j are CCSD indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
               intent(inout):: omega_aibj
!
   end subroutine omega_ccsd_a2_mlccsd
!
!
   module subroutine omega_ccsd_b2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega B2 term
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!      
!!    B2 =  sum_(cd)g_acbd * x_cidj = B2.1 + B.2.2
!!
!!    Structure: Batching over both a and b for B2.2.
!!                t^+_ci_dj = t_cidj + t_di_cj
!!                t^-_ci_dj = t_cidj - t_di_cj
!!                g^+_ac_bd = g_acbd + g_bc_ad
!!                g^-_ac_bd = g_acbd - g_bc_ad
!!
!!       omega_B2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!!       omega_B2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!!
!!
!!    INDEX RESTRICTIONS
!!
!!    a, b, i, j are CCSD indices
!!
!!    c, d, are CC2 indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
               intent(inout):: omega_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
   end subroutine omega_ccsd_b2_mlccsd
!
!
   module subroutine omega_ccsd_c2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    C2 = sum_(kl) x_ak_bl*(g_kilj + sum_(cd) x_cidj * g_kcld)
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD orbitals
!!
!!    k, l, c, d are CC2 + CCSD orbitals
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
               intent(inout):: omega_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
   end subroutine omega_ccsd_c2_mlccsd
!
!
   module subroutine omega_ccsd_d2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    D2 = -1/2 * sum_(ck) x_bk_cj*(g_kiac - 1/2 sum_(dl)x_al_di * g_kdlc)
!!                    - sum_(ck) x_bk_ci*(g_kj_ac - sum_(dl)x_al_dj * g_kdlc)
!!                    - 1/2 * sum_ck u_jk^bc g_acki
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, k, d, l are CCSD + CC2 indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout):: omega_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
   end subroutine omega_ccsd_d2_mlccsd
!
!
   module subroutine omega_ccsd_e2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega E2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    Calculates the D2 term,
!!
!!      E2: sum_ck u_jk^bc g_aikc
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!    where
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!    The first and second terms are referred to as E2.1 and E2.2.
!!    All terms are added to the omega vector of the wavefunction object wf.
!!
!!    The routine adds the terms in the following order: E2.2, E2.1
!!
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices 
!!
!!    c, d, k, l are CC2 + CCSD indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout):: omega_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
   end subroutine omega_ccsd_e2_mlccsd
!
!
   module subroutine omega_ccsd_f2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    Calculates the F2 term,
!!
!!      F2: sum_c x_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd)
!!        - sum_k x_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!    where
!!
!!        u_kl^bc = 2 * x_kl^bc - x_lk^bc.
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, d, k, l are CC2 + CCSD indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), & 
                  intent(inout):: omega_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
   end subroutine omega_ccsd_f2_mlccsd
