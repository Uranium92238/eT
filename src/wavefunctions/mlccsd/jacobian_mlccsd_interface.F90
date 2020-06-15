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
   module subroutine prepare_for_jacobian_mlccsd(wf)
!!
!!    Prepare for Jacobian
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none 
!
      class(mlccsd), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_mlccsd
!
!
   module subroutine jacobian_transformation_mlccsd(wf, c)
!!
!!    Jacobian transformation (MLCCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_ai = rho_ai,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine jacobian_transformation_mlccsd
!
!
   module subroutine jacobian_cc2_b2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    The doubles-doubles part of the CC2 Jacobian in non-canonical basis  
!!
!!       rho_aibj += F_bc c_cjai - F_kj c_aibk
!!
!!    INDEX RESTRICTIONS:
!!
!!    All indices are CC2 indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o,&
                     wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o), intent(inout) :: c_aibj
      real(dp), dimension(wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o,&
                     wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o), intent(inout) :: rho_aibj   
!
   end subroutine jacobian_cc2_b2_mlccsd
!
!
   module subroutine jacobian_ccsd_b2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian MLCCSD B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^B2 = - sum_kc (F_kc x_ij^ac c_bk + F_kc x_ik^ab c_cj)
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    k is full index in first term and CCSD + CC2 in second term 
!!    c is full index in second term and CCSD + CC2 in first term 
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                         :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_b2_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_1_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian MLCCSD C2-1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += sum_kcl g_ljkc x_ki^ac c_bl
!!
!!    INDEX RESTRICTIONS:
!!
!!     a, i, b, j are CCSD indices
!!
!!     c, k are CCSD + CC2, l is unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                                           intent(inout) :: rho_aibj
!
   end subroutine jacobian_ccsd_c2_1_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian MLCCSD C2-2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += sum_kcl g_ljkc x_li^bc c_ak 
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    l and c are CCSD + CC2, k is unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                              wf%n_ccsd_v, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
     end subroutine jacobian_ccsd_c2_2_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_3_mlccsd(wf, rho_aibj, c_ai, g_ljkc)
!!
!!    Jacobian MLCCSD C2-3
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += sum_kcl g_ljkc x_lk^ba c_ci
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    k and l are CCSD + CC2, c is unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  wf%n_ccsd_v, wf%n_ccsd_o), intent(inout) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o, &
                  wf%n_ccsd_o + wf%n_cc2_o, wf%n_v), intent(in) :: g_ljkc
!
     end subroutine jacobian_ccsd_c2_3_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_4_mlccsd(wf, rho_aibj, c_ai, L_ljck)
!!
!!    Jacobian MLCCSD C2-4
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += - sum_kcl L_ljkc x_il^ab c_ck  
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    l is CCSD + CC2, c and k are unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  wf%n_ccsd_v, wf%n_ccsd_o), intent(inout) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, &
                  wf%n_ccsd_o, wf%n_v, wf%n_o), intent(in) :: L_ljck
!
     end subroutine jacobian_ccsd_c2_4_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                                           intent(inout) :: rho_aibj
!
   end subroutine jacobian_ccsd_c2_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_1_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2-1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += - sum_kcd g_kcbd x_ij^cd c_ak
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout) :: rho_aibj
!
   end subroutine jacobian_ccsd_d2_1_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_2_mlccsd(wf, rho_aibj, c_ai, batch_b, g_dkbc)
!!
!!    Jacobian CCSD D2-2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += - sum_kcd g_kcbd x_kj^ad c_ci
!!
      implicit none
!
      class(mlccsd) :: wf
      type(batching_index), intent(in) :: batch_b
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o), &
                     intent(inout) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                           batch_b%length, wf%n_v), intent(in) :: g_dkbc
!
   end subroutine jacobian_ccsd_d2_2_mlccsd
!
!
   module subroutine jacobian_ccsd_d2_3_mlccsd(wf, rho_aibj, c_ai, batch_b, g_kcbd)
!!
!!    Jacobian CCSD D2-3
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += - sum_kcd g_kcbd x_ik^ca c_dj
!!
      implicit none
!
      class(mlccsd) :: wf
      type(batching_index), intent(in) :: batch_b
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                           batch_b%length, wf%n_ccsd_o), intent(inout) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, &
                  wf%n_ccsd_v + wf%n_cc2_v, batch_b%length, wf%n_v), intent(in) :: g_kcbd
!
   end subroutine jacobian_ccsd_d2_3_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_4_mlccsd(wf, rho_aibj, c_ai, batch_b, L_ckbd)
!!
!!    Jacobian CCSD D2-4
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += sum_kcd L_kcbd x_ik^ac c_dj
!!
      implicit none
!
      class(mlccsd) :: wf
      type(batching_index), intent(in) :: batch_b
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  batch_b%length, wf%n_ccsd_o), intent(inout) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, &
                  wf%n_ccsd_o + wf%n_cc2_o, batch_b%length, wf%n_v), intent(in) :: L_ckbd
!
   end subroutine jacobian_ccsd_d2_4_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_5_mlccsd(wf, rho_aibj, c_ai, batch_b, L_ckbd)
!!
!!    Jacobian CCSD D2-4
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 +=  sum_kcd L_kcbd x_ij^ad c_ck
!!
      implicit none
!
      class(mlccsd) :: wf
      type(batching_index), intent(in) :: batch_b
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  batch_b%length, wf%n_ccsd_o), intent(inout) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, batch_b%length, &
                  wf%n_ccsd_v + wf%n_cc2_v), intent(in) :: L_ckbd
!
   end subroutine jacobian_ccsd_d2_5_mlccsd
!
!
   module subroutine jacobian_ccsd_d2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                        intent(inout) :: rho_aibj
!
   end subroutine jacobian_ccsd_d2_mlccsd
!
!
    module subroutine jacobian_ccsd_e2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD E2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^E2 = 2 sum_dlck x_bj,dl * L_kc,ld * c_ai,ck
!!                  - sum_dlck x_bj,dl * L_kc,ld * c_ak,ci
!!
!!                = 2 Y_ck,bj * c_ai,ck - Y_ck,bj * c_ak,ci
!!
!!                = Y_ck,bj  (2 c_ai,ck - c_ak,ci)
!!
!!
!!    INDEX RESTRICTIONS:
!!
!!       a, i, b, j are CCSD indices
!!
!!       d, l, c, k are CC2 + CCSD indices
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_e2_mlccsd
!
!
   module subroutine jacobian_ccsd_f2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian MLCCSD F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!       rho_aibj^F2 = - sum_ckdl x_ai,dj * L_kc,ld * c_bl,ck
!!                     - sum_ckdl x_ai_bl * L_kc,ld * c_ck,dj
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
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
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_f2_mlccsd
!
!
   module subroutine jacobian_ccsd_g2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian MLCCSD G2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^G2 =  - sum_ckdl x_bl,dj * L_kc,ld * c_ai,ck
!!                   - sum_ckdl x_ck_bl * L_kc,ld * c_ai,dj
!!                   - sum_ckld x_ck,dj * L_kc,ld * c_ai,bl
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
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
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_g2_mlccsd
!
!
   module subroutine jacobian_ccsd_h2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD H2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!       rho_aibj^H2 =   sum_ckdl x_ci,ak * g_kc,ld * c_bl,dj
!!                     + sum_ckdl x_cj,al * g_kc,ld * c_bk,di
!!
!!    INDEX RESTRICTIONS
!!
!!    a, i, b, j are CCSD indices
!!    c, k, d, l are CC2 + CCSD indices
!!
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_h2_mlccsd
!
!
   module subroutine jacobian_ccsd_i2_1_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_kj * c_ai,bk
!!
!!    INDEX RESTRICTIONS:
!!
!!       a, i, b, j are CCSD indices
!!
!!       c is CC2 index 
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_i2_1_mlccsd
!
!
   module subroutine jacobian_ccsd_i2_2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^I2 =  sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj )
!!
!!    INDEX RESTRICTIONS:
!!
!!       a, i, b, j are CCSD indices
!!
!!       c, k is CC2 + CCSD index 
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_i2_2_mlccsd
!
!
   module subroutine jacobian_ccsd_i2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_kj * c_ai,bk
!!                   + sum_ck L_bj,kc * c_ai,ck
!!                   - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj )
!!
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_i2_mlccsd
!
!
   module subroutine jacobian_ccsd_j2_mlccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian MLCCSD J2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!       rho_abij^J2 =    sum_ckld x_ci,dj * g_kc,ld * c_ak,bl
!!                      + sum_ckdl x_ak,bl * g_kc,ld * c_ci,dj
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o) :: rho_abij
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_v + wf%n_cc2_v, &
                          wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_abij
!
   end subroutine jacobian_ccsd_j2_mlccsd
!
!
   module subroutine jacobian_ccsd_k2_mlccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian MLCCSD K2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_abij^K2 =    sum_kl g_kilj * c_akbl
!!                       + sum_cd g_acbd * c_cidj (in omega term)
!!
!!    For the last term we batch over a and b and
!!    add each batch to rho_aibj
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o) :: rho_abij
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_v + wf%n_cc2_v, &
                         wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_abij
!
   end subroutine jacobian_ccsd_k2_mlccsd
!
!
   module subroutine save_jacobian_c2_intermediates_mlccsd(wf)
!!
!!    Save jacobian c2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediates for Jacobian C2: 
!!
!!       X_ljai = sum_ck g_ljkc x_ki^ac - sum_kc L_ljkc x_ik^ac
!!
!!          Index restrictions: a, i, j are CCSD indices, 
!!          c and k are CCSD + CC2 indices, and l is unrestricted 
!!
!!       X_kjbi = g_ljkc x_li^bc
!!
!!          Index restrictions: b, i, j are CCSD indices, 
!!          c and l are CCSD + CC2 indices, and k is unrestricted 
!!
!!    used in the c2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_c2_intermediate_oovo_1
!!    and jacobian_c2_intermediate_oovo_2 which are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine save_jacobian_c2_intermediates_mlccsd
!
!
   module subroutine save_jacobian_d2_intermediate_mlccsd(wf)
!!
!!    Save jacobian d2 intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_kibj = sum_dl g_kcbd x_ij^cd 
!!
!!       Index restrictions: b, i, j are CCSD indices,  
!!       c and d are CC2 + CCSD indices, and k is unrestricted
!!
!!    used in the d2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_d2_intermediate
!!    which is a wf variable.
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine save_jacobian_d2_intermediate_mlccsd
!
!
   module subroutine save_jacobian_e2_intermediate_mlccsd(wf, L_ckdl)
!!
!!    Save jacobian e2 intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_ckbj = sum_dl L_kcld x_jl^bd 
!!
!!       Index restrictions: b, j are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    used in the e2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_e2_intermediate
!!    which is a wf variable.
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: L_ckdl
!
   end subroutine save_jacobian_e2_intermediate_mlccsd
!
!
   module subroutine save_jacobian_g2_intermediates_mlccsd(wf, L_ckdl)
!!
!!    Save jacobian g2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_ckbj = sum_dl L_kcld x_lj^bd 
!!
!!       X_db = sum_ckl L_kcld x_lk^bc 
!!
!!       X_lj = sum_cdl L_kcld x_kj^cd
!!
!!       Index restrictions: b, j are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    used in the g2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the files jacobian_G2_intermediate_vovo
!!    jacobian_G2_intermediate_vv, and jacobian_G2_intermediate_oo which
!!    are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: L_ckdl
!
   end subroutine save_jacobian_g2_intermediates_mlccsd
!
!
   module subroutine save_jacobian_h2_intermediates_mlccsd(wf, g_ckdl)
!!
!!    Save jacobian h2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_aidl = x^ca_ik g_kcld
!!
!!       X_ajdk = x^ca_jl g_kcld
!!
!!       Index restrictions: a and i are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    used in the h2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_h2_intermediate_ovov_1
!!    and jacobian_h2_intermediate_ovov_2 which are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: g_ckdl
!
   end subroutine save_jacobian_h2_intermediates_mlccsd
!
!
   module subroutine save_jacobian_j2_intermediates_mlccsd(wf, g_klcd)
!!
!!    Save jacobian j2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_klij = g_kcld x^cd_ij
!!
!!       Index restrictions: i and j are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    Also saves the integral g_kcld ordered as g_klcd
!!
!!    used in the j2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_j2_intermediate_oooo
!!    and jacobian_ji_intermediate_oovv which are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_v + wf%n_cc2_v), intent(in) :: g_klcd
!
   end subroutine save_jacobian_j2_intermediates_mlccsd
