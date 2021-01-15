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
   module subroutine prepare_for_jacobian_ccsd(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_ccsd
!
!
   module subroutine jacobian_transformation_ccsd(wf, c, rho)
!!
!!    Jacobian transformation 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
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
      class(ccsd), intent(inout) :: wf
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: c
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: rho
!
   end subroutine jacobian_transformation_ccsd
!
!
   module subroutine jacobian_ccsd_b2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^B2 = - sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)       :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)   :: rho_aibj
!
   end subroutine jacobian_ccsd_b2_ccsd
!
!
   module subroutine jacobian_ccsd_c2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^C2 = sum_kcl g_ljkc (t_ki^ac c_bl + t_li^bc c_ak + t_lk^ba c_ci)
!!                 - sum_kcl L_ljkc (t_il^ab c_ck + t_ik^ac c_bl)
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_c2_ccsd
!
!
  module subroutine jacobian_ccsd_d2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^D2 = - sum_kcd g_kcbd (t_ij^cd c_ak + t_kj^ad c_ci + t_ik^ca c_dj)
!!                       + sum_kcd L_kcbd (t_ik^ac c_dj + t_ij^ad c_ck)
!!
!!    Note: the code is structured so that we batch over the index b,
!!          where the integrals are made as g_kc_db = g_kcbd and held
!!          in some ordering or other throughout a given batch (i.e.,
!!          all five terms are constructed gradually in the batches).
!!
!!    The first term is constructed from the D2 intermediate: 
!!
!!       - sum_kcd g_kcbd t_ij^cd c_ak = - sum_d c_ak X_kijb  
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_d2_ccsd
!
!
    module subroutine jacobian_ccsd_e2_ccsd(wf, rho_aibj, c_aick)
!!
!!    Jacobian CCSD E2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^E2 = 2 sum_dlck t_bj,dl * L_kc,ld * c_ai,ck
!!                  - sum_ckld t_bj,dl * L_kc,ld * c_ak,ci
!!                = sum_ck Y_bjck (2c_aick - c_akci)
!!                = sum_ck Y_bjck v_aick
!!                        
!!    with
!!
!!       Y_bjck = t_bj,dl * L_kc,ld
!!       v_aick = (2c_aick - c_akci)
!!
!!    which is constructed in prepare_for_jacobian
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aick
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^F2 =   - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                       - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_f2_ccsd
!
!
   module subroutine jacobian_ccsd_g2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD G2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^G2 =  - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck
!!                       - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!!                       - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl
!!                = - sum_ck Y_bjck c_aick 
!!                  - sum_d Y_bd c_aidj
!!                  - sum_l Y_jl c_aibl
!!
!!    The intermediates are constructed once in prepare_for_jacobian
!!    in the routine save_jacobian_g2_intermediates_ccsd
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
   end subroutine jacobian_ccsd_g2_ccsd
!
!
   module subroutine jacobian_ccsd_h2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD H2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^H2 =  sum_ckdl t_ci,ak * g_kc,ld * c_bl,dj
!!                     + sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!!                   = sum_dl Y_aild c_bldj
!!                     + sum_dk Y_ajkd c_bkdi 
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
   end subroutine jacobian_ccsd_h2_ccsd
!
!
   module subroutine jacobian_ccsd_i2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_jk * c_ai,bk
!!                   + sum_ck L_bj,kc * c_ai,ck
!!                   - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj )
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
   end subroutine jacobian_ccsd_i2_ccsd
!
!
   module subroutine jacobian_ccsd_j2_ccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD J2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_abij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!                   =    sum_ckld Y_klij * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
   end subroutine jacobian_ccsd_j2_ccsd
!
!
   module subroutine jacobian_ccsd_k2_ccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD K2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_abij^K2 =    sum_kl g_kilj * c_akbl
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
   end subroutine jacobian_ccsd_k2_ccsd
!
!
   module subroutine save_jacobian_c2_intermediates_ccsd(wf)
!!
!!    Save jacobian c2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Constructs the intermediates for Jacobian C2: 
!!
!!       X_ljai = sum_ck g_ljkc t_ki^ac  
!!       X_kjbi = sum_lc g_ljkc t_li^bc 
!!       Y_ljai = sum_kc L_ljkc t_ik^ac 
!!
!!    used in the c2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_c2_intermediate
!!    which is a wf variable.
!!
      implicit none
!
      class(ccsd) :: wf
!
   end subroutine save_jacobian_c2_intermediates_ccsd
!
!
   module subroutine save_jacobian_d2_intermediate_ccsd(wf)
!!
!!    Save jacobian d2 intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Constructs the intermediate 
!!
!!       X_kbij = sum_dl g_kcbd t_ij^cd 
!!
!!    used in the d2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_d2_intermediate
!!    which is a wf variable.
!!
      implicit none
!
      class(ccsd) :: wf
!
   end subroutine save_jacobian_d2_intermediate_ccsd
!
!
   module subroutine save_jacobian_e2_intermediate_ccsd(wf)
!!
!!    Save jacobian e2 intermediate
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediate 
!!
!!       Y_bjck = sum_dl t_bjdl L_ldkc
!!
!!    used in the e2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_e2_intermediate
!!    which is a wf variable
!!
      implicit none
!
      class(ccsd) :: wf
!
   end subroutine save_jacobian_e2_intermediate_ccsd
!
!
   module subroutine save_jacobian_g2_intermediates_ccsd(wf)
!!
!!    Save jacobian g2 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediate 
!!
!!       Y_bjck = t_bldj L_kcdl 
!!
!!    The intermediates are stored in the file:
!!
!!       jacobian_g2_intermediate_vovo
!!
!!    G2 also needs intermediates formed for A1 term:
!!
!!       Y_jl   = t_ckdj L_kcld 
!!       Y_bd   = t_blck L_kcld 
!!
!!    Which are constructed in save_jacobian_a1_intermediates
!!    and stored on files
!!
!!       jacobian_a1_intermediate_oo
!!       jacobian_a1_intermediate_vv
!!
!!    which are wavefunction variables
!!
      implicit none
!
      class(ccsd) :: wf
!
   end subroutine save_jacobian_g2_intermediates_ccsd
!
!
   module subroutine save_jacobian_h2_intermediates_ccsd(wf)
!!
!!    Save jacobian h2 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediates 
!!
!!       Y_aild = t_ciak g_kcld
!!       Y_ajkd = t_cjal g_kcld
!!
!!    The intermediates are stored in the files:
!!
!!       jacobian_h2_intermediate_voov_1
!!       jacobian_h2_intermediate_voov_2
!!
!!    which are wavefunction variables
!!
      implicit none
!
      class(ccsd) :: wf
!
   end subroutine save_jacobian_h2_intermediates_ccsd
!
!
   module subroutine save_jacobian_j2_intermediate_ccsd(wf)
!!
!!    Save jacobian j2 intermediate
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediate
!!
!!       Y_klij = t_cidj g_kcld
!!
!!    The intermediates are stored in the files:
!!
!!       jacobian_j2_intermediate_oooo
!!
!!    which are wavefunction variables
!!
      implicit none
!
      class(ccsd) :: wf
!
   end subroutine save_jacobian_j2_intermediate_ccsd
