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
   module subroutine prep_cc3_integrals_t3_abc_batch_cc3(wf)
!!
!!    Prepares files containing the integrals needed 
!!    to construct t3-amplitudes in batches of a,b,c
!!
!!    (bd|ck) ordered as dk,bc
!!    (lj|ck) ordered as ljk,c
!!
!!    written by Alexander Paul, July 2019
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine prep_cc3_integrals_t3_abc_batch_cc3
!
!
   module subroutine prep_cc3_integrals_R3_abc_batch_cc3(wf, R_ai)
!!
!!    Prepares the files containing the R1-transformed 
!!    integrals needed to construct R3 in batches of a,b,c
!!    NB: The integrals (bd|ck) and (lj|ck) constructed in 
!!        prep_cc3_integrals_t3_abc_batch_cc3 are also needed
!!
!!    g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck')   ordered as dk,bc
!!    g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k)   ordered as ljk,c
!!
!!    written by Alexander Paul, July 2019
!!    Based on construct_c1_integrals_cc3 written by Rolf H. Myhre and A. Paul
!!
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine prep_cc3_integrals_R3_abc_batch_cc3
!
!
   module subroutine prep_cc3_integrals_L3_abc_batch_cc3(wf)
!!
!!    Prepares the files containing the integrals needed to 
!!    construct tbar3 and L3 in batches of a,b,c
!!
!!    (db|kc) ordered as dk,bc
!!    (jl|kc) ordered as ljk,c
!!    L_jbkc  ordered as jk,bc
!!
!!    written by Alexander Paul, July 2019
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine prep_cc3_integrals_L3_abc_batch_cc3
!
!
   module subroutine omega_cc3_W_calc_abc_batch_cc3(wf, a, b, c,           &
                                                   t_ijk, u_ijk, t_ijab,   &
                                                   g_ljak, g_ljbk, g_ljck, &
                                                   g_bdak, g_cdak, g_cdbk, &
                                                   g_adbk, g_adck, g_bdck, &
                                                   keep_t)
!!
!!    Calculate the the contributions to the t_3 amplitudes
!!    for virtual indices a,v,c
!!
!!    Contributions to W
!!     W^abc_ijk = P^abc_ijk(\sum_d t^ad_ij(bd|ck) - \sum_l t^ab_il(lj|ck))
!!
!!    Written by Rolf H. Myhre and Alexander Paul July 2019
!!    based on omega_cc3_W_calc written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_ijk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_ijk
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in)   :: t_ijab
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljak
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljbk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljck
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_bdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_cdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_cdbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_adbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_adck
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_bdck
!
      logical, optional, intent(in) :: keep_t ! If present and true, t_abc is not overwritten by first dgemm
!
   end subroutine omega_cc3_W_calc_abc_batch_cc3
!
!
   module subroutine omega_cc3_eps_abc_batch_cc3(wf, a, b, c, t_ijk, omega)
!!
!!    Divide W^abc_ijk by -epsilon^abc_ijk to obtain T^abc_ijk
!!    Optional argument omega for jacobian transformations
!!
!!    t^abc_ijk = -W^abc_ijk/epsilon^abc_ijk
!!
!!    Written by Alexander Paul, July 2019
!!    based on omega_cc3_eps_cc3 by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(inout) :: t_ijk
!
      real(dp), optional :: omega
!
   end subroutine omega_cc3_eps_abc_batch_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_calc_abc_batch_cc3(wf, a, b ,c, c_ia, c_ijab, &
                                                                  c_ijk, u_ijk, v_ijk, F_kc, &
                                                                  L_jakb, L_jakc, L_jbkc,    &
                                                                  g_jlka, g_jlkb, g_jlkc,    &
                                                                  g_dbka, g_dcka, g_dckb,    &
                                                                  g_dakb, g_dakc, g_dbkc)
!!
!!    Calculate the L3 amplitudes for fixed indices a,b,c
!!
!!    C^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc)
!!    + sum_l (C_ablk g_iljc - C_abil L_jlkc) - sum_d (C_adjk g_ibdc - C_adij L_dbkc)
!!
!!    Contibutions from outer products:
!!    P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc)
!!
!!    Contibutions from matrix multiplication:
!!      sum_l P^abc_ijk (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d P^abc_ijk (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
!!    Written by Alexander Paul, July 2019
!!    based on ojacobian_transpose_cc3_c3_calc written by A.Paul and Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: c_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: c_ijab
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: c_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: u_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: v_ijk
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: F_kc
!
!     L_ibjc ordered ij,bc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jakb
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jakc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jbkc
!
!     g_jlkc ordered jkl,c
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlka
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlkb
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlkc
!
!     g_dbkc ordered kd,bc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_dbka
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_dcka
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_dckb
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_dakb
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_dakc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_dbkc
!
   end subroutine jacobian_transpose_cc3_c3_calc_abc_batch_cc3