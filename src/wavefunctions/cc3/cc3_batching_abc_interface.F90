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
   module subroutine prepare_cc3_integrals_t3_abc_batch_cc3(wf)
!!
!!    Prepare integral files t3 amplitudes in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    (bd|ck) ordered as dk,bc
!!    (lj|ck) ordered as ljk,c
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_integrals_t3_abc_batch_cc3
!
!
   module subroutine prepare_cc3_integrals_R3_abc_batch_cc3(wf, R_ai)
!!
!!    Prepare integral files R3 amplitudes in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck')   ordered as dk,bc
!!    g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k)   ordered as ljk,c
!!
!!    NB: The integrals (bd|ck) and (lj|ck) constructed in 
!!        prepare_cc3_integrals_t3_abc_batch_cc3 are also needed
!!
!!    Based on construct_c1_integrals_cc3 
!!    written by Rolf H. Myhre and Alexander C. Paul
!!
      implicit none
!
      class(cc3) :: wf 
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine prepare_cc3_integrals_R3_abc_batch_cc3
!
!
   module subroutine prepare_cc3_integrals_L3_abc_batch_cc3(wf)
!!
!!    Prepare integral files for L3 amplitudes in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    (db|kc) ordered as dk,bc
!!    (jl|kc) ordered as ljk,c
!!    L_jbkc  ordered as jk,bc
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_integrals_L3_abc_batch_cc3
!
!
   module subroutine omega_cc3_W_calc_abc_batch_cc3(wf, a, b, c,           &
                                                   t_ijk, u_ijk, t_ijab,   &
                                                   g_ljak, g_ljbk, g_ljck, &
                                                   g_bdak, g_cdak, g_cdbk, &
                                                   g_adbk, g_adck, g_bdck, &
                                                   overwrite)
!!
!!    Omega CC3 intermediate W_ijk for fixed a,b,c
!!    Written by Rolf H. Myhre and Alexander C. Paul July 2019
!!
!!    Contributions to W
!!     W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il(lj|ck))
!!
!!    based on omega_cc3_W_calc written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: a, b, c
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_ijk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in)   :: t_ijab
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljak
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljbk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljck
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_bdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_cdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_cdbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_adbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_adck
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_bdck
      logical, optional, intent(in) :: overwrite
!
   end subroutine omega_cc3_W_calc_abc_batch_cc3
!
!
   module subroutine omega_cc3_eps_abc_batch_cc3(wf, a, b, c, t_ijk, omega)
!!
!!    Omega CC3 epsilon denominator in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019    
!!
!!    Divide W^abc_ijk by -epsilon^abc_ijk to obtain T^abc_ijk
!!    Optional argument omega for jacobian transformations
!!
!!    t^abc_ijk = -W^abc_ijk/epsilon^abc_ijk
!!
!!    based on omega_cc3_eps_cc3 by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: a, b, c
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(inout) :: t_ijk
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
!!    Written by Alexander C. Paul, July 2019
!!
!!    C^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + C_abij*F_kc - L_abik*F_jc)
!!    + sum_l (C_ablk g_iljc - C_abil L_jlkc) - sum_d (C_adjk g_ibdc - C_adij L_dbkc)
!!
!!    Contibutions from outer products:
!!    P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + C_abij*F_kc - C_abik*F_jc)
!!
!!    Contibutions from matrix multiplication:
!!      sum_l P^abc_ijk (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d P^abc_ijk (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
!!    based on jacobian_transpose_cc3_c3_calc written by A.Paul and Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: a, b, c
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: c_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: c_ijab
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: c_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: u_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: v_ijk
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: F_kc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jakb
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jakc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jbkc
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlka
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlkb
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlkc
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dbka
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dcka
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dckb
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dakb
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dakc
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dbkc
!
   end subroutine jacobian_transpose_cc3_c3_calc_abc_batch_cc3
!
!
   module subroutine get_triples_cvs_projector_abc_batch_cc3(wf, projector_ijk)
!!
!!    Get triples cvs projector for fixed a,b,c
!!    Written by Alexander C. Paul and Rolf H. Myhre, September 2019
!!
!!    Set up projector for cvs for the triples amplitudes 
!!    if we batch over the virtual indices
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: projector_ijk
!
   end subroutine get_triples_cvs_projector_abc_batch_cc3
