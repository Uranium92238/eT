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
   module subroutine construct_omega_cc3(wf, omega)
!!
!!    Construct omega (CC3)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Directs the construction of the projection vector < mu| exp(-T) H exp(T) |R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2)
!!
!!    CC3 Omega terms
!!    Written by Rolf H. Myhre, January 2019
!!
!!    t_mu3 = -< mu3|{U,T2}|HF > (epsilon_mu3)^-1
!!
!!    omega_mu1 += < mu1|[H,T3]|HF >
!!
!!    omega_mu2 += < mu2|[H,T3]|HF >
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Omega and store on disk
!!    (bd|ck) ordered as dbc,k
!!    (db|kc) ordered as bcd,k
!!    (lj|ck) ordered as lc,jk
!!    (jl|kc) ordered as cl,jk
!!    (jb|kc) stored as L_jbkc = 2g_jbkc - g_jckb ordered as bc,jk
!!
!!    Written by Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine omega_cc3_integrals_cc3
!
!
   module subroutine omega_cc3_a_n6_cc3(wf, i, j, k, t_abc, u_abc, &
                                        omega1, omega2, F_ov_ck,   &
                                        L_jbic, L_kbic, L_kbjc,    &
                                        L_ibjc, L_ibkc, L_jbkc)
!!
!!    omega_cc3_a_n6
!!
!!    Calculate the triples contribution to omega1 and
!!    the Fock contribution to omega2 scaling as n^6
!!
!!    Written by Rolf H. Myhre, January 2019
!!
!!    omega^a_i += sum_bcjk (t^abc_ijk - t^cba_ijk)*L_jbkc
!!    
!!    omega^ab_ij += P^{ab}_{ij}sum_ck (t^abc_ijk - t^cba_ijk)*F_kc
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_ov_ck
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!
   end subroutine omega_cc3_a_n6_cc3
!
!
   module subroutine omega_cc3_a_n7_cc3(wf, i, j, k, t_abc, u_abc, v_abc, omega2, &
                                        g_dbic, g_dbjc, g_dbkc, &
                                        g_jlic, g_klic, g_kljc, g_iljc, g_ilkc, g_jlkc)
!!
!!    omega_cc3_a_n7
!!
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Calculate the triples contribution to omega2. Scaling as n^7
!!
!!    omega_abli -= P^ab_li sum_cjk(2t^bac_ijk - t^bca_ijk - t^cab_ijk)*g_jlkc
!!
!!    omega_adij -= P^ad_ij sum_bck(2t^abc_ijk - t^cba_ijk - t^acb_ijk)*g_dbkc
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_jlkc
!
   end subroutine omega_cc3_a_n7_cc3
