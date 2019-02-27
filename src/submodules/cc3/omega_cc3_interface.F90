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
   module subroutine construct_omega_cc3(wf, omega)
!!
!!   Construct omega (CC3)
!!   Alex C. Paul and Rolf H. Myhre 2018
!!
     implicit none
!
     class(cc3), intent(inout) :: wf
!
     real(dp), dimension(wf%n_gs_amplitudes, 1), intent(inout) :: omega
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2)
!!
!!    CC3 Omega terms
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Omega and store on disk
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine omega_cc3_integrals_cc3
!
!
   module subroutine omega_cc3_vvv_reader_cc3(wf,batch_x,g_bdcx,g_dbxc)
!!
!!    Read in the intgrals needed in the current batches
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_bdcx
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_dbxc
!
   end subroutine omega_cc3_vvv_reader_cc3
!
!
   module subroutine omega_cc3_ov_vv_reader_cc3(wf,batch_y,batch_x,g_lycx,g_ylxc,L_ybxc)
!!
!!    Read the ljck, jlkc and jbkc integrals needed in the current batches
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      real(dp), dimension(wf%n_o,wf%n_v,batch_y%length,batch_x%length), intent(out) :: g_lycx
      real(dp), dimension(wf%n_v,wf%n_o,batch_y%length,batch_x%length), intent(out) :: g_ylxc
      real(dp), dimension(wf%n_v,wf%n_v,batch_y%length,batch_x%length), intent(out) :: L_ybxc
!
   end subroutine omega_cc3_ov_vv_reader_cc3
!
!
   module subroutine omega_cc3_W_calc_cc3(wf, i, j, k, t_abc, u_abc, t_abji, &
                                          g_bdci, g_bdcj, g_bdck, &
                                          g_ljci, g_lkci, g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Read the ljck, jlkc and jbkc integrals needed in the current batches
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abji
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci 
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj 
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck 
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljci 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkci 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkcj 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_licj 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lick 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljck 
!
!
   end subroutine omega_cc3_W_calc_cc3
!
!
   module subroutine omega_cc3_eps_cc3(wf, i, j, k, t_abc)
!!
!!    Divide W^abc_ijk with -epsilon^abc_ijk to obtain T^abc_ijk
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v), intent(inout) :: t_abc
!
!
   end subroutine omega_cc3_eps_cc3
!
!
   module subroutine omega_cc3_omega1_cc3(wf, i, j, k, t_abc, u_abc, omega1, omega2, F_kc, &
                                          L_jbic, L_kbic, L_kbjc, L_ibjc, L_ibkc, L_jbkc)
!!
!!    Calculate the triples contribution to omega1 and
!!    the Fock contribution to omega2
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_kc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!
   end subroutine omega_cc3_omega1_cc3
!
!
   module subroutine omega_cc3_omega2_cc3(wf, i, j, k, t_abc, u_abc, v_abc, omega2, &
                                          g_dbic, g_dbjc, g_dbkc, &
                                          g_jlic, g_klic, g_kljc, g_iljc, g_ilkc, g_jlkc)
!!
!!    Calculate the triples contribution to omega1 and
!!    the Fock contribution to omega2
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbic 
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbjc 
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbkc 
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jlic 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_klic 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kljc 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_iljc 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ilkc 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jlkc 
!
   end subroutine omega_cc3_omega2_cc3
!
!
