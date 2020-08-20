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
   module subroutine prepare_for_jacobian_cc3(wf)
!!
!!    Prepare for jacobian
!!    Written by Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_cc3
!
!
   module subroutine prepare_for_jacobian_transpose_cc3(wf)
!!
!!    Prepare for jacobian transpose transformation
!!    Written by Rolf H. Myhre, April 2019
!!
!!    Modified by Tor S. Haugland
!!
!!    Also prepares CCSD jacobian_transpose.
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_transpose_cc3
!
!
   module subroutine prepare_cc3_g_lbkc_t_file_cc3(wf)
!!
!!    Prepare ovov-integral 
!!    written by Rolf H. Myhre and Alexander C. Paul, April 2019
!!
!!    (lb|kc) ordered as bcl,k
!!
!!    only needed in the construction of the intermediates
!!    for the CC3 jacobian transformations and store on disk
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_g_lbkc_t_file_cc3
!
!
   module subroutine prepare_cc3_jacobian_intermediates_cc3(wf)
!!
!!
!!    Prepare intermediates for jacobian CC3 transformations
!!    written by Rolf H. Myhre and Alexander C. Paul, April 2019
!!
!!    Construct X_abdi and Y_akil needed in CC3 jacobian transpose and store on disk
!!    For that: construct t^abc_ijk in single batches of ijk 
!!    and contract with the respective integrals
!!
!!    t^abc_ijk = - (ε^abc_ijk)^-1 P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il(lj|ck))
!!
!!    X_abid = - sum_jck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!    X_ajil = - sum_bck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_jacobian_intermediates_cc3
!
!
   module subroutine construct_x_intermediates_cc3(wf, i, j, k, t_abc, u_abc, v_abc, x_alij,      &
                                                   X_abdi, X_abdj, X_abdk, g_lbic, g_lbjc, g_lbkc)
!!
!!    Construct X intermediates
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Constructs the intermediates X_vvvo and X_vooo
!!    used to compute the contributions to sigma_ai
!!
!!    X_abdi = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_kcjd
!!    X_ajil = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_lbkc
!!
!!    g_lbic, g_lbjc, g_lbkc can be used for g_pcqd as well: 
!!    The p(i,j,k) can be set in dgemm and q(i,j,k) is defined by the array used
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: v_abc
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out)  :: X_alij ! ordered alij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: X_abdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: X_abdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: X_abdk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)           :: g_lbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)           :: g_lbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)           :: g_lbkc
!
   end subroutine construct_x_intermediates_cc3
!
!
   module subroutine sort_x_to_abid_and_write_cc3(wf)
!!
!!    Read in intermediate X_abdi from file, resort to X_abid and write to file again
!!
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine sort_x_to_abid_and_write_cc3
