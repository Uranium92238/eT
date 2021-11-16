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
submodule (ccsd_class) jacobian_transpose_ccsd
!
!!
!!    Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * b_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    Transferred to the current eT program from the first version
!!    of eT by Andreas Skeidsvoll and Sarai D. Folkestad, 2018.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_transpose_ccsd(wf)
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
      type(timings) :: prep_timer
!
      real(dp), dimension(:,:,:,:), allocatable :: t_vovo
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ovov, g_ooov
      real(dp), dimension(:,:,:,:), allocatable :: L_ooov
      real(dp), dimension(:,:,:,:), allocatable :: L_vovo
!
      prep_timer = timings("Time preparing for CCSD Jacobian transpose", pl='normal')
      call prep_timer%turn_on()
!
!     Form t_vovo
!
      call mem%alloc(t_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_vovo, wf%n_t1)
!
      call wf%save_jacobian_transpose_d1_intermediates(t_vovo)
      call wf%save_jacobian_transpose_g1_intermediates(t_vovo)
!
!     Prepare integrals for f1 and e1 intermediates
!     g_ooov used in f1 and L_ooov used in e1
!
      call mem%alloc(g_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ooov', g_ooov)
!
      call wf%save_jacobian_transpose_f1_intermediates(t_vovo, g_ooov)
!
      call mem%alloc(L_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call copy_and_scale(two, g_ooov, L_ooov, wf%n_o**3 * wf%n_v)
      call add_3214_to_1234(-one, g_ooov, L_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%save_jacobian_transpose_e1_intermediates(t_vovo, L_ooov)
!
      call mem%dealloc(L_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%save_jacobian_transpose_a1_intermediates(wf%u_aibj)
!
!     Construct g_ovov and make g2 and i2 intermediates
!
      call mem%alloc(g_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_ovov)
!
      call wf%save_jacobian_transpose_g2_intermediates(t_vovo, g_ovov)
      call wf%save_jacobian_transpose_i2_intermediates(t_vovo, g_ovov)
!
!     Construct L_ovov ordered as L_vovo and make d2, e2, and f2 intermediates
!
      call mem%alloc(L_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_vovo, wf%n_o**2 * wf%n_v**2)
!
      call add_2143_to_1234(two, g_ovov, L_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o )
      call add_2341_to_1234(-one, g_ovov, L_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%save_jacobian_transpose_d2_intermediates(wf%u_aibj, L_vovo)
      call wf%save_jacobian_transpose_e2_oo_intermediate(t_vovo, L_vovo)
      call wf%save_jacobian_transpose_e2_vv_intermediate(t_vovo, L_vovo)
      call wf%save_jacobian_transpose_f2_intermediates(t_vovo, L_vovo)
!
      call mem%dealloc(t_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call prep_timer%turn_off()
!
   end subroutine prepare_for_jacobian_transpose_ccsd
!
!
   module subroutine jacobian_transpose_transformation_ccsd(wf, b, sigma)
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
!
!     Incoming vector b
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: b
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: sigma
!
!     Local unpacked and reordered vectors
!
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj ! Unpacked b_aibj
      real(dp), dimension(:,:,:,:), allocatable :: b_abij ! b_aibj, reordered
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_abij     ! sigma_aibj, reordered
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD', pl='normal')
      call timer%turn_on()
!
      call zero_array(sigma, wf%n_t1 + wf%n_t2)
!
!     Calculate and add the CCS contributions to the
!     singles transformed vector
!
      call wf%ccs%jacobian_transpose_transformation(b(1 : wf%n_t1), sigma(1 : wf%n_t1))
!
!     Calculate and add the CCSD contributions to the
!     singles transformed vector
!
      call wf%jacobian_transpose_doubles_a1(sigma(1 : wf%n_t1), b(1 : wf%n_t1), wf%u_aibj)
!
      call mem%alloc(b_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(b(wf%n_t1+1:wf%n_t1 + wf%n_t2), b_aibj, wf%n_v*wf%n_o)
!
      call wf%jacobian_transpose_doubles_b1(sigma(1 : wf%n_t1), b_aibj)
!
      call wf%jacobian_transpose_ccsd_d1(sigma(1 : wf%n_t1), b_aibj)
      call wf%jacobian_transpose_ccsd_e1_o3v(sigma(1 : wf%n_t1), b_aibj)
      call wf%jacobian_transpose_ccsd_e1_v3o(sigma(1 : wf%n_t1), b_aibj)
      call wf%jacobian_transpose_ccsd_f1(sigma(1 : wf%n_t1), b_aibj)
      call wf%jacobian_transpose_ccsd_g1(sigma(1 : wf%n_t1), b_aibj)
!
!     Add the CCSD contributions to the doubles vector arising from
!     the incoming singles vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(sigma_aibj, (wf%n_o*wf%n_v)**2)
!
      call wf%jacobian_transpose_doubles_a2(sigma_aibj, b(1 : wf%n_t1))
!
      call wf%jacobian_transpose_ccsd_b2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_c2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_d2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_e2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_f2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_g2(sigma_aibj, b_aibj)
!
!     Last two terms are already symmetric (h2 and i2). Perform the symmetrization
!     sigma_aibj = P_ij^ab sigma_aibj now, for convenience
!
      call symmetric_sum(sigma_aibj, wf%n_v*wf%n_o)
!
!     In preparation for last two terms, reorder b_aibj to b_abij
!
      call mem%alloc(b_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(b_aibj, b_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(b_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(sigma_aibj, sigma_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add the last two terms
!
      call wf%jacobian_transpose_ccsd_i2(sigma_abij, b_abij)
!
      call mem%dealloc(b_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call packin(sigma(wf%n_t1+1:), sigma_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%omega_ccsd_a2(sigma(wf%n_t1+1:), b(wf%n_t1+1:), right=.false., diagonal_factor=two)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_transformation_ccsd
!
!
   module subroutine save_jacobian_transpose_d1_intermediates(wf, t_aibj)
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
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ilck
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD D1 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Form intermediate X_ilck = sum_d F_i_d t_d_lck
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  (wf%n_o)**2 * wf%n_v, &
                  wf%n_v,               &
                  one,                  &
                  wf%fock_ia,           & ! F_i_d
                  wf%n_o,               &
                  t_aibj,               & ! t_d_lck
                  wf%n_v,               &
                  zero,                 &
                  X_ilck,               & ! X_i_lck
                  wf%n_o)
!
!     Write X_ilck to file
!
      wf%jacobian_transpose_d1_intermediate = sequential_file('jacobian_transpose_d1_intermediate')
      call wf%jacobian_transpose_d1_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_d1_intermediate%write_(X_ilck, wf%n_o**3 * wf%n_v)
!
      call wf%jacobian_transpose_d1_intermediate%close_()
!
      call mem%dealloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_d1_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_d1(wf, sigma_ai, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ilck ! intermediate, term 1
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl
      real(dp), dimension(:,:), allocatable     :: X_li   ! intermediate, term 2
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD D1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_ckdl b_ckal F_id t_kl^cd ::
!
!     Form intermediate X_ilck = sum_d F_id t_dlck from file
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_d1_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_d1_intermediate%read_(X_ilck, wf%n_o**3 * wf%n_v)
      call wf%jacobian_transpose_d1_intermediate%close_()
!
!     Add intermediate: sigma_ai += - sum_lck b_a_lck X_i_lck
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o)**2 * wf%n_v, &
                  -one,                 &
                  b_aibj,               & ! b_a_lck
                  wf%n_v,               &
                  X_ilck,               & ! X_i_lck
                  wf%n_o,               &
                  one,                  &
                  sigma_ai,             & ! sigma_a_i
                  wf%n_v)
!
      call mem%dealloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. - sum_ckdl b_ckdi F_la t_kl^cd
!
!     Form the intermediate X_li = sum_ckd t_ckd_l b_ckd_i
!
      call mem%alloc(X_li, wf%n_o, wf%n_o)
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckdl, wf%n_v * wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_ckdl,               & ! t_ckd_l
                  (wf%n_o)*(wf%n_v)**2, &
                  b_aibj,               & ! b_ckd_i
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_li,                 & ! X_l_i
                  wf%n_o)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add intermediate: sigma_ai += - sum_l F_la X_li = - sum_l F_i_a^T(a,l) X_li(l,i)
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%fock_ia, & ! F_i_a
                  wf%n_o,     &
                  X_li,       & ! X_l_i
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   & ! sigma_a_i
                  wf%n_v)
!
      call mem%dealloc(X_li, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_d1
!
!
   module subroutine save_jacobian_transpose_e1_intermediates(wf, t_aibj, L_ilmd)
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
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_v), intent(in) :: L_ilmd
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ildm
      real(dp), dimension(:,:,:,:), allocatable :: X_ilck
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD E1 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Sort L_ilmd into L_ildm
!
      call mem%alloc(L_ildm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1243(L_ilmd, L_ildm, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Form the intermediate X_ilck = sum_md L_ilmd t_mk^dc
!                                  = sum_md L_il_dm t_dm_ck
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ildm,            & ! L_il_dm
                  (wf%n_o)**2,       &
                  t_aibj,            & ! t_dm_ck
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ilck,            & ! X_il_ck
                  (wf%n_o)**2)
!
      call mem%dealloc(L_ildm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Write X_ilck to file
!
      wf%jacobian_transpose_e1_intermediate = sequential_file('jacobian_transpose_e1_intermediate')
      call wf%jacobian_transpose_e1_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_e1_intermediate%write_(X_ilck, wf%n_o**3 * wf%n_v)
!
      call wf%jacobian_transpose_e1_intermediate%close_()
!
      call mem%dealloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_e1_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_e1_o3v(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD E1 o3v
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E1 term,
!!
!!      -sum_ckdlm (b_ckal L_ilmd t_km^cd + b_ckdl L_mlia t_km^cd)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
!!    Modified by Tor S. Haugland, Nov 2019
!!    Reads intermediate for term 3 from the file 'jacobian_transpose_e1_intermediate'.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdm, g_mlia
!
      real(dp), dimension(:,:), allocatable :: X_ml
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ilck, L_aiml
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD E1 o3v', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_ckdlm b_ckal L_ilmd t_km^cd ::
!
!     Read the intermediate X_ilck = sum_md L_ilmd t_mk^dc = sum_md L_il_dm t_dm_ck
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_e1_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_e1_intermediate%read_(X_ilck, wf%n_o**3 * wf%n_v)
      call wf%jacobian_transpose_e1_intermediate%close_()
!
!     Add = - sum_ckl b_a_lck X_i_lck^T
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_aibj,               & ! b_a_lck
                  wf%n_v,               &
                  X_ilck,               & ! X_i_lck
                  wf%n_o,               &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. - sum_ckdlm b_ckdl L_mlia t_km^cd ::
!
!     Form the intermediate X_ml = sum_ckd t_km^cd b_ckdl
!                                = sum_ckd t_ckd_m^T b_ckd_l
!
      call mem%alloc(t_ckdm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckdm, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(X_ml, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_ckdm,               & ! t_ckd_m
                  (wf%n_o)*(wf%n_v)**2, &
                  b_aibj,               & ! b_aib_j
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_ml,                 &
                  wf%n_o)
!
      call mem%dealloc(t_ckdm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form g_mlia = g_mlia
!
      call mem%alloc(g_mlia, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ooov', g_mlia)
!
!     Form L_aiml = L_mlia = 2 * g_mlia - g_mail
!                          = 2 * g_mlia - g_ilma
!
      call mem%alloc(L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(L_aiml, wf%n_v*wf%n_o**3)
!
      call add_3421_to_1234(two, g_mlia, L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call add_2431_to_1234(-one, g_mlia, L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_mlia, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Add - sum_ckdlm b_ckdl L_mlia t_km^cd = - sum_lm L_ai_ml X_ml
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)**2,       &
                  -one,              &
                  L_aiml,            & ! L_ai_ml
                  (wf%n_o)*(wf%n_v), &
                  X_ml,              & ! X_ml
                  (wf%n_o)**2,       &
                  one,               &
                  sigma_ai,          & ! sigma_ai
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_ml, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_e1_o3v
!
!
   module subroutine jacobian_transpose_ccsd_e1_v3o(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD E1 v3o
!!    Written by Regina Matveeva and Alexander C. Paul, Sep 2021
!!    Based on jacobian_transpose_ccsd_e1 by Eirik F. Kjønstad, Sarai D. Folkestad
!!
!!    Calculates the E1 term
!!    sum_ckdle (b_ckdi L_dale t_kl^ce + b_ckdl L_deia t_kl^ce)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_elck
!
      real(dp), dimension(:,:,:,:), allocatable :: X_eldi, X_lide
!
      real(dp), dimension(:,:,:), allocatable :: L_Jvv, Y_Jdi, Y_Jvo
      real(dp), dimension(:,:,:), allocatable :: L_Jov, Y_Jid
      real(dp), dimension(:,:,:), allocatable :: Y_Jli, Y_liJ
!
      real(dp), dimension(:,:), allocatable :: X_de, X_ia
      real(dp), dimension(:), allocatable ::   Y_J
!
      integer :: req0, req1, batch
      type(batching_index) :: batch_v
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD E1 v3o', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(t_elck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_elck, wf%n_o*wf%n_v)
!
!     Form the intermediates from t2 and b2
!     X_eldi = sum_ck t_lk^ec b_ckdi = sum_ck t_el_ck b_ck_di
!     X_de = sum_ckl b_ckdl t_kl^ce = sum_ckl b_d_lck t_e_lck^T
!
      call mem%alloc(X_de, wf%n_v, wf%n_v)
!
      call mem%alloc(X_eldi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','T',           &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  b_aibj,           & ! b_d_lck = b_dlck = b_ckdl
                  wf%n_v,           &
                  t_elck,           & ! t_e_lck
                  wf%n_v,           &
                  zero,             &
                  X_de,             &
                  wf%n_v)
!
      call dgemm('N','N',        &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  t_elck,        & ! t_el_ck
                  wf%n_o*wf%n_v, &
                  b_aibj,        & ! b_ck_di
                  wf%n_o*wf%n_v, &
                  zero,          &
                  X_eldi,        & ! X_el_di
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(t_elck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add coulomb contribution from term 1 and exchange from term 2:
!     Y_J_di = sum_el 2 X_eldi L_J_le - sum_e X_de L_J_ie
!
!     Needs to be contracted with L_J_da
!     sigma_ai += Y_J_di L_J_da
!
      call mem%alloc(L_Jov, wf%eri%n_J, wf%n_o, wf%n_v)
      call wf%eri%get_cholesky_t1(L_Jov, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(Y_Jid, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call dgemm('N','T',               &
                  wf%eri%n_J*wf%n_o,    &
                  wf%n_v,               &
                  wf%n_v,               &
                  -one,                 &
                  L_Jov,                & ! L_Ji_e
                  wf%eri%n_J*wf%n_o,    &
                  X_de,                 & ! X_d_e
                  wf%n_v,               &
                  zero,                 &
                  Y_Jid,                &
                  wf%eri%n_J*wf%n_o)
!
      call mem%alloc(Y_Jdi, wf%eri%n_J, wf%n_v, wf%n_o)
      call sort_123_to_132(Y_Jid, Y_Jdi, wf%eri%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(Y_Jid, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call mem%alloc(Y_Jvo, wf%eri%n_J, wf%n_v, wf%n_o)
      call sort_123_to_132(L_Jov, Y_Jvo, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call dgemm('N','N',        &
                  wf%eri%n_J,    &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  two,           &
                  Y_Jvo,         & ! L_J_el
                  wf%eri%n_J,    &
                  X_eldi,        & ! X_el_di
                  wf%n_v*wf%n_o, &
                  one,           &
                  Y_Jdi,         & ! Y_J_di
                  wf%eri%n_J)
!
      call mem%dealloc(Y_Jvo, wf%eri%n_J, wf%n_v, wf%n_o)
!
!     Resort X_eldi for contraction with L_J_de in batches over e
!
      call mem%alloc(X_lide, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2431(X_eldi, X_lide, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_eldi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_liJ, wf%n_o, wf%n_o, wf%eri%n_J)
      call zero_array(Y_liJ, wf%n_o**2*wf%eri%n_J)
!
      call mem%alloc(X_ia, wf%n_o, wf%n_v)
      call zero_array(X_ia, wf%n_o*wf%n_v)
!
      call mem%alloc(Y_J, wf%eri%n_J)
      call zero_array(Y_J, wf%eri%n_J)
!
!     Contractions with L_J_vv
!
      req0 = 0
      req1 = wf%n_v*wf%eri%n_J
!
      batch_v = batching_index(wf%n_v)
      call mem%batch_setup(batch_v, req0, req1, 'jacobian_transpose_ccsd_e1_v3o')
!
      call mem%alloc(L_Jvv, wf%eri%n_J, wf%n_v, batch_v%max_length)
!
      do batch = 1, batch_v%num_batches
!
         call batch_v%determine_limits(batch)
!
         call wf%eri%get_cholesky_t1(L_Jvv, wf%n_o + 1, wf%n_mo, &
                                     wf%n_o + batch_v%first,     &
                                     wf%n_o + batch_v%get_last())
!
         call dgemm('T','N',                &
                    wf%n_o,                 &
                    batch_v%length,         &
                    wf%eri%n_J*wf%n_v,      &
                    one,                    &
                    Y_Jdi,                  & ! Y_Jd_i
                    wf%eri%n_J*wf%n_v,      &
                    L_Jvv,                  & ! L_Jd_a
                    wf%eri%n_J*wf%n_v,      &
                    one,                    &
                    X_ia(:,batch_v%first:), &
                    wf%n_o)
!
         call dgemm('N','T',                       &
                     wf%n_o**2,                    &
                     wf%eri%n_J,                   &
                     wf%n_v*batch_v%length,        &
                     one,                          &
                     X_lide(:,:,:,batch_v%first:), & ! X_li_de
                     wf%n_o**2,                    &
                     L_Jvv,                        & ! L_J_de
                     wf%eri%n_J,                   &
                     one,                          &
                     Y_liJ,                        & ! Y_li_J
                     wf%n_o**2)
!
         call dgemv('N',                    &
                    wf%eri%n_J,             &
                    wf%n_v*batch_v%length,  &
                    two,                    &
                    L_Jvv,                  & ! L_J_de
                    wf%eri%n_J,             &
                    X_de(:,batch_v%first:), & ! X_de
                    1,                      &
                    one,                    &
                    Y_J, 1)                   ! Y_J
!
      end do
!
      call mem%dealloc(L_Jvv, wf%eri%n_J, wf%n_v, batch_v%max_length)
!
      call mem%batch_finalize
!
      call mem%dealloc(X_lide, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(X_de, wf%n_v, wf%n_v)
      call mem%dealloc(Y_Jdi, wf%eri%n_J, wf%n_v, wf%n_o)
!
!     Term 2: Coulomb: sigma_ai += 2 L_J_de X_de L_J_ia = Y_J L_J_ia
!
      call dgemv('T',            &
                  wf%eri%n_J,    &
                  wf%n_o*wf%n_v, &
                  one,           &
                  L_Jov,         & ! L_J_ia
                  wf%eri%n_J,    &
                  Y_J, 1,        & ! Y_J
                  one,           &
                  X_ia, 1)         ! X_ia
!
      call mem%dealloc(Y_J, wf%eri%n_J)
!
      call add_21_to_12(one, X_ia, sigma_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_ia, wf%n_o, wf%n_v)
!
!     Term 1 Exchange: sigma_ai -= X_eldi L_J_de L_J_la = Y_J_li L_J_la
!
      call mem%alloc(Y_Jli, wf%eri%n_J, wf%n_o, wf%n_o)
      call sort_123_to_312(Y_liJ, Y_Jli, wf%n_o, wf%n_o, wf%eri%n_J)
      call mem%dealloc(Y_liJ, wf%n_o, wf%n_o, wf%eri%n_J)
!
      call dgemm('T','N',              &
                 wf%n_v,               &
                 wf%n_o,               &
                 wf%eri%n_J*wf%n_o,    &
                 -one,                 &
                 L_Jov,                & ! L_Jl_a
                 wf%eri%n_J*wf%n_o,    &
                 Y_Jli,                & ! Y_Jl_i
                 wf%eri%n_J*wf%n_o,    &
                 one,                  &
                 sigma_ai,             &
                 wf%n_v)
!
      call mem%dealloc(L_Jov, wf%eri%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(Y_Jli, wf%eri%n_J, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_e1_v3o
!
!
   module subroutine save_jacobian_transpose_f1_intermediates(wf, t_aibj, g_ikmc)
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
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_v), intent(in) :: g_ikmc
!
      real(dp), dimension(:,:,:,:), allocatable :: t_mcdl
      real(dp), dimension(:,:,:,:), allocatable :: g_limc
      real(dp), dimension(:,:,:,:), allocatable :: X_ikdl
      real(dp), dimension(:,:,:,:), allocatable :: X_kdli
      real(dp), dimension(:,:,:,:), allocatable :: X_lidk
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD F1 intermediate', pl='verbose')
      call timer%turn_on()
!
!     X_ikdl = sum_mc t_lm^cd g_ikmc = t_mcdl g_ikmc
!
!     Order amplitudes as t_mcdl = t_lm^cd
!
      call mem%alloc(t_mcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call sort_1234_to_4132(t_aibj, t_mcdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ikdl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  g_ikmc,            & ! g_ik_mc
                  (wf%n_o)**2,       &
                  t_mcdl,            & ! t_mc_dl
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ikdl,            & ! X_ik_dl
                  (wf%n_o)**2)
!
!     X_lidk = sum_mc t_mk^dc g_mlic = t_mcdk g_limc
!
!     Pretend that g_ikmc is g_mlic and reorder to g_limc
!
      call mem%alloc(g_limc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_2314(g_ikmc, g_limc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Pretend that t_mcdl is t_mcdk (t_lm^cd is t_km^cd)
!
      call mem%alloc(X_lidk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',       &
                  wf%n_o**2,     &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  g_limc,        & ! g_li_mc
                  wf%n_o**2,     &
                  t_mcdl,        & ! t_mc_dk
                  wf%n_o*wf%n_v, &
                  zero,          &
                  X_lidk,        & ! X_li_dk
                  wf%n_o**2)
!
      call mem%dealloc(t_mcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(g_limc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!      Reorder and add the intermediates to X_kdli
!
      call mem%alloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(X_kdli, wf%n_o**3 * wf%n_v)
!
      call add_4123_to_1234(one, X_ikdl, X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call add_3421_to_1234(one, X_lidk, X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_lidk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ikdl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Write X_kdli to file
!
      wf%jacobian_transpose_f1_intermediate = sequential_file('jacobian_transpose_f1_intermediate')
      call wf%jacobian_transpose_f1_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_f1_intermediate%write_(X_kdli, wf%n_o**3 * wf%n_v)
!
      call wf%jacobian_transpose_f1_intermediate%close_()
!
      call mem%dealloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_f1_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_f1(wf, sigma_ai, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cdml
!
      real(dp), dimension(:,:,:,:), allocatable :: b_kicd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_amkl
      real(dp), dimension(:,:,:,:), allocatable :: g_mkla
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kdli
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kiml
      real(dp), dimension(:,:,:,:), allocatable :: X_mkli
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD F1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. and Term 2.
!
!        sum_ckdlm b_akdl t_lm^cd g_ikmc + sum_ckdlm b_akdl t_mk^dc g_mlic
!        = sum_kdl b_akdl X_kdli
!
!     Read the intermediate X_kdli
!
      call mem%alloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%jacobian_transpose_f1_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_f1_intermediate%read_(X_kdli, wf%n_o**3 * wf%n_v)
      call wf%jacobian_transpose_f1_intermediate%close_()
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_aibj,               & ! b_a_kdl
                  wf%n_v,               &
                  X_kdli,               & ! X_kdl_i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 3. sum_ckdlm b_ckdi t_ml^cd g_mkla ::
!
!     X_kiml = sum_cd b_ckdi t_ml^cd
!
      call mem%alloc(t_cdml, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call squareup_and_sort_1234_to_1324(wf%t2, t_cdml, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(b_kicd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call sort_1234_to_2413(b_aibj, b_kicd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form intermediate X_kiml = sum_cd b_ckdi t_ml^cd = sum_cd b_ki_cd t_cd_ml
!
      call mem%alloc(X_kiml, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)**2,       &
                  (wf%n_v)**2,       &
                  one,               &
                  b_kicd,            & ! b_ki_cd
                  (wf%n_o)**2,       &
                  t_cdml,            & ! t_cd_ml
                  (wf%n_v)**2,       &
                  zero,              &
                  X_kiml,            & ! X_ki_ml
                  (wf%n_o)**2)
!
      call mem%dealloc(t_cdml, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(b_kicd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_mkla X_ki_ml
!
!     Reorder to X_mkli
!
      call mem%alloc(X_mkli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_3142(X_kiml, X_mkli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_kiml, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     We have g_kdm_i(kdm,i) = g_mkid
!     Reorder to g_amkl(a,m,k,l) = g_mkla = g_kdm_i(kam,l)
!
      call mem%alloc(g_mkla, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ooov', g_mkla)
!
      call mem%alloc(g_amkl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_4123(g_mkla, g_amkl, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(g_mkla, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Add sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_a_mkl X_mkl_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  (wf%n_o)**3, &
                  one,         &
                  g_amkl,      & ! g_a_mkl
                  wf%n_v,      &
                  X_mkli,      & ! X_mkl_i
                  (wf%n_o)**3, &
                  one,         &
                  sigma_ai,    &
                  wf%n_v)
!
      call mem%dealloc(g_amkl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_mkli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_f1
!
!
   module subroutine save_jacobian_transpose_g1_intermediates(wf, t_aibj)
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
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_idkl ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_kdli ! Reordered intermediate, term 1
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cekl
      real(dp), dimension(:,:,:,:), allocatable :: g_icde, g_idce
!
      type(timings), allocatable :: timer
!
!     Batching variables
!
      integer :: current_d_batch = 0
!
      type(batching_index) :: batch_d
!
      integer :: rec1, rec0
!
      timer = timings('Jacobian transpose CCSD G1 intermediate', pl='verbose')
      call timer%turn_on()
!
!     X_idkl = sum_ce t_kl^ce g_icde = sum_ce g_id_ce t_ce_kl
!
!     Reorder to t_cekl = t_kl^ce
!
      call mem%alloc(t_cekl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(t_aibj, t_cekl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_idkl, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(X_idkl, wf%n_v*wf%n_o**3)
!
!     Prepare for batching over d
!
      rec0 = wf%n_v*wf%n_o*wf%eri%n_J
      rec1 = 2*(wf%n_v**2)*(wf%n_o) + wf%n_v*wf%eri%n_J
!
      batch_d = batching_index(wf%n_v)
      call mem%batch_setup(batch_d, rec0, rec1, 'save_jacobian_transpose_g1_intermediates')
!
!     Loop over the d-batches
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_icde, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv', g_icde, first_r=batch_d%first, last_r=batch_d%get_last())
!
!        Reorder to g_id_ce = g_ic_de
!
         call mem%alloc(g_idce, wf%n_o, batch_d%length, wf%n_v, wf%n_v)
!
         call sort_1234_to_1324(g_icde, g_idce, wf%n_o, wf%n_v, (batch_d%length), wf%n_v)
!
         call mem%dealloc(g_icde, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call dgemm('N','N',                      &
                     (wf%n_o)*(batch_d%length),   &
                     (wf%n_o)**2,                 &
                     (wf%n_v)**2,                 &
                     one,                         &
                     g_idce,                      &
                     (wf%n_o)*(batch_d%length),   &
                     t_cekl,                      &
                     (wf%n_v)**2,                 &
                     one,                         &
                     X_idkl(1,batch_d%first,1,1), &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(g_idce, wf%n_o, batch_d%length, wf%n_v, wf%n_v)
!
      enddo ! End of batches over d
!
      call mem%batch_finalize()
!
      call mem%dealloc(t_cekl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     - sum_ckdle b_akdl t_kl^ce g_icde = sum_kdl b_akdl X_id_kl
!
      call mem%alloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_3241(X_idkl, X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_idkl, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Save X_kdli to file
!
      wf%jacobian_transpose_g1_intermediate = sequential_file('jacobian_transpose_g1_intermediate')
      call wf%jacobian_transpose_g1_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_g1_intermediate%write_(X_kdli, wf%n_o**3 * wf%n_v)
!
      call wf%jacobian_transpose_g1_intermediate%close_()
!
      call mem%dealloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_g1_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_g1(wf, sigma_ai, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: b_dicl, t_kecl
!
      real(dp), dimension(:,:,:,:), allocatable :: X_dike, X_kide, X_kdli
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ov, L_J_vv
      real(dp), dimension(:,:,:), allocatable :: Y_J_cl, Y_J_di, Y_J_oo
      real(dp), dimension(:,:),   allocatable :: s_ia
!
      type(batching_index) :: batch_v
      integer :: batch, req0, req1
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD G1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_ckdle b_akdl t_kl^ce g_icde ::
!
      call mem%alloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%jacobian_transpose_g1_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_g1_intermediate%read_(X_kdli, wf%n_o**3 * wf%n_v)
      call wf%jacobian_transpose_g1_intermediate%close_()
!
!     Add - sum_ckdle b_akdl t_kl^ce g_icde = - sum_dkl b_a_kdl X_kdl_i
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v*wf%n_o**2, &
                  -one,             &
                  b_aibj,           & ! b_a_kdl
                  wf%n_v,           &
                  X_kdli,           & ! X_kdl_i
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  sigma_ai,         &
                  wf%n_v)
!
      call mem%dealloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2 sigma^a_i -= sum_ckdle b^cd_il t^ce_kl g_kade
!     :: Term 3 sigma^a_i -= sum_ckdle b^cd_li t^ce_kl g_keda
!
      call mem%alloc(t_kecl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call squareup_and_sort_1234_to_2314(wf%t2, t_kecl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2: X_dike = sum_cl b_cidl t_kl^ce
!             Y_J_ki = X_dike L_J_de
!             sigma_ai += L_J_ka Y_J_ki
!
      call mem%alloc(b_dicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(b_aibj, b_dicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_dike, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','T',       &
                 wf%n_o*wf%n_v, &
                 wf%n_o*wf%n_v, &
                 wf%n_o*wf%n_v, &
                 one,           &
                 b_dicl,        & ! b_di_cl
                 wf%n_o*wf%n_v, &
                 t_kecl,        & ! t_ke_cl
                 wf%n_o*wf%n_v, &
                 zero,          &
                 X_dike,        &
                 wf%n_o*wf%n_v)
!
      call mem%dealloc(b_dicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_kide, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3214(X_dike, X_kide, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_dike, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Term 3: Y_J_cl = sum_ek  t^ce_kl L_J_ke
!             Y_J_di = sum_cdl b^cd_li Y_J_cl
!             sigma_ai += sum_Jd L_J_da Y_J_di
!
      call mem%alloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
      call wf%eri%get_cholesky_t1(L_J_ov, 1, wf%n_o, 1+wf%n_o, wf%n_mo)
!
      call mem%alloc(Y_J_cl, wf%eri%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N','N',        &
                  wf%eri%n_J,    &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  L_J_ov,        & ! L_J_ke
                  wf%eri%n_J,    &
                  t_kecl,        & ! t_ke_cl
                  wf%n_v*wf%n_o, &
                  zero,          &
                  Y_J_cl,        & ! Y_J_cl
                  wf%eri%n_J)
!
      call mem%dealloc(t_kecl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(Y_J_di, wf%eri%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N','N',        &
                  wf%eri%n_J,    &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  Y_J_cl,        & ! Y_J_cl
                  wf%eri%n_J,    &
                  b_aibj,        & ! b_cl_di
                  wf%n_v*wf%n_o, &
                  zero,          &
                  Y_J_di,        & ! Y_J_di
                  wf%eri%n_J)
!
      call mem%dealloc(Y_J_cl, wf%eri%n_J, wf%n_v, wf%n_o)
!
!     Contractions with L_J_vv
!
      call mem%alloc(s_ia, wf%n_o, wf%n_v)
      call zero_array(s_ia, wf%n_o*wf%n_v)
!
      call mem%alloc(Y_J_oo, wf%eri%n_J, wf%n_o, wf%n_o)
      call zero_array(Y_J_oo, wf%eri%n_J*wf%n_o**2)
!
      req0 = 0
      req1 = wf%eri%n_J*wf%n_v
!
      batch_v = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_v, req0, req1, 'jacobian_transpose_g1')
!
      call mem%alloc(L_J_vv, wf%eri%n_J, wf%n_v, batch_v%max_length)
!
      do batch = 1, batch_v%num_batches
!
         call batch_v%determine_limits(batch)
!
         call wf%eri%get_cholesky_t1(L_J_vv, wf%n_o + 1, wf%n_mo, &
                                     wf%n_o + batch_v%first,      &
                                     wf%n_o + batch_v%get_last())
!
         call dgemm('N','T',                       &
                     wf%eri%n_J,                   &
                     wf%n_o**2,                    &
                     wf%n_v*batch_v%length,        &
                     one,                          &
                     L_J_vv,                       & ! L_J_de
                     wf%eri%n_J,                   &
                     X_kide(:,:,:,batch_v%first:), & ! X_ki_de
                     wf%n_o**2,                    &
                     one,                          &
                     Y_J_oo,                       & ! Y_J_ki
                     wf%eri%n_J)
!
         call dgemm('T','N',                &
                    wf%n_o,                 &
                    batch_v%length,         &
                    wf%eri%n_J*wf%n_v,      &
                    one,                    &
                    Y_J_di,                 & ! Y_Jd_i
                    wf%eri%n_J*wf%n_v,      &
                    L_J_vv,                 & ! L_Jd_a
                    wf%eri%n_J*wf%n_v,      &
                    one,                    &
                    s_ia(:,batch_v%first:), &
                    wf%n_o)
!
      end do
!
      call mem%dealloc(L_J_vv, wf%eri%n_J, wf%n_v, batch_v%max_length)
      call mem%batch_finalize()
!
      call mem%dealloc(Y_J_di, wf%eri%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(X_kide, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call add_21_to_12(-one, s_ia, sigma_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(s_ia, wf%n_o, wf%n_v)
!
!     Final contraction for Term2
!
      call dgemm('T','N',              &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%eri%n_J*wf%n_o,   &
                  -one,                &
                  L_J_ov,              & ! L_Jk_a
                  wf%eri%n_J*wf%n_o,   &
                  Y_J_oo,              & ! L_Jk_i
                  wf%eri%n_J*wf%n_o,   &
                  one,                 &
                  sigma_ai,            & ! sigma_ai
                  wf%n_v)
!
      call mem%dealloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(Y_J_oo, wf%eri%n_J, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_g1
!
!
   module subroutine jacobian_transpose_ccsd_b2(wf, sigma_aibj, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: b_aijc ! b_aicj
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aijb ! sigma_aibj contribution
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbj ! g_ckjb & g_cbjk
!
      real(dp), dimension(:,:,:,:), allocatable :: g_cbjk_restricted ! g_cbjk, batch over b
!
      real(dp), dimension(:,:,:), allocatable :: L_J_vo, L_J_ov, W_vo_J
!
      integer :: k, j, b, c
!
!     Batching variables
!
      integer :: rec0, rec1
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_b
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD B2', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. sum_c b_aicj F_cb ::
!
!     Reorder to b_aijc = b_aicj
!
      call mem%alloc(b_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(b_aijc, (wf%n_o*wf%n_v)**2)
!
      call sort_1234_to_1243(b_aibj, b_aijc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate and add sum_c b_aicj F_cb = sum_c b_aij_c F_c_b
!
      call mem%alloc(sigma_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_v,               &
                  one,                  &
                  b_aijc,               & ! b_aij_c
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%fock_ab,           & ! F_c_b
                  wf%n_v,               &
                  zero,                 &
                  sigma_aijb,           &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(b_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call add_1243_to_1234(one, sigma_aijb, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 2. - sum_k b_aibk F_jk ::
!
!     - sum_k b_aibk F_jk = - sum_k b_aib_k F_jk^T(k,j)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  b_aibj,               & ! b_aib_k
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%fock_ij,           & ! F_j_k
                  wf%n_o,               &
                  one,                  &
                  sigma_aibj,           & ! sigma_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
!     :: Term 3. sum_ck b_aick L_ckjb ::
!
!     2 sum_ck b_aick g_ckjb
!
      call mem%alloc(L_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
!
      call wf%eri%get_cholesky_t1(L_J_vo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call mem%alloc(W_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
      call dgemm('N', 'T',          &
                  wf%n_v * wf%n_o,  &
                  wf%eri%n_J,       &
                  wf%n_v * wf%n_o,  &
                  one,              &
                  b_aibj,           & ! b_ai,ck
                  wf%n_v * wf%n_o,  &
                  L_J_vo,           & ! L_J,ck
                  wf%eri%n_J,       &
                  zero,             &
                  W_vo_J,           &
                  wf%n_v * wf%n_o)
!
!
!     Add 2 * sum_ck b_aick g_ckjb = 2 * W_ai_J L_jb^J
!
      call mem%alloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call wf%eri%get_cholesky_t1(L_J_ov, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call sort_123_to_132(L_J_ov, L_J_vo, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call mem%dealloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',          &
                 wf%n_v * wf%n_o,   &
                 wf%n_v * wf%n_o,   &
                 wf%eri%n_J,        &
                 two,               &
                 W_vo_J,            &
                 wf%n_v * wf%n_o,   &
                 L_J_vo,            &
                 wf%eri%n_J,        &
                 one,               &
                 sigma_aibj,        &
                 wf%n_v * wf%n_o)
!
      call mem%dealloc(L_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(W_vo_J, wf%n_o, wf%n_v, wf%eri%n_J)
!
!     - sum_ck b_aick g_cbjk
!
!     Prepare to batch over b to make g_cb_jk = g_cbjk successively
!
      call mem%alloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(g_ckbj, (wf%n_o*wf%n_v)**2) ! g_cbjk reordered
!
      rec0 = wf%n_o**2*wf%eri%n_J
      rec1 = wf%n_v*wf%eri%n_J  + (wf%n_o**2)*(wf%n_v)
!
!     Initialize batching variable
!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, rec0, rec1, 'jacobian_transpose_ccsd_b2')
!
!     Loop over the number of b batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cb_jk = g_cbjk
!
         call mem%alloc(g_cbjk_restricted, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
         call wf%eri%get_eri_t1('vvoo', g_cbjk_restricted, &
                                first_q=batch_b%first, last_q=batch_b%get_last())
!
!        Place in reordered full space vector and deallocate restricted vector
!
!$omp parallel do schedule(static) private(k,j,b,c)
         do j = 1, wf%n_o
            do b = batch_b%first, batch_b%get_last()
               do k = 1, wf%n_o
                  do c = 1, wf%n_v
!
                     g_ckbj(c,k,b,j) = g_cbjk_restricted(c,b-batch_b%first+1,j,k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(g_cbjk_restricted, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
      enddo ! End of batches over b
!
      call mem%batch_finalize()
!
!     Add  - sum_ck b_aick g_cbjk = - sum_ck b_ai_ck g_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_aibj,            & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            & ! g_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_b2
!

!
   module subroutine jacobian_transpose_ccsd_c2(wf, sigma_aibj, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_cbik
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbi ! g_cbik
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajib ! sigma_aibj contribution
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi ! sigma_aibj contribution
!
      real(dp), dimension(:,:,:,:), allocatable :: b_ajck ! b_akcj
!
      real(dp), dimension(:,:,:), allocatable :: L_J_vo, L_J_ov, W_vo_J
!
      integer :: k, i, b, c
!
!     Batching variables
!
      integer :: rec1, rec0
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_b
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD C2', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_ck b_ajck g_ibck ::
!
      call mem%alloc(L_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
      call wf%eri%get_cholesky_t1(L_J_vo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call mem%alloc(W_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
      call dgemm('N', 'T',          &
                  wf%n_o * wf%n_v,  &
                  wf%eri%n_J,       &
                  wf%n_o * wf%n_v,  &
                  one,              &
                  b_aibj,           & ! b_aj,ck
                  wf%n_o * wf%n_v,  &
                  L_J_vo,           & ! L_J,ck
                  wf%eri%n_J,       &
                  zero,             &
                  W_vo_J,           & ! W_aj,J
                  wf%n_o * wf%n_v)
!
      call mem%dealloc(L_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
!
      call mem%alloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
      call wf%eri%get_cholesky_t1(L_J_ov, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',          &
                  wf%n_o * wf%n_v,  &
                  wf%n_o * wf%n_v,  &
                  wf%eri%n_J,       &
                  -one,             &
                  W_vo_J,           & ! W_aj,J
                  wf%n_o * wf%n_v,  &
                  L_J_ov,           & ! L_J,ib
                  wf%eri%n_J,       &
                  zero,             &
                  sigma_ajib,       &
                  wf%n_o * wf%n_v)
!
      call mem%dealloc(W_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
      call mem%dealloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call add_1423_to_1234(one, sigma_ajib, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 2. - sum_ck b_akcj g_ikcb = - sum_ck b_akcj g_cbik ::
!
!     Make g_ckbi = g_cbik in batches over b
!
      call mem%alloc(g_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(g_ckbi, (wf%n_o*wf%n_v)**2)
!
      rec0 = wf%n_o**2*wf%eri%n_J
      rec1 = wf%n_v*wf%eri%n_J  + (wf%n_o**2)*(wf%n_v)
!
!     Initialize batching variable
!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, rec0, rec1, 'jacobian_transpose_ccsd_c2')
!
!     Loop over the b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cbik
!
         call mem%alloc(g_cbik, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
         call wf%eri%get_eri_t1('vvoo', g_cbik, first_q=batch_b%first, last_q=batch_b%get_last())
!
!        Place in reordered integral g_ckbi = g_cbik
!
!$omp parallel do schedule(static) private(k,i,b,c)
         do k = 1, wf%n_o
            do i = 1, wf%n_o
               do b = batch_b%first, batch_b%get_last()
                  do c = 1, wf%n_v
!
                     g_ckbi(c,k,b,i) = g_cbik(c,b-batch_b%first+1,i,k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(g_cbik, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
      enddo ! End of batches over b
!
      call mem%batch_finalize()
!
!     Reorder to b_ajck = b_akcj
!
      call mem%alloc(b_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(b_aibj, b_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form and add - sum_ck b_akcj g_cbik = - sum_ck b_ajck g_ck_bi
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_ajck,            & ! b_aj_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ckbi,            & ! g_ck_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_ajbi,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(b_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_c2
!
!
   module subroutine save_jacobian_transpose_d2_intermediates(wf, u_ckdl, L_dlbj)
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
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: L_dlbj ! Reordered L_jbld
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj ! An intermediate
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD D2 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Form the intermediate X_ckbj = sum_dl u_ck_dl L_dl_bj
!
      call mem%alloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ckdl,            & ! u_ck_dl
                  (wf%n_o)*(wf%n_v), &
                  L_dlbj,            & ! L_dl_bj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ckbj,            &
                  (wf%n_o)*(wf%n_v))
!
!     Save X_ckbj to file
!
      wf%jacobian_transpose_d2_intermediate = sequential_file('jacobian_transpose_d2_intermediate')
      call wf%jacobian_transpose_d2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_d2_intermediate%write_(X_ckbj, wf%n_o**2 * wf%n_v**2)
!
      call wf%jacobian_transpose_d2_intermediate%close_()
!
      call mem%dealloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_d2_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_d2(wf, sigma_aibj, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj ! An intermediate
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD D2', pl='verbose')
      call timer%turn_on()
!
!     Read the intermediate X_ckbj = sum_dl t_ck_dl L_dl_bj
!
      call mem%alloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_d2_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_d2_intermediate%read_(X_ckbj, wf%n_o**2 * wf%n_v**2)
      call wf%jacobian_transpose_d2_intermediate%close_()
!
!     Add sum_ckdl b_aick L_jbld t_kl^cd = sum_ck b_ai_ck X_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  b_aibj,            & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_d2
!
!
   module subroutine jacobian_transpose_ccsd_e2(wf, sigma_aibj, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:), allocatable     :: X_jl   ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable     :: X_cb   ! An intermediate, term 3
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD E2', pl='verbose')
      call timer%turn_on()
!
!     Term 1. - sum_ckdl b_aibl t_kl^cd L_kcjd
!             - sum_l X_jl b_aibl
!
!     Read X_jl to file
!
      call mem%alloc(X_jl, wf%n_o, wf%n_o)
!
      wf%jacobian_transpose_e2_oo_intermediate = sequential_file('jacobian_transpose_e2_oo_intermediate')
      call wf%jacobian_transpose_e2_oo_intermediate%open_('read', 'rewind')
!
      call wf%jacobian_transpose_e2_oo_intermediate%read_(X_jl, wf%n_o**2)
!
      call wf%jacobian_transpose_e2_oo_intermediate%close_()
!
!     Add - sum_ckdl b_aibl t_kl^cd L_kcjd
!         = - sum_l b_aib_l X_jl^T(l,j)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  b_aibj,               & ! "b_aib_l"
                  (wf%n_o)*(wf%n_v)**2, &
                  X_jl,                 &
                  wf%n_o,               &
                  one,                  &
                  sigma_aibj,           & ! "sigma_aib_j"
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_jl, wf%n_o, wf%n_o)
!
!     Term 2: - sum_ckdl b_aicj t_kl^cd L_ldkb  = - sum_c b_cjai X_cb
!
!     Read X_cb from file
!
      call mem%alloc(X_cb, wf%n_v, wf%n_v)
!
      wf%jacobian_transpose_e2_vv_intermediate = sequential_file('jacobian_transpose_e2_vv_intermediate')
      call wf%jacobian_transpose_e2_vv_intermediate%open_('read', 'rewind')
!
      call wf%jacobian_transpose_e2_vv_intermediate%read_(X_cb, wf%n_v**2)
!
      call wf%jacobian_transpose_e2_vv_intermediate%close_()
!
      call dgemm('T','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  X_cb,                 &
                  wf%n_v,               &
                  b_aibj,               & ! b_c_jai
                  wf%n_v,               &
                  one,                  &
                  sigma_aibj,           & ! sigma_b_jai, but we will symmetrize later on
                  (wf%n_v))
!
      call mem%dealloc(X_cb, wf%n_v, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_e2
!
!
   module subroutine save_jacobian_transpose_f2_intermediates(wf, t_ckdl, L_dlbi)
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
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_dlbi ! Reordered L_ldib
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbi ! An intermediate, term 2
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD F2 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Form the intermediate X_ckbi = sum_dl t_ck_dl L_dl_bi
!
      call mem%alloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  t_ckdl,            &
                  (wf%n_v)*(wf%n_o), &
                  L_dlbi,            &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ckbi,            &
                  (wf%n_v)*(wf%n_o))
!
!     Write X_ckbi to file
!
      wf%jacobian_transpose_f2_intermediate = sequential_file('jacobian_transpose_f2_intermediate')
      call wf%jacobian_transpose_f2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_f2_intermediate%write_(X_ckbi, wf%n_v**2 * wf%n_o**2)
!
      call wf%jacobian_transpose_f2_intermediate%close_()
!
      call mem%dealloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_f2_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_f2(wf, sigma_aibj, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbid
      real(dp), dimension(:,:,:,:), allocatable :: L_dibj ! L_jbid
      real(dp), dimension(:,:,:,:), allocatable :: L_dlbi ! L_ldib
!
      real(dp), dimension(:,:), allocatable :: X_ad   ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbi ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_lj   ! An intermediate, term 3
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi ! sigma_aibj contribution
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD F2', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_ckdl b_alck t_kl^cd L_jbid ::
!
!     X_ad = b_a_lck t_d_lck
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_ckdl, wf%n_t1)
!
!     Form the intermediate X_ad = sum_lck b_a_lck t_lck_d
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
!
      call dgemm('N','T',           &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  b_aibj,           & ! "b_a_lck"
                  wf%n_v,           &
                  t_ckdl,           & ! t_d_lck
                  wf%n_v,           &
                  zero,             &
                  X_ad,             &
                  wf%n_v)
!
!     Form g_jbid
!
      call mem%alloc(g_jbid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ovov', g_jbid)
!
!     Form L_dibj = L_jbid = 2 * g_jbid - g_jdib
!                          = 2 * g_jbid(j,b,i,d) - g_jbid(j,d,i,b)
!
      call mem%alloc(L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_dibj, (wf%n_o*wf%n_v)**2)
!
      call add_4321_to_1234(two, g_jbid, L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_jbid, L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jbid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Add - sum_ckdl b_alck t_kl^cd L_jbid
!         = - sum_d X_ad L_d_ibj
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  X_ad,                 &
                  wf%n_v,               &
                  L_dibj,               & ! L_d_ibj
                  wf%n_v,               &
                  one,                  &
                  sigma_aibj,           & ! "sigma_aibj"
                  wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
!
!     :: Term 2. - sum_ckdl b_ajck t_kl^cd L_ldib ::
!
!     X_ckbi = sum_dl t_ck_dl L_dl_bi
!
!     We have L_dibj = L_jbid => L_dibj(b,i,d,l) = L_ldib
!
!     Form L_dlbi = L_ldib = L_dibj(b,i,d,l)
!
      call mem%alloc(L_dlbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_dlbi, (wf%n_o*wf%n_v)**2)
!
      call sort_1234_to_3412(L_dibj, L_dlbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ckbi = sum_dl t_ck_dl L_dl_bi
!
      call mem%alloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_f2_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_f2_intermediate%read_(X_ckbi, wf%n_o**2 * wf%n_v**2)
      call wf%jacobian_transpose_f2_intermediate%close_()
!
!     Form and add - sum_ckdl b_ajck t_kl^cd L_ldib = - sum_ck b_aj_ck X_ck_bi
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  -one,              &
                  b_aibj,            & ! b_aj_ck
                  (wf%n_v)*(wf%n_o), &
                  X_ckbi,            & ! X_ck_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  sigma_ajbi,        &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3. - sum_ckdl b_djck t_kl^cd L_ialb = - sum_ckdl b_ckdj t_kl^cd L_ialb ::
!
!     X_lj = sum_ckd t_l_ckd b_ckd_j
!
!     Form the intermediate X_lj = sum_ckd t_kl^cd b_ckdj
!                                 = sum_ckd t_ckd_l b_ckd_j
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemm('T','N',           &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_o*wf%n_v**2, &
                  one,              &
                  t_ckdl,           & ! t_ckd_l
                  wf%n_o*wf%n_v**2, &
                  b_aibj,           & ! b_ckd_j
                  wf%n_o*wf%n_v**2, &
                  zero,             &
                  X_lj,             &
                  wf%n_o)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     - sum_ckdl b_ckdj t_kl^cd L_ialb = - sum_l L_ialb X_lj
!
!     We have L_dlbi = (L_ldib) = L_aib_l
!     Add - sum_ckdl b_ckdj t_kl^cd L_ialb = - sum_l L_aib_l X_lj
!
      call dgemm('N','N',                 &
                  (wf%n_o)*(wf%n_v)**2,   &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  L_dlbi,                 & ! L_dl_bi
                  (wf%n_o)*(wf%n_v)**2,   &
                  X_lj,                   &
                  wf%n_o,                 &
                  one,                    &
                  sigma_aibj,             & ! sigma_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(L_dlbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_f2
!
!
   module subroutine save_jacobian_transpose_g2_intermediates(wf, t_aibj, g_kdib)
!!
!!    Save Jacobian transpose g2 intermediates
!!    Written by Eirik F. Kjønstad and S. D. Folkestad, Nov 2019
!!
!!    Constructs intermediate
!!
!!       X_clbi = sum_dk t_ckdl g_kbid
!!
!!    and saves them to files 'jacobian_transpose_g2_intermediate'
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_kdib
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cldk ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_dkbi ! g_kbid
!
      real(dp), dimension(:,:,:,:), allocatable :: X_clbi ! An intermediate, term 1
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD G2 intermediates', pl='verbose')
      call timer%turn_on()
!
!     Form t_cldk = t_kl^cd
!
      call mem%alloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(t_aibj, t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder g_kbid to g_dkbi
!
      call mem%alloc(g_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_4123(g_kdib, g_dkbi, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     X_clbi = sum_dk t_cl_dk g_dk_bi
!
      call mem%alloc(X_clbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  t_cldk,            & ! t_cl_dk
                  (wf%n_v)*(wf%n_o), &
                  g_dkbi,            & ! g_dk_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_clbi,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      wf%jacobian_transpose_g2_intermediate = sequential_file('jacobian_transpose_g2_intermediate')
      call wf%jacobian_transpose_g2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_g2_intermediate%write_(X_clbi, wf%n_v**2 * wf%n_o**2)
!
      call wf%jacobian_transpose_g2_intermediate%close_()
!
      call mem%dealloc(X_clbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_g2_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_g2(wf, sigma_aibj, b_aibj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: b_ajcl     ! b_alcj
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi ! sigma_aibj contribution
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajib ! sigma_aibj contribution
!
      real(dp), dimension(:,:,:,:), allocatable :: X_clbi ! An intermediate, term 1
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj, t_clkd
!
      real(dp), dimension(:,:,:), allocatable :: W_vo_J, Z_vo_J, L_J_ov
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD G2', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. sum_ckdl b_alcj t_kl^cd g_kbid ::
!
!     Read pre-generated intermediate X_clbi = sum_dk t_ckdl g_kbid
!
      call mem%alloc(X_clbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_g2_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_g2_intermediate%read_(X_clbi, wf%n_v**2 * wf%n_o**2)
      call wf%jacobian_transpose_g2_intermediate%close_()
!
!     Reorder to b_ajcl = b_alcj
!
      call mem%alloc(b_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(b_aibj, b_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate and add sum_ckdl b_alcj t_kl^cd g_kbid
!                       = sum_cl b_aj_cl X_cl_bi
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_ajcl,            & ! b_aj_cl
                  (wf%n_v)*(wf%n_o), &
                  X_clbi,            & ! X_cl_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  sigma_ajbi,        &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_clbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(b_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. sum_ckdl b_ajcl t_kl^cd g_kdib ::
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_t1)
!
      call mem%alloc(t_clkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_1234_to_1423(t_aibj, t_clkd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate object W_cl,J = sum_kd t_clkd L_J,kd
!
      call mem%alloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
      call wf%eri%get_cholesky_t1(L_J_ov, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(W_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
      call dgemm('N', 'T',          &
                  wf%n_o * wf%n_v,  &
                  wf%eri%n_J,       &
                  wf%n_o * wf%n_v,  &
                  one,              &
                  t_clkd,           & ! t_cl,kd
                  wf%n_o * wf%n_v,  &
                  L_J_ov,           & ! L_J,kd
                  wf%eri%n_J,       &
                  zero,             &
                  W_vo_J,           & ! W_cl,J
                  wf%n_o * wf%n_v)
!
      call mem%dealloc(t_clkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     sigma_aibj += b_ajcl W_cl,J L_J,ib = Z_aj,J L_J,ib
!
      call mem%alloc(Z_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
      call dgemm('N', 'N',          &
                  wf%n_o * wf%n_v,  &
                  wf%eri%n_J,       &
                  wf%n_o * wf%n_v,  &
                  one,              &
                  b_aibj,           & ! b_aj,cl
                  wf%n_o * wf%n_v,  &
                  W_vo_J,           & ! W_cl,J
                  wf%n_o * wf%n_v,  &
                  zero,             &
                  Z_vo_J,           &
                  wf%n_o * wf%n_v)
!
      call mem%dealloc(W_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
!     sigma_ajib = Z_aj,J L_J,ib
!
      call mem%alloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',          &
                  wf%n_o * wf%n_v,  &
                  wf%n_o * wf%n_v,  &
                  wf%eri%n_J,       &
                  one,              &
                  Z_vo_J,           &
                  wf%n_o * wf%n_v,  &
                  L_J_ov,           &
                  wf%eri%n_J,       &
                  zero,             &
                  sigma_ajib,       &
                  wf%n_o * wf%n_v)
!
      call mem%dealloc(L_J_ov, wf%eri%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(Z_vo_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
      call add_1423_to_1234(one, sigma_ajib, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_g2
!
!
   module subroutine save_jacobian_transpose_i2_intermediates(wf, t_aibj, g_ovov)
!!
!!    Save Jacobian transpose i2 intermediates
!!    Written by Tor S. Haugland, Eirik F. Kjønstad and
!!    Sarai D. Folkestad, Nov 2019
!!
!!    Construct intermediate
!!
!!       X_klij = sum_cd t_kl^cd g_icjd + g_ikjl
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
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_ovov
!
      real(dp), dimension(:,:,:,:), allocatable :: t_klcd ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_cdij, g_ikjl, g_klij
!
      real(dp), dimension(:,:,:,:), allocatable :: X_klij ! An intermediate, terms 1 & 2
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD I2 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Reorder g_icjd to g_cdij
!
      call mem%alloc(g_cdij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_2413(g_ovov, g_cdij, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form t_ckdl and sort to t_klcd
!
      call mem%alloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2413(t_aibj, t_klcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_kl_ij = sum_cd t_kl_cd g_cd_ij
!
      call mem%alloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  t_klcd,      & ! t_kl_cd
                  (wf%n_o)**2, &
                  g_cdij,      & ! g_cd_ij
                  (wf%n_v)**2, &
                  zero,        &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2)
!
      call mem%dealloc(g_cdij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Add g_ikjl
!
      call mem%alloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%eri%get_eri_t1('oooo', g_ikjl)
      call sort_1234_to_2413(g_ikjl, g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call daxpy(wf%n_o**4, one, g_klij, 1, X_klij, 1)
!
      call mem%dealloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Write X_klij to file
!
      wf%jacobian_transpose_i2_intermediate = sequential_file('jacobian_transpose_i2_intermediate')
      call wf%jacobian_transpose_i2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_i2_intermediate%write_(X_klij, wf%n_o**4)
!
      call wf%jacobian_transpose_i2_intermediate%close_()
!
      call mem%dealloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_i2_intermediates
!
!
   module subroutine jacobian_transpose_ccsd_i2(wf, sigma_abij, b_abij)
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
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: t_klcd ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kalb ! g_kalb
      real(dp), dimension(:,:,:,:), allocatable :: g_abkl ! g_kalb
!
      real(dp), dimension(:,:,:,:), allocatable :: X_klij ! An intermediate, terms 1 & 2
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD I2', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. sum_ckdl b_cidj t_kl^cd g_kalb ::
!
!     sum_ckdl t_kl_cd b_cd_ij
!
!     Form t_kl_cd = t_kl^cd
!
      call mem%alloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call squareup_and_sort_1234_to_2413(wf%t2, t_klcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_klij = sum_cd t_kl_cd b_cd_ij
!
      call mem%alloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  t_klcd,      & ! t_kl_cd
                  (wf%n_o)**2, &
                  b_abij,      & ! "b_cd_ij"
                  (wf%n_v)**2, &
                  zero,        &
                  X_klij,      &
                  (wf%n_o)**2)
!
      call mem%dealloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!
!     Form g_kalb
!
      call mem%alloc(g_kalb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ovov', g_kalb)
!
!     Reorder to g_abkl = g_kalb
!
      call mem%alloc(g_abkl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_2413(g_kalb, g_abkl, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kalb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Add sum_ckdl b_cidj t_kl^cd g_kalb
!         = sum_kl g_ab_kl X_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  g_abkl,      &
                  (wf%n_v)**2, &
                  X_klij,      &
                  (wf%n_o)**2, &
                  one,         &
                  sigma_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(g_abkl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2. sum_ckdl b_akbl t_kl^cd g_icjd ::
!
!     Repurpose X_kl_ij to make sum_cd t_kl^cd g_icjd
!                               = sum_cd t_kl_cd g_cd_ij
!                               = sum_cd t_kl_cd g_ab_kl(cd,ij)
!
      call wf%jacobian_transpose_i2_intermediate%open_('read', 'rewind')
      call wf%jacobian_transpose_i2_intermediate%read_(X_klij, wf%n_o**4)
      call wf%jacobian_transpose_i2_intermediate%close_()
!
!     Add sum_ckdl b_akbl t_kl^cd g_icjd
!         = sum_kl b_ab_kl X_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  b_abij,      & ! "b_ab_kl"
                  (wf%n_v)**2, &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  sigma_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccsd_i2
!
!
   module subroutine save_jacobian_transpose_e2_oo_intermediate(wf, t_ckdl, L_ckdj)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_ckdj ! L_kcjd
!
      real(dp), dimension(:,:), allocatable :: X_jl
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD E2 oo intermediate', pl='verbose')
      call timer%turn_on()
!
!     X_jl = sum_kcd L_kcjd t_kl^cd = sum_kcd L_j_ckd t_ckd_l
!
      call mem%alloc(X_jl, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  L_ckdj,               & ! L_ckd_j
                  (wf%n_o)*(wf%n_v)**2, &
                  t_ckdl,               & ! t_ckd_l
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_jl,                 &
                  wf%n_o)
!
!     Write X_jl to file
!
      wf%jacobian_transpose_e2_oo_intermediate = sequential_file('jacobian_transpose_e2_oo_intermediate')
      call wf%jacobian_transpose_e2_oo_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_e2_oo_intermediate%write_(X_jl, wf%n_o**2)
!
      call wf%jacobian_transpose_e2_oo_intermediate%close_()
!
      call mem%dealloc(X_jl, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_e2_oo_intermediate
!
!
   module subroutine save_jacobian_transpose_e2_vv_intermediate(wf, t_ckdl, L_bkdl)
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
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_bkdl ! L_kbld
!
      real(dp), dimension(:,:), allocatable :: X_cb
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCSD E2 vv intermediate', pl='verbose')
      call timer%turn_on()
!
!     X_cd  sum_kdl t_kl^cd L_ldkb
!
      call mem%alloc(X_cb, wf%n_v, wf%n_v)
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  t_ckdl,               & ! t_c_kdl
                  wf%n_v,               &
                  L_bkdl,               & ! L_b_kdl
                  (wf%n_v),             &
                  zero,                 &
                  X_cb,                 &
                  wf%n_v)
!
!     Write X_cb to file
!
      wf%jacobian_transpose_e2_vv_intermediate = sequential_file('jacobian_transpose_e2_vv_intermediate')
      call wf%jacobian_transpose_e2_vv_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_transpose_e2_vv_intermediate%write_(X_cb, wf%n_v**2)
!
      call wf%jacobian_transpose_e2_vv_intermediate%close_()
!
      call mem%dealloc(X_cb, wf%n_v, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_e2_vv_intermediate
!
!
end submodule jacobian_transpose_ccsd
