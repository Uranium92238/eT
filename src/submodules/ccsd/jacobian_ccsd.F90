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
submodule (ccsd_class) jacobian_ccsd
!
!!
!!    Jacobian submodule (CCSD)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2017-2018
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    Transfered to the current eT program from the first version
!!    of eT by Andreas Skeidsvoll and Sarai D. Folkestad, 2018.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine jacobian_transform_trial_vector_ccsd(wf, c_i)
!!
!!    Jacobian transform trial vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2018
!!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
      call wf%jacobian_ccsd_transformation(c_i)
!
   end subroutine jacobian_transform_trial_vector_ccsd
!
!
   module subroutine jacobian_ccsd_transformation_ccsd(wf, c)
!!
!!    Jacobian transformation (CCSD)
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
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_abij
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      rho_ai = zero
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai)
!
         enddo
      enddo
!$omp end parallel do
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_ccsd_a1(rho_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c_aibj(a,i,b,j) = c(wf%n_o*wf%n_v + aibj)
                     c_aibj(b,j,a,i) = c(wf%n_o*wf%n_v + aibj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!$omp parallel do schedule(static) private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
         c_aibj(a,i,a,i) = two*c_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
!
      call wf%jacobian_ccsd_b1(rho_ai, c_aibj)
      call wf%jacobian_ccsd_c1(rho_ai, c_aibj)
      call wf%jacobian_ccsd_d1(rho_ai, c_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai) = rho_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do

      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     :: CCSD contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      rho_aibj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_ccsd_a2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_b2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_c2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_d2(rho_aibj, c_ai)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_f2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_g2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_h2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_i2(rho_aibj, c_aibj)
!
!     Last two terms are already symmetric (J2 and K2). Perform the symmetrization
!     rho_aibj = P_ij^ab rho_aibj now, for convenience
!
      call symmetric_sum(rho_aibj, (wf%n_v)*(wf%n_o))
!
!     In preparation for last two terms, reorder
!     rho_aibj to rho_abij, and c_aibj to c_abij
!
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(rho_aibj, rho_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_ccsd_j2(rho_abij, c_abij)
      call wf%jacobian_ccsd_k2(rho_abij, c_abij)
!
!     Done with reordered doubles c; deallocate
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Order rho_abij back into rho_aibj & divide by
!     the biorthonormal factor 1 + delta_ai,bj
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1324(rho_abij, rho_aibj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!$omp parallel do schedule(static) private(a,i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
         rho_aibj(a,i,a,i) = half*rho_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Overwrite the incoming doubles c vector & pack in
!
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj) = rho_aibj(a,i,b,j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Remaining deallocations
!
      call mem%dealloc(rho_aibj,wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_ccsd_transformation_ccsd
!
!
  module subroutine jacobian_ccsd_a1_ccsd(wf, rho_ai, c_ai)
!!
!!    Jacobian CCSD A1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^A1 = sum_ckdl L_lckd (u_li^ca c_dk  - t_li^cd c_ak - t_lk^ad c_ci)
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o)                :: rho_ai
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_lckd
      real(dp), dimension(:,:,:,:), allocatable :: L_lcdk
      real(dp), dimension(:,:,:,:), allocatable :: L_lkdc
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: X_lc
      real(dp), dimension(:,:), allocatable :: X_ik
      real(dp), dimension(:,:), allocatable :: X_ac
!
!     Amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ilcd
      real(dp), dimension(:,:,:,:), allocatable :: t_alkd
      real(dp), dimension(:,:,:,:), allocatable :: t_clai
      real(dp), dimension(:,:,:,:), allocatable :: u_ailc
!
      type(timings) :: jacobian_ccsd_a1_timer
!
      call jacobian_ccsd_a1_timer%init('jacobian ccsd a1')
      call jacobian_ccsd_a1_timer%start()
!
!     :: Term 1: sum_ckdl L_lckd u_li^ca c_dk ::
!
!     g_lc_kd = g_lckd
!
      call mem%alloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_lckd)
!
      call mem%alloc(L_lcdk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      L_lcdk = zero
!
!     L_lc_dk(lc,dk) = L_lckd = 2*g_lckd - g_ldkc = 2*g_lc_kd(lc,kd) - g_lc_kd(ld,kc)
!
      call add_1243_to_1234(two, g_lckd, L_lcdk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call add_1342_to_1234(-one, g_lckd, L_lcdk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     X_lc = sum_kd L_lckd c_dk = sum_kd L_lc_dk c_dk
!
      call mem%alloc(X_lc, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lcdk,            & ! L_lc,dk
                  (wf%n_o)*(wf%n_v), &
                  c_ai,              & ! c_dk
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lc,              &
                  (wf%n_o)*(wf%n_v))
!
!     Form u_ai_lc = u_li^ca = 2 * t_li^ca - t_il^ca = 2 * t2(clai,1) - t2(cial,1)
!
      call mem%alloc(t_clai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_clai, (wf%n_o)*(wf%n_v))
!
!     t_cl_ai(cl,ai) = t_li^ca
!
!     u_ai_lc(ai, lc) = 2 * t_li^ca - t_il^ca = 2 * t_cl_ai(cl, ai) - t_cl_ai(ci, al)
!
      call mem%alloc(u_ailc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      u_ailc = zero
!
      call add_4312_to_1234(two, t_clai, u_ailc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_4213_to_1234(-one, t_clai, u_ailc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_ai =+ sum_lc u_ai_lc X_lc
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  u_ailc,            & ! u_ai,lc
                  (wf%n_v)*(wf%n_o), &
                  X_lc,              & ! X_lc
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  rho_ai,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(u_ailc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(X_lc, wf%n_o, wf%n_v)
!
!     :: Term 2. - sum_ckdl L_lckd t_li^cd c_ak ::
!
!     t_il_cd(il, cd) = t_li^cd = t_cl_ai(cl, di)
!
      call mem%alloc(t_ilcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      t_ilcd = zero
!
      call sort_1234_to_4213(t_clai, t_ilcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_clai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate X_i_k = sum_cdl L_lcd_k t_i_lcd
!
      call mem%alloc(X_ik, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_ilcd,               & ! t_i,lcd
                  (wf%n_o),             &
                  L_lcdk,               & ! L_lcd,k
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_ik,                 &
                  wf%n_o)
!
!     Calculate rho_ai =+ - sum_k c_ai(a,k) X_i_k(i, k)
!
      call dgemm('N', 'T', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  c_ai,    & ! c_a,k
                  wf%n_v,  &
                  X_ik,    & ! X_i,k
                  wf%n_o,  &
                  one,     &
                  rho_ai,  &
                  wf%n_v)
!
      call mem%dealloc(X_ik, wf%n_o, wf%n_o)
!
!     :: Term 3: - sum_ckdl L_lckd t_lk^ad c_ci ::
!
!     Reorder to L_lk_dc(lk,dc) = L_lckd = L_lc_dk(lc,dk)
!
      call mem%alloc(L_lkdc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_1432(L_lcdk, L_lkdc, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_lcdk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Reorder amplitudes to t_al_kd(al,kd) = t_lk^ad = t_il_cd(kl, ad)
!
      call mem%alloc(t_alkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_3214(t_ilcd, t_alkd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_ilcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Calculate X_a_c = sum_kdl t_a,lkd L_lkd,c
!
      call mem%alloc(X_ac, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  t_alkd,               & ! t_a,lkd
                  wf%n_v,               &
                  L_lkdc,               & ! L_lkd,c
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_ac,                 & ! X_a,c
                  wf%n_v)
!
      call mem%dealloc(L_lkdc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(t_alkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Calculate rho_ai =+ - sum_c X_a_c(a,c) c_ai(c,i)
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  X_ac,    & ! X_a,c
                  wf%n_v,  &
                  c_ai,    & ! c_c,i
                  wf%n_v,  &
                  one,     &
                  rho_ai,  & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_ac, wf%n_v, wf%n_v)
!
      call jacobian_ccsd_a1_timer%freeze()
      call jacobian_ccsd_a1_timer%switch_off()
!
   end subroutine jacobian_ccsd_a1_ccsd
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, rho_ai, c_aibj)
!!
!!    Jacobian CCSD B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    rho_ai^B1 = sum_bj F_jb (2*c_aibj - c_ajbi)
!!              = sum_bj F_jb v_aijb
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: v_aijb
!
      type(timings) :: jacobian_ccsd_b1_timer
!
      call jacobian_ccsd_b1_timer%init('jacobian ccsd b1')
      call jacobian_ccsd_b1_timer%start()
!
!     Construct v_aibj = 2*c_aibj - c_ajbi ordered as
!
!        v_ai_jb(a,i,j,b) = 2*c_aibj(a,i,b,j) - c_aibj(a,j,b,i)
!
!     and do the matrix multiplication with F_jb
!
      call mem%alloc(v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      v_aijb = zero
!
      call add_1243_to_1234(two, c_aibj, v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one, c_aibj, v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  v_aijb,            & ! v_ai,jb
                  (wf%n_o)*(wf%n_v), &
                  wf%fock_ia,        & ! F_jb
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai,            & ! rho_ai
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call jacobian_ccsd_b1_timer%freeze()
      call jacobian_ccsd_b1_timer%switch_off()
!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, rho_ai, c_aibj)
!!
!!    Jacobian CCSD C1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^C1 = - sum_bjk L_jikb c_ajbk
!!              = - sum_bjk (2*g_jikb - g_kijb) c_ajbk
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jikb
      real(dp), dimension(:,:,:,:), allocatable :: L_jbki
!
      type(timings) :: jacobian_ccsd_c1_timer
!
      call jacobian_ccsd_c1_timer%init('jacobian ccsd c1')
      call jacobian_ccsd_c1_timer%start()
!
!     Construct L_jikb = 2*g_jikb - g_kijb as
!
!        L_jb_ki(jb,ki) = 2*g_ji_kb(ji,kb) - g_ji_kb(ki,jb)
!
!     and then contract with c_ajbk = c_aibj(aj,bk).
!
      call mem%alloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_ooov(g_jikb)
!
      call mem%alloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      L_jbki = zero
!
      call add_1432_to_1234(two, g_jikb, L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call add_3412_to_1234(-one, g_jikb, L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  c_aibj,                 & ! c_a,jbk
                  wf%n_v,                 &
                  L_jbki,                 & ! L_jbk,i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  rho_ai,                 & ! rho_ai
                  wf%n_v)
!
      call mem%dealloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call jacobian_ccsd_c1_timer%freeze()
      call jacobian_ccsd_c1_timer%switch_off()
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
   module subroutine jacobian_ccsd_d1_ccsd(wf, rho_ai, c_bicj)
!!
!!    Jacobian CCSD D1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bicj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
      integer :: current_a_batch
!
      type(batching_index) :: batch_a
!
      real(dp), dimension(:,:,:,:), allocatable :: c_bjci
      real(dp), dimension(:,:,:,:), allocatable :: g_abjc
      real(dp), dimension(:,:,:,:), allocatable :: L_abjc
!
      type(timings) :: jacobian_ccsd_d1_timer
!
      integer :: rec0, rec1
!
      call jacobian_ccsd_d1_timer%init('jacobian ccsd d1')
      call jacobian_ccsd_d1_timer%start()
!
!     Prepare for batching over index a
!
      rec0 = wf%n_o*wf%integrals%n_J*wf%n_v
!
      rec1 = wf%n_v*wf%integrals%n_J + (wf%n_v**2)*(wf%n_o)
!
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_a, rec0, rec1)
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Determine the limits for the current a-batch
!
         call batch_a%determine_limits(current_a_batch)
!
!        Construct L_abjc = 2 g_abjc - g_acjb
!
         call mem%alloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_abjc,                        &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
         call mem%alloc(L_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
         L_abjc = zero
!
         call daxpy((batch_a%length)*(wf%n_v)**2*(wf%n_o), two, g_abjc, 1, L_abjc, 1)
         call add_1432_to_1234(-one, g_abjc, L_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
!        Reorder c_bicj to c_bjci
!
         call mem%alloc(c_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
         call sort_1234_to_1432(c_bicj, c_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call mem%dealloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('N', 'N',                   &
                     batch_a%length,            &
                     wf%n_o,                    &
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     L_abjc,                    & ! L_a,bjc
                     batch_a%length,            &
                     c_bjci,                    & ! c_bjc,i
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     rho_ai(batch_a%first, 1),  &
                     wf%n_v)
!
         call mem%dealloc(L_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
         call mem%dealloc(c_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      enddo ! End batching over a
!
      call jacobian_ccsd_d1_timer%freeze()
      call jacobian_ccsd_d1_timer%switch_off()
!
   end subroutine jacobian_ccsd_d1_ccsd
!
!
   module subroutine jacobian_ccsd_a2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^A2 = sum_c g_aibc c_cj - sum_k g_aikj c_bk
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)          :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)      :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aikj
      real(dp), dimension(:,:,:,:), allocatable :: g_aijk
      real(dp), dimension(:,:,:,:), allocatable :: g_aibc
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij
!
      integer :: current_b_batch
!
      type(batching_index) :: batch_b
!
      type(timings) :: jacobian_ccsd_a2_timer
!
      integer :: rec0, rec1
!
      call jacobian_ccsd_a2_timer%init('jacobian ccsd a2')
      call jacobian_ccsd_a2_timer%start()
!
!     :: Term 1. - sum_k g_aikj c_bk ::
!
!     Calculate g_aikj
!
      call mem%alloc(g_aikj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%get_vooo(g_aikj)
!
!     Reorder g_aikj to g_aijk
!
      call mem%alloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1243(g_aikj, g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_aikj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%alloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b,k
                  wf%n_v,               &
                  g_aijk,               & ! g_aij,k
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  rho_baij,             & ! rho_b,aij
                  wf%n_v)
!
      call mem%dealloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call add_3124_to_1234(one, rho_baij, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2. rho_aibj =+ sum_c g_aibc c_cj ::
!
!     We do the matrix multiplication as g_aib_c c_cj, batching over b.
!
      rec0 = wf%n_o*wf%integrals%n_J*wf%n_v
!
      rec1 = wf%n_v*wf%integrals%n_J + (wf%n_v**2)*(wf%n_o)
!
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_b, rec0, rec1)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
!        Calculate the batch contribution to rho_aib_j = sum_c g_aibc c_cj
!
         call mem%alloc(g_aibc, wf%n_v, wf%n_o, wf%n_v, batch_b%length)
!
         call wf%get_vovv(g_aibc,        &
                          1,             &
                          wf%n_v,        &
                          1,             &
                          wf%n_o,        &
                          batch_b%first, &
                          batch_b%last,  &
                          1,             &
                          wf%n_v)
!
         call dgemm('N', 'N',                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_aibc,                             & ! g_aib,c
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_c,j
                     wf%n_v,                             &
                     one,                                &
                     rho_aibj(1,1, batch_b%first,1),     &
                     (wf%n_o)*(wf%n_v)**2)
!
         call mem%dealloc(g_aibc, wf%n_v, wf%n_o, wf%n_v, batch_b%length)
!
      enddo ! End of batches over b
!
      call jacobian_ccsd_a2_timer%freeze()
      call jacobian_ccsd_a2_timer%switch_off()
!
   end subroutine jacobian_ccsd_a2_ccsd
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
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)       :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)   :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aicj   ! t_ij^ac
      real(dp), dimension(:,:,:,:), allocatable :: t_caij   ! t_ij^ac
      real(dp), dimension(:,:,:,:), allocatable :: X_kaij   ! An intermediate
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij ! rho_aibj, reordered, term 1
!
      real(dp), dimension(:,:), allocatable :: X_kj         ! An intermediate
!
!
      type(timings) :: jacobian_ccsd_b2_timer
!
      call jacobian_ccsd_b2_timer%init('jacobian ccsd b2')
      call jacobian_ccsd_b2_timer%start()
!
!     :: Term 1. - sum_kc F_kc t_ij^ac c_bk ::
!
!     Order the amplitudes as t_ca_ij = t_ij^ac
!
      call mem%alloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_aicj, (wf%n_o)*(wf%n_v))
!
!     t_ai_cj to t_ca_ij
!
      call mem%alloc(t_caij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_3124(t_aicj, t_caij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_k_aij = sum_c F_k_c t_c_aij
!
      call mem%alloc(X_kaij, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  wf%fock_ia,           & ! F_k,c
                  wf%n_o,               &
                  t_caij,               & ! t_c,aij
                  wf%n_v,               &
                  zero,                 &
                  X_kaij,               &
                  wf%n_o)
!
      call mem%dealloc(t_caij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Form rho_b_aij = sum_k c_ai(b,k) X_k_aij(k,aij)
!
      call mem%alloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b,k
                  wf%n_v,               &
                  X_kaij,               & ! X_k_aij
                  wf%n_o,               &
                  zero,                 &
                  rho_baij,             & ! rho_b,aij
                  wf%n_v)
!
      call mem%dealloc(X_kaij, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Add rho_ba_ij(ba,ij) to rho_aibj(ai,bj)
!
      call add_3124_to_1234(one, rho_baij, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2. - sum_kc F_kc t_ik^ab c_cj ::
!
!     Form X_kj = sum_c F_kc c_cj = sum_c fock_ia(k,c) c_ai(c,j)
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
!
      call dgemm('N','N',     &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ia, & ! F_k,c
                  wf%n_o,     &
                  c_ai,       & ! c_c,j
                  wf%n_v,     &
                  zero,       &
                  X_kj,       &
                  wf%n_o)
!
!     Form rho_aib_j = - sum_k t_aib_k X_k_j
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_aicj,               & ! t_aib,k
                  (wf%n_o)*(wf%n_v)**2, &
                  X_kj,                 &
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_aib,j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call jacobian_ccsd_b2_timer%freeze()
      call jacobian_ccsd_b2_timer%switch_off()
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
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljkc ! g_ljkc
      real(dp), dimension(:,:,:,:), allocatable :: g_kjlc ! g_ljkc
      real(dp), dimension(:,:,:,:), allocatable :: L_ljck ! L_ljkc
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ljai   ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_kjbi   ! An intermediate, term 2
      real(dp), dimension(:,:,:,:), allocatable :: X_kjli   ! An intermediate, term 3
      real(dp), dimension(:,:,:,:), allocatable :: X_klij   ! X_kj_li reordered
      real(dp), dimension(:,:), allocatable     :: X_lj     ! An intermediate, term 4
      real(dp), dimension(:,:,:,:), allocatable :: Y_ljai  ! An intermediate, term 5
!
      real(dp), dimension(:,:,:,:), allocatable :: t_akci ! t_ki^ac
      real(dp), dimension(:,:,:,:), allocatable :: t_kcai ! t_ki^ac
      real(dp), dimension(:,:,:,:), allocatable :: t_bakl ! t_lk^ba
      real(dp), dimension(:,:,:,:), allocatable :: t_aibl ! t_il^ab
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_bjai ! rho_aibj, term 1 & 5
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi ! rho_aibj, term 2
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij ! rho_aibj, term 3
!
      type(timings) :: jacobian_ccsd_c2_timer
!
      call jacobian_ccsd_c2_timer%init('jacobian ccsd c2')
      call jacobian_ccsd_c2_timer%start()
!
!     :: Term 1. sum_kcl g_ljkc t_ki^ac c_bl ::
!
!     Form g_ljkc
!
      call mem%alloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_ljkc)
!
!     Read the amplitudes from disk
!
      !call wf%read_double_amplitudes
!
!     Square up (t_ak_ci = t_ki^ac)
!
      call mem%alloc(t_akci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_akci, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Order as t_kc_ai = t_ki^ac
!
      call mem%alloc(t_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Reorder t_ak_ci to t_kc_ai
!               1234       2314
!
      call sort_1234_to_2314(t_akci, t_kcai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_akci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form X_lj_ai = sum_ck g_lj_kc t_kc_ai
!
      call mem%alloc(X_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_ljkc,            & ! g_lj_kc
                  (wf%n_o)**2,       &
                  t_kcai,            & ! t_kc_ai
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ljai,            & ! X_lj_ai
                  (wf%n_o)**2)
!
!     Calculate rho_b_jai = sum_l c_bl X_ljai
!     (Interpret the X array as an X_l_jai object in the matrix multiplication)
!
      call mem%alloc(rho_bjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_ai,                 & ! c_b_l
                  wf%n_v,               &
                  X_ljai,               & ! X_l_jai
                  wf%n_o,               &
                  zero,                 &
                  rho_bjai,             & ! rho_b_jai
                  wf%n_v)
!
      call mem%dealloc(X_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Add rho_bj_ai (3412) to rho_aibj (1234)
!
      call add_3412_to_1234(one, rho_bjai, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Final deallocations for term 1 (keep g_ljkc for later use)
!
      call mem%dealloc(rho_bjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. sum_kcl g_ljkc t_li^bc c_ak ::
!
!     Reorder to g_kj_lc = g_lj_kc = g_ljkc
!                  3214      1234
!
      call mem%alloc(g_kjlc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_3214(g_ljkc, g_kjlc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Form the intermediate X_kjbi = sum_lc g_kjlc t_lcbi
!
      call mem%alloc(X_kjbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_kjlc,            & ! g_kj_lc
                  (wf%n_o)**2,       &
                  t_kcai,            & ! t_lc_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_kjbi,            & ! X_kj_bi
                  (wf%n_o)**2)
!
!     Calculate rho_a_jbi = sum_k c_ak X_kjbi
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_ai,                 & ! c_a_k
                  wf%n_v,               &
                  X_kjbi,               & ! X_k_jbi
                  wf%n_o,               &
                  zero,                 &
                  rho_ajbi,             & ! rho_a_jbi
                  wf%n_v)
!
      call mem%dealloc(X_kjbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Add rho_ajbi to rho_aibj
!             1432        1234
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Deallocations for term 2 (keep g_kj_lc = g_ljkc)
!
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3. sum_kcl g_ljkc t_lk^ba c_ci ::
!
!     Form the intermediate X_kjl_i = sum_c g_ljkc c_ci = sum_c g_kjlc c_ci
!
      call mem%alloc(X_kjli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**3, &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  g_kjlc,      & ! g_kjl_c
                  (wf%n_o)**3, &
                  c_ai,        & ! c_c_i
                  wf%n_v,      &
                  zero,        &
                  X_kjli,      & ! X_kjl_i
                  (wf%n_o)**3)
!
!     Reorder to X_kl_ij = X_kj_li
!                  1342      1234
!
      call mem%alloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1342(X_kjli, X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_kjli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Order amplitudes as t_ba_kl = t_lk^ba
!
      call mem%alloc(t_bakl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     t_kc_ai = t_ki^ac => t_kc_ai(la,bk) = t_lk^ba = t_ba_kl(ba,kl)
!                                  1234                       3241
!
      call sort_1234_to_3241(t_kcai, t_bakl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Calculate rho_ba_ij = sum_kcl g_ljkc t_lk^ba c_ci
!                         = sum kl ( sum_c g_ljkc c_ci ) t_lk^ba
!                         = sum_kl t_ba_kl X_kl_ij
!
      call mem%alloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_bakl,      & ! t_ba_kl
                  (wf%n_v)**2, &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  zero,        &
                  rho_baij,    & ! rho_ba_ij
                  (wf%n_v)**2)
!
      call mem%dealloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Add rho_ba_ij into rho_aibj
!             3124           1234
!
      call add_3124_to_1234(one, rho_baij, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Deallocations for term 3 (keep g_kj_lc = g_ljkc)
!
      call mem%dealloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 4. - sum_kcl L_ljkc t_il^ab c_ck ::
!
!     Form L_lj_ck(lj,ck) = L_ljkc = 2 * g_ljkc - g_lckj
!                  1234   = 2 * g_ljkc - g_kjlc = 2* g_kj_lc(kj,lc) - g_kj_lc(lj,kc)
!                                                            4213             1243
!
      call mem%alloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      L_ljck = zero
!
      call add_4213_to_1234(two, g_kjlc, L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_1243_to_1234(-one, g_kjlc, L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kjlc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Calculate the intermediate X_lj = sum_ck L_ljck c_ck
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)**2,       &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ljck,            & ! L_lj_ck
                  (wf%n_o)**2,       &
                  c_ai,              & ! c_ck
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lj,              & ! X_lj
                  (wf%n_o)**2)
!
!     Order the amplitudes as t_ai_bl = t_il^ab
!
      call mem%alloc(t_aibl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     t_lk^ba = t_ba_kl(ba,kl) => t_il^ab = t_bakl(ab,li) = t_aibl(ai,bl)
!                                                   1234             1423
!
      call sort_1234_to_1423(t_bakl, t_aibl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(t_bakl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Form rho_aibj =+ - sum_l t_il^ab X_lj = - sum_l t_aib_l X_lj
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_aibl,               & ! t_aib_l
                  (wf%n_o)*(wf%n_v)**2, &
                  X_lj,                 & ! X_l_j
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
!     :: Term 5. - sum_kcl L_ljkc t_ik^ac c_bl ::
!
!     t_il^ab = t_ai_bl(ai,bl) => t_ai_bl(ck,ai) = t_ki^ca = t_ik^ac
!
!     Form the intermediate Y_lj_ai = sum_kc L_ljkc t_ik^ac = sum_kc L_lj_ck t_ck_ai
!
      call mem%alloc(Y_ljai, wf%n_o, wf%n_o, wf%n_v,wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ljck,            & ! L_lj_ck
                  (wf%n_o)**2,       &
                  t_aibl,            & ! t_ck_ai
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ljai,            & ! Y_lj_ai
                  (wf%n_o)**2)
!
      call mem%dealloc(t_aibl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate rho_b_jai =+ - sum_l c_bl Y_lj_ai
!
!     Note: interpret Y_lj_ai as Y_l_jai in the matrix multiplication
!
      call mem%alloc(rho_bjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b_l
                  wf%n_v,               &
                  Y_ljai,               & ! Y_l_jai
                  wf%n_o,               &
                  zero,                 &
                  rho_bjai,             & ! rho_b_jai
                  wf%n_v)
!
!     Add rho_bjai to rho_aibj
!
      call add_3412_to_1234(one, rho_bjai, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Deallocations and cleanup
!
      call mem%dealloc(Y_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_bjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call jacobian_ccsd_c2_timer%freeze()
      call jacobian_ccsd_c2_timer%switch_off()
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
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bdkc ! g_kcbd
      real(dp), dimension(:,:,:,:), allocatable :: g_cdkb ! g_kcbd reordered
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbd ! g_kcbd reordered
      real(dp), dimension(:,:,:,:), allocatable :: L_ckbd ! L_kcbd = 2 g_kcbd - g-kdbc
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ijcd ! t_ij^cd
      real(dp), dimension(:,:,:,:), allocatable :: t_dkaj ! t_kj^ad
      real(dp), dimension(:,:,:,:), allocatable :: t_aick ! t_ik^ca
      real(dp), dimension(:,:,:,:), allocatable :: t_aijd ! t_ij^ad
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ij_kb ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_k_ijb ! The above intermediate, reordered
      real(dp), dimension(:,:,:,:), allocatable :: X_id_kb ! An intermediate, term 2
      real(dp), dimension(:,:,:,:), allocatable :: X_ib_dk ! The above intermediate, reordered
      real(dp), dimension(:,:,:,:), allocatable :: X_ckb_j ! An intermediate, term 3
      real(dp), dimension(:,:,:,:), allocatable :: Y_ckb_j ! An intermediate, term 4
!      
      real(dp), dimension(:,:), allocatable :: X_bd    ! An intermediate, term 5
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aijb ! rho_aibj, batching over b
      real(dp), dimension(:,:,:,:), allocatable :: rho_ibaj ! rho_aibj, batching over b
      real(dp), dimension(:,:,:,:), allocatable :: rho_aib_j ! rho_aibj, batching over b
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij ! rho_aibj, batching over b
!
      integer :: rec1, rec0
      integer :: current_b_batch
!
      type(batching_index) :: batch_b

      integer :: b, i, j, a
!
      type(timings) :: jacobian_ccsd_d2_timer
!
      call jacobian_ccsd_d2_timer%init('jacobian ccsd d2')
      call jacobian_ccsd_d2_timer%start()
!
!     Initialize batching variable
!
      rec0 = (wf%n_o**2)*(wf%n_v**2) + wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = wf%n_v**2*wf%n_o + wf%n_v*wf%integrals%n_J&
            + max((wf%n_o)*(wf%n_v**2), 2*(wf%n_o**3), (wf%n_o**3) + (wf%n_o**2)*(wf%n_v),&
               2*(wf%n_o)*(wf%n_v**2), 2*(wf%n_o**2)*(wf%n_v) )

!
      call batch_b%init(wf%n_v)
      call mem%batch_setup(batch_b, rec0, rec1)
!
!     Start looping over b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Get batching limits for current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
!        :: Term 1. - sum_kcd g_kcbd t_ij^cd c_ak ::
!
!        Form g_kc_db = g_kcbd
!
         call mem%alloc(g_bdkc, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_bdkc,                        &
                           batch_b%first, batch_b%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
!        Reorder g_bd_kc to g_cd_kb (= g_kcbd), i.e. 1234 to 4231
!
         call mem%alloc(g_cdkb, wf%n_v, wf%n_v, wf%n_o, batch_b%length)
!
         call sort_1234_to_4231(g_bdkc, g_cdkb, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_bdkc, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
!        Order amplitudes as t_ij_cd = t_ij^cd = t_ci_dj
!                              2413                1234
!
         call mem%alloc(t_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
         call squareup_and_sort_1234_to_2413(wf%t2, t_ijcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form intermediate X_ij_kb = sum_cd g_kcdb t_ij^cd
!                                  = sum_cd t_ij_cd g_cd_kb
!
         call mem%alloc(X_ij_kb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
!
         call dgemm('N', 'N',                   &
                     (wf%n_o)**2,               &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_v)**2,               &
                     one,                       &
                     t_ijcd,                    & ! t_ij_cd
                     (wf%n_o)**2,               &
                     g_cdkb,                    & ! g_cd_kb
                     (wf%n_v)**2,               &
                     zero,                      &
                     X_ij_kb,                   & ! X_ij_kb
                     (wf%n_o)**2)
!
         call mem%dealloc(t_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!        sum_kcd g_kcbd t_ij^cd c_ak = sum_k X_ij_kb c_ak
!        Reorder to X_k_ijb = X_ij_kb
!
         call mem%alloc(X_k_ijb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
         call sort_1234_to_3124(X_ij_kb, X_k_ijb, wf%n_o, wf%n_o, wf%n_o, (batch_b%length))
         call mem%dealloc(X_ij_kb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
!
!        Form rho_aijb = - sum_k c_ak X_k_ijb = - sum_k c_ai(a,k) X_k_ijb(k, ijb)
!
         call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
!
         call dgemm('N', 'N',                      &
                     wf%n_v,                       &
                     (batch_b%length)*(wf%n_o)**2, &
                     wf%n_o,                       &
                     -one,                         &
                     c_ai,                         & ! c_a_k
                     wf%n_v,                       &
                     X_k_ijb,                      & ! X_k_ijb
                     wf%n_o,                       &
                     zero,                         &
                     rho_aijb,                     & ! rho_aijb
                     wf%n_v)
!
         call mem%dealloc(X_k_ijb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
!
!        Add rho_aijb (batch over b) to rho_aibj (full space)
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do b = 1, batch_b%length
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_aijb(a,i,j,b)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Deallocations for term 1 (keep g_cd_kb = g_kcbd)
!
         call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
!
!        :: Term 2. - sum_kcd g_kcbd t_kj^ad c_ci ::
!
!        Form the intermediate X_i_dkb = sum_c g_kcbd c_ci
!                                      = sum_c c_ci g_cd_kb
!                                      = sum_c c_ai^T(i,c) g_cd_kb(c, dkb)
!
         call mem%alloc(X_id_kb, wf%n_o, wf%n_v, wf%n_o, batch_b%length)
!
         call dgemm('T','N',                             &
                     wf%n_o,                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_v,                             &
                     one,                                &
                     c_ai,                               & ! c_i_c
                     wf%n_v,                             &
                     g_cdkb,                             & ! g_c_dkb
                     wf%n_v,                             &
                     zero,                               &
                     X_id_kb,                            & ! X_i_dkb
                     wf%n_o)
!
!        sum_kcd g_kcbd t_kj^ad c_ci = sum_kd (sum_c c_ci g_kcbd) t_kj^ad
!                                    = sum_kd X_idkb t_kj^ad
!                                    = sum_kd X_ib_dk t_dk_aj
!
!        Reorder to X_ib_dk = X_id_kb
!
         call mem%alloc(X_ib_dk, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
!
         call sort_1234_to_1423(X_id_kb, X_ib_dk, wf%n_o, wf%n_v, wf%n_o, batch_b%length)
!
         call mem%dealloc(X_id_kb, wf%n_o, wf%n_v, wf%n_o, batch_b%length)
!
!        Order the amplitudes as t_dk_aj = t_kj^ad = t_jk^da (t_dj, ak)
!
         call mem%alloc(t_dkaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call squareup_and_sort_1234_to_1432(wf%t2, t_dkaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Calculate rho_ib_aj = - sum_kcd g_kcbd t_kj^ad c_ci
!                            = - sum_dk X_ib_dk t_dk_aj
!
         call mem%alloc(rho_ibaj, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_o)*(wf%n_v),         &
                     (wf%n_o)*(wf%n_v),         &
                     -one,                      &
                     X_ib_dk,                   & ! X_ib_dk
                     (wf%n_o)*(batch_b%length), &
                     t_dkaj,                    & ! t_dk_aj
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     rho_ibaj,                  & ! rho_ib_aj
                     (wf%n_o)*(batch_b%length))
!
         call mem%dealloc(X_ib_dk, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
         call mem%dealloc(t_dkaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Add rho_ibaj (batch over b) ro rho_aibj (full space)
!
!$omp parallel do private(b, j, i, a)
         do a = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, batch_b%length
                  do i = 1, wf%n_o
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_ibaj(i,b,a,j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Deallocations for term 2 (keep g_cd_kb = g_kcbd)
!
         call mem%dealloc(rho_ibaj, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
!
!        :: Term 3. - sum_kcd g_kcbd t_ik^ca c_dj ::
!
!        sum_d g_kcbd c_dj = sum_d g_cd_kb c_dj
!
!        Reorder integrals to g_cd_kb to g_ck_bd
!
         call mem%alloc(g_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
!
         call sort_1234_to_1342(g_cdkb, g_ckbd, wf%n_v, wf%n_v, wf%n_o, batch_b%length)
!
         call mem%dealloc(g_cdkb, wf%n_v, wf%n_v, wf%n_o, batch_b%length)
!
!        Form the intermediate X_ckb_j = sum_d g_kcbd c_dj = sum_d g_ckb_d c_d_j
!
         call mem%alloc(X_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_ckbd,                             & ! g_ckb_d
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_d_j
                     wf%n_v,                             &
                     zero,                               &
                     X_ckb_j,                            & ! X_ckb_j
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ca = t_ki^ac (t_ak,ci)
!
         call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call squareup_and_sort_1234_to_1432(wf%t2, t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form rho_aib_j = -sum_kcd g_kcbd t_ik^ca c_dj = sum_ck t_ai_ck X_ckb_j
!
         call mem%alloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_v)*(wf%n_o),         &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_v)*(wf%n_o),         &
                     -one,                      &
                     t_aick,                    & ! t_ai_ck
                     (wf%n_v)*(wf%n_o),         &
                     X_ckb_j,                   & ! X_ck_bj
                     (wf%n_v)*(wf%n_o),         &
                     zero,                      &
                     rho_aib_j,                  & ! rho_aibj
                     (wf%n_v)*(wf%n_o))
!
         call mem%dealloc(X_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
         call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Add rho_aibj to rho_aibj
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do b = 1, batch_b%length
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_aib_j(a,i,b,j)
!
                 enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Deallocations for term 3 (keep g_ckb_d = g_kcbd)
!
         call mem%dealloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
!        :: Term 4.  sum_kcd L_kcbd t_ik^ac c_dj ::
!
!        sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj
!
!        Form L_ckb_d = L_kcbd = 2 * g_kcbd - g_kdbc = 2 * g_ckb_d(ckb, d) - g_ckb_d(dkb, c)
!
         call mem%alloc(L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         L_ckbd = zero
!
!        Note: exchange g_ck_bd(dk,bc) -> 4231
!
         call add_4231_to_1234(-one, g_ckbd, L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         call daxpy((wf%n_v)**2*(wf%n_o)*(batch_b%length), two, g_ckbd, 1, L_ckbd, 1)
!
         call mem%dealloc(g_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
!
!        Form the intermediate Y_ckb_j = sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj
!
         call mem%alloc(Y_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     L_ckbd,                             & ! L_ckb_d
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_d_j
                     wf%n_v,                             &
                     zero,                               &
                     Y_ckb_j,                            & ! Y_ckb_j
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ac = t2(aick)
!
         call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call squareup(wf%t2, t_aick, (wf%n_o)*(wf%n_v))
!
!        Form rho_aib_j =  sum_ck t_ai_ck Y_ckb_j
!
!        Note: we interpret Y_ckb_j as Y_ck_bj in the matrix multiplication
!        Note: we interpret rho_aib_j as rho_aibj in the matrix multiplication
!
         call mem%alloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(wf%n_v),         &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_o)*(wf%n_v),         &
                     one,                       &
                     t_aick,                    & ! t_ai_ck
                     (wf%n_o)*(wf%n_v),         &
                     Y_ckb_j,                   & ! Y_ck_bj
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     rho_aib_j,                 & ! rho_ai_bj
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(Y_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
         call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Add rho_aib_j to rho_aibj
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do b = 1, batch_b%length
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_aib_j(a,i,b,j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Deallocations for term 4 (keep L_ckb_d = L_kcbd)
!
         call mem%dealloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
!        :: Term 5.  sum_kcd L_kcbd t_ij^ad c_ck ::
!
!        Form the intermediate X_1,bd = sum_ck c_ck L_kcbd = sum_ck c_1,ck L_ckb_d
!
!        Note: c_ai is interpreted as c_1,ai in the matrix multiplication
!
         call mem%alloc(X_bd, 1, (batch_b%length)*(wf%n_v))
!
         call dgemm('N','N',                    &
                     1,                         &
                     (batch_b%length)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v),         &
                     one,                       &
                     c_ai,                      & ! c_ck
                     1,                         &
                     L_ckbd,                    & ! L_ck_bd
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     X_bd,                      & ! X_bd
                     1)
!
!        Order amplitudes as t_ai_jd = t_ij^ad
!
         call mem%alloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
         call squareup_and_sort_1234_to_1243(wf%t2, t_aijd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form rho_b_aij =  sum_kcd L_kcbd t_ij^ad c_ck
!                       =  sum_d X_bd t_d_aij
!
!        Note: X_bd is interpreted as X_b_d in the matrix multiplication
!
         call mem%alloc(rho_baij, batch_b%length, wf%n_v, wf%n_o, wf%n_o)
!
         call dgemm('N','T',               &
                     batch_b%length,       &
                     (wf%n_v)*(wf%n_o)**2, &
                     wf%n_v,               &
                     one,                  &
                     X_bd,                 & ! X_b_d
                     batch_b%length,       &
                     t_aijd,               & ! t_aij_d
                     (wf%n_v)*(wf%n_o)**2, &
                     zero,                 &
                     rho_baij,             & ! rho_b_aij
                     batch_b%length)
!
         call mem%dealloc(X_bd, 1, (batch_b%length)*(wf%n_v))
         call mem%dealloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!        Add rho_baij to rho_aibj
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
                  do b = 1, batch_b%length
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_baij(b,a,i,j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Final deallocations in batching loop
!
         call mem%dealloc(rho_baij, batch_b%length, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(L_ckbd, wf%n_v, wf%n_o,batch_b%length, wf%n_v)
!
      enddo ! End of batches over b
!
      call jacobian_ccsd_d2_timer%freeze()
      call jacobian_ccsd_d2_timer%switch_off()
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
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aick
!
      real(dp), dimension(:,:,:,:), allocatable :: t_dlbj
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: L_ckdl
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj
!
      type(timings) :: jacobian_ccsd_e2_timer
!
      call jacobian_ccsd_e2_timer%init('jacobian ccsd e2')
      call jacobian_ccsd_e2_timer%start()
!
      call mem%alloc(t_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_dlbj, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Construct g_kcld
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
!     Construct L_kc,ld ordered as L_ck_dl
!
      call mem%alloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_ckdl = zero
!
!     L_ck_dl(ck, dl) = two*g_kcld(kc, ld) - g_kcld(kd, lc)
!             1234                  2143              2341
!
      call add_2341_to_1234(-one, g_kcld, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2143_to_1234(two, g_kcld, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Intermediate X_ck_bj = sum_dl L_ck_dl * t_dl_bj
!
      call mem%alloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ckdl,            & ! L_ck_dl
                  (wf%n_o)*(wf%n_v), &
                  t_dlbj,            & ! t_dl_bj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj = 2 * sum_ck c_ai_ck * X_ck_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  c_aick,            & ! c_ai_ck
                  (wf%n_o)*(wf%n_v), &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_ai_bj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call jacobian_ccsd_e2_timer%freeze()
      call jacobian_ccsd_e2_timer%switch_off()
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^F2 =   - sum_ckld t_ai,ck * L_kc,ld * c_bl,dj
!!                        - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                        - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ckdl
      real(dp), dimension(:,:,:,:), allocatable :: L_dlck
      real(dp), dimension(:,:,:,:), allocatable :: L_lckd
!
      real(dp), dimension(:,:,:,:), allocatable :: c_dlbj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aick
      real(dp), dimension(:,:,:,:), allocatable :: t_aijd
      real(dp), dimension(:,:,:,:), allocatable :: t_aibl
 !
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_l_j
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aijb
!
      type(timings) :: jacobian_ccsd_f2_timer
!
      call jacobian_ccsd_f2_timer%init('jacobian ccsd f2')
      call jacobian_ccsd_f2_timer%start()
!
!     :: Construct L_kc,ld ::
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
      call mem%alloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_ckdl = zero
!
!     Construct L_kc,ld ordered as L_ck_dl
!
!     L_ck_dl(ck, dl) = L_kcld = 2 * g_kcld(kc, ld) - g_kcld(kd, lc)
!                                            2143              2341
!
      call add_2341_to_1234(-one, g_kcld, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2143_to_1234(two, g_kcld, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 1: - sum_ckld t_aick * L_kcld * c_bldj ::
!
!     Reorder c_bl_dj as c_dl_bj
!
      call mem%alloc(c_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_3214(c_aibj, c_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_ckbj = sum_dl L_ckdl*c_dlbj = sum_dl L_kcld*c_bldj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ckdl,            & ! L_ck_dl
                  (wf%n_o)*(wf%n_v), &
                  c_dlbj,            & ! c_dl_bj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_aick, (wf%n_o)*(wf%n_v))
!
!     rho_aibj = sum_ck t_ai_ck*X_ckbj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  t_aick,            & ! t_ai_ck
                  (wf%n_o)*(wf%n_v), &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_ai_bj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2: - sum_ckdl t_aidj * L_kcld * c_blck
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
!     Construct L_ckdl reordered as L_dlck
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_dlck = zero
!
!     L_dlck(d,l,c,k) =- g_kcld(k,d,l,c) (4123)
!
      call add_4123_to_1234(-one, g_kcld, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     L_dlck(d,l,c,k) =+ 2*g_kcld(k,c,l,d) (4321)
!
      call add_4321_to_1234(two, g_kcld, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Y_d_b = sum_clk L_d_lck * c_b_lck
!     Here dgemm is tricked to believe that c_bl_ck is c_b_lck
!
      call mem%alloc(Y_d_b, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  wf%n_v,               &
                  ((wf%n_o)**2)*wf%n_v, &
                  one,                  &
                  L_dlck,               & ! L_d_lck
                  wf%n_v,               &
                  c_aibj,               & ! c_b_lck
                  wf%n_v,               &
                  zero,                 &
                  Y_d_b,                & ! Y_d_b
                  wf%n_v)
!
      call mem%dealloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Reorder T2 amplitudes
!
      call squareup_and_sort_1234_to_1243(wf%t2, t_aijd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_aij_b = sum_d t_aijd*Y_d_b
!
      call dgemm('N','N',                 &
                  ((wf%n_o)**2)*(wf%n_v), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  t_aijd,                 & ! t_aij_d
                  ((wf%n_o)**2)*(wf%n_v), &
                  Y_d_b,                  &
                  wf%n_v,                 &
                  zero,                   &
                  rho_aijb,               & ! rho_aij_b
                  ((wf%n_o)**2)*(wf%n_v))
!
      call mem%dealloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     Adding term 2 to rho_aibj
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 3: - sum_ckdl t_aibl * L_kcld * c_ckdj ::
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
      call mem%alloc(L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      L_lckd = zero
!
!     Construct L_kcld ordered as L_lckd
!
!     L_lckd = L_kcld = 2 * g_kcld - g_kdlc
!
      call add_3412_to_1234(-one, g_kcld, L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     L_lckd(lc, kd) =+ two * g_kcld (3214)
!
      call add_3214_to_1234(two, g_kcld, L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(Z_l_j, wf%n_o, wf%n_o)
!
!     Z_l_j = sum_ckd L_l_ckd * c_ckd_l
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  one,                  &
                  L_lckd,               & ! L_l_ckd
                  wf%n_o,               &
                  c_aibj,               & ! c_aibj(ck,dl)= c_ckd_l
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_l_j,                &
                  wf%n_o)
!
      call mem%dealloc(L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(t_aibl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibl, wf%n_o*(wf%n_v))
!
!     rho_aibj = sum_l t_aibl * Z_l_j
!
      call dgemm('N','N',                 &
                  ((wf%n_v)**2)*(wf%n_o), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  t_aibl,                 & ! t_aib_l
                  ((wf%n_v)**2)*(wf%n_o), &
                  Z_l_j,                  & ! Z_l_j
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               & ! rho_ai_bj
                  ((wf%n_v)**2)*(wf%n_o))
!
      call mem%dealloc(t_aibl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(Z_l_j, wf%n_o, wf%n_o)
!
      call jacobian_ccsd_f2_timer%freeze()
      call jacobian_ccsd_f2_timer%switch_off()
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
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ckdl
      real(dp), dimension(:,:,:,:), allocatable :: L_dclk
      real(dp), dimension(:,:,:,:), allocatable :: L_lckd
!
      real(dp), dimension(:,:,:,:), allocatable :: c_aijd
!
      real(dp), dimension(:,:,:,:), allocatable :: t_dlbj
      real(dp), dimension(:,:,:,:), allocatable :: t_clkb
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdj
 !
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_l_j
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aijb
!
      type(timings) :: jacobian_ccsd_g2_timer
!
      call jacobian_ccsd_g2_timer%init('jacobian ccsd g2')
      call jacobian_ccsd_g2_timer%start()
!
!     :: Term 1: - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck  ::
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
      call mem%alloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_ckdl = zero
!
!     Construct L_kcld = two*g_kcld - g_kdlc ordered as L_ckdl
!
!     L_ckdl(c,k,d,l) =- g_kcld(k,d,l,c) (2341)
!
      call add_2341_to_1234(-one, g_kcld, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     L_ck_dl(c,k,d,l) =+ two*g_kcld(k,c,l,d) (2143)
!
      call add_2143_to_1234(two, g_kcld, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Reorder t_bldj as t_dlbj
!
      call mem%alloc(t_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup_and_sort_1234_to_3214(wf%t2, t_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_ckbj = sum_dl t_bldj * L_kcld = sum_dl L_ckdl t_dlbj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ckdl,            & ! L_ck_dl
                  (wf%n_o)*(wf%n_v), &
                  t_dlbj,            & ! t_dl_bj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj =+ - sum_ck c_aick X_ckbj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_aibj,            & ! c_ai_bj
                  (wf%n_o)*(wf%n_v), &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_aibj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2: - sum_ckdl t_ckbl * L_kcld * c_aidj
!
!     Reorder L_ckdl to L_dclk
!
      call mem%alloc(L_dclk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_3142(L_ckdl, L_dclk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder t_ckbl as t_clkb
!
      call mem%alloc(t_clkb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call squareup_and_sort_1234_to_1423(wf%t2, t_clkb, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Y_d_b = sum_clk L_dclk * c_clkb
!
      call mem%alloc(Y_d_b, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  wf%n_v,               &
                  ((wf%n_o)**2)*wf%n_v, &
                  one,                  &
                  L_dclk,               & ! L_d_clk
                  wf%n_v,               &
                  t_clkb,               & ! t_clk_b
                  ((wf%n_o)**2)*wf%n_v, &
                  zero,                 &
                  Y_d_b,                & ! Y_d_b
                  wf%n_v)
!
      call mem%dealloc(t_clkb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(L_dclk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(c_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Reorder c_aidj to c_aijd
!
      call sort_1234_to_1243(c_aibj, c_aijd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_aijb = sum_d c_aijd * Y_d_b
!
      call dgemm('N','N',                 &
                  ((wf%n_o)**2)*(wf%n_v), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  c_aijd,                 & ! c_aij_d
                  ((wf%n_o)**2)*(wf%n_v), &
                  Y_d_b,                  & ! Y_d_b
                  wf%n_v,                 &
                  zero,                   &
                  rho_aijb,               & ! rho_aij_b
                  ((wf%n_o)**2)*(wf%n_v))
!
      call mem%dealloc(c_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     Adding term 2 to rho_aibj
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 3: - sum_ckld t_ckdj * L_kcld * c_aibl ::
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
      call mem%alloc(L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      L_lckd = zero
!
!     Construct L_kcld ordered as L_lckd
!
!     L_lckd =- g_kcld(k,d,l,c) (3412)
!
      call add_3412_to_1234(-one, g_kcld, L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     L_lckd =+ two*g_kcld(k,c,l,d) (3214)
!
      call add_3214_to_1234(two, g_kcld, L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Reorder t_ck,dj to t_ckd_j
!
      call mem%alloc(t_ckdj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckdj, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(Z_l_j, wf%n_o, wf%n_o)
!
!     Z_l_j = sum_ckd L_lckd*t_ckdj
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  one,                  &
                  L_lckd,               & ! L_l_ckd
                  wf%n_o,               &
                  t_ckdj,               & ! t_ckd_j
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_l_j,                & ! Z_l_j
                  wf%n_o)
!
      call mem%dealloc(L_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(t_ckdj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj = sum_l c_aib_l*Z_l_j
!
      call dgemm('N','N',                 &
                  ((wf%n_v)**2)*(wf%n_o), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  c_aibj,                 & ! c_ai_bj
                  ((wf%n_v)**2)*(wf%n_o), &
                  Z_l_j,                  & ! Z_l_j
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               & ! rho_aib_j
                  ((wf%n_v)**2)*(wf%n_o))
!
      call mem%dealloc(Z_l_j, wf%n_o, wf%n_o)
!
      call jacobian_ccsd_g2_timer%freeze()
      call jacobian_ccsd_g2_timer%switch_off()
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
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_lckd
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aikc
      real(dp), dimension(:,:,:,:), allocatable :: t_ajlc
!
      real(dp), dimension(:,:,:,:), allocatable :: c_ldbj
      real(dp), dimension(:,:,:,:), allocatable :: c_kdbi
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aild
      real(dp), dimension(:,:,:,:), allocatable :: Y_ajkd
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      type(timings) :: jacobian_ccsd_h2_timer
!
      call jacobian_ccsd_h2_timer%init('jacobian ccsd h2')
      call jacobian_ccsd_h2_timer%start()
!
!     :: Term 1: sum_ckld t_ci,ak * g_kc,ld * c_bl,dj ::
!
!     Construct g_kcld
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
!     t_ak,ci ordered as t_aikc
!
      call mem%alloc(t_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call squareup_and_sort_1234_to_1423(wf%t2, t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     X_ai_ld = sum_ck t_aikc*g_kcld
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_aikc,            & ! t_ai_kc
                  (wf%n_o)*(wf%n_v), &
                  g_kcld,            & ! g_kc_ld
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_aild,            & ! X_ai_ld
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(c_ldbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Reorder c_bldj as c_ldbj
!
      call sort_1234_to_2314(c_aibj, c_ldbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += sum_ld X_aild*c_ldbj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  X_aild,            & ! X_ai_ld
                  (wf%n_o)*(wf%n_v), &
                  c_ldbj,            & ! c_ld_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_ai_bj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_ldbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(X_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 2: sum_ckdl t_cjal * g_kcld * c_bkdi
!
!     Construct g_kcld
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
!        Reorder g_kcld to g_lckd
!
      call mem%alloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call sort_1234_to_3214(g_kcld, g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     t_alcj ordered as t_ajlc
!
      call mem%alloc(t_ajlc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call squareup_and_sort_1234_to_1423(wf%t2, t_ajlc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Y_ajkd = sum_lc t_ajlc * g_lckd
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_ajlc,            & ! t_aj_lc
                  (wf%n_o)*(wf%n_v), &
                  g_lckd,            & ! g_lckd
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ajkd,            & ! Y_aj_kd
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(t_ajlc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Reorder c_bkdi as c_kdbi
!
      call mem%alloc(c_kdbi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call sort_1234_to_2314(c_aibj, c_kdbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ajbi = sum_kd  Y_ajkd * c_kdbi
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  Y_ajkd,            & ! Y_aj_kd
                  (wf%n_o)*(wf%n_v), &
                  c_kdbi,            & ! c_kd_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  rho_ajbi,          & ! rho_aj_bi
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_kdbi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Reorder into rho_aibj
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call jacobian_ccsd_h2_timer%freeze()
      call jacobian_ccsd_h2_timer%switch_off()
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
!!    Batch over c to construct  g_ki_bc
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: c_aijc
      real(dp), dimension(:,:,:,:), allocatable :: c_aick
      real(dp), dimension(:,:,:,:), allocatable :: c_ajck
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aijb
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bjkc
      real(dp), dimension(:,:,:,:), allocatable :: g_bckj
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbj ! reordering of g_bj_kc and g_bc_kj
!
      type(timings) :: jacobian_ccsd_i2_timer
!
      call jacobian_ccsd_i2_timer%init('jacobian ccsd i2')
      call jacobian_ccsd_i2_timer%start()
!
!     :: sum_c F_bc * c_ai,cj ::
!
!     Reorder c_ai,cj to c_ai_jc
!
      call mem%alloc(c_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_1243(c_aibj, c_aijc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_aibj += sum_c F_bc * c_aicj = sum_c c_aijc F_ab(b,c)
!
      call dgemm('N','T',                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  one,                    &
                  c_aijc,                 & ! c_aij_c
                  (wf%n_v)*((wf%n_o)**2), &
                  wf%fock_ab,             & ! F_c_b
                  wf%n_v,                 &
                  zero,                   &
                  rho_aijb,               & ! rho_aij_b
                  (wf%n_v)*((wf%n_o)**2))
!
      call mem%dealloc(c_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Reorder rho_aijb into rho_aibj
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: - sum_k F_jk * c_aibk  ::
!
!     rho_aibj += - sum_k F_jk * c_aibk = - sum_k c_aibk F_ij(k,j)^T
!
      call dgemm('N', 'N',                &
                  (wf%n_o)*((wf%n_v)**2), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  c_aibj,                 & ! c_aib_k
                  (wf%n_o)*((wf%n_v)**2), &
                  wf%fock_ij,             & ! F_k_j
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               & ! rho_aib_j
                  (wf%n_o)*((wf%n_v)**2))
!
!     ::   sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj ) ::
!
!     sum_ck ( g_bj,kc*(2*c_ai,ck - c_ak,ci) - g_bc,kj*c_ai,ck - g_ki,bc*c_ak,cj )
!
!     Construct g_bj,kc
!
      call mem%alloc(g_bjkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_voov(g_bjkc)
!
!     Reordering g_bj_kc to g_ck_bj
!
      call mem%alloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_4312(g_bjkc, g_ckbj, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_bjkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     rho_aibj += sum_ck 2*c_aick * g_ckbj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  c_aibj,            & ! c_ai_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            & ! g_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_ai_bj
                  (wf%n_o)*(wf%n_v))
!
!     Reorder c_akci to c_aick
!
      call mem%alloc(c_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(c_aibj, c_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += - sum_ck g_ckbj*c_aick
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_aick,            & ! c_ai_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            & ! g_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_aibj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder g_bckj to g_ckbj
!
      call mem%alloc(g_bckj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_vvoo(g_bckj)
!
      call mem%alloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2314(g_bckj, g_ckbj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_bckj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     rho_aibj += - sum_ck c_aick * g_ckbj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_aibj,            & ! c_ai_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            & ! g_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_aibj,          & ! rho_ai_bj
                  (wf%n_o)*(wf%n_v))
!
!     Reorder  c_ak,cj to c_aj_ck
!
      call mem%alloc(c_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(c_aibj, c_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_ajck,            & ! c_aj_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            &  ! g_ck_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  rho_ajbi,          & ! rho_aj_bi
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder rho_ajbi into rho_aibj
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call jacobian_ccsd_i2_timer%freeze()
      call jacobian_ccsd_i2_timer%switch_off()
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
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: X_klij
!
      type(timings) :: jacobian_ccsd_j2_timer
!
      call jacobian_ccsd_j2_timer%init('jacobian ccsd j2')
      call jacobian_ccsd_j2_timer%start()
!
!     Constructing g_kcld
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
      call mem%alloc(g_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Reorder g_kcld to g_kl_cd
!
      call sort_1234_to_1324(g_kcld, g_klcd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Reordered T2 amplitudes
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_kl_ij = g_kl_cd * t_cd_ij
!
      call mem%alloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',     &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_klcd,      & ! g_kl_cd
                  (wf%n_o)**2, &
                  t_abij,      & ! t_cd_ij
                  (wf%n_v)**2, &
                  zero,        &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2)
!
!     rho_abij += c_ab_kl * X_kl_ij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  c_abij,      & ! c_ab_kl
                  (wf%n_v)**2, &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  rho_abij,    & ! rho_ab_ij
                  (wf%n_v)**2)
!
!     X_kl_ij = g_kl_cd * c_cd_ij
!
      call dgemm('N', 'N',     &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_klcd,      & ! g_kl_cd
                  (wf%n_o)**2, &
                  c_abij,      & ! c_cd_ij
                  (wf%n_v)**2, &
                  zero,        &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2)
!
!     rho_abij += t_abkl * X_klij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_abij,      & ! t_ab_kl
                  (wf%n_v)**2, &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  rho_abij,    & ! rho_ab_ij
                  (wf%n_v)**2)
!
      call mem%dealloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(g_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call jacobian_ccsd_j2_timer%freeze()
      call jacobian_ccsd_j2_timer%switch_off()
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
!!                       + sum_cd g_acbd * c_cidj
!!
!!    For the last term we batch over a and b and
!!    add each batch to rho_aibj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kilj
      real(dp), dimension(:,:,:,:), allocatable :: g_klij
      real(dp), dimension(:,:,:,:), allocatable :: g_acbd
      real(dp), dimension(:,:,:,:), allocatable :: g_abcd
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_batch_abij
!
      integer :: a = 0, b = 0, i = 0, j = 0
!
!     Batching and memory handling variables
!
      integer :: current_a_batch = 0
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
      integer :: rec0, rec1_a, rec1_b, rec2
!
      type(timings) :: jacobian_ccsd_k2_timer
!
      call jacobian_ccsd_k2_timer%init('jacobian ccsd k2')
      call jacobian_ccsd_k2_timer%start()
!
      call mem%alloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%get_oooo(g_kilj)
!
!     Reorder g_kilj to g_klij
!
      call mem%alloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(g_kilj, g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     rho_abij += sum_kl g_kilj * c_akbl = sum_kl c_abij(a,b,k,l) g_klij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  c_abij,      & ! c_ab_kl
                  (wf%n_v)**2, &
                  g_klij,      & ! g_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  rho_abij,    & ! rho_ab_ij
                  (wf%n_v)**2)
!
      call mem%dealloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Prepare for batching over a and b
!
!     ::  sum_cd g_ac,bd * c_ci,dj ::
!
      rec0 = 0
!
      rec1_a = wf%integrals%n_J*wf%n_v
      rec1_b = wf%integrals%n_J*wf%n_v
!
      rec2 = wf%n_v**2
!
!     Initialize batching variables
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_b, rec0, rec1_a, rec1_b, rec2)
!
!     Start looping over a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_b_batch = 1, batch_b%num_batches
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(g_acbd, (wf%n_v),(batch_a%length), (wf%n_v),(batch_b%length))
!
!           g_ca_db = sum_J L_ca_J*L_db_J
!
            call wf%get_vvvv(g_acbd,                     &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v,                    &
                           batch_b%first, batch_b%last,  &
                           1, wf%n_v)
!
!           sum_cd g_ac,bd * c_ci,dj = sum_cd g_ac,bd c_cd,ij = sum_cd g_abcd c_cd_ij
!
!           Reorder g_acbd into g_abcd (i.e., 1234 to 1324)
!
            call mem%alloc(g_abcd, (batch_a%length),(batch_b%length), (wf%n_v), (wf%n_v))
!
            call sort_1234_to_1324(g_acbd, g_abcd, batch_a%length, wf%n_v, batch_b%length, wf%n_v)
!
            call mem%dealloc(g_acbd, (wf%n_v),(batch_a%length), (wf%n_v),(batch_b%length))
!
            call mem%alloc(rho_batch_abij, batch_a%length, batch_b%length, wf%n_o, wf%n_o)
!
!           rho_abij += sum_cd g_acbd * c_cidj = sum_cd g_abcd(a,b,c,d) c_abij(c,d,i,j)
!
            call dgemm('N', 'N',                            &
                        (batch_a%length)*(batch_b%length),  &
                        (wf%n_o)**2,                        &
                        (wf%n_v)**2,                        &
                        one,                                &
                        g_abcd,                             & ! g_ab_cd
                        (batch_a%length)*(batch_b%length),  &
                        c_abij,                             & ! c_cd_ij
                        (wf%n_v)**2,                        &
                        zero,                               &
                        rho_batch_abij,                     & ! rho_ab_ij
                        (batch_a%length)*(batch_b%length))
!
            call mem%dealloc(g_abcd, (batch_a%length),(batch_b%length), (wf%n_v), (wf%n_v))
!
!           Reorder into rho_abij
!
!$omp parallel do private(b, a, i, j)
            do j = 1, wf%n_o
               do i = 1, wf%n_o
                  do b = 1, batch_b%length
                     do a = 1, batch_a%length
!
                        rho_abij(batch_a%first + a - 1, batch_b%first + b - 1, i, j) =                &
                                          rho_abij(batch_a%first + a - 1, batch_b%first + b - 1, i, j)&
                                          + rho_batch_abij(a,b,i,j)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(rho_batch_abij, batch_a%length, batch_b%length, wf%n_o, wf%n_o)
!
         enddo ! End batches of b
      enddo ! End batches of a
!
      call jacobian_ccsd_k2_timer%freeze()
      call jacobian_ccsd_k2_timer%switch_off()
!
   end subroutine jacobian_ccsd_k2_ccsd
!
!
end submodule jacobian_ccsd
