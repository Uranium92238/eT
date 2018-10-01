submodule (ccsd_class) jacobian_transpose_ccsd
!
!!
!!    Jacobian transpose submodule (CCSD)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    jacobian_transpose_ccsd_transformation: performs the transformation by the CCSD
!!                                            Jacobian transpose matrix A^T, placing the result in the
!!                                            incoming vector. 
!!
!!    jacobian_transpose_ccsd_x1:             adds the X1 term to the transformed singles vector; 
!!                                            x = a, b, c, ..., g 
!!    jacobian_transpose_ccsd_x2:             adds the X2 term to the transformed doubles vector; 
!!                                            x = a, b, ..., i
!! 
!
   implicit none 
!
!
contains
!
!
   module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD A1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the A1 term,
!!
!!       sum_ckdl b_ck L_iald u_kl^cd,
!! 
!!    abd adds it to the transformed vector sigma_a_i.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i 
!
      real(dp), dimension(:,:), allocatable :: u_ld_ck
!
      real(dp), dimension(:,:), allocatable :: t_ld_ck
!
      real(dp), dimension(:,:), allocatable :: X_ld ! An intermediate, see below
!
      real(dp), dimension(:,:), allocatable :: g_ia_ld ! g_iald 
      real(dp), dimension(:,:), allocatable :: L_ai_ld ! L_iald 
!
      integer(i15) :: k = 0, c = 0, d = 0, l = 0, ck = 0, dk = 0, dl = 0
      integer(i15) :: ld = 0, cl = 0, ckdl = 0, cldk = 0, i = 0, a = 0
      integer(i15) :: id = 0, la = 0, ia = 0, ai = 0
!
!     Read the amplitudes from disk 
! 
      !call wf%read_double_amplitudes 
!
!     Form u_ld_ck = u_kl^cd = 2 * t_kl^cd - t_lk^cd = 2 * t_ld_ck(ld, ck) - t_ld_ck(kd, cl)
!
      call mem%alloc(t_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call squareup_and_sort_1234_to_4312(wf%t2, t_ld_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      u_ld_ck = zero
!
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_ld_ck, 1, u_ld_ck, 1)
      call add_4231_to_1234(-one, t_ld_ck, u_ld_ck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_ld = sum_ck u_ld_ck b_ck  
!
      call mem%alloc(X_ld, (wf%n_v)*(wf%n_o), 1)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  u_ld_ck,           &
                  (wf%n_v)*(wf%n_o), &
                  b_a_i,             & ! "b_ck"
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ld,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form L_ai_ld = L_iald = 2 * g_iald - g_idla 
!                           = 2 * g_ia_ld(ia,ld) - g_ia_ld(id,la)
!
      call mem%alloc(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_ia_ld)
!
      call mem%alloc(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ai_ld = zero
!
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, g_ia_ld, 1, L_ai_ld, 1)
      call add_1432_to_1234(-one, g_ia_ld, L_ai_ld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate and add sum_ckdl b_ck L_iald u_kl^cd 
!                       = sum_ld L_ai_ld X_ld
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  X_ld,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_a_i,         & ! "sigma_ai"
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(X_ld, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD B1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the B1 term,
!!
!!       - sum_ckdl (b_al L_kcid t_kl^cd + b_ci L_ldka t_kl^cd),
!! 
!!    abd adds it to the transformed vector sigma_a_i.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i 
!
      real(dp), dimension(:,:), allocatable :: g_kc_id ! g_kcid 
      real(dp), dimension(:,:), allocatable :: L_kcd_i ! L_kcid 
!
      real(dp), dimension(:,:), allocatable :: t_l_kcd ! t_kl^cd 
!
      real(dp), dimension(:,:), allocatable :: L_a_ldk ! L_ldka
      real(dp), dimension(:,:), allocatable :: t_ldk_c ! t_kl^cd 
!
      real(dp), dimension(:,:), allocatable :: X_l_i ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_a_c ! An intermediate, term 2
!
!     Indices 
!
      integer(i15) :: i = 0, d = 0, c = 0, k = 0, id = 0, ic = 0, kd = 0, kc = 0, kcd = 0
      integer(i15) :: ck = 0, dl = 0, ckdl = 0, l = 0, ldk = 0, lda = 0, a = 0
!
!     :: Term 1. - sum_ckdl b_al L_kcid t_kl^cd ::
!
!     Form L_kcd_i = L_kcid = 2 * g_kcid - g_kdic
!                           = 2 * g_kc_id(kc,id) - g_kc_id(kd,ic)
!
      call mem%alloc(g_kc_id, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_id)
!
      call mem%alloc(L_kcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      L_kcd_i = zero
!
!     COMMENT: Maybe we should do the shifting of underscore in dgemm?
!$omp parallel do schedule(static) private(i,d,id,c,ic,k,kd,kc,kcd)
      do i = 1, wf%n_o
         do d = 1, wf%n_v
!
            id = index_two(i, d, wf%n_o)
!
            do c = 1, wf%n_v
!
               ic = index_two(i, c, wf%n_o)
!
               do k = 1, wf%n_o
!
                  kd = index_two(k, d, wf%n_o)
                  kc = index_two(k, c, wf%n_o)
!
                  kcd = index_three(k, c, d, wf%n_o, wf%n_v)
!
                  L_kcd_i(kcd, i) = two*g_kc_id(kc, id) - g_kc_id(kd, ic)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_kc_id, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form t_l_kcd = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_l_kcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      t_l_kcd = zero
!
!     COMMENT: Do an unpacking with sort, and shift underscore?
!$omp parallel do schedule(static) private(d,c,k,ck,kcd,l,dl,ckdl)
      do d = 1, wf%n_v
         do c = 1, wf%n_v
            do k = 1, wf%n_o
!
               ck = index_two(c, k, wf%n_v)
!
               kcd = index_three(k, c, d, wf%n_o, wf%n_v)
!
               do l = 1, wf%n_o
!
                  dl = index_two(d, l, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
!
                  t_l_kcd(l, kcd) = wf%t2(ckdl, 1) ! t_kl^cd 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      !call wf%destruct_double_amplitudes
!
!     Calculate the intermediate X_l_i = sum_kcd t_l_kcd L_kcd_i
!
      call mem%alloc(X_l_i, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_l_kcd,              &
                  wf%n_o,               &
                  L_kcd_i,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_l_i,                &
                  wf%n_o)
!
!     Add - sum_ckdl b_al L_kcid t_kl^cd = sum_l b_al X_l_i
!
      call dgemm('N','N',    &
                  wf%n_v,    &
                  wf%n_o,    &
                  wf%n_o,    &
                  -one,      &
                  b_a_i,     & ! b_al 
                  wf%n_v,    &
                  X_l_i,     &
                  wf%n_o,    &
                  one,       &
                  sigma_a_i, &
                  wf%n_v)
!
      call mem%dealloc(X_l_i, wf%n_o, wf%n_o)
!
!     :: Term 2. - sum_ckdl b_ci L_ldka t_kl^cd ::
!
!     Form L_a_ldk = L_ldka = L_kcd_i(lda,k)
!
      call mem%alloc(L_a_ldk, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      L_a_ldk = zero
!
!$omp parallel do schedule(static) private(k,d,l,ldk,a,lda)
      do k = 1, wf%n_o
         do d = 1, wf%n_v
            do l = 1, wf%n_o
!
               ldk = index_three(l, d, k, wf%n_o, wf%n_v)
!
               do a = 1, wf%n_v
!
                  lda = index_three(l, d, a, wf%n_o, wf%n_v)
!
                  L_a_ldk(a, ldk) = L_kcd_i(lda, k) ! L_ldka
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(L_kcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Form t_ldk_c = t_kl^cd = t_l_kcd(l, kcd)
!
      call mem%alloc(t_ldk_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_ldk_c = zero
!
!$omp parallel do schedule(static) private(c,k,d,kcd,l,ldk)
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do d = 1, wf%n_v
!
               kcd = index_three(k, c, d, wf%n_o, wf%n_v)
!
               do l = 1, wf%n_o
!
                  ldk = index_three(l, d, k, wf%n_o, wf%n_v)
!
                  t_ldk_c(ldk, c) = t_l_kcd(l, kcd) ! t_kl^cd 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(t_l_kcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Calculate the intermediate X_a_c = sum_ldk L_a_ldk t_ldk_c
!
      call mem%alloc(X_a_c, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  L_a_ldk,              &
                  wf%n_v,               &
                  t_ldk_c,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_a_c,                &
                  wf%n_v)
!
      call mem%dealloc(L_a_ldk, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call mem%dealloc(t_ldk_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     Add - sum_ckdl b_ci L_ldka t_kl^cd = - sum_c X_a_c b_ci 
!
      call dgemm('N','N',    &
                  wf%n_v,    &
                  wf%n_o,    &
                  wf%n_v,    &
                  -one,      &
                  X_a_c,     &
                  wf%n_v,    &
                  b_a_i,     & ! b_ci
                  wf%n_v,    &
                  one,       &
                  sigma_a_i, &
                  wf%n_v)
!
      call mem%dealloc(X_a_c, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_b1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_c1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD C1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the C1 term,
!!
!!       sum_cdl b_cidl g_dlca - sum_kdl b_akdl g_dlik,
!! 
!!    and adds it to the transformed vector sigma_a_i.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: g_dl_ca ! g_dlca 
      real(dp), dimension(:,:), allocatable :: g_a_dlc ! g_dlca 
!
      real(dp), dimension(:,:), allocatable :: b_dlc_i ! b_cidl 
!
      real(dp), dimension(:,:), allocatable :: g_ik_dl ! g_dlik = g_ikdl 
      real(dp), dimension(:,:), allocatable :: g_kdl_i ! g_dlik
!
      integer(i15) :: c = 0, l = 0, d = 0, i = 0, dl = 0, dlc = 0, ci = 0
      integer(i15) :: a = 0, ca = 0, k = 0, kdl = 0, ik = 0
!
!     Batching variables 
!
      integer(i15) :: required 
      integer(i15) :: current_a_batch 
!
      type(batching_index) :: batch_a 
!
!     :: Term 1. sum_cdl b_cidl g_dlca :: 
!
!     Reorder b_ci_dl = b_cidl to b_dlc_i
!
      call mem%alloc(b_dlc_i, (wf%n_o)*(wf%n_v)**2,  wf%n_o)
      b_dlc_i = zero
!
!$omp parallel do schedule(static) private(c,l,d,dl,dlc,i,ci)
      do c = 1, wf%n_v
         do l = 1, wf%n_o
            do d = 1, wf%n_v
!
               dl = index_two(d, l, wf%n_v)
!
               dlc = index_three(d, l, c, wf%n_v, wf%n_o)
!
               do i = 1, wf%n_o
!
                  ci = index_two(c, i, wf%n_v)
!
                  b_dlc_i(dlc, i) = b_ai_bj(ci, dl) ! b_cidl 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Prepare batching over index a 
!
      required = wf%integrals%get_required_vvvo()
!     
      call batch_a%init(wf%n_v)
      call mem%num_batch(batch_a, required)        
!
!     Loop over the a-batches 
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index 
!
         call batch_a%determine_limits(current_a_batch)
!
!        Form g_dl_ca
!
         call mem%alloc(g_dl_ca, (wf%n_v)*(wf%n_o), (wf%n_v)*batch_a%length)
!
         call wf%get_vovv(g_dl_ca,        &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last)
!
!        Add sum_dlc g_dlc_a^T b_dlc_i
!
         call dgemm('T','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     g_dl_ca,                    & ! "g_dlc_a"
                     (wf%n_o)*(wf%n_v)**2,       &
                     b_dlc_i,                    &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     sigma_a_i(batch_a%first,1), &
                     wf%n_v)
!
         call mem%dealloc(g_dl_ca, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_a%length))
!
      enddo ! End of batches over a 
!
      call mem%dealloc(b_dlc_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     :: Term 2. - sum_kdl b_akdl g_dlik = - sum_kdl b_akdl g_ikdl ::
!
!     Form g_ik_dl
!
      call mem%alloc(g_ik_dl, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_oovo(g_ik_dl)
!
!     Add - sum_kdl b_akdl g_dlik = - sum_kdl b_akdl g_kdl_i
!
!     Note: we interpret b_ai_bj as b_a_ibj, such that b_ai_bj(a,kdl) = b_akdl 
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_ai_bj,              & ! "b_a_kdl"
                  wf%n_v,               &
                  g_ik_dl,              & ! "g_i_kdl"
                  (wf%n_o),             &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call mem%dealloc(g_ik_dl, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_c1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD D1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the D1 term,
!!
!!       - sum_ckdl (b_ckal F_id t_kl^cd + b_ckdi F_la t_kl^cd),
!! 
!!    and adds it to the transformed vector sigma_a_i.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: t_lck_d ! t_kl^cd 
!
      real(dp), dimension(:,:), allocatable :: X_a_d ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_l_i ! An intermediate, term 2 
!
      integer(i15) :: l = 0, k = 0, c = 0, ck = 0, a = 0, al = 0
      integer(i15) :: d = 0, dl = 0, ckdl = 0, ckd = 0, lck = 0
!
!     :: Term 1. - sum_ckdl b_ckal F_id t_kl^cd ::
!
!     Read amplitudes and order as t_lck_d = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_lck_d = zero 
!
!$omp parallel do schedule(static) private(d,l,dl,k,c,ck,ckdl,lck)
      do d = 1, wf%n_v
         do l = 1, wf%n_o
!
            dl = index_two(d, l, wf%n_v)
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
!
                  lck = index_three(l, c, k, wf%n_o, wf%n_v)
!
                  t_lck_d(lck, d) = wf%t2(ckdl, 1) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_a_d = sum_ckl b_a_lck t_lck_d = sum_ckl b_ckal t_kl^cd
!  
      call mem%alloc(X_a_d, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_ai_bj,              &  ! b_a_lck = b_alck = b_ckal
                  wf%n_v,               &
                  t_lck_d,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_a_d,                &
                  wf%n_v)
!
!     Add - sum_ckdl b_ckal F_id t_kl^cd
!           = - sum_d X_a_d F_id 
!           = - sum_d X_a_d F_i_a^T(d,i)
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  X_a_d,      &
                  wf%n_v,     &
                  wf%fock_ia, & ! F_i_a
                  wf%n_o,     &
                  one,        &
                  sigma_a_i,  &
                  wf%n_v)
!
      call mem%dealloc(X_a_d, wf%n_v, wf%n_v)
!
!     :: Term 2. - sum_ckdl b_ckdi F_la t_kl^cd
!
!     Form the intermediate X_l_i = sum_ckd t_l_ckd b_ckd_i  = sum_ckd b_ckdi t_kl^cd 
!
!     Note: we interpret b_ai_bj as b_aib_j, such that b_aib_j(ckd, i) = b_ckdi
!           we interpret t_lck_d as t_l_ckd, such that t_l_ckd(l,ckd) = t_kl^cd 
!
      call mem%alloc(X_l_i, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_lck_d,              & ! "t_l_ckd"
                  wf%n_o,               &
                  b_ai_bj,              & ! "b_ckd_i"
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_l_i,                &
                  wf%n_o)
!
      call mem%dealloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     Add - sum_ckdl b_ckdi F_la t_kl^cd = - sum_l F_la X_l_i = - sum_l F_i_a^T(a,l) X_l_i(l,i)
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%fock_ia, &
                  wf%n_o,     &
                  X_l_i,      &
                  wf%n_o,     &
                  one,        &
                  sigma_a_i,  &
                  wf%n_v)
!
      call mem%dealloc(X_l_i, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd
end submodule jacobian_transpose_ccsd