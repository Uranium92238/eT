submodule (ccsd_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    jacobian_transpose_ccsd_transformation: performs the transposed Jacobian transformation (A^T)
!! 
!
   implicit none 
!
!
contains
!
!
   module subroutine jacobian_transpose_ccsd_transformation_ccsd(wf, b_a_i, b_aibj)
!!
!!    Jacobian transpose transformation (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
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
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
      real(dp), dimension(wf%n_t2am, 1)   :: b_aibj 
!
      real(dp), dimension(:,:), allocatable :: sigma_a_i
      real(dp), dimension(:,:), allocatable :: sigma_ai_bj
!
      real(dp), dimension(:,:), allocatable :: b_ai_bj ! b_aibj unpacked 
!
!     Allocate the transformed singles vector 
!
      call allocator(sigma_a_i, wf%n_v, wf%n_o)
      sigma_a_i = zero 
!
!     Calculate and add the CCS contributions to the 
!     singles transformed vector 
!
      call wf%jacobian_transpose_ccs_a1(sigma_a_i, b_a_i) 
      call wf%jacobian_transpose_ccs_b1(sigma_a_i, b_a_i) 
!
!     Calculate and add the CCSD contributions to the
!     singles transformed vector 
!
      call wf%jacobian_transpose_ccsd_a1(sigma_a_i, b_a_i) 
      call wf%jacobian_transpose_ccsd_b1(sigma_a_i, b_a_i) 
!
      call allocator(b_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      b_ai_bj = zero 
!
      call squareup(b_aibj, b_ai_bj, (wf%n_v)*(wf%n_o))
!
      call wf%jacobian_transpose_ccsd_c1(sigma_a_i, b_ai_bj) 
      call wf%jacobian_transpose_ccsd_d1(sigma_a_i, b_ai_bj) 
      call wf%jacobian_transpose_ccsd_e1(sigma_a_i, b_ai_bj) 
      call wf%jacobian_transpose_ccsd_f1(sigma_a_i, b_ai_bj)
      call wf%jacobian_transpose_ccsd_g1(sigma_a_i, b_ai_bj) 
!
!     Add the CCSD contributions to the doubles vector arising from 
!     the incoming singles vector  
!
      call allocator(sigma_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      sigma_ai_bj = zero 
!
      call wf%jacobian_transpose_ccsd_a2(sigma_ai_bj, b_a_i)
!
!     Copy the transformed singles over into the incoming singles vector
!
      b_a_i = sigma_a_i
!
!     End of routines - deallocations 
! 
      call deallocator(sigma_a_i, wf%n_v, wf%n_o)
      call deallocator(b_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
   end subroutine jacobian_transpose_ccsd_transformation_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD A1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
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
      real(dp), dimension(:,:), allocatable :: X_ld ! An intermediate, see below
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J 
      real(dp), dimension(:,:), allocatable :: g_ia_ld ! g_iald 
      real(dp), dimension(:,:), allocatable :: L_ai_ld ! L_iald 
!
      integer(i15) :: k = 0, c = 0, d = 0, l = 0, ck = 0, dk = 0, dl = 0
      integer(i15) :: ld = 0, cl = 0, ckdl = 0, cldk = 0, i = 0, a = 0, ld = 0
      integer(i15) :: id = 0, la = 0, ia = 0, ai = 0
!
!     Read the amplitudes from disk 
!
      call wf%initialize_amplitudes 
      call wf%read_double_amplitudes 
!
!     Form u_ld_ck = u_kl^cd = 2 * t_kl^cd - t_lk^cd 
!
      call allocator(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do d = 1, wf%n_v
!
               dk = index_two(d, k, wf%n_v)
!
               do l = 1, wf%n_o
!
                  dl = index_two(d, l, wf%n_v)
                  ld = index_two(l, d, wf%n_o)
                  cl = index_two(c, l, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
                  cldk = index_packed(cl, dk)
!
                  u_ld_ck(ld, ck) = two*(wf%t2am(ckdl,1)) - wf%t2am(cldk,1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_amplitudes
!
!     Form the intermediate X_ld = sum_ck u_ld_ck b_ck  
!
      call allocator(X_ld, (wf%n_v)*(wf%n_o), 1)
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
      call deallocator(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form L_ai_ld = L_iald = 2 * g_iald - g_idla 
!                           = 2 * g_ia_ld(ia,ld) - g_ia_ld(id,la)
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ia_J)
!
      call allocator(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_ld,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call allocator(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ai_ld = zero 
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o
!
            ld = index_two(l, d, wf%n_o)
!
            do i = 1, wf%n_o
!
               id = index_two(i, d, wf%n_o)
!
               do a = 1, wf%n_v
!
                  la = index_two(l, a, wf%n_o)
                  ia = index_two(i, a, wf%n_o)
                  ai = index_two(a, i, wf%n_v)
!
                  L_ai_ld(ai, ld) = two*g_ia_ld(ia, ld) - g_ia_ld(id, la)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call deallocator(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(X_ld, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD B1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
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
      real(dp), dimension(:,:), allocatable :: L_kc_J  ! L_kc^J 
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
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_kc_J)
!
      call allocator(g_kc_id, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_id,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call allocator(L_kcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      L_kcd_i = zero
!
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
!
      call deallocator(g_kc_id, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form t_l_kcd = t_kl^cd 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_l_kcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      t_l_kcd = zero
!
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
                  t_l_kcd(l, kcd) = wf%t2am(ckdl, 1) ! t_kl^cd 
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_amplitudes
!
!     Calculate the intermediate X_l_i = sum_kcd t_l_kcd L_kcd_i
!
      call allocator(X_l_i, wf%n_o, wf%n_o)
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
      call deallocator(X_l_i, wf%n_o, wf%n_o)
!
!     :: Term 2. - sum_ckdl b_ci L_ldka t_kl^cd ::
!
!     Form L_a_ldk = L_ldka = L_kcd_i(lda,k)
!
      call allocator(L_a_ldk, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      L_a_ldk = zero
!
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
!
      call deallocator(L_kcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Form t_ldk_c = t_kl^cd = t_l_kcd(l, kcd)
!
      call allocator(t_ldk_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_ldk_c = zero
!
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
!
      call deallocator(t_l_kcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Calculate the intermediate X_a_c = sum_ldk L_a_ldk t_ldk_c
!
      call allocator(X_a_c, wf%n_v, wf%n_v)
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
      call deallocator(L_a_ldk, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call deallocator(t_ldk_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
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
      call deallocator(X_a_c, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_b1_ccsd 
!
!
   module subroutine jacobian_transpose_ccsd_c1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD C1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
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
      real(dp), dimension(:,:), allocatable :: L_dl_J ! L_dl^J 
      real(dp), dimension(:,:), allocatable :: L_ca_J ! L_ca^J 
!
      real(dp), dimension(:,:), allocatable :: g_dl_ca ! g_dlca 
      real(dp), dimension(:,:), allocatable :: g_a_dlc ! g_dlca 
!
      real(dp), dimension(:,:), allocatable :: b_dlc_i ! b_cidl 
!
      real(dp), dimension(:,:), allocatable :: g_dl_ik ! g_dlik 
      real(dp), dimension(:,:), allocatable :: g_kdl_i ! g_dlik
!
      real(dp), dimension(:,:), allocatable :: L_ik_J ! L_ik^J
!
      integer(i15) :: c = 0, l = 0, d = 0, i = 0, dl = 0, dlc = 0, ci = 0
      integer(i15) :: a = 0, ca = 0, k = 0, kdl = 0, ik = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0 
      integer(i15) :: batch_dimension = 0, max_batch_length = 0, n_batch = 0
!
      integer(i15) :: a_batch = 0, a_length = 0, a_first = 0, a_last = 0
!
!     :: Term 1. sum_cdl b_cidl g_dlca :: 
!
!     Reorder b_ci_dl = b_cidl to b_dlc_i
!
      call allocator(b_dlc_i, (wf%n_o)*(wf%n_v)**2,  wf%n_o)
      b_dlc_i = zero
!
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
!
!     Prepare batching over index a 
!
      call allocator(L_dl_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call wf%get_cholesky_ai(L_dl_J)
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
!        Form g_dl_ca
!
         call allocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_ca_J, 1, wf%n_v, a_first, a_last)
!
         call allocator(g_dl_ca, (wf%n_v)*(wf%n_o), (wf%n_v)*a_length)
!
         call dgemm('N','T',            &
                     (wf%n_v)*(wf%n_o), & 
                     (wf%n_v)*a_length, &
                     wf%n_J,            &
                     one,               &
                     L_dl_J,            &
                     (wf%n_v)*(wf%n_o), &
                     L_ca_J,            &
                     (wf%n_v)*a_length, &
                     zero,              &
                     g_dl_ca,           &
                     (wf%n_v)*(wf%n_o))
!
         call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
!        Reorder g_dl_ca to g_a_dlc 
!
         call allocator(g_a_dlc, a_length, (wf%n_o)*(wf%n_v)**2)
!
         do c = 1, wf%n_v
            do l = 1, wf%n_o
               do d = 1, wf%n_v
!
                  dl = index_two(d, l, wf%n_v)
!
                  dlc = index_three(d, l, c, wf%n_v, wf%n_o)
!
                  do a = 1, a_length
!
                     ca = index_two(c, a, wf%n_v)
!
                     g_a_dlc(a, dlc) = g_dl_ca(dl, ca) ! g_dlca 
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_dl_ca, (wf%n_v)*(wf%n_o), (wf%n_v)*a_length)
!
!        Add sum_dlc g_a_dlc b_dlc_i
!
         call dgemm('N','N',               &
                     a_length,             &
                     wf%n_o,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     g_a_dlc,              &
                     a_length,             &
                     b_dlc_i,              &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     sigma_a_i(a_first,1), &
                     wf%n_v)
!
         call deallocator(g_a_dlc, a_length, (wf%n_o)*(wf%n_v)**2)
!
      enddo ! End of batches over a 
!
      call deallocator(b_dlc_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     :: Term 2. - sum_kdl b_akdl g_dlik ::
!
!     Form g_dl_ik and reorder to g_kdl_i
!
      call allocator(L_ik_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ij(L_ik_J)
!
      call allocator(g_dl_ik, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)**2,       &
                  wf%n_J,            &
                  one,               &
                  L_dl_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ik_J,            &
                  (wf%n_o)**2,       &
                  zero,              &
                  g_dl_ik,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ik_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_dl_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Reorder to g_kdl_i 
!
      call allocator(g_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      g_kdl_i = zero
!
      do i = 1, wf%n_o
         do l = 1, wf%n_o
            do d = 1, wf%n_v
!
               dl = index_two(d, l, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ik = index_two(i, k, wf%n_o)
!  
                  kdl = index_three(k, d, l, wf%n_o, wf%n_v)
!
                  g_kdl_i(kdl, i) = g_dl_ik(dl, ik) ! g_dlik 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_dl_ik, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
!     Add - sum_kdl b_akdl g_dlik = - sum_kdl b_akdl g_kdl_i
!
!     Note: we interpret b_ai_bj as b_a_ibj, such that b_ai_bj(a,kdl) = b_akdl 
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_ai_bj,              & ! "b_a_ibj"
                  wf%n_v,               &
                  g_kdl_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call deallocator(g_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_c1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD D1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
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
      real(dp), dimension(:,:), allocatable :: b_a_ckl ! b_ckal 
      real(dp), dimension(:,:), allocatable :: t_ckl_d ! t_kl^cd 
!
      real(dp), dimension(:,:), allocatable :: X_a_d ! An intermediate, term 1
!
      real(dp), dimension(:,:), allocatable :: t_l_ckd ! t_kl^cd 
!
      real(dp), dimension(:,:), allocatable :: X_l_i ! An intermediate, term 2 
!
      integer(i15) :: l = 0, k = 0, c = 0, ck = 0, ckl = 0, a = 0, al = 0
      integer(i15) :: d = 0, dl = 0, ckdl = 0, ckd = 0
!
!
!     :: Term 1. - sum_ckdl b_ckal F_id t_kl^cd ::
!
!     Reorder b_ai_bj to b_a_ckl(a,ckl) = b_ckal 
!
      call allocator(b_a_ckl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      b_a_ckl = zero
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               ckl = index_three(c, k, l, wf%n_v, wf%n_o)
!
               do a = 1, wf%n_v
!
                  al = index_two(a, l, wf%n_v)
!
                  b_a_ckl(a, ckl) = b_ai_bj(ck, al) ! b_ckal
!
               enddo
            enddo
         enddo
      enddo
!
!     Read amplitudes and order as t_ckl_d = t_kl^cd 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_ckl_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_ckl_d = zero 
!
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
                  ckl = index_three(c, k, l, wf%n_v, wf%n_o)
!
                  t_ckl_d(ckl, d) = wf%t2am(ckdl, 1) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_amplitudes
!
!     Form the intermediate X_a_d = sum_ckl b_a_ckl t_ckl_d 
!  
      call allocator(X_a_d, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_a_ckl,              &
                  wf%n_v,               &
                  t_ckl_d,              &
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
      call deallocator(X_a_d, wf%n_v, wf%n_v)
      call deallocator(b_a_ckl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 2. - sum_ckdl b_ckdi F_la t_kl^cd
!
!     Order amplitudes as t_l_ckd = t_kl^cd = t_ckl_d
!
      call allocator(t_l_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      t_l_ckd = zero 
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
               do l = 1, wf%n_o
!
                  ckl = index_three(c, k, l, wf%n_v, wf%n_o)
!
                  t_l_ckd(l, ckd) = t_ckl_d(ckl, d)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(t_ckl_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     Form the intermediate X_l_i = sum_ckd t_l_ckd b_ckd_i  
!
!     Note: we interpret b_ai_bj as b_aib_j, such that b_aib_j(ckd, i) = b_ckdi
!
      call allocator(X_l_i, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_l_ckd,              &
                  wf%n_o,               &
                  b_ai_bj,              & ! "b_ckd_i"
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_l_i,                &
                  wf%n_o)
!
      call deallocator(t_l_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
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
      call deallocator(X_l_i, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD E1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the E1 term,
!!
!!       sum_ckdle (b_ckdi L_dale t_kl^ce + b_ckdl L_deia t_kl^ce)
!!      -sum_ckdlm (b_ckal L_ilmd t_km^cd + b_ckdl L_mlia t_km^cd)
!! 
!!    and adds it to the transformed vector sigma_a_i.
!!
!!    The routine adds the third and forth terms first.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: t_dm_ck ! t_km^cd 
      real(dp), dimension(:,:), allocatable :: t_m_ckd ! t_km^cd
      real(dp), dimension(:,:), allocatable :: t_el_ck ! t_lk^ec
      real(dp), dimension(:,:), allocatable :: t_ckl_e ! t_lk^ec 
!
      real(dp), dimension(:,:), allocatable :: L_il_J ! L_il^J 
      real(dp), dimension(:,:), allocatable :: L_md_J ! L_md^J
      real(dp), dimension(:,:), allocatable :: L_ml_J ! L_ml^J 
      real(dp), dimension(:,:), allocatable :: L_ia_J ! L_ia^J 
      real(dp), dimension(:,:), allocatable :: L_de_J ! L_de^J 
!
      real(dp), dimension(:,:), allocatable :: g_il_md ! g_ilmd
      real(dp), dimension(:,:), allocatable :: g_ml_ia ! g_mlia
!
      real(dp), dimension(:,:), allocatable :: L_il_dm ! L_ilmd 
      real(dp), dimension(:,:), allocatable :: L_ai_ml ! L_mlia
!
      real(dp), dimension(:,:), allocatable :: X_il_ck ! An intermediate, term 3
      real(dp), dimension(:,:), allocatable :: X_lck_i ! Reordered intermediate, term 3
!
      real(dp), dimension(:,:), allocatable :: X_m_l ! An intermediate, term 4
!
      real(dp), dimension(:,:), allocatable :: b_a_lck ! b_ckal
      real(dp), dimension(:,:), allocatable :: b_d_ckl ! b_ckdl
!
      real(dp), dimension(:,:), allocatable :: L_le_J ! L_le^J
      real(dp), dimension(:,:), allocatable :: L_da_J ! L_da^J 
!
      real(dp), dimension(:,:), allocatable :: g_da_le ! g_dale
      real(dp), dimension(:,:), allocatable :: L_a_eld ! L_dale
!
      real(dp), dimension(:,:), allocatable :: g_de_ia ! g_deia 
      real(dp), dimension(:,:), allocatable :: L_ai_de ! L_deia 
!
      real(dp), dimension(:,:), allocatable :: X_el_di ! An intermediate, term 1
!
      real(dp), dimension(:,:), allocatable :: X_d_e ! An intermediate, term 2
!
      integer(i15) :: ml = 0, md = 0, ma = 0, m = 0, lck = 0, l = 0, il = 0, ck = 0, i = 0
      integer(i15) :: dm = 0, ckd = 0, k = 0, id = 0, ia = 0, c = 0, d = 0, al = 0, ai = 0
      integer(i15) :: a = 0, le = 0, la = 0, eld = 0, e = 0, de = 0, da = 0, ie = 0, el = 0
      integer(i15) :: dl = 0, ckl = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0 
      integer(i15) :: batch_dimension = 0, max_batch_length = 0, n_batch = 0
!
      integer(i15) :: a_batch = 0, a_length = 0, a_first = 0, a_last = 0
!
      integer(i15) :: e_batch = 0, e_length = 0, e_first = 0, e_last = 0
!
!     :: Term 3. - sum_ckdlm b_ckal L_ilmd t_km^cd ::
!
!     Read the amplitudes from disk 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_dm_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_dm_ck = zero
!
      call squareup(wf%t2am, t_dm_ck, (wf%n_o)*(wf%n_v)) ! t_dm_ck(dm,ck) = t_mk^dc = t_km^cd 
!
      call wf%destruct_amplitudes
!
!     Form g_il_md = g_ilmd 
!
      call allocator(L_il_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ij(L_il_J)
!
      call allocator(L_md_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_md_J)
!
      call allocator(g_il_md, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_il_J,            &
                  (wf%n_o)**2,       &
                  L_md_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_il_md,           &
                  (wf%n_o)**2)
!
      call deallocator(L_il_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_md_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Form L_il_dm = L_ilmd = 2 * g_ilmd - g_idml 
!                           = 2 * g_ilmd - g_mlid
!                           = 2 * g_il_md(il,md) - g_il_md(ml,id)
!
      call allocator(L_il_dm, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      L_il_dm = zero
!
      do m = 1, wf%n_o
         do d = 1, wf%n_v
!
            dm = index_two(d, m, wf%n_v)
            md = index_two(m, d, wf%n_o)
!
            do l = 1, wf%n_o
!
               ml = index_two(m, l, wf%n_o)
!
               do i = 1, wf%n_o
!
                  id = index_two(i, d, wf%n_o)
                  il = index_two(i, l, wf%n_o)
!
                  L_il_dm(il, dm) = two*g_il_md(il, md) - g_il_md(ml, id) 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_il_md, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_il_ck = sum_md L_ilmd t_km^cd 
!                                   = sum_md L_ilmd t_mk^dc 
!                                   = sum_md L_il_dm t_dm_ck
!
      call allocator(X_il_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_il_dm,           &
                  (wf%n_o)**2,       &
                  t_dm_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_il_ck,           &
                  (wf%n_o)**2)
!
      call deallocator(L_il_dm, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Order X_il_ck as X_lck_i 
!
      call allocator(X_lck_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      X_lck_i = zero 
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do l = 1, wf%n_o
!
               lck = index_three(l, c, k, wf%n_o, wf%n_v)
!
               do i = 1, wf%n_o
!
                  il = index_two(i, l, wf%n_o)
!
                  X_lck_i(lck, i) = X_il_ck(il, ck)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(X_il_ck, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Reorder to b_a_lck = b_ckal
!
      call allocator(b_a_lck, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      b_a_lck = zero
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do l = 1, wf%n_o
!
               lck = index_three(l, c, k, wf%n_o, wf%n_v)
!
               do a = 1, wf%n_v
!
                  al = index_two(a, l, wf%n_v)
!
                  b_a_lck(a, lck) = b_ai_bj(ck, al)
!
               enddo
            enddo
         enddo
      enddo
!
!     Add - sum_ckdlm b_ckal L_ilmd t_km^cd
!         = - sum_ckl b_a_lck X_lck_i
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_a_lck,              &
                  wf%n_v,               &
                  X_lck_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call deallocator(X_lck_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      call deallocator(b_a_lck, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 4. - sum_ckdlm b_ckdl L_mlia t_km^cd ::
!
!     sum_ckd t_km^cd b_ckd_l = sum_ckd t_dm_ck b_ckd_l
!
!     Reorder t_dm_ck to t_m_ckd
!
      call allocator(t_m_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      t_m_ckd = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
               do m = 1, wf%n_o
!
                  dm = index_two(d, m, wf%n_v)
!
                  t_m_ckd(m, ckd) = t_dm_ck(dm, ck) ! t_mk^dc = t_km^cd 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(t_dm_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_m_l = sum_ckd t_km^cd b_ckdl
!                                 = sum_ckd t_m_ckd b_ckd_l
!
!     Note: we interpret b_ai_bj as b_aib_j, such that b_aib_j(ckd,l) = b_ckdl
!
      call allocator(X_m_l, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_m_ckd,              &
                  wf%n_o,               &
                  b_ai_bj,              & ! "b_aib_j"
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_m_l,                &
                  wf%n_o)
!
      call deallocator(t_m_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Form g_ml_ia = g_mlia 
!
      call allocator(L_ml_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ij(L_ml_J)
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ia_J)
!
      call allocator(g_ml_ia, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)**2,       & 
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ml_J,            &
                  (wf%n_o)**2,       &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ml_ia,           &
                  (wf%n_o)**2)
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call deallocator(L_ml_J, (wf%n_o)**2, wf%n_J)
!
!     Form L_ai_ml = L_mlia = 2 * g_mlia - g_mail
!                           = 2 * g_mlia - g_ilma 
!
      call allocator(L_ai_ml, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      L_ai_ml = zero
!
      do l = 1, wf%n_o
         do m = 1, wf%n_o
!
            ml = index_two(m, l, wf%n_o)
!
            do i = 1, wf%n_o
!
               il = index_two(i, l, wf%n_o)
!
               do a = 1, wf%n_v
!
                  ma = index_two(m, a, wf%n_o)
                  ia = index_two(i, a, wf%n_o)
                  ai = index_two(a, i, wf%n_v)
!
                  L_ai_ml(ai, ml) = two*g_ml_ia(ml, ia) - g_ml_ia(il, ma)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ml_ia, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Add - sum_ckdlm b_ckdl L_mlia t_km^cd = - sum_lm L_ai_ml X_m_l 
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)**2,       &
                  -one,              &
                  L_ai_ml,           &
                  (wf%n_o)*(wf%n_v), &
                  X_m_l,             & ! "X_ml"
                  (wf%n_o)**2,       &
                  one,               &
                  sigma_a_i,         & ! "sigma_ai"
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ai_ml, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      call deallocator(X_m_l, wf%n_o, wf%n_o)
!
!     :: Term 1. sum_ckdle b_ckdi L_dale t_kl^ce :: 
!
!     Read amplitudes from disk 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_el_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_el_ck = zero 
!
      call squareup(wf%t2am, t_el_ck, (wf%n_o)*(wf%n_v))
!
      call wf%destruct_amplitudes
!
!     Form the intermediate X_el_di = sum_ck t_kl^ce b_ckdi = sum_ck t_lk^ec b_ckdi
!                                   = sum_ck t_el_ck b_ck_di  
!
      call allocator(X_el_di, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_el_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  b_ai_bj,           & ! "b_ck_di"
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_el_di,           &
                  (wf%n_o)*(wf%n_v))
!
!     sum_dle L_dale X_el_di 
! ... reorder L_dale to L_a_eld & interpret X_el_di as X_eld_i 
!
      call allocator(L_le_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_le_J)
!
!     Prepare batching over index a 
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1
!
!        Form g_da_le 
!
         call allocator(L_da_J, (wf%n_v)*a_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_da_J, 1, wf%n_v, a_first, a_last)
!
         call allocator(g_da_le, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
!
         call dgemm('N','T',            &
                     (wf%n_v)*a_length, & 
                     (wf%n_o)*(wf%n_v), &
                     wf%n_J,            &
                     one,               &
                     L_da_J,            &
                     (wf%n_v)*a_length, &
                     L_le_J,            &
                     (wf%n_v)*(wf%n_o), &
                     zero,              &
                     g_da_le,           &
                     (wf%n_v)*a_length)
!
         call deallocator(L_da_J, (wf%n_v)*a_length, wf%n_J)
!
!        Form  L_a_eld = L_dale = 2 * g_dale - g_dela 
!                               = 2 * g_da_le(da, le) - g_da_le(de, la)      
!
         call allocator(L_a_eld, a_length, (wf%n_o)*(wf%n_v)**2)
         L_a_eld = zero 
! 
         do d = 1, wf%n_v
            do l = 1, wf%n_o
               do e = 1, wf%n_v
!     
                  de = index_two(d, e, wf%n_v)
                  le = index_two(l, e, wf%n_o)
!
                  eld = index_three(e, l, d, wf%n_v, wf%n_o)
!
                  do a = 1, a_length
!
                     la = index_two(l, a, wf%n_o)
                     da = index_two(d, a, wf%n_v)
!
                     L_a_eld(a, eld) = two*g_da_le(da, le) - g_da_le(de, la) ! L_dale
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_da_le, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
!
!        Add sum_ckdle b_ckdi L_dale t_kl^ce
!            = sum_eld L_a_eld X_el_di
!
!        Note: we interpret X_el_di as X_eld_i 
!
         call dgemm('N','N',               &
                     a_length,             &
                     wf%n_o,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     L_a_eld,              &
                     a_length,             &
                     X_el_di,              & ! "X_eld_i"
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     sigma_a_i(a_first,1), &
                     wf%n_v)
!
         call deallocator(L_a_eld, a_length, (wf%n_o)*(wf%n_v)**2)
!
      enddo ! End of batches over a 
!
      call deallocator(L_le_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     :: Term 2. sum_ckdle b_ckdl L_deia t_kl^ce ::
!
!     sum_ckl b_ckdl t_kl^ce = sum_ckl b_d_ckl t_ckl_e
!
!     Reorder t_el_ck to t_ckl_e (= t_kl^ce)
!
      call allocator(t_ckl_e, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_ckl_e = zero 
!
      do e = 1, wf%n_v
         do l = 1, wf%n_o
!
            el = index_two(e, l, wf%n_v)
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
!  
                  ckl = index_three(c, k, l, wf%n_v, wf%n_o)
!
                  t_ckl_e(ckl, e) = t_el_ck(el, ck) ! t_kl^ce 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(t_el_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Order b_ck_dl as b_d_ckl
!
      call allocator(b_d_ckl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      b_d_ckl = zero 
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               ckl = index_three(c, k, l, wf%n_v, wf%n_o)
!
               do d = 1, wf%n_v
!
                  dl = index_two(d, l, wf%n_v)
!
                  b_d_ckl(d, ckl) = b_ai_bj(ck, dl) ! b_ckdl
!
               enddo
            enddo
         enddo
      enddo
!
!     Form the intermediate X_d_e = sum_ckl b_ckdl t_kl^ce = sum_ckl b_d_ckl t_ckl_e
!
      call allocator(X_d_e, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_d_ckl,              &
                  wf%n_v,               &
                  t_ckl_e,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_d_e,                &
                  wf%n_v)
!
      call deallocator(b_d_ckl, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      call deallocator(t_ckl_e, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     sum_ckdle b_ckdl L_deia t_kl^ce = sum_de L_deia X_d_e
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Prepare batching over index e
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do e_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(e_first, e_last, e_batch, max_batch_length, batch_dimension)
         e_length = e_last - e_first + 1    
!
!        Form g_de_ia = g_deia 
!
         call allocator(L_de_J, (wf%n_v)*e_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_de_J, 1, wf%n_v, e_first, e_last)
!
         call allocator(g_de_ia, (wf%n_v)*e_length, (wf%n_v)*(wf%n_o))
!
         call dgemm('N','T',            &
                     (wf%n_v)*e_length, & 
                     (wf%n_v)*(wf%n_o), &
                     wf%n_J,            &
                     one,               &
                     L_de_J,            &
                     (wf%n_v)*e_length, &
                     L_ia_J,            &
                     (wf%n_v)*(wf%n_o), &
                     zero,              &
                     g_de_ia,           &
                     (wf%n_v)*e_length)
!
         call deallocator(L_de_J, (wf%n_v)*e_length, wf%n_J)
!
!        Form L_ai_de = L_deia = 2 * g_deia - g_daie
!                              = 2 * g_de_ia(de,ia) - g_de_ia(da,ie)
!
!        E: This will not work when batching... a is not restricted; e is restricted...!!
!
         call allocator(L_ai_de, (wf%n_o)*(wf%n_v), (wf%n_v)*e_length)
         L_ai_de = zero
!
         do e = 1, e_length
            do d = 1, wf%n_v
!
               de = index_two(d, e, wf%n_v)
!
               do i = 1, wf%n_o
!
                  ie = index_two(i, e, wf%n_o)
!
                  do a = 1, wf%n_v
!
                     da = index_two(d, a, wf%n_v)
                     ia = index_two(i, a, wf%n_o)
                     ai = index_two(a, i, wf%n_v)
!
                     L_ai_de(ai, de) = two*g_de_ia(de, ia) - g_de_ia(da, ie) ! L_deia
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_de_ia, (wf%n_v)*e_length, (wf%n_o)*(wf%n_v))
!
!        Calculate the contribution to the sum, 
!
!           sum_de L_ai_de X_d_e 
!
!        for the given batch of e:
!
         call dgemm('N','N',            &
                     (wf%n_v)*(wf%n_o), &
                     1,                 &
                     (wf%n_v)*e_length, &
                     one,               &
                     L_ai_de,           &
                     (wf%n_v)*(wf%n_o), &
                     X_d_e(1,e_first),  & ! Trick dgemm into thinking this is an X_de array,
                     (wf%n_v)*e_length, & ! with e restricted. 
                     one,               &
                     sigma_a_i,         &
                     (wf%n_v)*(wf%n_o))
!
         call deallocator(L_ai_de, (wf%n_o)*(wf%n_v), (wf%n_v)*e_length)
!
      enddo ! End of batches over e 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call deallocator(X_d_e, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD F1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the F1 term,
!!
!!       sum_ckdlm (b_akdl t_lm^cd g_ikmc + b_ckal t_ml^cd g_mkid + b_ckdi t_ml^cd g_mkla)
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
      real(dp), dimension(:,:), allocatable :: t_mc_dl ! t_lm^cd 
      real(dp), dimension(:,:), allocatable :: t_cl_dm ! t_ml^cd 
      real(dp), dimension(:,:), allocatable :: t_cd_ml ! t_ml^cd
!
      real(dp), dimension(:,:), allocatable :: b_ak_cl ! b_ckal 
      real(dp), dimension(:,:), allocatable :: b_ki_cd ! b_ckdi
!
      real(dp), dimension(:,:), allocatable :: L_ik_J ! L_ik^J 
      real(dp), dimension(:,:), allocatable :: L_mc_J ! L_mc^J
!
      real(dp), dimension(:,:), allocatable :: g_ik_mc ! g_ikmc 
      real(dp), dimension(:,:), allocatable :: g_kdm_i ! g_mkid
      real(dp), dimension(:,:), allocatable :: g_a_mkl ! g_mkla
!
      real(dp), dimension(:,:), allocatable :: X_ik_dl ! An intermediate, term 1 
      real(dp), dimension(:,:), allocatable :: X_kdl_i ! Reordered intermediate, term 1
!
      real(dp), dimension(:,:), allocatable :: X_ak_dm ! An intermediate, term 2
!
      real(dp), dimension(:,:), allocatable :: X_ki_ml ! An intermediate, term 3 
      real(dp), dimension(:,:), allocatable :: X_mkl_i ! Reordered intermediate, term 3
!
      integer(i15) :: m = 0, c = 0, mc = 0, k = 0, l = 0, kdl = 0, d = 0, i = 0, ik = 0
      integer(i15) :: dl = 0, dm = 0, cldm = 0, cl = 0, mkl = 0, ml = 0, mk = 0, lc = 0
      integer(i15) :: ki = 0, kdm = 0, kam = 0, id = 0, di = 0, ck = 0, cd = 0, al = 0
      integer(i15) :: ak = 0, a = 0
!
!     :: Term 1. sum_ckdlm b_akdl t_lm^cd g_ikmc :: 
!
!     X_ik_dl = sum_mc t_lm^cd g_ikmc = sum_mc g_ik_mc t_mc_dl
!
!     Order amplitudes as t_mc_dl = t_lm^cd 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
! 
      call allocator(t_mc_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_mc_dl = zero 
!
      do l = 1, wf%n_o
         do d = 1, wf%n_v
!
            dl = index_two(d, l, wf%n_v)
!
            do c = 1, wf%n_v
!
               cl = index_two(c, l, wf%n_v)
!
               do m = 1, wf%n_o
!
                  dm = index_two(d, m, wf%n_v)
                  mc = index_two(m, c, wf%n_o)
!
                  cldm = index_packed(cl, dm)
!
                  t_mc_dl(mc, dl) = wf%t2am(cldm, 1) ! t_lm^cd 
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_amplitudes
!
!     Form the integral g_ik_mc 
!
      call allocator(L_ik_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ij(L_ik_J)
!
      call allocator(L_mc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_mc_J)
!
      call allocator(g_ik_mc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)**2,       & 
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ik_J,            &
                  (wf%n_o)**2,       &
                  L_mc_J,            & 
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ik_mc,           &
                  (wf%n_o)**2)
!
      call deallocator(L_ik_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_mc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Form the intermediate X_ik_dl = sum_mc t_lm^cd g_ikmc = sum_mc g_ik_mc t_mc_dl
!
      call allocator(X_ik_dl, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       & 
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  g_ik_mc,           &
                  (wf%n_o)**2,       &
                  t_mc_dl,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ik_dl,           &
                  (wf%n_o)**2)
!
!     Reorder X_ik_dl to X_kdl_i
!
      call allocator(X_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      X_kdl_i = zero 
!
      do i = 1, wf%n_o
         do l = 1, wf%n_o
            do d = 1, wf%n_v
!
               dl = index_two(d, l, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ik = index_two(i, k, wf%n_o)
!
                  kdl = index_three(k, d, l, wf%n_o, wf%n_v)
!
                  X_kdl_i(kdl, i) = X_ik_dl(ik, dl)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(X_ik_dl, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Add sum_ckdlm b_akdl t_lm^cd g_ikmc
!         = sum_kdl b_a_kdl X_kdl_i
!
!     Note: we interpret b_ai_bj as b_a_ibj, such that b_a_ibj(a, kdl) = b_akdl
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_ai_bj,              & ! "b_a_ibj"
                  wf%n_v,               &
                  X_kdl_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call deallocator(X_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
!     :: Term 2. sum_ckdlm b_ckal t_ml^cd g_mkid ::
!
!     X_ak_dm = sum_cl b_ckal t_ml^cd
!             = sum_cl b_ak_cl t_cl_dm
!
!     We have t_mc_dl(mc,dl) = t_lm^cd 
!     Reorder t_cl_dm(cl,dm) = t_mc_dl(lc,dm) = t_ml^cd  
!
      call allocator(t_cl_dm, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      t_cl_dm = zero 
!
      do m = 1, wf%n_o
         do d = 1, wf%n_v
!
            dm = index_two(d, m, wf%n_v)
!
            do l = 1, wf%n_o
               do c = 1, wf%n_v
!
                  lc = index_two(l, c, wf%n_o)
                  cl = index_two(c, l, wf%n_v)
!
                  t_cl_dm(cl, dm) = t_mc_dl(lc, dm)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(t_mc_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder to b_ak_cl = b_ckal 
!
      call allocator(b_ak_cl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) 
      b_ak_cl = zero 
!
      do l = 1, wf%n_o
         do c = 1, wf%n_v
!
            cl = index_two(c, l, wf%n_v)
!
            do k = 1, wf%n_o
!
               ck = index_two(c, k, wf%n_v)
!
               do a = 1, wf%n_v
!
                  al = index_two(a, l, wf%n_v)
                  ak = index_two(a, k, wf%n_v)
!
                  b_ak_cl(ak, cl) = b_ai_bj(ck, al) ! b_ckal
!
               enddo
            enddo 
         enddo
      enddo  
!
!     Form the intermediate X_ak_dm = sum_cl b_ak_cl t_cl_dm
!
      call allocator(X_ak_dm, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), & 
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_ak_cl,           &
                  (wf%n_v)*(wf%n_o), &
                  t_cl_dm,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ak_dm,           &
                  (wf%n_v)*(wf%n_o))
!
      call deallocator(b_ak_cl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     sum_ckdlm b_ckal t_ml^cd g_mkid = sum_kdm X_ak_dm g_mkid
!
!     We have g_ik_mc(ik,mc) = g_ikmc 
!     Reorder to g_kdm_i(kdm,i) = g_mkid = g_ik_mc(mk, id)
!
      call allocator(g_kdm_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      g_kdm_i = zero 
!
      do i = 1, wf%n_o
         do m = 1, wf%n_o
            do d = 1, wf%n_v
!
               id = index_two(i, d, wf%n_o)
!
               do k = 1, wf%n_o
!
                  mk = index_two(m, k, wf%n_o)
!
                  kdm = index_three(k, d, m, wf%n_o, wf%n_v)
!
                  g_kdm_i(kdm, i) = g_ik_mc(mk, id) ! g_mkid
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ik_mc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Add sum_ckdlm b_ckal t_ml^cd g_mkid = sum_kdm X_ak_dm g_mkid
!                                         = sum_kdm X_ak_dm g_kdm_i
!
!     Note: we interpret X_ak_dm as X_a_kdm 
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_ak_dm,              & ! "X_a_kdm"
                  wf%n_v,               &
                  g_kdm_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call deallocator(X_ak_dm, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 3. sum_ckdlm b_ckdi t_ml^cd g_mkla ::
!
!     X_ki_ml = sum_cd b_ckdi t_ml^cd 
!
!     We have t_cl_dm(cl,dm) = t_ml^cd
!     Reorder to t_cd_ml(cd,ml) = t_cl_dm(cl,dm) = t_ml^cd 
!
      call allocator(t_cd_ml, (wf%n_v)**2, (wf%n_o)**2)
      t_cd_ml = zero 
!
      do l = 1, wf%n_o
         do m = 1, wf%n_o
!
            ml = index_two(m, l, wf%n_o)
!
            do d = 1, wf%n_v
!
               dm = index_two(d, m, wf%n_v)
!
               do c = 1, wf%n_v
!
                  cl = index_two(c, l, wf%n_v)
                  cd = index_two(c, d, wf%n_v)
!
                  t_cd_ml(cd, ml) = t_cl_dm(cl, dm) ! t_ml^cd 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(t_cl_dm, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder to b_ki_cd = b_ckdi 
!
      call allocator(b_ki_cd, (wf%n_o)**2, (wf%n_v)**2)
      b_ki_cd = zero 
!
      do d = 1, wf%n_v
         do c = 1, wf%n_v
!
            cd = index_two(c, d, wf%n_v)
!
            do i = 1, wf%n_o
!
               di = index_two(d, i, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ck = index_two(c, k, wf%n_v)
                  ki = index_two(k, i, wf%n_o)
!
                  b_ki_cd(ki, cd) = b_ai_bj(ck, di) ! b_ckdi
!
               enddo
            enddo
         enddo
      enddo
!
!     Form intermediate X_ki_ml = sum_cd b_ckdi t_ml^cd = sum_cd b_ki_cd t_cd_ml
!
      call allocator(X_ki_ml, (wf%n_o)**2, (wf%n_o)**2)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       & 
                  (wf%n_o)**2,       &
                  (wf%n_v)**2,       &
                  one,               &
                  b_ki_cd,           &
                  (wf%n_o)**2,       &
                  t_cd_ml,           &
                  (wf%n_v)**2,       &
                  zero,              &
                  X_ki_ml,           &
                  (wf%n_o)**2)
!
      call deallocator(t_cd_ml, (wf%n_v)**2, (wf%n_o)**2)
      call deallocator(b_ki_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!     sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_mkla X_ki_ml 
!
!     Reorder to X_mkl_i
!
      call allocator(X_mkl_i, (wf%n_o)**3, wf%n_o)
      X_mkl_i = zero 
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do m = 1, wf%n_o
!
               ml = index_two(m, l, wf%n_o)
!  
               mkl = index_three(m, k, l, wf%n_o, wf%n_o)
!
               do i = 1, wf%n_o
!
                  ki = index_two(k, i, wf%n_o)
!
                  X_mkl_i(mkl, i) = X_ki_ml(ki, ml)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(X_ki_ml, (wf%n_o)**2, (wf%n_o)**2)
!
!     We have g_kdm_i(kdm,i) = g_mkid
!     Reorder to g_a_mkl(a,mkl) = g_mkla = g_kdm_i(kam,l)
!
      call allocator(g_a_mkl, wf%n_v, (wf%n_o)**3)
      g_a_mkl = zero 
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do m = 1, wf%n_o
!
               mkl = index_three(m, k, l, wf%n_o, wf%n_o)
!
               do a = 1, wf%n_v
!
                  kam = index_three(k, a, m, wf%n_o, wf%n_v)
!
                  g_a_mkl(a, mkl) = g_kdm_i(kam, l) ! g_mkla
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_kdm_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
!     Add sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_a_mkl X_mkl_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  (wf%n_o)**3, &
                  one,         &
                  g_a_mkl,     &
                  wf%n_v,      &
                  X_mkl_i,     &
                  (wf%n_o)**3, &
                  one,         & 
                  sigma_a_i,   &
                  wf%n_v)
!
      call deallocator(g_a_mkl, wf%n_v, (wf%n_o)**3)
      call deallocator(X_mkl_i, (wf%n_o)**3, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD G1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the G1 term,
!!
!!       - sum_ckdle (b_akdl t_kl^ce g_icde + b_cidl t_kl^ce g_kade + b_cldi t_kl^ce g_keda)
!! 
!!    and adds it to the transformed vector sigma_a_i.
!!
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
      real(dp), dimension(:,:), allocatable :: L_ka_J ! L_ka^J 
      real(dp), dimension(:,:), allocatable :: L_de_J ! L_de^J 
      real(dp), dimension(:,:), allocatable :: L_ke_J ! L_ke^J
      real(dp), dimension(:,:), allocatable :: L_ic_J ! L_ic^J
      real(dp), dimension(:,:), allocatable :: L_da_J ! L_da^J 
!
      real(dp), dimension(:,:), allocatable :: b_di_cl ! b_cidl 
!
      real(dp), dimension(:,:), allocatable :: X_di_ek ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_kde_i ! Reordered intermediate, term 2
!
      real(dp), dimension(:,:), allocatable :: X_id_kl ! An intermediate, term 1 
      real(dp), dimension(:,:), allocatable :: X_kdl_i ! Reordered intermediate, term 1
!
      real(dp), dimension(:,:), allocatable :: t_cl_ek ! t_kl^ce 
      real(dp), dimension(:,:), allocatable :: t_ce_kl ! t_kl^ce 
!
      real(dp), dimension(:,:), allocatable :: g_ka_de ! g_kade 
      real(dp), dimension(:,:), allocatable :: g_a_kde ! g_kade
      real(dp), dimension(:,:), allocatable :: g_ke_da ! g_keda  
      real(dp), dimension(:,:), allocatable :: g_ic_de ! g_icde
      real(dp), dimension(:,:), allocatable :: g_id_ce ! g_icde
!
      integer(i15) :: l = 0, kde = 0, ka = 0, k = 0, i = 0, el = 0, ek = 0, e = 0
      integer(i15) :: dl = 0, di = 0, de = 0, d = 0, cl = 0, ckel = 0, ci = 0, ck = 0
      integer(i15) :: c = 0, a = 0, kl = 0, ke = 0, id = 0, ic = 0, kdl = 0, da = 0, ce = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0 
      integer(i15) :: batch_dimension = 0, max_batch_length = 0, n_batch = 0
!
      integer(i15) :: a_batch = 0, a_length = 0, a_first = 0, a_last = 0
      integer(i15) :: e_batch = 0, e_length = 0, e_first = 0, e_last = 0, offset_kde = 0
      integer(i15) :: d_batch = 0, d_length = 0, d_first = 0, d_last = 0, offset_id = 0
!
!     :: Term 2. - sum_ckdle b_cidl t_kl^ce g_kade ::
!
!     X_di_ek = sum_cl b_cidl t_kl^ce = sum_cl b_di_cl t_cl_ek
!
!     Reorder to b_di_cl = b_cidl
!
      call allocator(b_di_cl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      b_di_cl = zero 
!
      do l = 1, wf%n_o
         do c = 1, wf%n_v
!
            cl = index_two(c, l, wf%n_v)
!
            do i = 1, wf%n_o
!
               ci = index_two(c, i, wf%n_v)
!
               do d = 1, wf%n_v
!
                  di = index_two(d, i, wf%n_v)
                  dl = index_two(d, l, wf%n_v)
!
                  b_di_cl(di, cl) = b_ai_bj(ci, dl) ! b_cidl
!
               enddo
            enddo
         enddo
      enddo
!
!     Order amplitudes as t_cl_ek = t_kl^ce
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_cl_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_cl_ek = zero 
!
      do k = 1, wf%n_o
         do e = 1, wf%n_v
!
            ek = index_two(e, k, wf%n_v)
!
            do l = 1, wf%n_o
!
               el = index_two(e, l, wf%n_v)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  cl = index_two(c, l, wf%n_v)
!
                  ckel = index_packed(ck, el) 
!
                  t_cl_ek(cl, ek) = wf%t2am(ckel, 1) ! t_kl^ce
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_amplitudes
!
!     Form the intermediate X_di_ek = sum_cl b_di_cl t_cl_ek
!
      call allocator(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  b_di_cl,           &
                  (wf%n_o)*(wf%n_v), &
                  t_cl_ek,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_di_ek,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(b_di_cl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     - sum_ckdle b_cidl t_kl^ce g_kade = sum_kde g_kade X_di_ek
!
!     Reorder X_di_ek to X_kde_i
!
      call allocator(X_kde_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      X_kde_i = zero 
!
      do i = 1, wf%n_o
         do e = 1, wf%n_v
            do d = 1, wf%n_v
!
               di = index_two(d, i, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ek = index_two(e, k, wf%n_v)
!
                  kde = index_three(k, d, e, wf%n_o, wf%n_v)
!
                  X_kde_i(kde, i) = X_di_ek(di, ek)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Get L_ka^J 
!
      call allocator(L_ka_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ka_J)
!
!     Prepare batching over index e
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do e_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(e_first, e_last, e_batch, max_batch_length, batch_dimension)
         e_length = e_last - e_first + 1    
!
!        Form g_ka_de = g_kade 
!
         call allocator(L_de_J, (wf%n_v)*e_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_de_J, 1, wf%n_v, e_first, e_last)
!
         call allocator(g_ka_de, (wf%n_o)*(wf%n_v), (wf%n_v)*e_length)
!
         call dgemm('N','T',            &
                     (wf%n_o)*(wf%n_v), & 
                     (wf%n_v)*e_length, &
                     wf%n_J,            &
                     one,               &
                     L_ka_J,            &
                     (wf%n_o)*(wf%n_v), &
                     L_de_J,            &
                     (wf%n_v)*e_length, &
                     zero,              &
                     g_ka_de,           &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(L_de_J, (wf%n_v)*e_length, wf%n_J)
!
!        Reorder to g_a_kde = g_ka_de = g_kade 
!
         call allocator(g_a_kde, wf%n_v, (wf%n_o)*(wf%n_v)*e_length)
         g_a_kde = zero 
!
         do e = 1, e_length
            do d = 1, wf%n_v
!
               de = index_two(d, e, wf%n_v)
!  
               do k = 1, wf%n_o
!
                  kde = index_three(k, d, e, wf%n_o, wf%n_v)
!
                  do a = 1, wf%n_v
!
                     ka = index_two(k, a, wf%n_o)
!
                     g_a_kde(a, kde) = g_ka_de(ka, de)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ka_de, (wf%n_v)*(wf%n_o), (wf%n_v)*e_length)
!
!        Add the contribution of
!
!           - sum_ckdle b_cidl t_kl^ce g_kade = sum_kde g_a_kde X_kde_i
!
!        arising from the present batch over e 
!
         offset_kde = index_three(1, 1, e_first, wf%n_o, wf%n_v)
!
         call dgemm('N','N',                     &
                     wf%n_v,                     &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)*e_length, &
                     -one,                       &
                     g_a_kde,                    &
                     wf%n_v,                     &
                     X_kde_i(offset_kde,1),      & ! First element to use 
                     (wf%n_o)*(wf%n_v)**2,       & ! Full space dimension of X_kde_i 
                     one,                        &
                     sigma_a_i,                  &
                     wf%n_v)
!
         call deallocator(g_a_kde, wf%n_v, (wf%n_o)*(wf%n_v)*e_length)
!
      enddo ! End of batches over e 
!
      call deallocator(L_ka_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call deallocator(X_kde_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     :: Term 3. - sum_ckdle b_cldi t_kl^ce g_keda ::
!
!     X_di_ek = sum_cl b_cldi t_kl^ce = sum_cl b_di_cl t_cl_ek
!
!     We have t_cl_ek = t_kl^ce, so this can be used unaltered
!     b_ai_bj(cl,di) = b_cldi & hence b_ai^bj^T(di,cl) = b_cldi
!
      call allocator(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('T','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  b_ai_bj,           & ! "b_cl_di"
                  (wf%n_o)*(wf%n_v), &
                  t_cl_ek,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_di_ek,           &
                  (wf%n_o)*(wf%n_v))
!
!     - sum_kde X_di_ek g_keda
!
!     Reorder X_di_ek to X_kde_i
!
      call allocator(X_kde_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      X_kde_i = zero
!
      do i = 1, wf%n_o
         do e = 1, wf%n_v
            do d = 1, wf%n_v
!
               di = index_two(d, i, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ek = index_two(e, k, wf%n_v)
!
                  kde = index_three(k, d, e, wf%n_o, wf%n_v)
!
                  X_kde_i(kde, i) = X_di_ek(di, ek)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Get L_ke^J 
!
      call allocator(L_ke_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ke_J)
!
!     Prepare batching over a 
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1    
!
!        Form g_ke_da = g_keda
!
         call allocator(L_da_J, (wf%n_v)*a_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_da_J, 1, wf%n_v, a_first, a_last)
!
         call allocator(g_ke_da, (wf%n_o)*(wf%n_v), (wf%n_v)*a_length)
!
         call dgemm('N','T',            &
                     (wf%n_o)*(wf%n_v), & 
                     (wf%n_v)*a_length, &
                     wf%n_J,            &
                     one,               &
                     L_ke_J,            &
                     (wf%n_o)*(wf%n_v), &
                     L_da_J,            &
                     (wf%n_v)*a_length, &
                     zero,              &
                     g_ke_da,           &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(L_da_J, (wf%n_v)*a_length, wf%n_J)
!
!        Reorder to g_a_kde = g_ke_da 
!
         call allocator(g_a_kde, a_length, (wf%n_o)*(wf%n_v)**2)
         g_a_kde = zero 
!
         do e = 1, wf%n_v
            do d = 1, wf%n_v
               do k = 1, wf%n_o
!
                  ke = index_two(k, e, wf%n_o)
!
                  kde = index_three(k, d, e, wf%n_o, wf%n_v)
!     
                  do a = 1, a_length
!
                     da = index_two(d, a, wf%n_v)
!
                     g_a_kde(a, kde) = g_ke_da(ke, da) ! g_keda
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ke_da, (wf%n_o)*(wf%n_v), (wf%n_v)*a_length)
!
!        Add 
!
!            - sum_ckdle b_cldi t_kl^ce g_keda
!                 = -sum_kde g_a_kde X_kde_i 
!
!        for the current batch of a's
!
         call dgemm('N','N',               &
                     a_length,             &
                     wf%n_o,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     -one,                 &
                     g_a_kde,              &
                     a_length,             &
                     X_kde_i,              &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     sigma_a_i(a_first,1), &
                     wf%n_v)
!
         call deallocator(g_a_kde, a_length, (wf%n_o)*(wf%n_v)**2)
!
      enddo ! End of batches over a 
!
      call deallocator(X_kde_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call deallocator(L_ke_J, (wf%n_o)*(wf%n_v), wf%n_J)
!  
!     E: There's a symmetry bug in the following term: 
!
!     :: Term 1. - sum_ckdle b_akdl t_kl^ce g_icde :: 
! 
!     X_id_kl = sum_ce t_kl^ce g_icde = sum_ce g_id_ce t_ce_kl
!
!     We have t_cl_ek = t_kl^ce
!     Reorder to t_ce_kl = t_cl_ek = t_kl^ce
!
      call allocator(t_ce_kl, (wf%n_v)**2, (wf%n_o)**2)
      t_ce_kl = zero 
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
!
            kl = index_two(k, l, wf%n_o)
!
            do e = 1, wf%n_v
!
               ek = index_two(e, k, wf%n_v)
!
               do c = 1, wf%n_v
!
                  cl = index_two(c, l, wf%n_v)
                  ce = index_two(c, e, wf%n_v)
!
                  t_ce_kl(ce, kl) = t_cl_ek(cl, ek) ! t_kl^ce
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(t_cl_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Get L_ic^J
!
      call allocator(L_ic_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ic_J)
!
      call allocator(X_id_kl, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      X_id_kl = zero
!
!     Prepare for batching over d 
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index d
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do d_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(d_first, d_last, d_batch, max_batch_length, batch_dimension)
         d_length = d_last - d_first + 1    
!
!        Form g_ic_de = g_icde 
!
         call allocator(L_de_J, (wf%n_v)*d_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_de_J, d_first, d_last, 1, wf%n_v)
!
         call allocator(g_ic_de, (wf%n_o)*(wf%n_v), (wf%n_v)*d_length)
!
         call dgemm('N','T',            &
                     (wf%n_o)*(wf%n_v), & 
                     (wf%n_v)*d_length, &
                     wf%n_J,            &
                     one,               &
                     L_ic_J,            &
                     (wf%n_o)*(wf%n_v), &
                     L_de_J,            &
                     (wf%n_v)*d_length, & 
                     zero,              &
                     g_ic_de,           &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(L_de_J, d_length*(wf%n_v), wf%n_J)
!
!        Reorder to g_id_ce = g_ic_de
!
         call allocator(g_id_ce, (wf%n_o)*d_length, (wf%n_v)**2)
         g_id_ce = zero 
!
         do e = 1, wf%n_v
            do c = 1, wf%n_v
!
               ce = index_two(c, e, wf%n_v)
!
               do d = 1, d_length
!
                  de = index_two(d, e, d_length)
!
                  do i = 1, wf%n_o
!
                     ic = index_two(i, c, wf%n_o)
                     id = index_two(i, d, wf%n_o)
!
                     g_id_ce(id, ce) = g_ic_de(ic, de)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ic_de, (wf%n_o)*(wf%n_v), (wf%n_v)**2)
!
!        Add the contribution 
!
!           X_id_kl = sum_ce t_kl^ce g_icde = sum_ce g_id_ce t_ce_kl
!
!        from the current batch of d 
!
         offset_id = index_two(1, d_first, wf%n_o)
!
         call dgemm('N','N',               &
                     (wf%n_o)*d_length,    &
                     (wf%n_o)**2,          &
                     (wf%n_v)**2,          &
                     one,                  &
                     g_id_ce,              &
                     (wf%n_o)*d_length,    &
                     t_ce_kl,              &
                     (wf%n_v)**2,          &
                     one,                  &
                     X_id_kl(offset_id,1), &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(g_id_ce, (wf%n_o)*d_length, (wf%n_v)**2)
!
      enddo ! End of batches over d 
!
      call deallocator(t_ce_kl, (wf%n_v)**2, (wf%n_o)**2)
      call deallocator(L_ic_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     - sum_ckdle b_akdl t_kl^ce g_icde = sum_kdl b_akdl X_id_kl
!
!     Reorder to X_kdl_i = X_id_kl 
!
      call allocator(X_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      X_kdl_i = zero
!
      do i = 1, wf%n_o
         do l = 1, wf%n_o
            do k = 1, wf%n_o
!
               kl = index_two(k, l, wf%n_o)
!
               do d = 1, wf%n_v
!
                  id = index_two(i, d, wf%n_o)
!
                  kdl = index_three(k, d, l, wf%n_o, wf%n_v)
!
                  X_kdl_i(kdl, i) = X_id_kl(id, kl)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(X_id_kl, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
!     Add - sum_ckdle b_akdl t_kl^ce g_icde = - sum_dkl b_a_kdl X_kdl_i
!
!     Note: we interpret b_ai_bj as b_a_ibj
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_ai_bj,              & ! "b_a_ibj"
                  wf%n_v,               &
                  X_kdl_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call deallocator(X_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_g1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a2_ccsd(wf, sigma_ai_bj, b_a_i)
!!
!!    Jacobian transpose CCSD A2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A2 term,
!!
!!       2 F_jb b_ai - F_ib b_aj - sum_k L_ikjb b_ak + sum_c L_cajb b_ci 
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: b_a_i  
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: L_ik_J
      real(dp), dimension(:,:), allocatable :: L_jb_J
      real(dp), dimension(:,:), allocatable :: L_ca_J
      real(dp), dimension(:,:), allocatable :: L_cb_J
!
      real(dp), dimension(:,:), allocatable :: g_ik_jb ! g_ikjb 
      real(dp), dimension(:,:), allocatable :: L_k_ibj ! L_ikjb
!
      real(dp), dimension(:,:), allocatable :: g_ca_jb ! g_cajb 
      real(dp), dimension(:,:), allocatable :: g_cb_ja ! g_cbja 
!
      real(dp), dimension(:,:), allocatable :: sigma_i_ajb ! sigma_ai_bj contribution 
      real(dp), dimension(:,:), allocatable :: sigma_i_bja ! sigma_ai_bj contribution 
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, ai = 0, bj = 0, k = 0, jk = 0, jb = 0
      integer(i15) :: ik = 0, ibj = 0, ib = 0, bja = 0, ajb = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0 
      integer(i15) :: batch_dimension = 0, max_batch_length = 0, n_batch = 0
!
      integer(i15) :: a_batch = 0, a_length = 0, a_first = 0, a_last = 0
      integer(i15) :: b_batch = 0, b_length = 0, b_first = 0, b_last = 0
!
!     :: Term 1 & 2. 2 F_jb b_ai - F_ib b_aj :: 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!  
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) &
                                      + two*(wf%fock_ia(j, b))*b_a_i(a, i) &
                                      - (wf%fock_ia(i, b))*b_a_i(a, j)
!
               enddo
            enddo
         enddo
      enddo
!
!     :: Term 3. - sum_k L_ikjb b_ak ::
!
!     Form g_ik_jb = g_ikjb
!
      call allocator(L_ik_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ij(L_ik_J)
!
      call allocator(L_jb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_jb_J)
!
      call allocator(g_ik_jb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ik_J,            &
                  (wf%n_o)**2,       &
                  L_jb_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ik_jb,           &
                  (wf%n_o)**2)
!
!     Deallocate Cholesky (keep L_jb_J for the next term)
!
      call deallocator(L_ik_J, (wf%n_o)**2, wf%n_J)
!
!     Form L_k_ibj = L_ikjb = 2 * g_ikjb - g_ibjk
!                           = 2 * g_ikjb - g_jkib
!                           = 2 * g_ik_jb(ik,jb) - g_ik_jb(jk,ib)
!
      call allocator(L_k_ibj, wf%n_o, (wf%n_v)*(wf%n_o)**2)
      L_k_ibj = zero 
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            jb = index_two(j, b, wf%n_o)
!
            do i = 1, wf%n_o
!
               ib = index_two(i, b, wf%n_o)
!
               ibj = index_three(i, b, j, wf%n_o, wf%n_v)
!
               do k = 1, wf%n_o
!
                  jk = index_two(j, k, wf%n_o)
                  ik = index_two(i, k, wf%n_o)
!
                  L_k_ibj(k, ibj) = two*g_ik_jb(ik, jb) - g_ik_jb(jk, ib) ! L_ikjb
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ik_jb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Add - sum_k L_ikjb b_ak = - sum_k b_ak L_k_ibj
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  b_a_i,                & ! "b_a_k"
                  wf%n_v,               &
                  L_k_ibj,              &
                  wf%n_o,               &
                  one,                  &
                  sigma_ai_bj,          & ! "sigma_a_ibj"
                  wf%n_v)
!
      call deallocator(L_k_ibj, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 4. 2 sum_c g_cajb b_ci - sum_c g_cbja b_ci :: 
!
!     2 sum_c g_cajb
!
!     Prepare for batching over a 
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
!        Form g_ca_jb
!
         call allocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_ca_J, 1, wf%n_v, a_first, a_last)
!
         call allocator(g_ca_jb, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
!
         call dgemm('N','T',            &
                     (wf%n_v)*a_length, & 
                     (wf%n_o)*(wf%n_v), &
                     wf%n_J,            &
                     one,               &
                     L_ca_J,            &
                     (wf%n_v)*a_length, &
                     L_jb_J,            &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     g_ca_jb,           &
                     (wf%n_v)*a_length)
!
         call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
!        Add 2 sum_c g_cajb b_ci = 2 sum_c b_c_i^T(i,c) g_c_ajb  
!
         call allocator(sigma_i_ajb, wf%n_o, (wf%n_o)*(wf%n_v)*a_length)
!
         call dgemm('T','N',                     &
                     wf%n_o,                     & 
                     (wf%n_o)*(wf%n_v)*a_length, &
                     wf%n_v,                     &
                     two,                        &
                     b_a_i,                      & ! "b_c_i"
                     wf%n_v,                     &
                     g_ca_jb,                    & ! "g_c_ajb"
                     wf%n_v,                     &
                     zero,                       &
                     sigma_i_ajb,                &
                     wf%n_o)
!
         call deallocator(g_ca_jb, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
!
         do i = 1, wf%n_o
            do a = 1, a_length
!
               Ai = index_two(a + a_first - 1, i, wf%n_v)
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     ajb = index_three(a, j, b, a_length, wf%n_o)
!
                     sigma_ai_bj(Ai, bj) = sigma_ai_bj(Ai, bj) + sigma_i_ajb(i, ajb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(sigma_i_ajb, wf%n_o, (wf%n_o)*(wf%n_v)*a_length)
!
      enddo ! End of batches over a
!
!     - sum_c g_cbja b_ci
!
!     Prepare for batching over b 
!
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index b
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of b batches 
!
      do b_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(b_first, b_last, b_batch, max_batch_length, batch_dimension)
         b_length = b_last - b_first + 1 
!
!        Form g_cb_ja = g_cbja 
! 
         call allocator(L_cb_J, (wf%n_v)*b_length, wf%n_J)
!
         call wf%get_cholesky_ab(L_cb_J, 1, wf%n_v, b_first, b_last)
!
         call allocator(g_cb_ja, (wf%n_v)*b_length, (wf%n_o)*(wf%n_v))
!
         call dgemm('N','T',            &
                     (wf%n_v)*b_length, & 
                     (wf%n_o)*(wf%n_v), &
                     wf%n_J,            &
                     one,               &
                     L_cb_J,            &
                     (wf%n_v)*b_length, & 
                     L_jb_J,            & ! "L_ja_J" full space
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     g_cb_ja,           &
                     (wf%n_v)*b_length)
!
!         Form - sum_c g_cbja b_ci = - sum_c b_ci^T(i,c) g_c_bja 
!
         call allocator(sigma_i_bja, wf%n_o, (wf%n_o)*(wf%n_v)*b_length)
!
         call dgemm('T','N',                     &
                     wf%n_o,                     & 
                     (wf%n_o)*(wf%n_v)*b_length, &
                     wf%n_v,                     &
                     -one,                       &
                     b_a_i,                      & ! "b_c_i"
                     wf%n_v,                     &
                     g_cb_ja,                    & ! "g_c_bja"
                     wf%n_v,                     &
                     zero,                       &
                     sigma_i_bja,                &
                     wf%n_o)
!
         call deallocator(g_cb_ja, (wf%n_v)*b_length, (wf%n_o)*(wf%n_v))
!
!        Add it to sigma_ai_bj 
!
         do j = 1, wf%n_o
            do b = 1, b_length
!
               Bj = index_two(b + b_first - 1, j, wf%n_v)
!
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     ai = index_two(a, i, wf%n_v)
!
                     bja = index_three(b, j, a, b_length, wf%n_o)
!
                     sigma_ai_bj(ai, Bj) = sigma_ai_bj(ai, Bj) &
                                         + sigma_i_bja(i, bja)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(sigma_i_bja, wf%n_o, b_length*(wf%n_o)*(wf%n_v))
!
      enddo ! End of batches over b
!
      call deallocator(L_jb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
   end subroutine jacobian_transpose_ccsd_a2_ccsd
!
!
end submodule jacobian_transpose