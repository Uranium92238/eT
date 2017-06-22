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
!
!     End of routine ... deallocs:
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
!     :: Term 3. - sum_ckdlm b_ckal L_ilmd t_km^cd ::
!
!     X_il_ck = sum_md L_ilmd t_km^cd = sum_md L_ilmd t_mk^dc 
!
!     Read the amplitudes from disk 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_dm_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call squareup(wf%t2am, t_dm_ck, (wf%n_o)*(wf%n_v)) 
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
      call dgemm('N','T',&
                  (wf%n_o)**2,&
                  (wf%n_o)*(wf%n_v),&
                  wf%n_J,&
                  one,&
                  )
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
!
!
!
!
!
end submodule jacobian_transpose