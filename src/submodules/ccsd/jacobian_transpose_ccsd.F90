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
      call add_1243_to_1234(two, g_kc_id, L_kcd_i, wf%n_o, wf%n_v, wf%n_v, wf% n_o)
      call add_1342_to_1234(-one, g_kc_id, L_kcd_i, wf%n_o, wf%n_v, wf%n_v, wf% n_o)
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
      call squareup_and_sort_1234_to_4213(wf%t2, t_l_kcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
!
      call sort_1234_to_3124(L_kcd_i, L_a_ldk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_kcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Form t_ldk_c = t_kl^cd = t_l_kcd(l, kcd)
!
      call mem%alloc(t_ldk_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_ldk_c = zero
!
      call sort_1234_to_1423(t_l_kcd, t_ldk_c, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
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
                     g_dl_ca,                    & ! g_dlc_a
                     (wf%n_o)*(wf%n_v)**2,       &
                     b_ai_bj,                    & ! b_dlc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     sigma_a_i(batch_a%first,1), &
                     wf%n_v)
!
         call mem%dealloc(g_dl_ca, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_a%length))
!
      enddo ! End of batches over a 
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
!     :: Term 1. - sum_ckdl b_ckal F_id t_kl^cd ::
!
!     Read amplitudes and order as t_lck_d = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      call squareup_and_sort_1234_to_2341(wf%t2, t_lck_d, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
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
                  (wf%n_v),             &
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
                  t_lck_d,              & ! t_l_ckd
                  (wf%n_o),             &
                  b_ai_bj,              & ! b_ckd_i
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
!
!
   module subroutine jacobian_transpose_ccsd_e1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD E1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
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
      real(dp), dimension(:,:), allocatable :: g_da_le ! g_dale
      real(dp), dimension(:,:), allocatable :: L_a_eld ! L_dale
!
      real(dp), dimension(:,:), allocatable :: g_de_ia ! g_deia 
      real(dp), dimension(:,:), allocatable :: L_ai_ed ! L_deia 
!
      real(dp), dimension(:,:), allocatable :: X_el_di ! An intermediate, term 1
!
      real(dp), dimension(:,:), allocatable :: X_d_e ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_e_d ! Reordered intermediate, term 2
!
!     Batching variables 
!
      integer(i15) :: required = 0
!
      integer(i15) :: current_d_batch
      integer(i15) :: current_a_batch 
!  
      type(batching_index) :: batch_d 
      type(batching_index) :: batch_a
!
!     :: Term 3. - sum_ckdlm b_ckal L_ilmd t_km^cd ::
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_dm_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_dm_ck, (wf%n_o)*(wf%n_v)) 
!
      !call wf%destruct_double_amplitudes
!
      call mem%alloc(g_il_md, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_ooov(g_il_md)
!
!     Form L_il_dm = L_ilmd = 2 * g_ilmd - g_mlid
!                           = 2 * g_il_md(il,md) - g_il_md(ml,id)
!
      call mem%alloc(L_il_dm, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      L_il_dm = zero
!
      call daxpy((wf%n_o)**3*(wf%n_v), two, g_il_md, 1, L_il_dm,1)
      call add_3214_to_1234(-one, g_il_md, L_il_dm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_il_md, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_il_ck = sum_md L_ilmd t_mk^dc 
!                                   = sum_md L_il_dm t_dm_ck
!
      call mem%alloc(X_il_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
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
      call mem%dealloc(L_il_dm, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Add - sum_ckdlm b_ckal L_ilmd t_km^cd
!         = - sum_ckl b_a_lck X_i_lck^T
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_ai_bj,              & ! b_a_lck (= b_al_ck = b_ai_bj)
                  wf%n_v,               &
                  X_il_ck,              & ! X_i_lck
                  wf%n_o,               &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call mem%dealloc(X_il_ck, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     :: Term 4. - sum_ckdlm b_ckdl L_mlia t_km^cd ::
!
!     Form the intermediate X_m_l = sum_ckd t_km^cd b_ckdl
!                                 = sum_ckd t_ckd_m^T b_ckd_l
!
      call mem%alloc(X_m_l, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_dm_ck,              & ! t_ckd_m
                  (wf%n_o)*(wf%n_v)**2, &
                  b_ai_bj,              & ! b_aib_j
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_m_l,                &
                  wf%n_o)
!
      call mem%dealloc(t_dm_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form g_ml_ia = g_mlia 
!
      call mem%alloc(g_ml_ia, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_ooov(g_ml_ia)
!
!     Form L_ai_ml = L_mlia = 2 * g_mlia - g_mail
!                           = 2 * g_mlia - g_ilma 
!
      call mem%alloc(L_ai_ml, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      L_ai_ml = zero
!
      call daxpy((wf%n_o)**3*(wf%n_v), two, g_ml_ia, 1, L_ai_ml,1)
      call add_3214_to_1234(-one, g_ml_ia, L_ai_ml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ml_ia, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
                  X_m_l,             & ! X_ml
                  (wf%n_o)**2,       &
                  one,               &
                  sigma_a_i,         & ! sigma_ai
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ai_ml, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      call mem%dealloc(X_m_l, wf%n_o, wf%n_o)
!
!     :: Term 1. sum_ckdle b_ckdi L_dale t_kl^ce :: 
!
!     Read amplitudes from disk 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_el_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_el_ck, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_el_di = sum_ck t_lk^ec b_ckdi
!                                   = sum_ck t_el_ck b_ck_di  
!
      call mem%alloc(X_el_di, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_el_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  b_ai_bj,           & ! b_ck_di
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_el_di,           &
                  (wf%n_o)*(wf%n_v))
!
!     sum_dle L_dale X_el_di 
!
!     Prepare batching over index a 
!
      required =  wf%integrals%get_required_vvov() + (wf%n_v**3)*(wf%n_o)
!     
      call batch_a%init(wf%n_v)
      call mem%num_batch(batch_a, required)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_da_le, (wf%n_v)*(batch_a%length), (wf%n_o)*(wf%n_v))
!
         call wf%get_vvov(g_da_le,        &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Form  L_a_eld = L_dale = 2 * g_dale - g_dela 
!                               = 2 * g_da_le(da, le) - g_da_le(de, la)      
!
         call mem%alloc(L_a_eld, batch_a%length, (wf%n_o)*(wf%n_v)**2)
         L_a_eld = zero
!
         call add_4132_to_1234(two, g_da_le, L_a_eld, batch_a%length, (wf%n_v), (wf%n_o), (wf%n_v))
         call add_4231_to_1234(-one, g_da_le, L_a_eld, batch_a%length, (wf%n_v), (wf%n_o), (wf%n_v))
!
         call mem%dealloc(g_da_le, (wf%n_v)*(batch_a%length), (wf%n_o)*(wf%n_v))
!
!        Add sum_ckdle b_ckdi L_dale t_kl^ce
!            = sum_eld L_a_eld X_el_di
!
         call dgemm('N','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     L_a_eld,                    &
                     batch_a%length,             &
                     X_el_di,                    & ! X_eld_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     sigma_a_i(batch_a%first,1), &
                     wf%n_v)
!
         call mem%dealloc(L_a_eld, batch_a%length, (wf%n_o)*(wf%n_v)**2)
!
      enddo ! End of batches over a 
!
!     :: Term 2. sum_ckdle b_ckdl L_deia t_kl^ce ::
!
!     Form the intermediate X_d_e = sum_ckl b_ckdl t_kl^ce = sum_ckl b_d_lck t_e_lck^T
!
      call mem%alloc(X_d_e, wf%n_v, wf%n_v)
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_ai_bj,              & ! b_d_lck = b_dlck = b_ckdl
                  wf%n_v,               &
                  t_el_ck,              & ! t_e_lck
                  wf%n_v,               &
                  zero,                 &
                  X_d_e,                &
                  wf%n_v)
!
      call mem%dealloc(t_el_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call mem%alloc(X_e_d, wf%n_v, wf%n_v)
!
      call sort_12_to_21(X_d_e, X_e_d, wf%n_v, wf%n_v)
!
      call mem%dealloc(X_d_e, wf%n_v, wf%n_v)
!
!     sum_ckdle b_ckdl L_deia t_kl^ce = sum_de L_deia X_e_d
!
!     Prepare batching over index d 
!
      required = wf%integrals%get_required_vvov() + (wf%n_v**3)*(wf%n_o)
!
      call batch_d%init(wf%n_v)
      call mem%num_batch(batch_d, required)
!
      do current_d_batch = 1, batch_d%num_batches
!
!        For each batch, get the limits for the d index 
!
         call batch_d%determine_limits(current_d_batch)
!
!        Form g_de_ia = g_deia 
!
         call mem%alloc(g_de_ia, (wf%n_v)*(batch_d%length), (wf%n_v)*(wf%n_o))
!
         call wf%get_vvov(g_de_ia,        &
                           batch_d%first, &
                           batch_d%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Form L_ai_ed = L_deia = 2 * g_deia - g_daie
!                              = 2 * g_de_ia(de,ia) - g_de_ia(da,ie)
!
         call mem%alloc(L_ai_ed, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_d%length))
         L_ai_ed = zero
!
         call add_4321_to_1234(two, g_de_ia, L_ai_ed, (wf%n_v), (wf%n_o), (wf%n_v), (batch_d%length))
         call add_4123_to_1234(two, g_de_ia, L_ai_ed, (wf%n_v), (wf%n_o), (wf%n_v), (batch_d%length))
!
         call mem%dealloc(g_de_ia, (wf%n_v)*(batch_d%length), (wf%n_o)*(wf%n_v))
!
         call dgemm('N','N',                    &
                     (wf%n_v)*(wf%n_o),         &
                     1,                         &
                     (wf%n_v)*(batch_d%length), &
                     one,                       &
                     L_ai_ed,                   &
                     (wf%n_v)*(wf%n_o),         &
                     X_e_d(1,batch_d%first),    & 
                     (wf%n_v)*(batch_d%length), & 
                     one,                       &
                     sigma_a_i,                 &
                     (wf%n_v)*(wf%n_o))
!
         call mem%dealloc(L_ai_ed, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_d%length))
!
      enddo ! End of batches over d
!
      call mem%dealloc(X_e_d, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD F1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
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
      real(dp), dimension(:,:), allocatable :: g_ik_mc ! g_ikmc 
      real(dp), dimension(:,:), allocatable :: g_kdm_i ! g_mkid
      real(dp), dimension(:,:), allocatable :: g_a_mkl ! g_mkla
!
      real(dp), dimension(:,:), allocatable :: X_ik_dl ! An intermediate, term 1 
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
      !call wf%read_double_amplitudes
! 
      call mem%alloc(t_mc_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call squareup_and_sort_1234_to_4132(wf%t2, t_mc_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      !call wf%destruct_double_amplitudes
!
!     Form the integral g_ik_mc 
!
      call mem%alloc(g_ik_mc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_ooov(g_ik_mc)
!
!     Form the intermediate X_ik_dl = sum_mc t_lm^cd g_ikmc = sum_mc g_ik_mc t_mc_dl
!
      call mem%alloc(X_ik_dl, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
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
!     Add sum_ckdlm b_akdl t_lm^cd g_ikmc
!         = sum_kdl b_a_kdl X_i_kdl^T 
!
!     Note: we interpret b_ai_bj as b_a_ibj
!           we interpret X_ik_dl as X_i_kdl
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_ai_bj,              & ! "b_a_ibj"
                  wf%n_v,               &
                  X_ik_dl,              & ! "X_i_kdl"
                  wf%n_o,               &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
       call mem%dealloc(X_ik_dl, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     :: Term 2. sum_ckdlm b_ckal t_ml^cd g_mkid ::
!
!     X_ak_dm = sum_cl b_ckal t_ml^cd
!             = sum_cl b_ak_cl t_cl_dm
!
!     We have t_mc_dl(mc,dl) = t_lm^cd 
!     Reorder t_cl_dm(cl,dm) = t_mc_dl(lc,dm) = t_ml^cd  
!
      call mem%alloc(t_cl_dm, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call sort_1234_to_2134(t_mc_dl, t_cl_dm, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_mc_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder to b_ak_cl = b_ckal 
!
      call mem%alloc(b_ak_cl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) 
!
      call sort_1234_to_3214(b_ai_bj, b_ak_cl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ak_dm = sum_cl b_ak_cl t_cl_dm
!
      call mem%alloc(X_ak_dm, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
      call mem%dealloc(b_ak_cl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     sum_ckdlm b_ckal t_ml^cd g_mkid = sum_kdm X_ak_dm g_mkid
!
!     We have g_ik_mc(ik,mc) = g_ikmc 
!     Reorder to g_kdm_i(kdm,i) = g_mkid = g_ik_mc(mk, id)
!
      call mem%alloc(g_kdm_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
      call sort_1234_to_2413(g_ik_mc, g_kdm_i, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ik_mc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
      call mem%dealloc(X_ak_dm, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 3. sum_ckdlm b_ckdi t_ml^cd g_mkla ::
!
!     X_ki_ml = sum_cd b_ckdi t_ml^cd 
!
!     We have t_cl_dm(cl,dm) = t_ml^cd
!     Reorder to t_cd_ml(cd,ml) = t_cl_dm(cl,dm) = t_ml^cd 
!
      call mem%alloc(t_cd_ml, (wf%n_v)**2, (wf%n_o)**2)
!
      call sort_1234_to_1342(t_cl_dm, t_cd_ml, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_cl_dm, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder to b_ki_cd = b_ckdi 
!
      call mem%alloc(b_ki_cd, (wf%n_o)**2, (wf%n_v)**2)
!
      call sort_1234_to_2413(b_ai_bj, b_ki_cd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form intermediate X_ki_ml = sum_cd b_ckdi t_ml^cd = sum_cd b_ki_cd t_cd_ml
!
      call mem%alloc(X_ki_ml, (wf%n_o)**2, (wf%n_o)**2)
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
      call mem%dealloc(t_cd_ml, (wf%n_v)**2, (wf%n_o)**2)
      call mem%dealloc(b_ki_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!     sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_mkla X_ki_ml 
!
!     Reorder to X_mkl_i
!
      call mem%alloc(X_mkl_i, (wf%n_o)**3, wf%n_o)
      X_mkl_i = zero 
!
      call sort_1234_to_3142(X_ki_ml, X_mkl,i, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_ki_ml, (wf%n_o)**2, (wf%n_o)**2)
!
!     We have g_kdm_i(kdm,i) = g_mkid
!     Reorder to g_a_mkl(a,mkl) = g_mkla = g_kdm_i(kam,l)
!
      call mem%alloc(g_a_mkl, wf%n_v, (wf%n_o)**3)
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
      call mem%dealloc(g_kdm_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
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
      call mem%dealloc(g_a_mkl, wf%n_v, (wf%n_o)**3)
      call mem%dealloc(X_mkl_i, (wf%n_o)**3, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD G1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
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
      real(dp), dimension(:,:), allocatable :: b_di_cl ! b_cidl 
!
      real(dp), dimension(:,:), allocatable :: X_di_ek ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_kde_i ! Reordered intermediate, term 2
!
      real(dp), dimension(:,:), allocatable :: X_ked_i ! Reordered intermediate, term 3
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
!     Batching variables 
!
      integer(i15) :: required = 0 
!
      integer(i15) :: current_a_batch = 0
      integer(i15) :: current_d_batch = 0
      integer(i15) :: current_e_batch = 0
!
      integer(i15) :: offset_id = 0
      integer(i15) :: offset_kde = 0
!
      type(batching_index) :: batch_a 
      type(batching_index) :: batch_d
      type(batching_index) :: batch_e 
!
!     :: Term 2. - sum_ckdle b_cidl t_kl^ce g_kade ::
!
!     X_di_ek = sum_cl b_cidl t_kl^ce = sum_cl b_di_cl t_cl_ek
!
!     Reorder to b_di_cl = b_cidl
!
      call mem%alloc(b_di_cl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_3214(b_ai_bj, b_di_cl, wf%n_v, wf%n_o, wf%n_v, wf%n_o) ! b_ai_bj = b_ci_dl
!
!     Order amplitudes as t_cl_ek = t_kl^ce
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_cl_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup_and_sort_1234_to_1432(wf%t2, t_cl_ek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_di_ek = sum_cl b_di_cl t_cl_ek
!
      call mem%alloc(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call mem%dealloc(b_di_cl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     - sum_ckdle b_cidl t_kl^ce g_kade = sum_kde g_kade X_di_ek
!
!     Reorder X_di_ek to X_kde_i
!
      call mem%alloc(X_kde_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call sort_1234_to_4132(X_di_ek, X_kde_i, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Prepare batching over index e
!
      required = wf%integrals%get_required_vvov() + (wf%n_v**3)*(wf%n_o)
!     
!     Initialize batching variable 
!
      call batch_e%init(wf%n_v)
      call mem%num_batch(batch_e, required)          
!
!     Loop over the e-batches
!
      do current_e_batch = 1, batch_e%num_batches
!
!        For each batch, get the limits for the e index 
!
         call batch_e%determine_limits(current_e_batch)
!
!        Form g_ka_de = g_kade 
!
         call mem%alloc(g_ka_de, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_e%length))
!
         call wf%get_ovvv(g_ka_de,       &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_v,        &
                           batch_e%first, &
                           batch_e%last)
!
!        Reorder to g_a_kde = g_ka_de = g_kade 
!
         call mem%alloc(g_a_kde, wf%n_v, (wf%n_o)*(wf%n_v)*(batch_e%length))
         call sort_1234_to_2134(g_ka_de, g_a_kde, wf%n_o, wf%n_v, wf%n_v, batch_e%length)
!
         call mem%dealloc(g_ka_de, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_e%length))
!
         offset_kde = index_three(1, 1, batch_e%first, wf%n_o, wf%n_v)
!
         call dgemm('N','N',                             &
                     wf%n_v,                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_e%length), &
                     -one,                               &
                     g_a_kde,                            &
                     wf%n_v,                             &
                     X_kde_i(offset_kde,1),              & 
                     (wf%n_o)*(wf%n_v)**2,               & 
                     one,                                &
                     sigma_a_i,                          &
                     wf%n_v)
!
         call mem%dealloc(g_a_kde, wf%n_v, (wf%n_o)*(wf%n_v)*(batch_e%length))
!
      enddo ! End of batches over e 
!
      call mem%dealloc(X_kde_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     :: Term 3. - sum_ckdle b_cldi t_kl^ce g_keda ::
!
!     X_di_ek = sum_cl b_cldi t_kl^ce = sum_cl b_di_cl t_cl_ek
!
      call mem%alloc(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('T','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  b_ai_bj,           & ! b_cl_di
                  (wf%n_o)*(wf%n_v), &
                  t_cl_ek,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_di_ek,           &
                  (wf%n_o)*(wf%n_v))
!
!     - sum_kde X_di_ek g_keda
!
!     Reorder X_di_ek to X_ked_i
!
      call mem%alloc(X_ked_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call sort_1234_to_4312(X_di_ek, X_ked_i, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_di_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Prepare batching over a 
!
      required = wf%integrals%get_required_vvov()
!     
!     Initialize batching variable 
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
!        Form g_ke_da = g_keda
!
         call mem%alloc(g_ke_da, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_a%length))
!
         call wf%get_ovvv(g_ke_da,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last)
!
         call dgemm('T','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     -one,                       &
                     g_ke_da,                    & ! g_ked_a
                     (wf%n_o)*(wf%n_v)**2,       &
                     X_ked_i,                    &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     sigma_a_i(batch_a%first,1), &
                     wf%n_v)
!
         call mem%dealloc(g_ke_da, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_a%length))
!
      enddo ! End of batches over a 
!
      call mem%dealloc(X_ked_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     :: Term 1. - sum_ckdle b_akdl t_kl^ce g_icde :: 
! 
!     X_id_kl = sum_ce t_kl^ce g_icde = sum_ce g_id_ce t_ce_kl
!
!     Reorder to t_ce_kl = t_cl_ek = t_kl^ce
!
      call mem%alloc(t_ce_kl, (wf%n_v)**2, (wf%n_o)**2)
      call sort_1234_to_1342(t_cl_ek, t_ce_kl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_cl_ek, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call mem%alloc(X_id_kl, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      X_id_kl = zero
!
!     Prepare for batching over d 
!
      required = wf%integrals%get_required_vvov() + (wf%n_v**3)*(wf%n_o)
!
      call batch_d%init(wf%n_v)
      call mem%num_batch(batch_d, required)         
!
!     Loop over the d-batches
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)   
!
         call mem%alloc(g_ic_de, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_d%length))
!
         call wf%get_ovvv(g_ic_de,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           batch_d%first, &
                           batch_d%last,  &
                           1,             &
                           wf%n_v)
!
!        Reorder to g_id_ce = g_ic_de
!
         call mem%alloc(g_id_ce, (wf%n_o)*(batch_d%length), (wf%n_v)**2)
         call soret_1234_to_1324(g_ic_de, g_id_ce, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call mem%dealloc(g_ic_de, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_d%length))
!
         offset_id = index_two(1, batch_d%first, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(batch_d%length), &
                     (wf%n_o)**2,               &
                     (wf%n_v)**2,               &
                     one,                       &
                     g_id_ce,                   &
                     (wf%n_o)*(batch_d%length), &
                     t_ce_kl,                   &
                     (wf%n_v)**2,               &
                     one,                       &
                     X_id_kl(offset_id,1),      &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(g_id_ce, (wf%n_o)*(batch_d%length), (wf%n_v)**2)
!
      enddo ! End of batches over d 
!
      call mem%dealloc(t_ce_kl, (wf%n_v)**2, (wf%n_o)**2)
!
!     - sum_ckdle b_akdl t_kl^ce g_icde = sum_kdl b_akdl X_id_kl
!
      call mem%alloc(X_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
      call sort_1234_to_3241(X_id_kl, X_kdl_i, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_id_kl, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
!     Add - sum_ckdle b_akdl t_kl^ce g_icde = - sum_dkl b_a_kdl X_kdl_i
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_ai_bj,              & ! b_a_ibj
                  wf%n_v,               &
                  X_kdl_i,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_a_i,            &
                  wf%n_v)
!
      call mem%dealloc(X_kdl_i, (wf%n_v)*(wf%n_o)**2, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_g1_ccsd
!
!
end submodule jacobian_transpose_ccsd