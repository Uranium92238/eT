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
      call mem%alloc(g_ik_dl, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_oovo(g_ik_dl)
!
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_ai_bj,              & ! b_a_kdl
                  wf%n_v,               &
                  g_ik_dl,              & ! g_i_kdl
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
         call add_4123_to_1234(-one, g_de_ia, L_ai_ed, (wf%n_v), (wf%n_o), (wf%n_v), (batch_d%length))
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
                  b_ai_bj,              & ! b_a_ibj
                  wf%n_v,               &
                  X_ik_dl,              & ! X_i_kdl
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
!
      call sort_1234_to_3142(X_ki_ml, X_mkl_i, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_ki_ml, (wf%n_o)**2, (wf%n_o)**2)
!
!     We have g_kdm_i(kdm,i) = g_mkid
!     Reorder to g_a_mkl(a,mkl) = g_mkla = g_kdm_i(kam,l)
!
      call mem%alloc(g_a_mkl, wf%n_v, (wf%n_o)**3)
!
      call sort_1234_to_2314(g_kdm_i, g_a_mkl, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
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
!
         call sort_1234_to_1324(g_ic_de, g_id_ce, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
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
   module subroutine jacobian_transpose_ccsd_a2_ccsd(wf, sigma_ai_bj, b_a_i)
!!
!!    Jacobian transpose CCSD A2 
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
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
      integer(i15) :: required = 0 
!
      integer(i15) :: current_a_batch = 0
      integer(i15) :: current_b_batch = 0
!
      type(batching_index) :: batch_a 
      type(batching_index) :: batch_b 
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
!
      call mem%alloc(g_ik_jb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_ooov(g_ik_jb)
!
!     Form L_k_ibj = L_ikjb = 2 * g_ikjb - g_ibjk
!                           = 2 * g_ikjb - g_jkib
!                           = 2 * g_ik_jb(ik,jb) - g_ik_jb(jk,ib)
!
      call mem%alloc(L_k_ibj, wf%n_o, (wf%n_v)*(wf%n_o)**2)
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
      call mem%dealloc(g_ik_jb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
      call mem%dealloc(L_k_ibj, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 4. 2 sum_c g_cajb b_ci - sum_c g_cbja b_ci :: 
!
!     2 sum_c g_cajb
!
!     Prepare for batching over a 
!
      required = wf%integrals%get_required_vvov() + (wf%n_o**2)*(wf%n_v**2)
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
!        Form g_ca_jb
!
         call mem%alloc(g_ca_jb, (wf%n_v)*(batch_a%length), (wf%n_o)*(wf%n_v))
!
         call wf%get_vvov(g_ca_jb,        &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Add 2 sum_c g_cajb b_ci = 2 sum_c b_c_i^T(i,c) g_c_ajb  
!
         call mem%alloc(sigma_i_ajb, wf%n_o, (wf%n_o)*(wf%n_v)*(batch_a%length))
!
         call dgemm('T','N',                             &
                     wf%n_o,                             & 
                     (wf%n_o)*(wf%n_v)*(batch_a%length), &
                     wf%n_v,                             &
                     two,                                &
                     b_a_i,                              & ! "b_c_i"
                     wf%n_v,                             &
                     g_ca_jb,                            & ! "g_c_ajb"
                     wf%n_v,                             &
                     zero,                               &
                     sigma_i_ajb,                        &
                     wf%n_o)
!
         call mem%dealloc(g_ca_jb, (wf%n_v)*(batch_a%length), (wf%n_o)*(wf%n_v))
!
         do i = 1, wf%n_o
            do a = 1, batch_a%length
!
               Ai = index_two(a + batch_a%first - 1, i, wf%n_v)
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     ajb = index_three(a, j, b, batch_a%length, wf%n_o)
!
                     sigma_ai_bj(Ai, bj) = sigma_ai_bj(Ai, bj) + sigma_i_ajb(i, ajb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(sigma_i_ajb, wf%n_o, (wf%n_o)*(wf%n_v)*(batch_a%length))
!
      enddo ! End of batches over a
!
!     - sum_c g_cbja b_ci
!
!     Prepare for batching over b 
!
      required = wf%integrals%get_required_vvov() + (wf%n_v**2)*(wf%n_o**2)
!     
!     Initialize batching variable          
!
      call batch_b%init(wf%n_v)
      call mem%num_batch(batch_b, required)
!
!     Loop over the number of b batches 
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index 
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cb_ja = g_cbja 
!
         call mem%alloc(g_cb_ja, (wf%n_v)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
         call wf%get_vvov(g_cb_ja,        &
                           1,             &
                           wf%n_v,        &
                           batch_b%first, &
                           batch_b%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!         Form - sum_c g_cbja b_ci = - sum_c b_ci^T(i,c) g_c_bja 
!
         call mem%alloc(sigma_i_bja, wf%n_o, (wf%n_o)*(wf%n_v)*(batch_b%length))
!
         call dgemm('T','N',                             &
                     wf%n_o,                             & 
                     (wf%n_o)*(wf%n_v)*(batch_b%length), &
                     wf%n_v,                             &
                     -one,                               &
                     b_a_i,                              & ! "b_c_i"
                     wf%n_v,                             &
                     g_cb_ja,                            & ! "g_c_bja"
                     wf%n_v,                             &
                     zero,                               &
                     sigma_i_bja,                        &
                     wf%n_o)
!
         call mem%dealloc(g_cb_ja, (wf%n_v)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
!        Add it to sigma_ai_bj 
!
         do j = 1, wf%n_o
            do b = 1, batch_b%length
!
               Bj = index_two(b + batch_b%first - 1, j, wf%n_v)
!
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     ai = index_two(a, i, wf%n_v)
!
                     bja = index_three(b, j, a, batch_b%length, wf%n_o)
!
                     sigma_ai_bj(ai, Bj) = sigma_ai_bj(ai, Bj) &
                                         + sigma_i_bja(i, bja)
!
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(sigma_i_bja, wf%n_o, (batch_b%length)*(wf%n_o)*(wf%n_v))
!
      enddo ! End of batches over b
!
   end subroutine jacobian_transpose_ccsd_a2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD B2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the B2 term,
!!
!!       sum_c b_aicj F_cb - sum_k b_aibk F_jk + sum_ck b_aick L_ckjb 
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: b_aij_c ! b_aicj 
!
      real(dp), dimension(:,:), allocatable :: sigma_aij_b ! sigma_ai_bj contribution
!
      real(dp), dimension(:,:), allocatable :: g_ck_jb ! g_ckjb
      real(dp), dimension(:,:), allocatable :: g_ck_bj ! g_ckjb & g_cbjk
!
      real(dp), dimension(:,:), allocatable :: g_cb_jk ! g_cbjk 
      real(dp), dimension(:,:), allocatable :: g_cb_jk_restricted ! g_cbjk, batch over b 
!
      integer(i15) :: c = 0, j = 0, i = 0, a = 0, cj = 0, ai = 0, bj = 0, aij = 0, b = 0
      integer(i15) :: jb = 0, jk = 0, k = 0, ck = 0, cb = 0, cb_restricted = 0, bj_full = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0
      integer(i15) :: current_b_batch = 0
!
      type(batching_index) :: batch_b
!
      integer(i15) :: cb_offset = 0
!
!     :: Term 1. sum_c b_aicj F_cb ::
!
!     Reorder to b_aij_c = b_aicj 
!
      call mem%alloc(b_aij_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      b_aij_c = zero 
!
      do c = 1, wf%n_v
         do j = 1, wf%n_o
!
            cj = index_two(c, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  b_aij_c(aij, c) = b_ai_bj(ai, cj) ! b_aicj 
!
               enddo
            enddo
         enddo
      enddo
!
!     Calculate and add sum_c b_aicj F_cb = sum_c b_aij_c F_c_b
!
      call mem%alloc(sigma_aij_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, & 
                  wf%n_v,               &
                  wf%n_v,               &
                  one,                  &
                  b_aij_c,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%fock_ab,           & ! "F_c_b"
                  wf%n_v,               &
                  zero,                 &   
                  sigma_aij_b,          &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(b_aij_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aij_b(aij, b)
!  
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(sigma_aij_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
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
                  b_ai_bj,              & ! "b_aib_k"
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%fock_ij,           & ! "F_j_k"
                  wf%n_o,               &
                  one,                  &
                  sigma_ai_bj,          & ! "sigma_aib_j"
                  (wf%n_o)*(wf%n_v)**2)
!
!     :: Term 3. sum_ck b_aick L_ckjb :: 
!
!     Form g_ck_jb 
!
      call mem%alloc(g_ck_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_voov(g_ck_jb)
!
!     Reorder to g_ck_bj = g_ck_jb = g_ckjb 
!
      call mem%alloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ck_bj = zero 
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
            jb = index_two(j, b, wf%n_o)
!
            do I = 1, (wf%n_o)*(wf%n_v)
!
               g_ck_bj(I, bj) = g_ck_jb(I, jb) ! g_ckjb 
!
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ck_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add 2 * sum_ck b_aick g_ckjb = 2 * sum_ck b_aick g_ck_bj 
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  b_ai_bj,           & ! "b_ai_ck" 
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_ai_bj,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     - sum_ck b_aick g_cbjk
!
!     Prepare to batch over b to make g_cb_jk = g_cbjk successively
!
      call mem%alloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) ! g_cbjk reordered
      g_ck_bj = zero
!
      required = wf%integrals%get_required_vvoo()
!     
!     Initialize batching variable 
!
      call batch_b%init(wf%n_v)
      call mem%num_batch(batch_b, required)         
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
         call mem%alloc(g_cb_jk_restricted, (wf%n_v)*(batch_b%length), (wf%n_o)**2)
!
         call wf%get_vvoo(g_cb_jk_restricted, &
                           1,                  &
                           wf%n_v,             &
                           batch_b%first,      &
                           batch_b%last,       &
                           1,                  &
                           wf%n_o,             &
                           1,                  &
                           wf%n_o)
!
!        Place in reordered full space vector and deallocate restricted vector   
!
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do b = batch_b%first, batch_b%last
                  do c = 1, wf%n_v
!
                     ck = index_two(c, k, wf%n_v)
                     jk = index_two(j, k, wf%n_o)
!
                     cb_restricted = index_two(c, b - batch_b%first + 1, wf%n_v)
!
                     bj_full = index_two(b, j, wf%n_v)
!
                     g_ck_bj(ck, bj_full) = g_cb_jk_restricted(cb_restricted, jk)
!
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(g_cb_jk_restricted, (wf%n_v)*(batch_b%length), (wf%n_o)**2)
!
      enddo ! End of batches over b 
!
!     Add  - sum_ck b_aick g_cbjk = - sum_ck b_ai_ck g_ck_bj 
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_ai_bj,           & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_ai_bj,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_b2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_c2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD C2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll
!!
!!    Calculates the C2 term,
!!
!!       - sum_ck (b_ajck g_ibck + b_akcj g_ikcb) 
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: g_ib_ck ! g_ibck 
      real(dp), dimension(:,:), allocatable :: g_cb_ik ! g_cbik
      real(dp), dimension(:,:), allocatable :: g_ck_bi ! g_cbik
!
      real(dp), dimension(:,:), allocatable :: sigma_aj_ib ! sigma_ai_bj contribution
      real(dp), dimension(:,:), allocatable :: sigma_aj_bi ! sigma_ai_bj contribution
!
      real(dp), dimension(:,:), allocatable :: b_aj_ck ! b_akcj
!
      integer(i15) :: k = 0, j = 0, ik = 0, ib = 0, i = 0, ck = 0, cj = 0, cb = 0, c = 0
      integer(i15) :: bj = 0, bi = 0, b = 0, ak = 0, aj = 0, ai = 0, a = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0
      integer(i15) :: current_b_batch = 0
!
      type(batching_index) :: batch_b 
!
      integer(i15) :: cb_offset = 0
!
!     :: Term 1. - sum_ck b_ajck g_ibck ::
!
!     Form g_ib_ck
!
      call mem%alloc(g_ib_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovvo(g_ib_ck)
!
!     Calculate and add - sum_ck b_ajck g_ibck = - sum_ck b_aj_ck g_ib_ck^T(ck,ib)
!
      call mem%alloc(sigma_aj_ib, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_ai_bj,           & ! "b_aj_ck"
                  (wf%n_o)*(wf%n_v), &
                  g_ib_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_aj_ib,       &     
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ib_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
!
               ib = index_two(i, b, wf%n_o)
!
               do a = 1, wf%n_v
!
                  aj = index_two(a, j, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aj_ib(aj, ib)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(sigma_aj_ib, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 2. - sum_ck b_akcj g_ikcb = - sum_ck b_akcj g_cbik ::
!
!     Make g_ck_bi = g_cbik in batches over b 
!
      call mem%alloc(g_ck_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ck_bi = zero
!
      required = wf%integrals%get_required_vvoo()
!     
!     Initialize batching variable 
!
      call batch_b%init(wf%n_v)
      call mem%num_batch(batch_b, required)           
!
!     Loop over the b-batches 
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index 
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cb_ik = g_cbik 
!
         call mem%alloc(g_cb_ik, (wf%n_v)*(batch_b%length), (wf%n_o)**2)
!
         call wf%get_vvoo(g_cb_ik,        &
                           1,             &
                           wf%n_v,        &
                           batch_b%first, &
                           batch_b%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)
!
!        Place in reordered integral g_ck_bi = g_cbik  
!
         do k = 1, wf%n_o
            do i = 1, wf%n_o
               do b = batch_b%first, batch_b%last
                  do c = 1, wf%n_v
!
                     ck = index_two(c, k, wf%n_v)
                     bi = index_two(b, i, wf%n_v) ! Full space 
                     cb = index_two(c, b - batch_b%first + 1, wf%n_v) ! Restricted 
                     ik = index_two(i, k, wf%n_o)
!
                     g_ck_bi(ck, bi) = g_cb_ik(cb, ik) ! g_cbik
!
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(g_cb_ik, (wf%n_v)*(batch_b%length), (wf%n_o)**2)
!
      enddo ! End of batches over b 
!
!     Reorder to b_aj_ck = b_akcj 
!
      call mem%alloc(b_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      b_aj_ck = zero 
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do j = 1, wf%n_o
!
               cj = index_two(c, j, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ak = index_two(a, k, wf%n_v)
                  aj = index_two(a, j, wf%n_v)
!
                  b_aj_ck(aj, ck) = b_ai_bj(ak, cj) ! b_akcj
!  
               enddo
            enddo
         enddo
      enddo
!
!     Form and add - sum_ck b_akcj g_cbik = - sum_ck b_aj_ck g_ck_bi
!
      call mem%alloc(sigma_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_aj_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bi,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_aj_bi,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(b_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(g_ck_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
!
               bi = index_two(b, i, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  aj = index_two(a, j, wf%n_v)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aj_bi(aj, bi)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(sigma_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_c2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD D2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the D2 term,
!!
!!       2 * sum_ckdl b_aick L_jbld t_kl^cd 
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      real(dp), dimension(:,:), allocatable :: t_ck_dl ! t_kl^cd  
!
      real(dp), dimension(:,:), allocatable :: X_ck_bj ! An intermediate 
!
      real(dp), dimension(:,:), allocatable :: g_jb_ld ! g_jbld
      real(dp), dimension(:,:), allocatable :: L_dl_bj ! L_jbld
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, l = 0, lb = 0, jd = 0
      integer(i15) :: dl = 0, d = 0, ld = 0, bj = 0, jb = 0
!
!     Form g_jb_ld = g_jbld 
!
      call mem%alloc(g_jb_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_jb_ld)
!
!     Form L_dl_bj = L_jbld = 2 * g_jbld - g_jdlb
!                           = 2 * g_jb_ld(jb,ld) - g_jb_ld(jd,lb)
!
      call mem%alloc(L_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_dl_bj = zero
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
            jb = index_two(j, b, wf%n_o)
!
            do l = 1, wf%n_o
!
               lb = index_two(l, b, wf%n_o)
!
               do d = 1, wf%n_v
!
                  jd = index_two(j, d, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
                  dl = index_two(d, l, wf%n_v)
!
                  L_dl_bj(dl, bj) = two*g_jb_ld(jb, ld) - g_jb_ld(jd, lb) ! L_jbld
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jb_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form t_ck_dl = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ck_dl = zero 
!
      call squareup(wf%t2, t_ck_dl, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_ck_bj = sum_dl t_ck_dl L_dl_bj
!
      call mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_ck_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  L_dl_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(L_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add 2 * sum_ckdl b_aick L_jbld t_kl^cd = 2 * sum_ck b_ai_ck X_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  b_ai_bj,           & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_ai_bj,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_d2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD E2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the E2 term,
!!
!!       - sum_ckdl (b_aibl t_kl^cd L_kcjd + b_aicl t_kl^cd L_jbkd + b_aicj t_kl^cd L_ldkb)
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_kc_jd
      real(dp), dimension(:,:), allocatable :: L_j_ckd ! L_kcjd
      real(dp), dimension(:,:), allocatable :: L_dk_bj ! L_jbkd
      real(dp), dimension(:,:), allocatable :: L_ldk_b ! L_ldkb
!
      real(dp), dimension(:,:), allocatable :: t_ck_dl ! t_kl^cd
      real(dp), dimension(:,:), allocatable :: t_cl_dk ! t_kl^cd
!
      real(dp), dimension(:,:), allocatable :: X_j_l   ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_cl_bj ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_c_b   ! An intermediate, term 3
!
      real(dp), dimension(:,:), allocatable :: sigma_aij_b ! sigma_ai_bj contribution 
      real(dp), dimension(:,:), allocatable :: b_aij_c     ! b_aicj 
!
      integer(i15) :: kd = 0, kc = 0, k = 0, jd = 0, jc = 0, j = 0, d = 0
      integer(i15) :: ckd = 0, c = 0, l = 0, dk = 0, cl = 0, bjd = 0, bj = 0
      integer(i15) :: b = 0, ldk = 0, i = 0, dl = 0, cj = 0, bk = 0, aij = 0
      integer(i15) :: a = 0, ai = 0, ck = 0
!
!     :: Term 1. - sum_ckdl b_aibl t_kl^cd L_kcjd ::
!
!     X_j_l = sum_kcd L_kcjd t_kl^cd = sum_kcd L_j_ckd t_ck_dl 
!
!     Form g_kc_jd = g_kcjd
!
!
      call mem%alloc(g_kc_jd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_jd)
!
!     Form L_j_ckd = L_kcjd = 2 * g_kcjd - g_kdjc
!                           = 2 * g_kc_jd(kc,jd) - g_kc_jd(kd,jc)
!
      call mem%alloc(L_j_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      L_j_ckd = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
!
            kd = index_two(k, d, wf%n_o)
!
            do c = 1, wf%n_v
!
               kc = index_two(k, c, wf%n_o)
!
               ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
               do j = 1, wf%n_o
!
                  jd = index_two(j, d, wf%n_o)
                  jc = index_two(j, c, wf%n_o)
!
                  L_j_ckd(j, ckd) = two*g_kc_jd(kc, jd) - g_kc_jd(kd, jc) ! L_kcjd
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kc_jd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form t_ck_dl = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ck_dl = zero
!
      call squareup(wf%t2, t_ck_dl, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_j_l = sum_kcd L_kcjd t_kl^cd 
!                                 = sum_kcd L_j_ckd t_ckd_l
!
      call mem%alloc(X_j_l, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               & 
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  L_j_ckd,              &
                  wf%n_o,               &
                  t_ck_dl,              & ! "t_ckd_l"
                  (wf%n_o)*(wf%n_v)**2, & 
                  zero,                 &
                  X_j_l,                &
                  wf%n_o)
!
!
!     Add - sum_ckdl b_aibl t_kl^cd L_kcjd 
!         = - sum_l b_aib_l X_j_l^T(l,j)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  b_ai_bj,              & ! "b_aib_l"
                  (wf%n_o)*(wf%n_v)**2, &
                  X_j_l,                &
                  wf%n_o,               &
                  one,                  &
                  sigma_ai_bj,          & ! "sigma_aib_j"
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_j_l, wf%n_o, wf%n_o)
!
!     :: Term 2. -sum_ckdl b_aicl t_kl^cd L_jbkd ::
!
!     X_cl_bj = sum_kd t_kl^cd L_jbkd = sum_kd t_cl_dk L_dk_bj
!
!     We have L_j_ckd = L_kcjd   => L_j_ckd(k, bjd) = L_jbkd
!     We have t_ckd_l = t_kl^cd 
!
!     Reorder to L_dk_bj(dk,bj) = L_jbkd = L_j_ckd(k, bjd)
!
      call mem%alloc(L_dk_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_dk_bj = zero
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do k = 1, wf%n_o
               do d = 1, wf%n_v
!
                  dk = index_two(d, k, wf%n_v)
!
                  bjd = index_three(b, j, d, wf%n_v, wf%n_o)
!
                  L_dk_bj(dk, bj) = L_j_ckd(k, bjd) ! L_jbkd
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(L_j_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Reorder to t_cl_dk = t_kl^cd = t_ck_dl 
!
      call mem%alloc(t_cl_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_cl_dk = zero
!
      do k = 1, wf%n_o
         do d = 1, wf%n_v
!
            dk = index_two(d, k, wf%n_v)
!
            do l = 1, wf%n_o
!
               dl = index_two(d, l, wf%n_v)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  cl = index_two(c, l, wf%n_v)
!
                  t_cl_dk(cl, dk) = t_ck_dl(ck, dl) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_cl_bj = sum_kd t_kl^cd L_jbkd 
!                                   = sum_kd t_cl_dk L_dk_bj
!
      call mem%alloc(X_cl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_cl_dk,           &
                  (wf%n_o)*(wf%n_v), &
                  L_dk_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_cl_bj,           &
                  (wf%n_o)*(wf%n_v))
!
!     Add - sum_ckdl b_aicl t_kl^cd L_jbkd
!           = - sum_cl b_ai_cl X_cl_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_ai_bj,           & ! "b_ai_cl"
                  (wf%n_o)*(wf%n_v), &
                  X_cl_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_ai_bj,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_cl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 3. - sum_ckdl b_aicj t_kl^cd L_ldkb ::
!
!     - sum_c b_aij_c X_c_b,   X_c_b = sum_kdl t_kl^cd L_ldkb
!                                    = sum_kdl t_cl_dk L_ldk_b
!
!     We have L_dk_bj(dk,bj) = L_jbkd => L_dk_bj(bk,dl) = L_ldkb 
!
!     Reorder to L_ldk_b = L_ldkb = L_dk_bj(bk,dl)
!
      call mem%alloc(L_ldk_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      L_ldk_b = zero
!
      do b = 1, wf%n_v
         do k = 1, wf%n_o
!
            bk = index_two(b, k, wf%n_v)
!
            do d = 1, wf%n_v
               do l = 1, wf%n_o
!
                  dl = index_two(d, l, wf%n_v)
!
                  ldk = index_three(l, d, k, wf%n_o, wf%n_v)
!
                  L_ldk_b(ldk, b) = L_dk_bj(bk, dl) ! L_ldkb 
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(L_dk_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form the intermediate X_c_b = sum_kdl t_kl^cd L_ldkb
!                                 = sum_kdl t_cl_dk L_ldk_b
!
      call mem%alloc(X_c_b, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, & 
                  one,                  &
                  t_cl_dk,              & ! "t_c_ldk"
                  wf%n_v,               &
                  L_ldk_b,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_c_b,                &
                  wf%n_v)
!
      call mem%dealloc(t_cl_dk, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(L_ldk_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     Reorder to b_aij_c = b_aicj 
!
      call mem%alloc(b_aij_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      b_aij_c = zero
!
      do c = 1, wf%n_v
         do j = 1, wf%n_o
!
            cj = index_two(c, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  b_aij_c(aij, c) = b_ai_bj(ai, cj) ! b_aicj
!
               enddo
            enddo
         enddo
      enddo
!
!     Form and add - sum_ckdl b_aicj t_kl^cd L_ldkb
!                  = - sum_c b_aij_c X_c_b
!
      call mem%alloc(sigma_aij_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, & 
                  wf%n_v,               &
                  wf%n_v,               &
                  -one,                 &
                  b_aij_c,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  X_c_b,                &
                  wf%n_v,               &
                  zero,                 &
                  sigma_aij_b,          &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(X_c_b, wf%n_v, wf%n_v)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aij_b(aij, b)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(sigma_aij_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_e2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD F2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the F2 term,
!!
!!       - sum_ckdl (b_alck t_kl^cd L_jbid + b_ajck t_kl^cd L_ldib + b_djck t_kl^cd L_ialb)
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      real(dp), dimension(:,:), allocatable :: t_lck_d ! t_kl^cd
      real(dp), dimension(:,:), allocatable :: t_l_ckd ! t_kl^cd
      real(dp), dimension(:,:), allocatable :: t_ck_dl ! t_kl^cd
!
      real(dp), dimension(:,:), allocatable :: L_jb_J
!
      real(dp), dimension(:,:), allocatable :: g_jb_id ! g_jbid 
      real(dp), dimension(:,:), allocatable :: L_d_ibj ! L_jbid
      real(dp), dimension(:,:), allocatable :: L_aib_l ! L_ialb
      real(dp), dimension(:,:), allocatable :: L_dl_bi ! L_ldib
!
      real(dp), dimension(:,:), allocatable :: X_a_d   ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_ck_bi ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_l_j   ! An intermediate, term 3
!
      real(dp), dimension(:,:), allocatable :: sigma_aj_bi ! sigma_ai_bj contriution
!
      integer(i15) :: k = 0, l = 0, lck = 0, d = 0, c = 0, jd = 0, jb = 0
      integer(i15) :: id = 0, j = 0, i = 0, ibj = 0, ib = 0, dl = 0, ckdl = 0
      integer(i15) :: b = 0, ck = 0, idl = 0, ckd = 0, bl = 0, bj = 0, bi = 0
      integer(i15) :: aj = 0, ai = 0, aib = 0, a = 0
!
!     :: Term 1. - sum_ckdl b_alck t_kl^cd L_jbid :: 
!
!     X_a_d = b_a_lck t_lck_d
!
!     Form t_lck_d = t_kl^cd
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      t_lck_d = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!     
               do l = 1, wf%n_o
!
                  dl = index_two(d, l, wf%n_v)
!
                  lck = index_three(l, c, k, wf%n_o, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
!
                  t_lck_d(lck, d) = wf%t2(ckdl, 1) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_a_d = sum_lck b_a_lck t_lck_d
!
      call mem%alloc(X_a_d, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               & 
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_ai_bj,              & ! "b_a_lck"
                  wf%n_v,               &
                  t_lck_d,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_a_d,                &
                  wf%n_v)
!
!     Form g_jb_id = g_jbid 
!
      call mem%alloc(g_jb_id, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%get_ovov(g_jb_id)
!
!     Form L_d_ibj = L_jbid = 2 * g_jbid - g_jdib 
!                           = 2 * g_jb_id(jb,id) - g_jb_id(jd,ib)
!
      call mem%alloc(L_d_ibj, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      L_d_ibj = zero
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
               do d = 1, wf%n_v
!
                  id = index_two(i, d, wf%n_o)
                  jd = index_two(j, d, wf%n_o)
!
                  L_d_ibj(d, ibj) = two*g_jb_id(jb, id) - g_jb_id(jd, ib) ! L_jbid
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jb_id, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Add - sum_ckdl b_alck t_kl^cd L_jbid
!         = - sum_d X_a_d L_d_ibj
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  X_a_d,                &
                  wf%n_v,               &
                  L_d_ibj,              &
                  wf%n_v,               &
                  one,                  &
                  sigma_ai_bj,          & ! "sigma_a_ibj"
                  wf%n_v)
!
      call mem%dealloc(X_a_d, wf%n_v, wf%n_v)
!
!     :: Term 2. - sum_ckdl b_ajck t_kl^cd L_ldib ::
!
!     X_ck_bi = sum_dl t_ck_dl L_dl_bi
!
!     Reorder to t_ck_dl = t_lck_d = t_kl^cd
!
      call mem%alloc(t_ck_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      t_ck_dl = zero
!
      do l = 1, wf%n_o
         do d = 1, wf%n_v
!
            dl = index_two(d, l, wf%n_v)
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
!
                  lck = index_three(l, c, k, wf%n_o, wf%n_v)
!
                  t_ck_dl(ck, dl) = t_lck_d(lck, d) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     We have L_d_ibj = L_jbid => L_d_ibj(b,idl) = L_ldib
!  
!     Form L_dl_bi = L_ldib = L_d_ibj(b,idl)
!
      call mem%alloc(L_dl_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      L_dl_bi = zero
!
      do i = 1, wf%n_o
         do b = 1, wf%n_v
!
            bi = index_two(b, i, wf%n_v)
!
            do l = 1, wf%n_o
               do d = 1, wf%n_v
!
                  dl = index_two(d, l, wf%n_v)
!
                  idl = index_three(i, d, l, wf%n_o, wf%n_v)
!
                  L_dl_bi(dl, bi) = L_d_ibj(b, idl) ! L_ldib
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(L_d_ibj, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     Form the intermediate X_ck_bi = sum_dl t_ck_dl L_dl_bi
!
      call mem%alloc(X_ck_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), & 
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  t_ck_dl,           &
                  (wf%n_v)*(wf%n_o), &
                  L_dl_bi,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ck_bi,           &
                  (wf%n_v)*(wf%n_o))
!
!     Form and add - sum_ckdl b_ajck t_kl^cd L_ldib = - sum_ck b_aj_ck X_ck_bi
!
      call mem%alloc(sigma_aj_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            & 
                  (wf%n_v)*(wf%n_o), & 
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  -one,              &
                  b_ai_bj,           & ! "b_aj_ck"
                  (wf%n_v)*(wf%n_o), &
                  X_ck_bi,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  sigma_aj_bi,       &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_ck_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
!
               bi = index_two(b, i, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  aj = index_two(a, j, wf%n_v)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aj_bi(aj, bi)
!
               enddo
            enddo
         enddo
      enddo  
!
      call mem%dealloc(sigma_aj_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 3. - sum_ckdl b_djck t_kl^cd L_ialb = - sum_ckdl b_ckdj t_kl^cd L_ialb ::
!
!     X_l_j = sum_ckd t_l_ckd b_ckd_j
!
!     Reorder to t_l_ckd = t_ck_dl = t_kl^cd 
!
      call mem%alloc(t_l_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      t_l_ckd = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
               do l = 1, wf%n_o
!
                  dl = index_two(d, l, wf%n_v)
!
                  t_l_ckd(l, ckd) = t_ck_dl(ck, dl) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_l_j = sum_ckd t_kl^cd b_ckdj 
!                                 = sum_ckd t_l_ckd b_ckd_j 
!
      call mem%alloc(X_l_j, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_l_ckd,              &
                  wf%n_o,               &
                  b_ai_bj,              & ! "b_ckd_j"
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_l_j,                &
                  wf%n_o)
!
      call mem%dealloc(t_l_ckd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     - sum_ckdl b_ckdj t_kl^cd L_ialb = - sum_l L_ialb X_l_j
!
!     We have L_dl_bi = L_ldib 
!
!     Form L_aib_l = L_ialb = L_dl_bi(ai,bl)
!
      call mem%alloc(L_aib_l, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      L_aib_l = zero
!
      do l = 1, wf%n_o
         do b = 1, wf%n_v
!
            bl = index_two(b, l, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                  L_aib_l(aib, l) = L_dl_bi(ai, bl) ! L_ialb
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(L_dl_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add - sum_ckdl b_ckdj t_kl^cd L_ialb = - sum_l L_aib_l X_l_j
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  L_aib_l,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  X_l_j,                &
                  wf%n_o,               &
                  one,                  &
                  sigma_ai_bj,          & ! "sigma_aib_j"
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(L_aib_l, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call mem%dealloc(X_l_j, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_f2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD G2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the G2 term,
!!
!!       sum_ckdl (b_alcj t_kl^cd g_kbid + b_ajcl t_kl^cd g_kdib)
!! 
!!    and adds it to the transformed vector sigma_ai_bj.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      real(dp), dimension(:,:), allocatable :: t_cl_dk ! t_kl^cd
      real(dp), dimension(:,:), allocatable :: t_cl_kd ! t_kl^cd
!
      real(dp), dimension(:,:), allocatable :: g_kb_id ! g_kbid
      real(dp), dimension(:,:), allocatable :: g_dk_bi ! g_kbid
      real(dp), dimension(:,:), allocatable :: g_kd_ib ! g_kdib
!
      real(dp), dimension(:,:), allocatable :: b_aj_cl     ! b_alcj
      real(dp), dimension(:,:), allocatable :: sigma_aj_bi ! sigma_ai_bj contribution
      real(dp), dimension(:,:), allocatable :: sigma_aj_ib ! sigma_ai_bj contribution
!
      real(dp), dimension(:,:), allocatable :: X_cl_bi ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_aj_kd ! An intermediate, term 2
!
      integer(i15) :: k = 0, d = 0, l = 0, c = 0, dk = 0, dl = 0, ck = 0, cl = 0
      integer(i15) :: ckdl = 0, i = 0, b = 0, bi = 0, kb = 0, id = 0, aj = 0
      integer(i15) :: j = 0, cj = 0, bj = 0, al = 0, ai = 0, a = 0, kd = 0, ib = 0
!
!     :: Term 1. sum_ckdl b_alcj t_kl^cd g_kbid ::
!
!     Form t_cl_dk = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_cl_dk, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      t_cl_dk = zero
!
      do k = 1, wf%n_o
         do d = 1, wf%n_v
!
            dk = index_two(d, k, wf%n_v)
!
            do l = 1, wf%n_o
!
               dl = index_two(d, l, wf%n_v)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  cl = index_two(c, l, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
!
                  t_cl_dk(cl, dk) = wf%t2(ckdl, 1) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      !call wf%destruct_double_amplitudes
!
!     Form g_kb_id = g_kbid
!
      call mem%alloc(g_kb_id, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%get_ovov(g_kb_id)
!
!     Reorder to g_dk_bi = g_kb_id
!
      call mem%alloc(g_dk_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      g_dk_bi = zero
!
      do i = 1, wf%n_o
         do b = 1, wf%n_v
!
            bi = index_two(b, i, wf%n_v)
!
            do k = 1, wf%n_o
!
               kb = index_two(k, b, wf%n_o)
!
               do d = 1, wf%n_v
!
                  id = index_two(i, d, wf%n_o)
                  dk = index_two(d, k, wf%n_v)
!
                  g_dk_bi(dk, bi) = g_kb_id(kb, id) ! g_kbid
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kb_id, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form the intermediate X_cl_bi = sum_dk t_cl_dk g_dk_bi
!
      call mem%alloc(X_cl_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  t_cl_dk,           &
                  (wf%n_v)*(wf%n_o), &
                  g_dk_bi,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_cl_bi,           &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(t_cl_dk, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(g_dk_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Reorder to b_aj_cl = b_alcj
!
      call mem%alloc(b_aj_cl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      b_aj_cl = zero 
!
      do l = 1, wf%n_o
         do c = 1, wf%n_v
!
            cl = index_two(c, l, wf%n_v)
!
            do j = 1, wf%n_o
!
               cj = index_two(c, j, wf%n_v)
!
               do a = 1, wf%n_v
!
                  al = index_two(a, l, wf%n_v)
                  aj = index_two(a, j, wf%n_v)
!
                  b_aj_cl(aj, cl) = b_ai_bj(al, cj) ! b_alcj
!
               enddo
            enddo
         enddo
      enddo
!
!     Calculate and add sum_ckdl b_alcj t_kl^cd g_kbid
!                       = sum_cl b_aj_cl X_cl_bi
!
      call mem%alloc(sigma_aj_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), & 
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_aj_cl,           &
                  (wf%n_v)*(wf%n_o), &
                  X_cl_bi,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  sigma_aj_bi,       &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_cl_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(b_aj_cl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
!
               bi = index_two(b, i, wf%n_v)
!
               do a = 1, wf%n_v
!
                  aj = index_two(a, j, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aj_bi(aj, bi)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(sigma_aj_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 2. sum_ckdl b_ajcl t_kl^cd g_kdib ::
!
!     Form t_cl_kd = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_cl_kd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      t_cl_kd = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
!
            kd = index_two(k, d, wf%n_o)
!
            do l = 1, wf%n_o
!
               dl = index_two(d, l, wf%n_v)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  cl = index_two(c, l, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
!
                  t_cl_kd(cl, kd) = wf%t2(ckdl, 1) ! t_kl^cd
!
               enddo
            enddo
         enddo
      enddo
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_aj_kd = sum_cl b_aj_cl t_cl_kd 
!
      call mem%alloc(X_aj_kd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_ai_bj,           & ! "b_aj_cl"
                  (wf%n_v)*(wf%n_o), &
                  t_cl_kd,           &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_aj_kd,           &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(t_cl_kd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form g_kd_ib = g_kdib 
!
      call mem%alloc(g_kd_ib, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kd_ib)
!
!     Form and add sum_ckdl b_ajcl t_kl^cd g_kdib
!                  = sum_kd X_aj_kd g_kd_ib
!
      call mem%alloc(sigma_aj_ib, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  X_aj_kd,           &
                  (wf%n_o)*(wf%n_v), &
                  g_kd_ib,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_aj_ib,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_kd_ib, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(X_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
!
               ib = index_two(i, b, wf%n_o)
!
               do a = 1, wf%n_v
!
                  aj = index_two(a, j, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
!
                  sigma_ai_bj(ai, bj) = sigma_ai_bj(ai, bj) + sigma_aj_ib(aj, ib)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(sigma_aj_ib, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_g2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_h2_ccsd(wf, sigma_ab_ij, b_ab_ij)
!!
!!    Jacobian transpose CCSD H2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the H2 term,
!!
!!       sum_kl b_akbl g_ikjl + sum_cd b_cidj g_cadb 
!! 
!!    and adds it to the transformed vector sigma_ab_ij.
!!
!!    In this routine, the b and sigma vectors are ordered as
!!
!!       b_ab_ij = b_ai_bj 
!!       sigma_ab_ij = sigma_ab_ij
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: b_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: sigma_ab_ij
!
      real(dp), dimension(:,:), allocatable :: sigma_ab_ij_batch
!
      real(dp), dimension(:,:), allocatable :: g_ca_db ! g_cadb 
      real(dp), dimension(:,:), allocatable :: g_ab_cd ! g_cadb 
!
      real(dp), dimension(:,:), allocatable :: g_ik_jl ! g_ikjl
      real(dp), dimension(:,:), allocatable :: g_kl_ij ! g_ikjl
!
      integer(i15) :: l = 0, kl = 0, k = 0, jl = 0, j = 0, ik = 0, ij = 0
      integer(i15) :: i = 0, db = 0, d = 0, cd = 0, ca = 0, c = 0, b = 0, a = 0
!
!     Batching variables 
!
      integer(i15) :: required = 0 
!
      integer(i15) :: current_a_batch = 0
      integer(i15) :: current_b_batch = 0 
!
      type(batching_index) :: batch_a 
      type(batching_index) :: batch_b
!
      integer(i15) :: ab = 0, ab_full = 0
!
!     :: Term 1. sum_kl b_akbl g_ikjl ::
!
!     Form g_ik_jl 
!
      call mem%alloc(g_ik_jl, (wf%n_o)**2, (wf%n_o)**2)
!
      call wf%get_oooo(g_ik_jl)
!
!     Reorder to g_kl_ij = g_ik_jl = g_ikjl 
!
      call mem%alloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
      g_kl_ij = zero
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do l = 1, wf%n_o
!
               jl = index_two(j, l, wf%n_o)
!
               do k = 1, wf%n_o
!
                  ik = index_two(i, k, wf%n_o)
                  kl = index_two(k, l, wf%n_o)
!
                  g_kl_ij(kl, ij) = g_ik_jl(ik, jl) ! g_ikjl 
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ik_jl, (wf%n_o)**2, (wf%n_o)**2)
!
!     Add sum_kl b_akbl g_ikjl = sum_kl b_ab_kl g_kl_ij 
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  b_ab_ij,     & ! "b_ab_kl"
                  (wf%n_v)**2, &
                  g_kl_ij,     &
                  (wf%n_o)**2, &
                  one,         &
                  sigma_ab_ij, &
                  (wf%n_v)**2)
!
      call mem%dealloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!  
!     :: Term 2. sum_cd b_cidj g_cadb ::
!
!     sum_cd b_cidj g_cadb = sum_cd g_cadb b_cd_ij
!                          = sum_cd g_ab_cd b_cd_ij
!
!     Prepare batching over a and b 
!
      required = wf%integrals%get_required_vvvv() + (wf%n_v**4) + (wf%n_o**2)*(wf%n_v**2)
!     
!     Initialize batching indices 
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v) 
!
      call mem%num_two_batch(batch_a, batch_b, required)
!
!     Loop over a-batches 
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each a batch, get the limits for the a index 
!
         call batch_a%determine_limits(current_a_batch)
!
!        Loop over b-batches
!
         do current_b_batch = 1, batch_b%num_batches
!
!           For each b batch, get the limits for the b index 
!
            call batch_b%determine_limits(current_b_batch)
!
!           Form g_ca_db = g_cadb 
!
            call mem%alloc(g_ca_db, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
            call wf%get_vvvv(g_ca_db,        &
                              1,             &
                              wf%n_v,        &
                              batch_a%first, &
                              batch_a%last,  &
                              1,             &
                              wf%n_v,        &
                              batch_b%first, &
                              batch_b%last)
!
!           Reorder to g_ab_cd = g_ca_db = g_cadb 
!
            call mem%alloc(g_ab_cd, (batch_a%length)*(batch_b%length), (wf%n_v)**2)
            g_ab_cd = zero 
!
            do d = 1, wf%n_v
               do c = 1, wf%n_v
!
                  cd = index_two(c, d, wf%n_v)
!
                  do b = 1, batch_b%length
!
                     db = index_two(d, b, wf%n_v)
!
                     do a = 1, batch_a%length
!
                        ca = index_two(c, a, wf%n_v)
                        ab = index_two(a, b, batch_a%length)
!
                        g_ab_cd(ab, cd) = g_ca_db(ca, db) ! g_cadb
!
                     enddo
                  enddo
               enddo
            enddo
!  
            call mem%dealloc(g_ca_db, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!           Calculate sigma_ab_ij_batch = sum_cd g_ab_cd b_cd_ij
!           and add it to the full space sigma vector 
!
            call mem%alloc(sigma_ab_ij_batch, (batch_a%length)*(batch_b%length), (wf%n_o)**2)
!
            call dgemm('N','N',                            &
                        (batch_a%length)*(batch_b%length), & 
                        (wf%n_o)**2,                       &
                        (wf%n_v)**2,                       &
                        one,                               &
                        g_ab_cd,                           &
                        (batch_a%length)*(batch_b%length), &
                        b_ab_ij,                           & ! "b_cd_ij"
                        (wf%n_v)**2,                       &
                        zero,                              &
                        sigma_ab_ij_batch,                 &
                        (batch_a%length)*(batch_b%length))
!
            call mem%dealloc(g_ab_cd, (batch_a%length)*(batch_b%length), (wf%n_v)**2)
!
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  ij = index_two(i, j, wf%n_o)
!
                  do b = 1, batch_b%length
                     do a = 1, batch_a%length
!
                        ab = index_two(a, b, batch_a%length)
!
                        ab_full = index_two(a + batch_a%first - 1, b + batch_b%first - 1, wf%n_v)
!
                        sigma_ab_ij(ab_full, ij) = sigma_ab_ij(ab_full, ij) &
                                                 + sigma_ab_ij_batch(ab, ij)
!
                     enddo
                  enddo
               enddo
            enddo
!
            call mem%dealloc(sigma_ab_ij_batch, (batch_a%length)*(batch_b%length), (wf%n_o)**2)
!
         enddo ! End of batches over b 
      enddo ! End of batches over a
!
   end subroutine jacobian_transpose_ccsd_h2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_i2_ccsd(wf, sigma_ab_ij, b_ab_ij)
!!
!!    Jacobian transpose CCSD I2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the I2 term,
!!
!!       sum_ckdl b_cidj t_kl^cd g_kalb + sum_ckdl b_akbl t_kl^cd g_icjd 
!! 
!!    and adds it to the transformed vector sigma_ab_ij.
!!
!!    In this routine, the b and sigma vectors are ordered as
!!
!!       b_ab_ij = b_ai_bj 
!!       sigma_ab_ij = sigma_ab_ij
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: b_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: sigma_ab_ij
!
      real(dp), dimension(:,:), allocatable :: t_kl_cd ! t_kl^cd 
!
      real(dp), dimension(:,:), allocatable :: g_ka_lb ! g_kalb 
      real(dp), dimension(:,:), allocatable :: g_ab_kl ! g_kalb
!
      real(dp), dimension(:,:), allocatable :: X_kl_ij ! An intermediate, terms 1 & 2
!
      integer(i15) :: d = 0, c = 0, cd = 0, l = 0, dl = 0, k = 0, ck = 0, kl = 0
      integer(i15) :: ckdl = 0, lb = 0, ka = 0, b = 0, ab = 0, a = 0
!
!     :: Term 1. sum_ckdl b_cidj t_kl^cd g_kalb :: 
!
!     sum_ckdl t_kl_cd b_cd_ij 
!
!     Form t_kl_cd = t_kl^cd 
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
      t_kl_cd = zero 
!
      do d = 1, wf%n_v
         do c = 1, wf%n_v
!
            cd = index_two(c, d, wf%n_v)
!
            do l = 1, wf%n_o
!
               dl = index_two(d, l, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ck = index_two(c, k, wf%n_v)
                  kl = index_two(k, l, wf%n_o)
!
                  ckdl = index_packed(ck, dl)
!
                  t_kl_cd(kl, cd) = wf%t2(ckdl, 1) ! t_kl^cd 
!
               enddo
            enddo
         enddo
      enddo
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_kl_ij = sum_cd t_kl_cd b_cd_ij 
!
      call mem%alloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  t_kl_cd,     &
                  (wf%n_o)**2, &
                  b_ab_ij,     & ! "b_cd_ij"
                  (wf%n_v)**2, &
                  zero,        &
                  X_kl_ij,     &
                  (wf%n_o)**2)
!
!     Form g_ka_lb = g_kalb 
!
      call mem%alloc(g_ka_lb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_ka_lb)
!
!     Reorder to g_ab_kl = g_ka_lb = g_kalb 
!
      call mem%alloc(g_ab_kl, (wf%n_v)**2, (wf%n_o)**2)
      g_ab_kl = zero 
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
!
            kl = index_two(k, l, wf%n_o)
!
            do b = 1, wf%n_v
!
               lb = index_two(l, b, wf%n_o)
!
               do a = 1, wf%n_v
!
                  ka = index_two(k, a, wf%n_o)
                  ab = index_two(a, b, wf%n_v)
!
                  g_ab_kl(ab, kl) = g_ka_lb(ka, lb) ! g_kalb 
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ka_lb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add sum_ckdl b_cidj t_kl^cd g_kalb
!         = sum_kl g_ab_kl X_kl_ij 
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  g_ab_kl,     &
                  (wf%n_v)**2, &
                  X_kl_ij,     &
                  (wf%n_o)**2, &
                  one,         &
                  sigma_ab_ij, &
                  (wf%n_v)**2)
!
!     :: Term 2. sum_ckdl b_akbl t_kl^cd g_icjd :: 
!
!     Repurpose X_kl_ij to make sum_cd t_kl^cd g_icjd
!                               = sum_cd t_kl_cd g_cd_ij
!                               = sum_cd t_kl_cd g_ab_kl(cd,ij)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  t_kl_cd,     &
                  (wf%n_o)**2, &
                  g_ab_kl,     & ! "g_cd_ij"
                  (wf%n_v)**2, &
                  zero,        &
                  X_kl_ij,     &
                  (wf%n_o)**2)
!
      call mem%dealloc(g_ab_kl, (wf%n_v)**2, (wf%n_o)**2)
      call mem%dealloc(t_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!     Add sum_ckdl b_akbl t_kl^cd g_icjd
!         = sum_kl b_ab_kl X_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  b_ab_ij,     & ! "b_ab_kl"
                  (wf%n_v)**2, &
                  X_kl_ij,     &
                  (wf%n_o)**2, &
                  one,         &
                  sigma_ab_ij, &
                  (wf%n_v)**2)
!
      call mem%dealloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
   end subroutine jacobian_transpose_ccsd_i2_ccsd
!
!
end submodule jacobian_transpose_ccsd