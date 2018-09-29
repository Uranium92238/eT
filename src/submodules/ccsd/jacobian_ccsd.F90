submodule (ccsd_class) jacobian_ccsd
!
!!
!!    Jacobian submodule (ccsd)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Andreas Skeidsvoll, 2018
!!
!
   implicit none
!
!
contains
!
!
  module subroutine jacobian_ccsd_a1_ccsd(wf, rho_a_i, c_a_i)
!!
!!    Jacobian CCSD A1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai^A1 = sum_ckdl L_lckd (u_li^ca c_dk  - t_li^cd c_ak - t_lk^ad c_ci)
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^B1,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_lc_kd 
      real(dp), dimension(:,:), allocatable :: L_lc_dk 
      real(dp), dimension(:,:), allocatable :: L_lk_dc 
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: X_lc
      real(dp), dimension(:,:), allocatable :: X_i_k
      real(dp), dimension(:,:), allocatable :: X_a_c
!
!     Amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_il_cd 
      real(dp), dimension(:,:), allocatable :: t_al_kd 
      real(dp), dimension(:,:), allocatable :: t_cl_ai 
      real(dp), dimension(:,:), allocatable :: u_ai_lc 
!
!     :: Term 1: sum_ckdl L_lckd u_li^ca c_dk ::
!
!     g_lc_kd = g_lckd
!
      call mem%alloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_lc_kd)
!
      call mem%alloc(L_lc_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_lc_dk = zero
!
!     L_lc_dk(lc,dk) = L_lckd = 2*g_lckd - g_ldkc = 2*g_lc_kd(lc,kd) - g_lc_kd(ld,kc)
!                                                             1243             1342
!
      call add_1342_to_1234(-one, g_lc_kd, L_lc_dk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call add_1243_to_1234(two, g_lc_kd, L_lc_dk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Deallocate g_lc_kd
!
      call mem%dealloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     X_lc = sum_kd L_lckd c_dk = sum_kd L_lc_dk c_dk
!
      call mem%alloc(X_lc, (wf%n_o)*(wf%n_v), 1)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lc_dk,           &
                  (wf%n_o)*(wf%n_v), &
                  c_a_i,             &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lc,              &
                  (wf%n_o)*(wf%n_v))
!
!
!     Form u_ai_lc = u_li^ca = 2 * t_li^ca - t_il^ca = 2 * t2(clai,1) - t2(cial,1)
!
!     Allocate & read double amplitudes
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_cl_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_cl_ai, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     t_cl_ai(cl,ai) = t_li^ca
!
!     u_ai_lc(ai, lc) = 2 * t_li^ca - t_il^ca = 2 * t_cl_ai(cl, ai) - t_cl_ai(ci, al)
!                                                           4312              4213
!
      call mem%alloc(u_ai_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      u_ai_lc = zero
!
      call add_4312_to_1234(two, t_cl_ai, u_ai_lc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_4213_to_1234(-one, t_cl_ai, u_ai_lc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_a_i =+ sum_lc u_ai_lc X_lc
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  u_ai_lc,           &
                  (wf%n_v)*(wf%n_o), &
                  X_lc,              &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  rho_a_i,           &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocations (keep L_lc_dk = L_lckd)
!
      call mem%dealloc(u_ai_lc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(X_lc, (wf%n_v)*(wf%n_o), 1)
!
!
!     :: Term 2. - sum_ckdl L_lckd t_li^cd c_ak ::
!
!     Reorder amplitudes to t_il_cd = t_li^cd
!
      call mem%alloc(t_il_cd, (wf%n_o)**2, (wf%n_v)**2)
      t_il_cd = zero
!
!     t_il_cd(il, cd) = t_li^cd = t_cl_ai(cl, di)
!             4213                        1234
!
      call sort_1234_to_4213(t_cl_ai, t_il_cd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_cl_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate X_i_k = sum_cdl L_lcd_k t_i_lcd
!
      call mem%alloc(X_i_k, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_il_cd,              & ! t_i_lcd
                  (wf%n_o),             &
                  L_lc_dk,              & ! L_lcd_k
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_i_k,                &
                  wf%n_o)
!
!     Calculate rho_a_i =+ - sum_k c_a_i(a,k) X_i_k(i, k)
!
      call dgemm('N', 'T', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  c_a_i,   & ! c_a_k
                  wf%n_v,  &
                  X_i_k,   &
                  wf%n_o,  &
                  one,     &
                  rho_a_i, &
                  wf%n_v)
!
!     Deallocations (keep L_lcd_k = L_lc,kd)
!
      call mem%dealloc(X_i_k, wf%n_o, wf%n_o)
!
!     :: Term 3: - sum_ckdl L_lckd t_lk^ad c_ci ::
!
!     Reorder to L_lk_dc(lk,dc) = L_lckd = L_lc_dk(lc,dk)
!                        1432                      1234
!
      call mem%alloc(L_lk_dc, (wf%n_o)**2, (wf%n_v)**2)
      call sort_1234_to_1432(L_lc_dk, L_lk_dc, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_lc_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder amplitudes to t_al_kd(al,kd) = t_lk^ad = t_il_cd(kl, ad)
!                                   3214                       1234
!
!
      call mem%alloc(t_al_kd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call sort_1234_to_3214(t_il_cd, t_al_kd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_il_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!     Calculate X_a_c = sum_kdl t_a_lkd L_lkd_c
!
      call mem%alloc(X_a_c, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  t_al_kd,              & ! t_a_lkd
                  wf%n_v,               &
                  L_lk_dc,              & ! L_lkd_c
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_a_c,                &
                  wf%n_v)
!
      call mem%dealloc(L_lk_dc, (wf%n_o)**2, (wf%n_v)**2)
      call mem%dealloc(t_al_kd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Calculate rho_a_i =+ - sum_c X_a_c(a,c) c_a_i(c,i)
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  X_a_c,   &
                  wf%n_v,  &
                  c_a_i,   &
                  wf%n_v,  &
                  one,     &
                  rho_a_i, &
                  wf%n_v)
!
      call mem%dealloc(X_a_c, wf%n_v, wf%n_v)
!
   end subroutine jacobian_ccsd_a1_ccsd
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai^B1 = sum_bj F_jb (2*c_ai_bj  -  c_aj_bi)
!!              = sum_bj F_jb v_ai_jb
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^B1,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      real(dp), dimension(:,:), allocatable :: v_ai_jb
!
!     Construct v_aibj = 2*c_aibj - c_ajbi
!
!               v_ai_jb(ai,jb) = 2*c_ai_bj(ai,bj) - c_ai_bj(aj,bi)
!                        1234              1243             1342
!
      call mem%alloc(v_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      v_ai_jb = zero
!
      call add_1243_to_1234(two, c_ai_bj, v_ai_jb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one, c_ai_bj, v_ai_jb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',           &                  
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  v_ai_jb,           &
                  (wf%n_o)*(wf%n_v), &
                  wf%fock_ia,        & ! F_jb
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_a_i,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(v_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))

!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai^C1 = - sum_bjk L_jikb c_aj_bk
!!              = - sum_bjk (2*g_jikb - g_kijb) c_aj_bk
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^C1,
!!    where c_ai_bj(ai,bj) = c_aibj above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      real(dp), dimension(:,:), allocatable :: g_ji_kb
      real(dp), dimension(:,:), allocatable :: L_jb_ki
      real(dp), dimension(:,:), allocatable :: c_a_jbk
!
!     Construct the integral g_ji_kb = sum_J L_ji_J * L_kb_J
!
      call mem%alloc(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_ooov(g_ji_kb)
!
!     Constructing L_jikb = 2*g_jikb - g_kijb
!
!                  L_jb_ki(jb,ki) = 2*g_ji_kb(ji,kb) - g_ji_kb(ki,jb)
!                           1234              1432             3412
!
!
      call mem%alloc(L_jb_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
      L_jb_ki = zero
!
      call add_1432_to_1234(two, g_ji_kb, L_jb_ki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call add_3412_to_1234(-one, g_ji_kb, L_jb_ki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  c_ai_bj,                & ! c_a_jbk
                  wf%n_v,                 &
                  L_jb_ki,                & ! L_jkb_i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  rho_a_i,                &
                  wf%n_v)
!
      call mem%dealloc(L_jb_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
   module subroutine jacobian_ccsd_d1_ccsd(wf, rho_a_i, c_bi_cj)
!!
!!    Jacobian CCSD D1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^D1,
!!    where c_bi_cj(bi,cj) = c_bicj above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_bi_cj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      real(dp), dimension(:,:), allocatable :: c_bj_ci
!
!     Variables for batching
!
      integer(i15) :: required = 0
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ab_jc ! g_abjc
      real(dp), dimension(:,:), allocatable :: L_ab_jc ! L_abjc
!
!     Prepare for batching over index a
!
      required = wf%integrals%get_required_vvov() + (wf%n_v**3)*(wf%n_o)
!
!     Initialize batching variable
!
      call batch_a%init(wf%n_v)
      call mem%num_batch(batch_a, required)
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Determine the limits for the current a-batch
!
         call batch_a%determine_limits(current_a_batch)
!
!        Form g_ab_jc = g_abjc
!
         call mem%alloc(g_ab_jc, (batch_a%length)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
         call wf%get_vvov(g_ab_jc,        &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Construct L_ab_jc
!
         call mem%alloc(L_ab_jc, (batch_a%length)*(wf%n_v), (wf%n_o)*(wf%n_v))
         L_ab_jc = zero
!
         call add_1432_to_1234(-one, g_ab_jc, L_ab_jc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
         call daxpy((batch_a%length)*(wf%n_v)**2*(wf%n_o), two, g_ab_jc, 1, L_ab_jc, 1)
!
!        Reorder c_bi_cj to c_bj_ci
!
         call mem%alloc(c_bj_ci, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
         call sort_1234_to_1432(c_bi_cj, c_bj_ci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call mem%dealloc(g_ab_jc, (batch_a%length)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
         call dgemm('N', 'N',                   &
                     batch_a%length,            &
                     wf%n_o,                    &
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     L_ab_jc,                   & ! L_a_bjc
                     batch_a%length,            &
                     c_bj_ci,                   & ! c_bjc_i
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     rho_a_i(batch_a%first, 1), &
                     wf%n_v)
!
         call mem%dealloc(L_ab_jc, (batch_a%length)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call mem%dealloc(c_bj_ci, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      enddo ! End batching over a
!
   end subroutine jacobian_ccsd_d1_ccsd
!
!
   module subroutine jacobian_ccsd_a2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^A2 = sum_c g_aibc c_cj - sum_k g_aikj c_bk
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^A2,
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_ai_kj 
      real(dp), dimension(:,:), allocatable :: g_ai_jk 
      real(dp), dimension(:,:), allocatable :: g_ai_bc 
!
      real(dp), dimension(:,:), allocatable :: rho_ba_ij 
!
!     Batching variables
!
      integer(i15) :: required = 0
      integer(i15) :: current_b_batch = 0
      integer(i15) :: aib_offset = 0
!
      type(batching_index) :: batch_b
!
!     :: Term 1. - sum_k g_aikj c_bk ::
!
!     Calculate g_ai_kj
!
      call mem%alloc(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call wf%get_vooo(g_ai_kj)
!
!     Reorder to g_ai_jk = g_aikj = g_ai_kj
!
      call mem%alloc(g_ai_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call sort_1234_to_1243(g_ai_kj, g_ai_jk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call mem%alloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                &
                  wf%n_v,               &
                  g_ai_jk,              & ! g_aij_k
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  rho_ba_ij,            & ! rho_b_aij
                  wf%n_v)
!  
      call mem%dealloc(g_ai_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call add_3124_to_1234(one, rho_ba_ij, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     :: Term 2. rho_ai_bj =+ sum_c g_aibc c_cj ::
!
!     We do the matrix multiplication as g_aib_c c_cj,
!     batching over the b index.
!
      required = wf%integrals%get_required_vvvo() + (wf%n_v)*(wf%n_o)*(wf%n_v)*(batch_b%length)
!
!     Initialize batching variable
!
      call batch_b%init(wf%n_v)
      call mem%num_batch(batch_b, required)
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Determine limits for the current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
         call mem%alloc(g_ai_bc, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_b%length))
!
         call wf%get_vovv(g_ai_bc,       &
                          1,             &
                          wf%n_v,        &
                          1,             &
                          wf%n_o,        &
                          batch_b%first, &
                          batch_b%last,  &
                          1,             &
                          wf%n_v)
!
!        Calculate the contribution to rho_aib_j = sum_c g_aib_c c_cj
!
         aib_offset = index_three(1, 1, batch_b%first, wf%n_v, wf%n_o)
!
         call dgemm('N', 'N',                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_ai_bc,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_a_i,                              &
                     wf%n_v,                             &
                     one,                                &
                     rho_ai_bj(aib_offset,1),            &
                     (wf%n_o)*(wf%n_v)**2)
!
          call mem%dealloc(g_ai_bc, (wf%n_v)*(wf%n_o), wf%n_v*(batch_b%length))
!
      enddo ! End of batches over b
!
   end subroutine jacobian_ccsd_a2_ccsd
!
!
   module subroutine jacobian_ccsd_b2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD B2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai_bj^B2 = - sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^B2,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: t_ai_cj ! t_ij^ac
      real(dp), dimension(:,:), allocatable :: t_ca_ij ! t_ij^ac
!
      real(dp), dimension(:,:), allocatable :: X_k_aij ! An intermediate
      real(dp), dimension(:,:), allocatable :: X_k_j   ! An intermediate
!
      real(dp), dimension(:,:), allocatable :: rho_ba_ij ! rho_ai_bj, reordered, term 1
!
!     :: Term 1. - sum_kc F_kc t_ij^ac c_bk ::
!
!     Allocate & read the amplitudes from disk
!
      !call wf%read_double_amplitudes
!
!     Order the amplitudes as t_ca_ij = t_ij^ac
!
      call mem%alloc(t_ai_cj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ai_cj = zero
!
      call squareup(wf%t2, t_ai_cj, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     t_ai_cj to t_ca_ij
!
      call mem%alloc(t_ca_ij, (wf%n_v)**2, (wf%n_o)**2)
      t_ca_ij = zero
!
      call sort_1234_to_3124(t_ai_cj, t_ca_ij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_k_aij = sum_c F_k_c t_c_aij
!
      call mem%alloc(X_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  wf%fock_ia,           &
                  wf%n_o,               &
                  t_ca_ij,              & ! t_c_aij
                  wf%n_v,               &
                  zero,                 &
                  X_k_aij,              &
                  wf%n_o)
!
      call mem%dealloc(t_ca_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     Form rho_b_aij = sum_k c_a_i(b,k) X_k_aij(k,aij)
!
      call mem%alloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                &
                  wf%n_v,               &
                  X_k_aij,              &
                  wf%n_o,               &
                  zero,                 &
                  rho_ba_ij,            & ! rho_b_aij
                  wf%n_v)
!
      call mem%dealloc(X_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
!     Add rho_ba_ij(ba,ij) to rho_ai_bj(ai,bj)
!                   3124                1234
!
      call add_3124_to_1234(one, rho_ba_ij, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     :: Term 2. - sum_kc F_kc t_ik^ab c_cj ::
!
!     Form X_k_j = sum_c F_kc c_cj = sum_c fock_ia(k,c) c_a_i(c,j)
!
      call mem%alloc(X_k_j, wf%n_o, wf%n_o)
!
      call dgemm('N','N',     &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ia, &
                  wf%n_o,     &
                  c_a_i,      &
                  wf%n_v,     &
                  zero,       &
                  X_k_j,      &
                  wf%n_o)
!
!     Form rho_aib_j = - sum_k t_aib_k X_k_j
!     (Interpret rho_ai_bj as rho_aib_j)
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_ai_cj,              & ! t_aib_k
                  (wf%n_o)*(wf%n_v)**2, &
                  X_k_j,                &
                  wf%n_o,               &
                  one,                  &
                  rho_ai_bj,            & ! rho_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_k_j, wf%n_o, wf%n_o)
      call mem%dealloc(t_ai_cj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_b2_ccsd
!
!
   module subroutine jacobian_ccsd_c2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai_bj^C2 = sum_kcl g_ljkc (t_ki^ac c_bl + t_li^bc c_ak + t_lk^ba c_ci)
!!                 - sum_kcl L_ljkc (t_il^ab c_ck + t_ik^ac c_bl)
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_lj_kc ! g_ljkc
      real(dp), dimension(:,:), allocatable :: g_kj_lc ! g_ljkc
      real(dp), dimension(:,:), allocatable :: L_lj_ck ! L_ljkc
!
      real(dp), dimension(:,:), allocatable :: X_lj_ai ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_kj_bi ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_kj_li ! An intermediate, term 3
      real(dp), dimension(:,:), allocatable :: X_kl_ij ! X_kj_li reordered
      real(dp), dimension(:,:), allocatable :: X_lj    ! An intermediate, term 4
      real(dp), dimension(:,:), allocatable :: Y_lj_ai ! An intermediate, term 5
!
      real(dp), dimension(:,:), allocatable :: t_ak_ci ! t_ki^ac
      real(dp), dimension(:,:), allocatable :: t_kc_ai ! t_ki^ac
      real(dp), dimension(:,:), allocatable :: t_lc_bi ! t_li^bc
      real(dp), dimension(:,:), allocatable :: t_ba_kl ! t_lk^ba
      real(dp), dimension(:,:), allocatable :: t_ai_bl ! t_il^ab
      real(dp), dimension(:,:), allocatable :: t_ck_ai ! t_ik^ac
!
      real(dp), dimension(:,:), allocatable :: rho_bj_ai ! rho_ai_bj, term 1 & 5
      real(dp), dimension(:,:), allocatable :: rho_aj_bi ! rho_ai_bj, term 2
      real(dp), dimension(:,:), allocatable :: rho_ba_ij ! rho_ai_bj, term 3
!
!     :: Term 1. sum_kcl g_ljkc t_ki^ac c_bl ::
!
!     Form g_lj_kc = g_ljkc
!
      call mem%alloc(g_lj_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call wf%get_ooov(g_lj_kc)
!
!     Read the amplitudes from disk
!
      !call wf%read_double_amplitudes
!
!     Square up (t_ak_ci = t_ki^ac)
!
      call mem%alloc(t_ak_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call squareup(wf%t2, t_ak_ci, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Order as t_kc_ai = t_ki^ac
!
      call mem%alloc(t_kc_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ak_ci to t_kc_ai
!               1234       2314
!
      call sort_1234_to_2314(t_ak_ci, t_kc_ai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ak_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form X_lj_ai = sum_ck g_lj_kc t_kc_ai
!
      call mem%alloc(X_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_lj_kc,           &
                  (wf%n_o)**2,       &
                  t_kc_ai,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lj_ai,           &
                  (wf%n_o)**2)
!
!     Calculate rho_b_jai = sum_l c_bl X_lj_ai
!     (Interpret the X array as an X_l_jai object in the matrix multiplication)
!
      call mem%alloc(rho_bj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_a_i,                &
                  wf%n_v,               &
                  X_lj_ai,              & ! X_l_jai
                  wf%n_o,               &
                  zero,                 &
                  rho_bj_ai,            & ! rho_b_jai
                  wf%n_v)
!
      call mem%dealloc(X_lj_ai, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Add rho_bj_ai (3412) to rho_ai_bj (1234)
!
      call add_3412_to_1234(one, rho_bj_ai, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Final deallocations for term 1 (keep g_lj_kc for later use)
!
      call mem%dealloc(rho_bj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 2. sum_kcl g_ljkc t_li^bc c_ak ::
!
!     Reorder to g_kj_lc = g_lj_kc = g_ljkc
!                  3214      1234
!
      call mem%alloc(g_kj_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call sort_1234_to_3214(g_lj_kc, g_kj_lc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_lj_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_kj_bi = sum_lc g_kj_lc t_lc_bi
!
      call mem%alloc(X_kj_bi, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'N',           &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_kj_lc,           &
                  (wf%n_o)**2,       &
                  t_kc_ai,           & ! t_lc_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_kj_bi,           &
                  (wf%n_o)**2)
!
!     Calculate rho_a_jbi = sum_k c_ak X_kj_bi
!
      call mem%alloc(rho_aj_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_a_i,                &
                  wf%n_v,               &
                  X_kj_bi,              & ! X_k_jbi
                  wf%n_o,               &
                  zero,                 &
                  rho_aj_bi,            & ! rho_a_jbi
                  wf%n_v)
!
      call mem%dealloc(X_kj_bi, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Add rho_aj_bi to rho_ai_bj
!             1432         1234
!
      call add_1432_to_1234(one, rho_aj_bi, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Deallocations for term 2 (keep g_kj_lc = g_ljkc)
!
      call mem%dealloc(rho_aj_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 3. sum_kcl g_ljkc t_lk^ba c_ci ::
!
!     Form the intermediate X_kjl_i = sum_c g_ljkc c_ci = sum_c g_kj_lc c_c_i
!
!     Note: interpret g_kj_lc as g_kjl_c in matrix multiplication.
!
      call mem%alloc(X_kj_li, (wf%n_o)**2, (wf%n_o)**2)
!
      call dgemm('N','N',      &
                  (wf%n_o)**3, &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  g_kj_lc,     & ! g_kjl_c
                  (wf%n_o)**3, &
                  c_a_i,       &
                  wf%n_v,      &
                  zero,        &
                  X_kj_li,     & ! X_kjl_i
                  (wf%n_o)**3)
!
!     Reorder to X_kl_ij = X_kj_li
!                  1342      1234
!
      call mem%alloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
      call sort_1234_to_1342(X_kj_li, X_kl_ij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_kj_li, (wf%n_o)**2, (wf%n_o)**2)
!
!     Order amplitudes as t_ba_kl = t_lk^ba
!
      call mem%alloc(t_ba_kl, (wf%n_v)**2, (wf%n_o)**2)
!
!     t_kc_ai = t_ki^ac => t_kc_ai(la,bk) = t_lk^ba = t_ba_kl(ba,kl)
!                                  1234                       3241
!
      call sort_1234_to_3241(t_kc_ai, t_ba_kl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_kc_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate rho_ba_ij = sum_kcl g_ljkc t_lk^ba c_ci
!                         = sum kl ( sum_c g_ljkc c_ci ) t_lk^ba
!                         = sum_kl t_ba_kl X_kl_ij
!
      call mem%alloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_ba_kl,     &
                  (wf%n_v)**2, &
                  X_kl_ij,     &
                  (wf%n_o)**2, &
                  zero,        &
                  rho_ba_ij,   &
                  (wf%n_v)**2)
!
      call mem%dealloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Add rho_ba_ij into rho_ai_bj
!             3124           1234
!
      call add_3124_to_1234(one, rho_ba_ij, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Deallocations for term 3 (keep g_kj_lc = g_ljkc)
!
      call mem%dealloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     :: Term 4. - sum_kcl L_ljkc t_il^ab c_ck ::
!
!     Form L_lj_ck(lj,ck) = L_ljkc = 2 * g_ljkc - g_lckj
!                  1234   = 2 * g_ljkc - g_kjlc = 2* g_kj_lc(kj,lc) - g_kj_lc(lj,kc)
!                                                            4213             1243
!
      call mem%alloc(L_lj_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
      L_lj_ck = zero
!
      call add_4213_to_1234(two, g_kj_lc, L_lj_ck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_1243_to_1234(-one, g_kj_lc, L_lj_ck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kj_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Calculate the intermediate X_lj = sum_ck L_lj_ck c_ck
!
      call mem%alloc(X_lj, (wf%n_o)**2, 1)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)**2,       &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lj_ck,           &
                  (wf%n_o)**2,       &
                  c_a_i,             & ! c_ck
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lj,              &
                  (wf%n_o)**2)
!
!     Order the amplitudes as t_ai_bl = t_il^ab
!
      call mem%alloc(t_ai_bl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     t_lk^ba = t_ba_kl(ba,kl) => t_il^ab = t_ba_kl(ab,li) = t_ai_bl(ai,bl)
!                                                   1234             1423
!
      call sort_1234_to_1423(t_ba_kl, t_ai_bl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(t_ba_kl, (wf%n_v)**2, (wf%n_o)**2)
!
!     Form rho_ai_bj =+ - sum_l t_il^ab X_lj = - sum_l t_aib_l X_lj
!
!     Note: interpret rho_ai_bj as rho_aib_j in the
!           matrix multiplication, and X_lj as X_l_j)
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_ai_bl,              & ! t_aib_l
                  (wf%n_o)*(wf%n_v)**2, &
                  X_lj,                 & ! X_l_j
                  wf%n_o,               &
                  one,                  &
                  rho_ai_bj,            & ! rho_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_lj, (wf%n_o)**2, 1)
!
!     :: Term 5. - sum_kcl L_ljkc t_ik^ac c_bl ::
!
!     t_il^ab = t_ai_bl(ai,bl) => t_ai_bl(ck,ai) = t_ki^ca = t_ik^ac
!
!     Form the intermediate Y_lj_ai = sum_kc L_ljkc t_ik^ac = sum_kc L_lj_ck t_ck_ai
!
      call mem%alloc(Y_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lj_ck,           &
                  (wf%n_o)**2,       &
                  t_ai_bl,           & ! t_ck_ai
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_lj_ai,           &
                  (wf%n_o)**2)
!
      call mem%dealloc(t_ai_bl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate rho_b_jai =+ - sum_l c_bl Y_lj_ai
!
!     Note: interpret Y_lj_ai as Y_l_jai in the matrix multiplication
!
      call mem%alloc(rho_bj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                &
                  wf%n_v,               &
                  Y_lj_ai,              & ! Y_l_jai
                  wf%n_o,               &
                  zero,                 &
                  rho_bj_ai,            & ! rho_b_jai
                  wf%n_v)
!
!     Add rho_bj_ai to rho_ai_bj
!
      call add_3412_to_1234(one, rho_bj_ai, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Deallocations and cleanup
!
      call mem%dealloc(Y_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      call mem%dealloc(L_lj_ck, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      call mem%dealloc(rho_bj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
   end subroutine jacobian_ccsd_c2_ccsd
!
!
   end submodule jacobian_ccsd
