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
   module subroutine jacobi_transform_trial_vector_ccsd(wf, c_i)
!!
!!    Jacobi transform trial vector 
!!    Written by Sarai D. Folkestad, Sep 2018
!!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c_i
!
      call wf%jacobian_ccsd_transformation(c_i)
!
   end subroutine jacobi_transform_trial_vector_ccsd
!
!
   module subroutine jacobian_ccsd_transformation_ccsd(wf, c)
!!
!!    Jacobian transformation (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_a_i = rho_a_i,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c
!
      real(dp), dimension(:,:), allocatable :: c_a_i
      real(dp), dimension(:,:), allocatable :: c_ai_bj, c_ab_ij 
!
      real(dp), dimension(:,:), allocatable :: rho_a_i   
      real(dp), dimension(:,:), allocatable :: rho_ai_bj, rho_ab_ij 
!
      integer(i15) :: i, j, a, b, ai, bj, aibj ! Index
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
      call mem%alloc(c_a_i, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_o*wf%n_v, c, 1, c_a_i, 1)
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_ccsd_a1(rho_a_i, c_a_i)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
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
                     c_ai_bj(ai, bj) = c(wf%n_o*wf%n_v + aibj, 1)
                     c_ai_bj(bj, ai) = c(wf%n_o*wf%n_v + aibj, 1)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!!$omp parallel do schedule(static) private(ai) 
      do ai = 1, (wf%n_o)*(wf%n_v)
!
         c_ai_bj(ai,ai) = two*c_ai_bj(ai,ai)
!
      enddo
!!$omp end parallel do
!
      call wf%jacobian_ccsd_b1(rho_a_i, c_ai_bj)
      call wf%jacobian_ccsd_c1(rho_a_i, c_ai_bj)
      call wf%jacobian_ccsd_d1(rho_a_i, c_ai_bj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
      call dcopy((wf%n_o)*(wf%n_v), rho_a_i, 1, c, 1)

      call mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
!     :: CCSD contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      rho_ai_bj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_ccsd_a2(rho_ai_bj, c_a_i)
      call wf%jacobian_ccsd_b2(rho_ai_bj, c_a_i)
      call wf%jacobian_ccsd_c2(rho_ai_bj, c_a_i)
      call wf%jacobian_ccsd_d2(rho_ai_bj, c_a_i)
!
      call mem%dealloc(c_a_i, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_f2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_g2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_h2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_i2(rho_ai_bj, c_ai_bj)
!
!     Last two terms are already symmetric (J2 and K2). Perform the symmetrization
!     rho_ai_bj = P_ij^ab rho_ai_bj now, for convenience
!
      call symmetric_sum(rho_ai_bj, (wf%n_v)*(wf%n_o))
!
!     In preparation for last two terms, reorder
!     rho_ai_bj to rho_ab_ij, and c_ai_bj to c_ab_ij
!
      call mem%alloc(rho_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
      call mem%alloc(c_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call sort_1234_to_1324(c_ai_bj, c_ab_ij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(rho_ai_bj, rho_ab_ij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%jacobian_ccsd_j2(rho_ab_ij, c_ab_ij)
      call wf%jacobian_ccsd_k2(rho_ab_ij, c_ab_ij)
!
!     Done with reordered doubles c; deallocate
!
      call mem%dealloc(c_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     Order rho_ab_ij back into rho_ai_bj & divide by
!     the biorthonormal factor 1 + delta_ai,bj
!
      call mem%alloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call sort_1234_to_1324(rho_ab_ij, rho_ai_bj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!!$omp parallel do schedule(static) private(ai) 
      do ai = 1, (wf%n_o)*(wf%n_v)
!
         rho_ai_bj(ai,ai) = half*rho_ai_bj(ai,ai)
!
      enddo
!!$omp end parallel do
!
      call mem%dealloc(rho_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     Overwrite the incoming doubles c vector & pack in
!
!
!!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
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
                     c(wf%n_o*wf%n_v + aibj, 1) = rho_ai_bj(ai, bj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
!     Remaining deallocations
!
      call mem%dealloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_transformation_ccsd
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
  module subroutine jacobian_ccsd_d2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai_bj^D2 = - sum_kcd g_kcbd (t_ij^cd c_ak + t_kj^ad c_ci + t_ik^ca c_dj)
!!                       + sum_kcd L_kcbd (t_ik^ac c_dj + t_ij^ad c_ck)
!!
!!    Note: the code is structured so that we batch over the index b,
!!          where the integrals are made as g_kc_db = g_kcbd and held
!!          in some ordering or other throughout a given batch (i.e.,
!!          all five terms are constructed gradually in the batches).
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^D2,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_bd_kc ! g_kcbd
      real(dp), dimension(:,:), allocatable :: g_cd_kb ! g_kcbd reordered
      real(dp), dimension(:,:), allocatable :: g_ck_bd ! g_kcbd reordered
      real(dp), dimension(:,:), allocatable :: L_ck_bd ! L_kcbd = 2 g_kcbd - g-kdbc
!
      real(dp), dimension(:,:), allocatable :: t_ij_cd ! t_ij^cd
      real(dp), dimension(:,:), allocatable :: t_dk_aj ! t_kj^ad
      real(dp), dimension(:,:), allocatable :: t_ai_ck ! t_ik^ca
      real(dp), dimension(:,:), allocatable :: t_ai_jd ! t_ij^ad
!
      real(dp), dimension(:,:), allocatable :: X_ij_kb ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_k_ijb ! The above intermediate, reordered
      real(dp), dimension(:,:), allocatable :: X_id_kb ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_ib_dk ! The above intermediate, reordered
      real(dp), dimension(:,:), allocatable :: X_ckb_j ! An intermediate, term 3
      real(dp), dimension(:,:), allocatable :: Y_ckb_j ! An intermediate, term 4
      real(dp), dimension(:,:), allocatable :: X_bd    ! An intermediate, term 5
!
      real(dp), dimension(:,:), allocatable :: rho_a_ijb ! rho_ai_bj, batching over b
      real(dp), dimension(:,:), allocatable :: rho_ib_aj ! rho_ai_bj, batching over b
      real(dp), dimension(:,:), allocatable :: rho_aib_j ! rho_ai_bj, batching over b
      real(dp), dimension(:,:), allocatable :: rho_b_aij ! rho_ai_bj, batching over b
!
      logical :: reorder
!
!     Batching variables
!
      integer(i15) :: required = 0
      integer(i15) :: current_b_batch = 0
!
      type(batching_index) :: batch_b
!
!     Indices
!
      integer(i15) :: b = 0, c = 0, cd = 0, ci = 0, dj = 0, cidj = 0, d = 0, db = 0, bd = 0
      integer(i15) :: k = 0, j = 0, kb = 0, kc = 0, i = 0, ij = 0, ijb = 0
      integer(i15) :: a = 0, ai = 0, bj = 0, ib = 0, dkb = 0, dk = 0, akdj = 0, ak = 0
      integer(i15) :: aj = 0, ck = 0, ckb = 0, ciak = 0, aib = 0, aidj = 0, aij = 0
!
!     Read amplitudes from disk
!
      !call wf%read_double_amplitudes
!
!     Determine batch size, etc.
!     (Redo estimate once loop is done)
!
      required = wf%integrals%get_required_vvov() + &
                  (max((wf%n_o)*(wf%n_v**3), &
                  (wf%n_o**2)*(wf%n_v**2) + (wf%n_o**3)*(wf%n_v), &
                  2*(wf%n_o**3)*(wf%n_v), &
                  3*(wf%n_o**2)*(wf%n_v**2)))
!
!     Initialize batching variable
!
      call batch_b%init(wf%n_v)
      call mem%num_batch(batch_b, required)
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
         call mem%alloc(g_bd_kc, (wf%n_v)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
         call wf%get_vvov(g_bd_kc,        &
                           batch_b%first, &
                           batch_b%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Reorder g_bd_kc to g_cd_kb (= g_kcbd), i.e. 1234 to 4231
!
         call mem%alloc(g_cd_kb, (wf%n_v)**2, (wf%n_o)*(batch_b%length))
!
         call sort_1234_to_4231(g_bd_kc, g_cd_kb, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_bd_kc, (wf%n_v)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
!        Order amplitudes as t_ij_cd = t_ij^cd = t_ci_dj
!                              2413                1234
!
         call mem%alloc(t_ij_cd, (wf%n_o)**2, (wf%n_v)**2)
         t_ij_cd = zero
!
         call squareup_and_sort_1234_to_2413(wf%t2, t_ij_cd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form intermediate X_ij_kb = sum_cd g_kcdb t_ij^cd
!                                  = sum_cd t_ij_cd g_cd_kb
!
         call mem%alloc(X_ij_kb, (wf%n_o)**2, (wf%n_o)*(batch_b%length))
!
         call dgemm('N', 'N',                   &
                     (wf%n_o)**2,               &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_v)**2,               &
                     one,                       &
                     t_ij_cd,                   &
                     (wf%n_o)**2,               &
                     g_cd_kb,                   &
                     (wf%n_v)**2,               &
                     zero,                      &
                     X_ij_kb,                   &
                     (wf%n_o)**2)
!
         call mem%dealloc(t_ij_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!        sum_kcd g_kcbd t_ij^cd c_ak = sum_k X_ij_kb c_ak
!        Reorder to X_k_ijb = X_ij_kb
!
         call mem%alloc(X_k_ijb, (wf%n_o), (batch_b%length)*(wf%n_o)**2)
         X_k_ijb = zero
!
         do b = 1, batch_b%length
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  ij = index_two(i, j, wf%n_o)
!
                  ijb = index_three(i, j, b, wf%n_o, wf%n_o)
!
                  do k = 1, wf%n_o
!
                     kb = index_two(k, b, wf%n_o)
!
                     X_k_ijb(k, ijb) = X_ij_kb(ij, kb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(X_ij_kb, (wf%n_o)**2, (wf%n_o)*(batch_b%length))
!
!        Form rho_a_ijb = - sum_k c_ak X_k_ijb = - sum_k c_a_i(a,k) X_k_ijb(k, ijb)
!
         call mem%alloc(rho_a_ijb, wf%n_v, (batch_b%length)*(wf%n_o)**2)
!
         call dgemm('N', 'N',                      &
                     wf%n_v,                       &
                     (batch_b%length)*(wf%n_o)**2, &
                     wf%n_o,                       &
                     -one,                         &
                     c_a_i,                        &
                     wf%n_v,                       &
                     X_k_ijb,                      &
                     wf%n_o,                       &
                     zero,                         &
                     rho_a_ijb,                    &
                     wf%n_v)
!
         call mem%dealloc(X_k_ijb, wf%n_o, (batch_b%length)*(wf%n_o)**2)
!
!        Add rho_a_ijb (batch over b) to rho_ai_bj (full space)
!
         do b = 1, batch_b%length ! Loop over restricted space
            do j = 1, wf%n_o
!
               Bj = index_two(b + batch_b%first - 1, j, wf%n_v) ! b in full space
!
               do i = 1, wf%n_o
!
                  ijb = index_three(i, j, b, wf%n_o, wf%n_o)
!
                  do a = 1, wf%n_v
!
                     ai = index_two(a, i, wf%n_v)
!
                     rho_ai_bj(ai,Bj) = rho_ai_bj(ai,Bj) + rho_a_ijb(a, ijb)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocations for term 1 (keep g_cd_kb = g_kcbd)
!
         call mem%dealloc(rho_a_ijb, wf%n_v, (batch_b%length)*(wf%n_o)**2)
!
!
!        :: Term 2. - sum_kcd g_kcbd t_kj^ad c_ci ::
!
!        Form the intermediate X_i_dkb = sum_c g_kcbd c_ci
!                                      = sum_c c_ci g_cd_kb
!                                      = sum_c c_a_i^T(i,c) g_cd_kb(c, dkb)
!
!        Note: g_cd_kb is interpreted as g_c_dkb in the matrix multiplication.
!
         call mem%alloc(X_id_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(batch_b%length))
!
         call dgemm('T','N',                             &
                     wf%n_o,                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_v,                             &
                     one,                                &
                     c_a_i,                              &
                     wf%n_v,                             &
                     g_cd_kb,                            & ! "g_c_dkb"
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
         call mem%alloc(X_ib_dk, (wf%n_o)*(batch_b%length), (wf%n_v)*(wf%n_o))
!
         call sort_1234_to_1423(X_id_kb, X_ib_dk, wf%n_o, wf%n_v, wf%n_o, batch_b%length)
!
         call mem%dealloc(X_id_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(batch_b%length))
!
!        Order the amplitudes as t_dk_aj = t_kj^ad = t_jk^da (t_dj, ak)
!
         call mem%alloc(t_dk_aj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_dk_aj = zero
!
         call squareup_and_sort_1234_to_1432(wf%t2, t_dk_aj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Calculate rho_ib_aj = - sum_kcd g_kcbd t_kj^ad c_ci
!                            = - sum_dk X_ib_dk t_dk_aj
!
         call mem%alloc(rho_ib_aj, (wf%n_o)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_o)*(wf%n_v),         &
                     (wf%n_o)*(wf%n_v),         &
                     -one,                      &
                     X_ib_dk,                   &
                     (wf%n_o)*(batch_b%length), &
                     t_dk_aj,                   &
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     rho_ib_aj,                 &
                     (wf%n_o)*(batch_b%length))
!
         call mem%dealloc(X_ib_dk, (wf%n_o)*(batch_b%length), (wf%n_o)*(wf%n_v))
         call mem%dealloc(t_dk_aj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Add rho_ib_aj (batch over b) ro rho_ai_bj (full space)
!
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               aj = index_two(a, j, wf%n_v)
!
               do b = 1, batch_b%length
!
                  Bj = index_two(b + batch_b%first - 1, j, wf%n_v) ! b is full space index
!
                  do i = 1, wf%n_o
!
                     ai = index_two(a, i, wf%n_v)
                     ib = index_two(i, b, wf%n_o)
!
                     rho_ai_bj(ai, Bj) = rho_ai_bj(ai, Bj) + rho_ib_aj(ib, aj)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocations for term 2 (keep g_cd_kb = g_kcbd)
!
         call mem%dealloc(rho_ib_aj, (wf%n_o)*(batch_b%length), (wf%n_v)*(wf%n_o))
!
!
!        :: Term 3. - sum_kcd g_kcbd t_ik^ca c_dj ::
!
!        sum_d g_kcbd c_dj = sum_d g_cd_kb c_dj
!
!        Reorder integrals to g_cd_kb to g_ck_bd
!
         call mem%alloc(g_ck_bd, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_b%length))
!
         call sort_1234_to_1342(g_cd_kb, g_ck_bd, wf%n_v, wf%n_v, wf%n_o, batch_b%length)
!
         call mem%dealloc(g_cd_kb, (wf%n_v)**2, (wf%n_o)*(batch_b%length))
!
!        Form the intermediate X_ckb_j = sum_d g_kcbd c_dj = sum_d g_ckb_d c_d_j
!
         call mem%alloc(X_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_ck_bd,                            & ! g_ckb_d
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_a_i,                              &
                     wf%n_v,                             &
                     zero,                               &
                     X_ckb_j,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ca = t_ki^ac (t_ak,ci)
!
         call mem%alloc(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
         t_ai_ck = zero
!
         call squareup_and_sort_1234_to_1432(wf%t2, t_ai_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form rho_aib_j = -sum_kcd g_kcbd t_ik^ca c_dj = sum_ck t_ai_ck X_ckb_j
!
!        Note: X_ckb_j is interpreted as X_ck_bj in the matrix multiplication.
!        Note: rho_aib_j is interpreted as rho_ai_bj in the matrix multiplication.
!
         call mem%alloc(rho_aib_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_v)*(wf%n_o),         &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_v)*(wf%n_o),         &
                     -one,                      &
                     t_ai_ck,                   &
                     (wf%n_v)*(wf%n_o),         &
                     X_ckb_j,                   & ! "X_ck_bj"
                     (wf%n_v)*(wf%n_o),         &
                     zero,                      &
                     rho_aib_j,                 & ! "rho_ai_bj"
                     (wf%n_v)*(wf%n_o))
!
         call mem%dealloc(X_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
         call mem%dealloc(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!        Add rho_aib_j to rho_ai_bj
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
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     rho_ai_bj(ai, Bj) = rho_ai_bj(ai, Bj) + rho_aib_j(aib, j)
!
                 enddo
               enddo
            enddo
         enddo
!
!        Deallocations for term 3 (keep g_ckb_d = g_kcbd)
!
         call mem%dealloc(rho_aib_j, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_o)
!
!        :: Term 4.  sum_kcd L_kcbd t_ik^ac c_dj ::
!
!        sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj
!
!        Form L_ckb_d = L_kcbd = 2 * g_kcbd - g_kdbc = 2 * g_ckb_d(ckb, d) - g_ckb_d(dkb, c)
!
         call mem%alloc(L_ck_bd, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_b%length))
         L_ck_bd = zero
!
!        Note: exchange g_ck_bd(dk,bc) -> 4231
!
         call add_4231_to_1234(-one, g_ck_bd, L_ck_bd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         call daxpy((wf%n_v)**2*(wf%n_o)*(batch_b%length), two, g_ck_bd, 1, L_ck_bd, 1)
!
         call mem%dealloc(g_ck_bd, (wf%n_o)*(wf%n_v), (wf%n_v)*(batch_b%length))
!
!        Form the intermediate Y_ckb_j = sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj
!
         call mem%alloc(Y_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     L_ck_bd,                            & ! L_ckb_d
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_a_i,                              &
                     wf%n_v,                             &
                     zero,                               &
                     Y_ckb_j,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ac = t2(aick, 1)
!
         call mem%alloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_ai_ck = zero
!
         call squareup(wf%t2, t_ai_ck, (wf%n_o)*(wf%n_v))
!
!        Form rho_aib_j =  sum_ck t_ai_ck Y_ckb_j
!
!        Note: we interpret Y_ckb_j as Y_ck_bj in the matrix multiplication
!        Note: we interpret rho_aib_j as rho_ai_bj in the matrix multiplication
!
         call mem%alloc(rho_aib_j, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(wf%n_v),         &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_o)*(wf%n_v),         &
                     one,                       &
                     t_ai_ck,                   &
                     (wf%n_o)*(wf%n_v),         &
                     Y_ckb_j,                   & ! "Y_ck_bj"
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     rho_aib_j,                 & ! "rho_ai_bj"
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(Y_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
         call mem%dealloc(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!        Add rho_aib_j to rho_ai_bj
!
         do j = 1, wf%n_o
            do b = 1, batch_b%length
!
               Bj = index_two(b + batch_b%first - 1, j, wf%n_v)
!
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     ai = index_two(a, i, wf%n_v)
!
                     rho_ai_bj(ai, Bj) = rho_ai_bj(ai, Bj) + rho_aib_j(aib, j)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocations for term 4 (keep L_ckb_d = L_kcbd)
!
         call mem%dealloc(rho_aib_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
!
!        :: Term 5.  sum_kcd L_kcbd t_ij^ad c_ck ::
!
!        Form the intermediate X_1,bd = sum_ck c_ck L_kcbd = sum_ck c_1,ck L_ckb_d
!
!        Note: c_a_i is interpreted as c_1,ai in the matrix multiplication
!
         call mem%alloc(X_bd, 1, (batch_b%length)*(wf%n_v))
!
         call dgemm('N','N',                    &
                     1,                         &
                     (batch_b%length)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v),         &
                     one,                       &
                     c_a_i,                     & ! "c_1,ai"
                     1,                         &
                     L_ck_bd,                   &
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     X_bd,                      &
                     1)
!
!        Order amplitudes as t_ai_jd = t_ij^ad
!
         call mem%alloc(t_ai_jd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
         call squareup_and_sort_1234_to_1243(wf%t2, t_ai_jd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form rho_b_aij =  sum_kcd L_kcbd t_ij^ad c_ck
!                       =  sum_d X_bd t_d_aij
!
!        Note: X_bd is interpreted as X_b_d in the matrix multiplication
!
         call mem%alloc(rho_b_aij, batch_b%length, (wf%n_v)*(wf%n_o)**2)
!
         call dgemm('N','T',               &
                     batch_b%length,       &
                     (wf%n_v)*(wf%n_o)**2, &
                     wf%n_v,               &
                     one,                  &
                     X_bd,                 & ! X_b_d
                     batch_b%length,       &
                     t_ai_jd,              & ! t_aij_d
                     (wf%n_v)*(wf%n_o)**2, &
                     zero,                 &
                     rho_b_aij,            &
                     batch_b%length)
!
         call mem%dealloc(X_bd, 1, (batch_b%length)*(wf%n_v))
         call mem%dealloc(t_ai_jd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!        Add rho_b_aij to rho_ai_bj
!
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  do b = 1, batch_b%length
!
                     Bj = index_two(b + batch_b%first - 1, j, wf%n_v)
!
                     rho_ai_bj(ai, Bj) = rho_ai_bj(ai, Bj) + rho_b_aij(b, aij)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Final deallocations in batching loop
!
         call mem%dealloc(rho_b_aij, batch_b%length, (wf%n_v)*(wf%n_o)**2)
         call mem%dealloc(L_ck_bd, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_b%length))
!
     enddo ! End of batches over b
!
!    Destroy amplitudes from memory
!
     !call wf%destruct_double_amplitudes
!
   end subroutine jacobian_ccsd_d2_ccsd
!

    module subroutine jacobian_ccsd_e2_ccsd(wf, rho_ai_bj, c_ai_ck)
!!
!!    Jacobian CCSD E2
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai_bj^E2 = 2 sum_dlck t_bj,dl * L_kc,ld * c_ai,ck
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_ck
!
      real(dp), dimension(:,:), allocatable :: t_dl_bj
      real(dp), dimension(:,:), allocatable :: g_kc_ld
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: X_ck_bj
!
!     Read T2 amplitudes from disk
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_dl_bj = zero
!
      call squareup(wf%t2, t_dl_bj, (wf%n_o)*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
!     Construct g_kcld
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
!     Construct L_kc,ld ordered as L_ck_dl
!
      call mem%alloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ck_dl = zero
!
!     L_ck_dl(ck, dl) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!             1234                  2143              2341
!
      call add_2341_to_1234(-one, g_kc_ld, L_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2143_to_1234(two, g_kc_ld, L_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Intermediate X_ck_bj = sum_dl L_ck_dl * t_dl_bj
!
      call mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ck_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  t_dl_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     rho_ai_bj = 2 * sum_ck c_ai_ck * X_ck_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  c_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai_bj,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD F2
!!       Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!       and Andreas Skeidsvoll, 2018
!!
!!       rho_ai_bj^F2 =   - sum_ckld t_ai,ck * L_kc,ld * c_bl,dj
!!                        - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                        - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj
!!
!!       L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kc_ld(kc,ld) - 2*g_kc_ld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_kc_ld
!
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: L_dl_ck
      real(dp), dimension(:,:), allocatable :: L_lc_kd
!
      real(dp), dimension(:,:), allocatable :: c_dl_bj
      real(dp), dimension(:,:), allocatable :: c_clk_b
      real(dp), dimension(:,:), allocatable :: c_ckd_j
!
      real(dp), dimension(:,:), allocatable :: t_ai_ck
      real(dp), dimension(:,:), allocatable :: t_ai_jd
      real(dp), dimension(:,:), allocatable :: t_aib_l
 !
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_l_j
      real(dp), dimension(:,:), allocatable :: X_ck_bj
!
      real(dp), dimension(:,:), allocatable :: rho_ai_jb
      real(dp), dimension(:,:), allocatable :: rho_aib_j
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ai = 0, bj = 0, bk = 0, bl= 0, ck = 0, cl = 0, dj = 0, dl = 0
      integer(i15) :: kc = 0, kd = 0, ld = 0, lc = 0
!
      integer(i15) :: aij = 0, aib = 0, lck = 0, ckd = 0
!
      integer(i15) :: bldj = 0, aidj = 0, bkcl = 0, aibl = 0
!
!     :: Construct L_kc,ld ::

      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
      call mem%alloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ck_dl = zero
!
!     Construct L_kc,ld ordered as L_ck_dl
!
!     L_ck_dl(ck, dl) = L_kcld = 2 * g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!                                            2143              2341
!
      call add_2341_to_1234(-one, g_kc_ld, L_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2143_to_1234(two, g_kc_ld, L_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 1: - sum_ckld t_ai,ck * L_kc,ld * c_bl,dj ::
!
!     Reorder c_bl_dj as c_dl_bj
!
      call mem%alloc(c_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      c_dl_bj = zero
!
      call sort_1234_to_3214(c_ai_bj, c_dl_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     X_ck_bj = sum_dl L_ck_dl*c_dl_bj = sum_dl L_kc,ld*c_bl,dj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ck_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  c_dl_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ai_ck = zero
!
      call squareup(wf%t2, t_ai_ck, (wf%n_o)*(wf%n_v))
      !call wf%destruct_double_amplitudes
!
!     rho_ai_bj = sum_ck t_ai_ck*X_ck_bj
!
      call dgemm('N', 'N',          &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  t_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai_bj,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 2: - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
!     Construct L_ck,dl reordered as L_dl_ck
!
      call mem%alloc(L_dl_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      L_dl_ck = zero
!
!     L_dl_ck(dl,ck) =- g_kc_ld(kd, lc) (4123)
!
      call add_4123_to_1234(-one, g_kc_ld, L_dl_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     L_dl_ck(dl,ck) =+ 2*g_kc_ld(kc, ld) (4321)
!
      call add_4321_to_1234(two, g_kc_ld, L_dl_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
                  L_dl_ck,              & ! L_d_lck
                  wf%n_v,               & ! c_b_lck
                  c_ai_bj,              &
                  wf%n_v,               &
                  zero,                 &
                  Y_d_b,                &
                  wf%n_v)
!
      call mem%dealloc(L_dl_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_ai_jd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      t_ai_jd = zero
!
!     Reorder T2 amplitudes
!
      call squareup_and_sort_1234_to_1243(wf%t2, t_ai_jd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      !call wf%destruct_double_amplitudes
!
      call mem%alloc(rho_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_aij_b = sum_d t_aij_d*Y_d_b
!
      call dgemm('N','N',                 &
                  ((wf%n_o)**2)*(wf%n_v), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  t_ai_jd,                & ! t_aij_d
                  ((wf%n_o)**2)*(wf%n_v), &
                  Y_d_b,                  &
                  wf%n_v,                 &
                  zero,                   &
                  rho_ai_jb,              & ! rho_aij_b
                  ((wf%n_o)**2)*(wf%n_v))
!
      call mem%dealloc(t_ai_jd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     Adding term 2 to rho_ai_bj
!
      call add_1243_to_1234(one, rho_ai_jb, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 3: - sum_ckdl t_ai,bl * L_kc,ld * c_ck,dj ::
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
      call mem%alloc(L_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_lc_kd = zero
!
!     Construct L_kc,ld ordered as L_lc_kd
!
!     L_lc_kd(lc,kd) = L_kcld = 2 * g_kcld - g_kdlc = 2 * g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
!     L_lc_kd(lc, kd) =- g_kc_ld(kd, lc) (3412)
!
      call add_3412_to_1234(-one, g_kc_ld, L_lc_kd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     L_lc_kd(lc, kd) =+ two * g_kc_ld(kc, ld) (3214)
!
      call add_3214_to_1234(two, g_kc_ld, L_lc_kd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
                  L_lc_kd,              & ! L_l_ckd
                  wf%n_o,               &
                  c_ai_bj,              & ! c_ai_bj(ck,dl)= c_ckd_l
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_l_j,                &
                  wf%n_o)
!
      call mem%dealloc(L_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_aib_l, (wf%n_o)*((wf%n_v)), wf%n_o*(wf%n_v))
      t_aib_l = zero
      call squareup(wf%t2, t_aib_l, wf%n_o*(wf%n_v))

!
      !call wf%destruct_double_amplitudes
!
!     rho_ai_bj = sum_l t_aib_l * Z_l_j
!
      call dgemm('N','N',                 &
                  ((wf%n_v)**2)*(wf%n_o), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  t_aib_l,                &
                  ((wf%n_v)**2)*(wf%n_o), &
                  Z_l_j,                  &
                  wf%n_o,                 &
                  one,                    &
                  rho_ai_bj,              &
                  ((wf%n_v)**2)*(wf%n_o))
!
      call mem%dealloc(t_aib_l, (wf%n_o)*((wf%n_v)), wf%n_o*(wf%n_v))
      call mem%dealloc(Z_l_j, wf%n_o, wf%n_o)
!
   end subroutine jacobian_ccsd_f2_ccsd
!
!
   module subroutine jacobian_ccsd_g2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!    Jacobian CCSD G2
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai_bj^G2 =  - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck
!!                       - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!!                       - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kc_ld(kc,ld) - 2*g_kc_ld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_kc_ld
!
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: L_d_clk
      real(dp), dimension(:,:), allocatable :: L_lc_kd
!
      real(dp), dimension(:,:), allocatable :: c_ai_ck
      real(dp), dimension(:,:), allocatable :: c_aib_l
      real(dp), dimension(:,:), allocatable :: c_ai_jd
!
      real(dp), dimension(:,:), allocatable :: t_dl_bj
      real(dp), dimension(:,:), allocatable :: t_clk_b
      real(dp), dimension(:,:), allocatable :: t_ckd_j
 !
      real(dp), dimension(:,:), allocatable :: X_ck_bj
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_l_j
!
      real(dp), dimension(:,:), allocatable :: rho_ai_jb
      real(dp), dimension(:,:), allocatable :: rho_aib_j
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ai = 0, bj = 0, bk = 0, bl = 0, ck = 0, cl = 0, dj = 0, dl = 0
      integer(i15) :: kc = 0, lc = 0, kd = 0, ld = 0
!
      integer(i15) :: aib = 0, aij = 0, ckd = 0, clk = 0
!
      integer(i15) :: ckbl = 0, ckdj = 0, bldj = 0
!
!     :: Term 1: - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck  ::
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
      call mem%alloc(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ck_dl = zero
!
!     Construct L_kc_ld = two*g_kcld - g_kdlc ordered as L_ck_dl
!
!     L_ck_dl(ck, dl) =- g_kc_ld(kd, lc) (2341)
!
      call add_2341_to_1234(-one, g_kc_ld, L_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     L_ck_dl(ck, dl) =+ two*g_kc_ld(kc, ld) (2143)
!
      call add_2143_to_1234(two, g_kc_ld, L_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_bl_dj as t_dl_bj
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_dl_bj = zero
!
      do l = 1, wf%n_o
         do j = 1, wf%n_o
            do d = 1, wf%n_v
!
               dj = index_two(d, j, wf%n_v)
               dl = index_two(d, l, wf%n_v)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
                  bl = index_two(b, l, wf%n_v)
!
                  bldj = index_packed(bl, dj)
!
                  t_dl_bj(dl, bj) = wf%t2(bldj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      !call wf%destruct_double_amplitudes
!
      call mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     X_ck_bj = sum_dl t_bl,dj * L_kc,ld = sum_dl L_ck_dl t_dl_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ck_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  t_dl_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     rho_ai_bj =+ - sum_ck c_ai,ck X_ck_bj
!
      call dgemm('N', 'N',          &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_ai_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  X_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai_bj,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 2: - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
!     Reorder L_ck_dl to L_d_clk
!
      call mem%alloc(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2))
      L_d_clk = zero
!
      do k = 1, wf%n_o
         do l = 1, wf%n_o
            do c = 1, wf%n_v
!
               clk = index_three(c, l, k, wf%n_v, wf%n_o)
!
               kc = index_two(k, c, wf%n_o)
               lc = index_two(l, c, wf%n_o)
!
               do d = 1, wf%n_v
!
                  ld = index_two(l, d, wf%n_o)
                  kd = index_two(k, d, wf%n_o)
!
                  L_d_clk(d, clk) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ck,bl as t_clk_b
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      t_clk_b = zero
!
      do k = 1, wf%n_o
         do l = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               clk = index_three(c, l, k, wf%n_v, wf%n_o)
!
               do b = 1, wf%n_v
!
                  bl = index_two(b, l, wf%n_v)
!
                  ckbl = index_packed(ck, bl)
!
                  t_clk_b(clk, b) = wf%t2(ckbl, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      !call wf%destruct_double_amplitudes
!
!     Y_d_b = sum_clk L_d_clk * c_clk_b
!
      call mem%alloc(Y_d_b, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  wf%n_v,               &
                  ((wf%n_o)**2)*wf%n_v, &
                  one,                  &
                  L_d_clk,              &
                  wf%n_v,               &
                  t_clk_b,              &
                  ((wf%n_o)**2)*wf%n_v, &
                  zero,                 &
                  Y_d_b,                &
                  wf%n_v)
!
      call mem%dealloc(t_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      call mem%dealloc(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2))
!
      call mem%alloc(c_ai_jd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      c_ai_jd = zero
!
!     Reorder c_ai_dj to c_aij_d
!
      call sort_1234_to_1243(c_ai_bj, c_ai_jd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      !call wf%destruct_double_amplitudes
!
      call mem%alloc(rho_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_aij_b = sum_d c_aij_d * Y_d_b
!
      call dgemm('N','N',                 &
                  ((wf%n_o)**2)*(wf%n_v), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  c_ai_jd,                & ! c_aij_d
                  ((wf%n_o)**2)*(wf%n_v), &
                  Y_d_b,                  &
                  wf%n_v,                 &
                  zero,                   &
                  rho_ai_jb,              & ! rho_aij_b
                  ((wf%n_o)**2)*(wf%n_v))
!
      call mem%dealloc(c_ai_jd, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     Adding term 2 to rho_ai_bj
!
      call add_1243_to_1234(one, rho_ai_jb, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: Term 3: - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl ::
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
      call mem%alloc(L_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_lc_kd = zero
!
!     Construct L_kc_ld ordered as L_lc_kd
!
!     L_lc_kd(lc,kd) =- g_kc_ld(kd,lc) (3412)
!
      call add_3412_to_1234(-one, g_kc_ld, L_lc_kd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     L_lc_kd(lc,kd) =+ two*g_kc_ld(kc,ld) (3214)
!
      call add_3214_to_1234(two, g_kc_ld, L_lc_kd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ck,dj to t_ckd_j
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_ckd_j, ((wf%n_v))*(wf%n_o), wf%n_o*(wf%n_v))
      t_ckd_j = zero
      call squareup(wf%t2, t_ckd_j, wf%n_o*(wf%n_v))
!
      !call wf%destruct_double_amplitudes
!
      call mem%alloc(Z_l_j, wf%n_o, wf%n_o)
!
!     Z_l_j = sum_ckd L_l_ckd*t_ckd_j
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  one,                  &
                  L_lc_kd,              & ! L_l_ckd
                  wf%n_o,               &
                  t_ckd_j,              &
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_l_j,                &
                  wf%n_o)
!
      call mem%dealloc(L_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(t_ckd_j, ((wf%n_v))*(wf%n_o), wf%n_o*(wf%n_v))
!
!     rho_aib_j = sum_l c_aib_l*Z_l_j
!
      call dgemm('N','N',                 &
                  ((wf%n_v)**2)*(wf%n_o), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai_bj,                &
                  ((wf%n_v)**2)*(wf%n_o), &
                  Z_l_j,                  &
                  wf%n_o,                 &
                  one,                    &
                  rho_ai_bj,              &
                  ((wf%n_v)**2)*(wf%n_o))
!
      call mem%dealloc(Z_l_j, wf%n_o, wf%n_o)
!
   end subroutine jacobian_ccsd_g2_ccsd
!
!
   module subroutine jacobian_ccsd_h2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD H2
!!       Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!       and Andreas Skeidsvoll
!!
!!       rho_ai_bj^H2 =  sum_ckdl t_ci,ak * g_kc,ld * c_bl,dj
!!                     + sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
         real(dp), dimension(:,:), allocatable :: L_ia_J
         real(dp), dimension(:,:), allocatable :: g_kc_ld
         real(dp), dimension(:,:), allocatable :: g_lc_kd
!
         real(dp), dimension(:,:), allocatable :: t_ai_kc
         real(dp), dimension(:,:), allocatable :: t_aj_lc
!
         real(dp), dimension(:,:), allocatable :: c_ld_bj
         real(dp), dimension(:,:), allocatable :: c_kd_bi
!
         real(dp), dimension(:,:), allocatable :: X_ai_ld
         real(dp), dimension(:,:), allocatable :: Y_aj_kd
!
         real(dp), dimension(:,:), allocatable :: rho_aj_bi
!
         integer(i15) :: a = 0, b = 0, c = 0, d = 0
         integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
         integer(i15) :: ak = 0, ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, bl = 0, ci = 0, cj = 0, di = 0, dj = 0
         integer(i15) :: kc = 0, kd = 0, ld = 0, lc = 0
         integer(i15) :: akci = 0, alcj = 0
!
!        :: Term 1: sum_ckld t_ci,ak * g_kc,ld * c_bl,dj ::
!
!        Construct g_kc_ld
!
         call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_ovov(g_kc_ld)
!
!        t_ak,ci ordered as t_ai_kc
!
         !call wf%read_double_amplitudes
!
         call mem%alloc(t_ai_kc, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         t_ai_kc = zero
!
         call squareup_and_sort_1234_to_1423(wf%t2, t_ai_kc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         !call wf%destruct_double_amplitudes
!
         call mem%alloc(X_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        X_ai_ld = sum_ck t_ai_kc*g_kc_ld
!
         call dgemm('N', 'N',           &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     t_ai_kc,           &
                     (wf%n_o)*(wf%n_v), &
                     g_kc_ld,           &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     X_ai_ld,           &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(t_ai_kc, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call mem%alloc(c_ld_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         c_ld_bj = zero
!
!        Reorder c_bl,dj as c_ld_bj
!
         call sort_1234_to_2314(c_ai_bj, c_ld_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        rho_ai_bj += sum_ld X_ai_ld*c_ld_bj
!
         call dgemm('N', 'N',           &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     X_ai_ld,           &
                     (wf%n_o)*(wf%n_v), &
                     c_ld_bj,           &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     rho_ai_bj,         &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(c_ld_bj, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         call mem%dealloc(X_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        :: Term 2: sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!
!        Construct g_kc_ld
!
         call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_ovov(g_kc_ld)
!
!        Reorder g_kc_ld to g_lc_kd
!
         call mem%alloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         g_lc_kd = zero
!
         call sort_1234_to_3214(g_kc_ld, g_lc_kd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        t_al,cj ordered as t_aj_lc
!
         !call wf%read_double_amplitudes
!
         call mem%alloc(t_aj_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_aj_lc = zero
!
         call squareup_and_sort_1234_to_1423(wf%t2, t_aj_lc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         !call wf%destruct_double_amplitudes
!
         call mem%alloc(Y_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Y_aj_kd = sum_lc t_aj_lc * g_lc_kd
!
         call dgemm('N', 'N',           &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     t_aj_lc,           &
                     (wf%n_o)*(wf%n_v), &
                     g_lc_kd,           &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     Y_aj_kd,           &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call mem%dealloc(t_aj_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder c_bk,di as c_kd_bi
!
         call mem%alloc(c_kd_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         c_kd_bi = zero
!
         call sort_1234_to_2314(c_ai_bj, c_kd_bi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call mem%alloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        rho_aj_bi = sum_kd  Y_aj_kd * c_kd_bi
!
         call dgemm('N', 'N',           &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     Y_aj_kd,           &
                     (wf%n_o)*(wf%n_v), &
                     c_kd_bi,           &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     rho_aj_bi,         &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(c_kd_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call mem%dealloc(Y_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder into rho_ai_bj
!
         call add_1432_to_1234(one, rho_aj_bi, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call mem%dealloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      end subroutine jacobian_ccsd_h2_ccsd
!
!
   module subroutine jacobian_ccsd_i2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!    Jacobian CCSD I2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai_bj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_jk * c_ai,bk
!!                   + sum_ck L_bj,kc * c_ai,ck
!!                   - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj )
!!
!!    Batch over c to construct  g_ki_bc
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      real(dp), dimension(:,:), allocatable :: c_ai_jc
      real(dp), dimension(:,:), allocatable :: c_aib_k
      real(dp), dimension(:,:), allocatable :: c_ai_ck
      real(dp), dimension(:,:), allocatable :: c_aj_ck
!
      real(dp), dimension(:,:), allocatable :: rho_ai_jb
      real(dp), dimension(:,:), allocatable :: rho_aib_j
      real(dp), dimension(:,:), allocatable :: rho_aj_bi
!
      real(dp), dimension(:,:), allocatable :: g_bj_kc
      real(dp), dimension(:,:), allocatable :: g_bc_kj
      real(dp), dimension(:,:), allocatable :: g_ck_bj ! reordering of g_bj_kc and g_bc_kj
!
      integer(i15) :: a = 0, b = 0, c = 0
      integer(i15) :: i = 0, j = 0, k = 0
!
      integer(i15) :: ai = 0, aj = 0, ak = 0, bi = 0, bj = 0, bk = 0, ci = 0, cj = 0, ck = 0
      integer(i15) :: bc = 0
      integer(i15) :: kc = 0
      integer(i15) :: kj = 0
!
      integer(i15) :: aij = 0, aib = 0
!
!     Batching variables
!
      integer(i15) :: required = 0
      integer(i15) :: current_c_batch = 0
      integer(i15) :: offset = 0
!
      type(batching_index) :: batch_c
!
!     :: sum_c F_bc * c_ai,cj ::
!
!     Reorder c_ai,cj to c_ai_jc
!
      call mem%alloc(c_ai_jc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      c_ai_jc = zero
!
      call sort_1234_to_1243(c_ai_bj, c_ai_jc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_ai_bj += sum_c F_bc * c_ai,cj = sum_c c_aij_c(aij,c) F_ab(b,c) = sum_c c_aij_c(aij,c) F_ab^T(c,b)
!
      call dgemm('N','T',                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  one,                    &
                  c_ai_jc,                & ! c_aij_c
                  (wf%n_v)*((wf%n_o)**2), &
                  wf%fock_ab,             &
                  wf%n_v,                 &
                  zero,                   &
                  rho_ai_jb,              & ! rho_aij_b
                  (wf%n_v)*((wf%n_o)**2))
!
      call mem%dealloc(c_ai_jc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Reorder rho_aij_b into rho_ai_bj
!
      call add_1243_to_1234(one, rho_ai_jb, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     :: - sum_k F_jk * c_ai,bk  ::
!
!     rho_ai_bj += - sum_k F_jk * c_ai,bk = - sum_k c_aib_k(aib,k) F_ij(k,j)^T
!
      call dgemm('N', 'N',                &
                  (wf%n_o)*((wf%n_v)**2), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai_bj,                &
                  (wf%n_o)*((wf%n_v)**2), &
                  wf%fock_ij,             &
                  wf%n_o,                 &
                  one,                    &
                  rho_ai_bj,              &
                  (wf%n_o)*((wf%n_v)**2))
!
!     ::   sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj ) ::
!
!     sum_ck ( g_bj,kc*(2*c_ai,ck - c_ak,ci) - g_bc,kj*c_ai,ck - g_ki,bc*c_ak,cj )
!
!     Construct g_bj,kc
!
      call mem%alloc(g_bj_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_voov(g_bj_kc)
!
!     Reordering g_bj_kc to g_ck_bj
!
      call mem%alloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ck_bj = zero
!
      call sort_1234_to_4312(g_bj_kc, g_ck_bj, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_bj_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     rho_ai_bj += sum_ck 2*c_ai_ck * g_ck_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  c_ai_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai_bj,         &
                  (wf%n_o)*(wf%n_v))
!
!     Reorder c_ak,ci to c_ai_ck
!
      call mem%alloc(c_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      c_ai_ck = zero
!
      call sort_1234_to_1432(c_ai_bj, c_ai_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai_bj += - sum_ck g_ck_bj*c_ai_ck
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai_bj,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call mem%alloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ck_bj = zero
!
!     Start batching over c
!
      required = wf%integrals%get_required_vvoo()
!
!     Initialize batching variable
!
      call batch_c%init(wf%n_v)
      call mem%num_batch(batch_c, required)
!
!     Loop over the number of c batches
!
      do current_c_batch = 1, batch_c%num_batches
!
!        For each batch, get the limits for the c index
!
         call batch_c%determine_limits(current_c_batch)
!
!        Construct g_bc_kj
!
         call mem%alloc(g_bc_kj, (wf%n_v)*(batch_c%length), (wf%n_o)**2)
         g_bc_kj = zero
!
         call wf%get_vvoo(g_bc_kj,           &
                           1,             &
                           wf%n_v,        &
                           batch_c%first, &
                           batch_c%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)
!
!
!        Reorder g_bc_kj
!
         do c = 1, batch_c%length
            do b = 1, wf%n_v
 !
               bc = index_two(b, c, wf%n_v)
 !
               do j = 1, wf%n_o
 !
                  bj = index_two(b, j, wf%n_v)
 !
                  do k = 1, wf%n_o
 !
                     kj = index_two(k, j, wf%n_o)
                     ck = index_two(c + batch_c%first - 1, k, wf%n_v)
 !
                     g_ck_bj(ck, bj) = g_bc_kj(bc, kj)
 !
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(g_bc_kj, (wf%n_v)*(batch_c%length), (wf%n_o)**2)
!
      enddo ! End of c-batches
!
!     rho_ai_bj += - sum_ck c_ai_ck * g_ck_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_ai_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai_bj,         &
                  (wf%n_o)*(wf%n_v))
!
!     Reorder  c_ak,cj to c_aj_ck
!
      call mem%alloc(c_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      c_aj_ck = zero
!
      call sort_1234_to_1432(c_ai_bj, c_aj_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_aj_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  g_ck_bj,           &  ! g_ck_bi(ck,bi) = g_bc,ki
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  rho_aj_bi,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(c_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder rho_aj_bi into rho_ai_bj
!
      call add_1432_to_1234(one, rho_aj_bi, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_i2_ccsd
!
!
   module subroutine jacobian_ccsd_j2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!     Jacobian CCSD J2
!!     Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!     and Andreas Skeidsvoll, 2018
!!
!!       rho_ab_ij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2), intent(in) :: c_ab_ij
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_kc_ld
      real(dp), dimension(:,:), allocatable :: g_kl_cd
!
      real(dp), dimension(:,:), allocatable :: t_ab_ij
!
      real(dp), dimension(:,:), allocatable :: X_kl_ij
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ab = 0, cd = 0
      integer(i15) :: ai = 0, bj = 0
      integer(i15) :: kl = 0, ij = 0
      integer(i15) :: kc = 0, ld = 0
!
      integer(i15) :: aibj = 0
!
!     Constructing g_kc_ld
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
      call mem%alloc(g_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
      g_kl_cd = zero
!
!     Reorder g_kc_ld to g_kl_cd
!
      call sort_1234_to_1324(g_kc_ld, g_kl_cd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reordered T2 amplitudes
!
      !call wf%read_double_amplitudes
!
      call mem%alloc(t_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
      t_ab_ij = zero
!
      call squareup_and_sort_1234_to_1324(wf%t2, t_ab_ij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      !call wf%destruct_double_amplitudes
!
!     X_kl_ij = g_kl_cd * t_cd_ij
!
      call mem%alloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
      call dgemm('N', 'N',     &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_cd,     &
                  (wf%n_o)**2, &
                  t_ab_ij,     &
                  (wf%n_v)**2, &
                  zero,        &
                  X_kl_ij,     &
                  (wf%n_o)**2)
!
!     rho_ab_ij += c_ab_kl * X_kl_ij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  c_ab_ij,     &
                  (wf%n_v)**2, &
                  X_kl_ij,     &
                  (wf%n_o)**2, &
                  one,         &
                  rho_ab_ij,   &
                  (wf%n_v)**2)
!
!     X_kl_ij = g_kl_cd * c_cd_ij
!
      call dgemm('N', 'N',     &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_cd,     &
                  (wf%n_o)**2, &
                  c_ab_ij,     &
                  (wf%n_v)**2, &
                  zero,        &
                  X_kl_ij,     &
                  (wf%n_o)**2)
!
!     rho_ab_ij += t_ab_kl * X_kl_ij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_ab_ij,     &
                  (wf%n_v)**2, &
                  X_kl_ij,     &
                  (wf%n_o)**2, &
                  one,         &
                  rho_ab_ij,   &
                  (wf%n_v)**2)
!
      call mem%dealloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
      call mem%dealloc(g_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
      call mem%dealloc(t_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
   end subroutine jacobian_ccsd_j2_ccsd
!
!
   module subroutine jacobian_ccsd_k2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!    Jacobian CCSD K2
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ab_ij^K2 =    sum_kl g_ki,lj * c_ak,bl
!!                       + sum_cd g_ac,bd * c_ci,dj
!!
!!    For the last term we batch over a and b and
!!    add each batch to rho_ai_bj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2), intent(in) :: c_ab_ij
!
      real(dp), dimension(:,:), allocatable :: g_ki_lj
      real(dp), dimension(:,:), allocatable :: g_kl_ij
      real(dp), dimension(:,:), allocatable :: g_ac_bd
      real(dp), dimension(:,:), allocatable :: g_ab_cd
!
      real(dp), dimension(:,:), allocatable :: rho_batch_ab_ij
!
      integer(i15) :: a = 0, b = 0, i = 0, j = 0, ab = 0, full_ab = 0, ij = 0
!
!     Batching and memory handling variables
!
      integer(i15) :: required = 0
!
      integer(i15) :: current_a_batch = 0
      integer(i15) :: current_b_batch = 0
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
      call mem%alloc(g_ki_lj, (wf%n_o)**2, (wf%n_o)**2)
!
      call wf%get_oooo(g_ki_lj)
!
!     Reorder g_ki_lj to g_kl_ij
!
      call mem%alloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
      g_kl_ij = zero
!
      call sort_1234_to_1324(g_ki_lj, g_kl_ij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ki_lj, (wf%n_o)**2, (wf%n_o)**2)
!
!     rho_ab_ij += sum_kl g_ki,lj * c_ak,bl = sum_kl c_ab_ij(ab,kl) g_kl_ij(kl,ij)
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  c_ab_ij,     &
                  (wf%n_v)**2, &
                  g_kl_ij,     &
                  (wf%n_o)**2, &
                  one,         &
                  rho_ab_ij,   &
                  (wf%n_v)**2)
!
      call mem%dealloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Prepare for batching over a and b
!
!     ::  sum_cd g_ac,bd * c_ci,dj ::
!
      required = wf%integrals%get_required_vvvv() + (wf%n_v**4)
!
!     Initialize batching variables
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%num_two_batch(batch_a, batch_b, required)
!
!     Start looping over a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Determine limits for current a-batch
!
         call batch_a%determine_limits(current_a_batch)
!
!        Start looping over b-batches
!
         do current_b_batch = 1, batch_b%num_batches
!
!           Determine limits for current b-batch
!
            call batch_b%determine_limits(current_b_batch)
!
!           Allocate g_ca_db = g_acbd
!
            call mem%alloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!           g_ca_db = sum_J L_ca_J*L_db_J
!
            call wf%get_vvvv(g_ac_bd,        &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_v,        &
                           batch_b%first, &
                           batch_b%last,  &
                           1,             &
                           wf%n_v)
!
!           sum_cd g_ac,bd * c_ci,dj = sum_cd g_ac,bd c_cd,ij = sum_cd g_ab_cd c_cd_ij
!
!           Reorder g_ac_bd into g_ab_cd (i.e., 1234 to 1324)
!
            call mem%alloc(g_ab_cd, (batch_a%length)*(batch_b%length), (wf%n_v)**2)
!
            call sort_1234_to_1324(g_ac_bd, g_ab_cd, batch_a%length, wf%n_v, batch_b%length, wf%n_v)
!
            call mem%dealloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
            call mem%alloc(rho_batch_ab_ij, (batch_a%length)*(batch_b%length), (wf%n_o)**2)
!
!           rho_ab_ij += sum_cd g_ac,bd * c_ci,dj = sum_cd g_ab_cd(ab, cd) c_ab_ij(cd, ij)
!
            call dgemm('N', 'N',                            &
                        (batch_a%length)*(batch_b%length), &
                        (wf%n_o)**2,                       &
                        (wf%n_v)**2,                       &
                        one,                               &
                        g_ab_cd,                           &
                        (batch_a%length)*(batch_b%length), &
                        c_ab_ij,                           &
                        (wf%n_v)**2,                       &
                        zero,                              &
                        rho_batch_ab_ij,                   &
                        (batch_a%length)*(batch_b%length))
!
            call mem%dealloc(g_ab_cd, (batch_a%length)*(batch_b%length), (wf%n_v)**2)
!
!           Reorder into rho_ab_ij
!
            do b = 1, batch_b%length
               do a = 1, batch_a%length
!
                  ab = index_two(a, b, batch_a%length)
!
                  full_ab = index_two(a + batch_a%first - 1, b + batch_b%first - 1, wf%n_v)
!
                  do i = 1, wf%n_o
                     do j = 1, wf%n_o
!
                        ij = index_two(i, j, wf%n_o)
!
                        rho_ab_ij(full_ab, ij) = rho_ab_ij(full_ab, ij) + rho_batch_ab_ij(ab, ij)
!
                     enddo
                  enddo
               enddo
            enddo
!
            call mem%dealloc(rho_batch_ab_ij,  (batch_a%length)*(batch_b%length), (wf%n_o)**2)
!
         enddo ! End batches of b
      enddo ! End batches of a
!
   end subroutine jacobian_ccsd_k2_ccsd
!
!
   end submodule jacobian_ccsd
