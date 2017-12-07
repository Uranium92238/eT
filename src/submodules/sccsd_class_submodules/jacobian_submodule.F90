submodule (sccsd_class) jacobian
!
!!
!!    Jacobian submodule (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Contains the following family of procedures of the SCCSD class:
!!
!!    jacobian_transformation: performs the transformation by the SCCSD
!!                             Jacobian matrix A, placing the result in the
!!                             incoming vector.
!!    jacobian_sccsd_x2_sccsd: adds the contributions to the doubles transformed 
!!                             vector arising from the triples; x = a, b, c 
!!
!
   implicit none 
!
!
contains 
!
!
   module subroutine jacobian_ccsd_transformation_sccsd(wf, c_a_i, c_aibj)
!!
!!    Jacobian CCSD transformation (SCCSD)
!!    Written by Eirik F. Kjønstad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix 
!!    (with the SCCSD cluster operator, see below),
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
!!    The cluster operator is T = T_1 + T_2 + t_IJK^ABC P_IJK^ABC tau_IJK^ABC,
!!    where P_IJK^ABC performs all six permutations of the paired indices AI,
!!    BJ, and CK. The name "triples_amplitude" is given to t_IJK^ABC.
!!
      implicit none 
!
      class(sccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: c_aibj ! c_aibj 
!
      real(dp), dimension(:,:), allocatable :: rho_a_i        ! rho_ai = (A c)_ai
      real(dp), dimension(:,:), allocatable :: rho_ai_bj      ! rho_aibj = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_corr ! Holds the SCCSD correction,
                                                              ! temporarily. 
!
      real(dp), dimension(:,:), allocatable :: c_a_i_copy ! Copy of c_ai 
!
      real(dp), dimension(:,:), allocatable :: L_kc_J  ! L_kc^J
!  
      real(dp), dimension(:,:), allocatable :: g_kc_ld ! g_kcld 
      real(dp), dimension(:,:), allocatable :: L_ld_ck ! L_kcld
!
      real(dp), dimension(:,:), allocatable :: X_l_d    ! An intermediate, see below 
      real(dp), dimension(:,:), allocatable :: Y_kc_lj ! An intermediate, see below 
      real(dp), dimension(:,:), allocatable :: Z_ac_ld ! An intermediate, see below
!
      real(dp) :: norm_correction, ddot
!
      integer(i15) :: k = 0, c = 0, ck = 0, kc = 0, kd = 0, d = 0, l = 0
      integer(i15) :: lc = 0, ld = 0, i = 0, a = 0, j = 0, b = 0, ai = 0
      integer(i15) :: bj = 0
!
!     :: The CCSD contributions ::
!     ::::::::::::::::::::::::::::
!
!     Make a copy of the singles contribution (which
!     is needed for the later, SCCSD terms)
!
      call allocator(c_a_i_copy, wf%n_v, wf%n_o)
      c_a_i_copy = c_a_i
!
!     Performs CCSD Jacobian transformation (non-T3 terms)
!
      call jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
!
!     Copy the transformed vector (c_a_i, c_aibj) into (rho_a_i, rho_ai_bj)
!
      call allocator(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
      call allocator(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      rho_ai_bj = zero
!
      rho_a_i = c_a_i
      call squareup(c_aibj, rho_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Place back the copy of the singles c 
!
      c_a_i = c_a_i_copy 
      call deallocator(c_a_i_copy, wf%n_v, wf%n_o)
!
!
!     :: The SCCSD contributions ::
!     :::::::::::::::::::::::::::::
!
!     :: Make the X, Y, and Z intermediates ::
!
!        X_l_d   = sum_ck L_kcld c_ck 
!        Y_kc_lj = sum_d  g_kcld c_dj
!        Z_ac_ld = -sum_k c_ak g_kcld
!
!     - Y intermediate
!
!     Form g_kc_ld = g_kcld 
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call wf%get_cholesky_ia(L_kc_J)
!
      call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Form Y_kc_lj = sum_d g_kc_ld c_dj
!
!     Note: we interpret g_kc_ld as g_kcl_d and 
!           Y_kc_lj as Y_kcl_j
!
      call allocator(Y_kc_lj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  (wf%n_o),             &
                  wf%n_v,               &
                  one,                  &
                  g_kc_ld,              & ! "g_kcl_d"
                  (wf%n_v)*(wf%n_o)**2, &
                  c_a_i,                & ! "c_d_j"
                  wf%n_v,               &
                  zero,                 &
                  Y_kc_lj,              & ! "Y_kcl_j"
                  (wf%n_v)*(wf%n_o)**2)
!
!     - Z intermediate 
!
!     Z_ac_ld = - sum_k c_ak g_kc_ld
!
!     Note: we interpret g_kc_ld as g_k_cld and
!           Z_ac_ld as Z_a_cld 
!
      call allocator(Z_ac_ld, (wf%n_v)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                & ! "c_a_k"
                  wf%n_v,               &
                  g_kc_ld,              & ! "g_k_cld"
                  wf%n_o,               &
                  zero,                 &
                  Z_ac_ld,              & ! "Z_a_cld"
                  wf%n_v)
!
!     - X intermediate 
!
!     Form L_ld_ck = L_kcld = 2 * g_kcld - g_kdlc
!                           = 2 * g_kc_ld(kc,ld) - g_kc_ld(kd,lc)
!
      call allocator(L_ld_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ld_ck = zero 
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
            kc = index_two(k, c, wf%n_o)
!
            do d = 1, wf%n_v
!
               kd = index_two(k, d, wf%n_o)
!
               do l = 1, wf%n_o
!
                  lc = index_two(l, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
!
                  L_ld_ck(ld, ck) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc) ! L_kcld 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form X_ld = sum_ck L_ld_ck c_ck
!
!     Note: we interpret X_l_d as X_ld 
!           and c_a_i as c_ai 
!
      call allocator(X_l_d, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ld_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  c_a_i,             & ! "c_ck"
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_l_d,             & ! "X_ld"
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ld_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the contributions to the Jacobi transformation
!     resulting from, respectively, the X, Y, and Z intermediates
!
      call allocator(rho_ai_bj_corr, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      rho_ai_bj_corr = zero
!
      call wf%jacobian_sccsd_a2(rho_ai_bj_corr, X_l_d)
      call wf%jacobian_sccsd_b2(rho_ai_bj_corr, Y_kc_lj)
      call wf%jacobian_sccsd_c2(rho_ai_bj_corr, Z_ac_ld)
!
      call deallocator(Y_kc_lj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      call deallocator(Z_ac_ld, (wf%n_v)**2, (wf%n_o)*(wf%n_v))
      call deallocator(X_l_d, wf%n_o, wf%n_v)
!
!     Permute ai <-> bj, multiply by 1/(1+delta_ai,bj),
!     and add the correction to the transformed vector 
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
                  if (a .eq. b .and. i .eq. j) then
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) &
                                       + half*(rho_ai_bj_corr(ai, bj) &
                                             + rho_ai_bj_corr(bj, ai))
!
                  else
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) &
                                       + rho_ai_bj_corr(ai, bj) &
                                       + rho_ai_bj_corr(bj, ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      ! norm_correction = ddot((wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v), rho_ai_bj_corr, 1, rho_ai_bj_corr, 1)
      ! write(unit_output,*) 'Norm of correction:',norm_correction; flush(unit_output)
!
      call deallocator(rho_ai_bj_corr, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Pack the transformed vector into the incoming doubles
!
      c_aibj = zero
      call packin(c_aibj, rho_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Set the singles transformed vector 
!
      c_a_i = rho_a_i
      call deallocator(rho_a_i, wf%n_v, wf%n_o)
!
!     Final deallocations 
!
      call deallocator(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_transformation_sccsd
!
!
   module subroutine jacobian_sccsd_a2_sccsd(wf, rho, X)
!!
!!    Jacobian SCCSD A2 
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A2 term,
!!
!!       t * P_IJK^ABC (delta_aibj,AIBJ X_KC - delta_aibj,AIBK X_JC),
!!
!!    and adds it to rho(ai,bj) = rho_aibj. Capital
!!    letters denote the indices of the triples amplitude 
!!    t_IJK^ABC, and t its value.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_o, wf%n_v) :: X                         ! X(l,d) = X_ld 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho ! rho(ai,bj) = rho_aibj 
!
      integer(i15) :: AI = 0, BJ = 0, CK = 0, BK = 0, CJ = 0
      integer(i15) :: AK = 0, CI = 0, AJ = 0, BI = 0
!
!     Scale X by the triples amplitude 
!
      call dscal((wf%n_o)*(wf%n_v), wf%triples, X, 1) 
!
!     Add the twelve contributions
!
      AI = index_two(wf%A, wf%I, wf%n_v)
      BJ = index_two(wf%B, wf%J, wf%n_v)
!
      rho(AI, BJ) = rho(AI, BJ) + X(wf%K, wf%C) ! 1
      rho(BJ, AI) = rho(BJ, AI) + X(wf%K, wf%C) ! 2
!
      CK = index_two(wf%C, wf%K, wf%n_v)
!
      rho(AI, CK) = rho(AI, CK) + X(wf%J, wf%B) ! 3 
      rho(CK, AI) = rho(CK, AI) + X(wf%J, wf%B) ! 4
!
      rho(BJ, CK) = rho(BJ, CK) + X(wf%I, wf%A) ! 5
      rho(CK, BJ) = rho(CK, BJ) + X(wf%I, wf%A) ! 6
!
      BK = index_two(wf%B, wf%K, wf%n_v)
      CJ = index_two(wf%C, wf%J, wf%n_v)
!
      rho(AI, BK) = rho(AI, BK) - X(wf%J, wf%C) ! 7
      rho(AI, CJ) = rho(AI, CJ) - X(wf%K, wf%B) ! 8
!
      AK = index_two(wf%A, wf%K, wf%n_v)
      CI = index_two(wf%C, wf%I, wf%n_v)
!
      rho(BJ, AK) = rho(BJ, AK) - X(wf%I, wf%C) ! 9
      rho(BJ, CI) = rho(BJ, CI) - X(wf%K, wf%A) ! 10
!
      AJ = index_two(wf%A, wf%J, wf%n_v)
      BI = index_two(wf%B, wf%I, wf%n_v)
!
      rho(CK, AJ) = rho(CK, AJ) - X(wf%I, wf%B) ! 11
      rho(CK, BI) = rho(CK, BI) - X(wf%J, wf%A) ! 12
!
   end subroutine jacobian_sccsd_a2_sccsd  
!
!
   module subroutine jacobian_sccsd_b2_sccsd(wf, rho, Y)
!!
!!    Jacobian SCCSD B2 
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the B2 term,
!!
!!       t * P_IJK^ABC (delta_baj,ABK Y_IC_Ji
!!                +     delta_bja,AIB Y_JC_Ki
!!                - 2 * delta_bja,AIB Y_KC_Ji),
!!
!!    and adds it to rho(ai,bj) = rho_aibj. Capital
!!    letters denote the indices of the triples amplitude 
!!    t_IJK^ABC, and t its value.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)**2) :: Y ! Y(kc,lj) = Y_kclj 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho ! rho(ai,bj) = rho_aibj
!
      integer(i15) :: i = 0
!
      integer(i15) :: Bi = 0, AK = 0, IC = 0, Ji = 0
      integer(i15) :: Ci = 0, AJ = 0, IB = 0, Ki = 0
      integer(i15) :: Ai = 0, BK = 0, JC = 0, Ii = 0
      integer(i15) :: KC = 0, KB = 0, KA = 0, JB = 0
      integer(i15) :: JA = 0, IA = 0, CK = 0, CJ = 0
      integer(i15) :: BJ = 0
!
!     Scale Y by the triples amplitude 
!
      call dscal((wf%n_v)*(wf%n_o)**3, (wf%triples), Y, 1) 
!
!     Add the 18 contributions to rho(ai,bj) = rho_aibj
!
      do i = 1, wf%n_o
!
         Bi = index_two(wf%B, i, wf%n_v)
         AK = index_two(wf%A, wf%K, wf%n_v)
         IC = index_two(wf%I, wf%C, wf%n_o)
         Ji = index_two(wf%J, i, wf%n_o)
!
         rho(Bi, AK) = rho(Bi, AK) + Y(IC, Ji) ! 1
!
         Ci = index_two(wf%C, i, wf%n_v)
         AJ = index_two(wf%A, wf%J, wf%n_v)
         IB = index_two(wf%I, wf%B, wf%n_o)
         Ki = index_two(wf%K, i, wf%n_o)
!
         rho(Ci, AJ) = rho(Ci, AJ) + Y(IB, Ki) ! 2
!
         Ai = index_two(wf%A, i, wf%n_v)
         BK = index_two(wf%B, wf%K, wf%n_v)
         JC = index_two(wf%J, wf%C, wf%n_o)
         Ii = index_two(wf%I, i, wf%n_o)
!
         rho(Ai, BK) = rho(Ai, BK) + Y(JC, Ii) ! 3
!
         Ci = index_two(wf%C, i, wf%n_v)
         BI = index_two(wf%B, wf%I, wf%n_v)
         JA = index_two(wf%J, wf%A, wf%n_o)
         Ki = index_two(wf%K, i, wf%n_o)
!
         rho(Ci, BI) = rho(Ci, BI) + Y(JA, Ki) ! 4
!
         Ai = index_two(wf%A, i, wf%n_v)
         CJ = index_two(wf%C, wf%J, wf%n_v)
         KB = index_two(wf%K, wf%B, wf%n_o)
         Ii = index_two(wf%I, i, wf%n_o)
!
         rho(Ai, CJ) = rho(Ai, CJ) + Y(KB, Ii) ! 5
!
         Bi = index_two(wf%B, i, wf%n_v)
         CI = index_two(wf%C, wf%I, wf%n_v)
         KA = index_two(wf%K, wf%A, wf%n_o)
         Ji = index_two(wf%J, i, wf%n_o)
!
         rho(Bi, CI) = rho(Bi, CI) + Y(KA, Ji) ! 6
!
         Bi = index_two(wf%B, i, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         JC = index_two(wf%J, wf%C, wf%n_o)
         Ki = index_two(wf%K, i, wf%n_o)
!
         rho(Bi, AI) = rho(Bi, AI) + Y(JC, Ki) ! 7
!
         Ci = index_two(wf%C, i, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         KB = index_two(wf%K, wf%B, wf%n_o)
         Ji = index_two(wf%J, i, wf%n_o)
!
         rho(Ci, AI) = rho(Ci, AI) + Y(KB, Ji) ! 8
!
         Ai = index_two(wf%A, i, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         IC = index_two(wf%I, wf%C, wf%n_o)
         Ki = index_two(wf%K, i, wf%n_o)
!
         rho(Ai, BJ) = rho(Ai, BJ) + Y(IC, Ki) ! 9
!
         Ci = index_two(wf%C, i, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         KA = index_two(wf%K, wf%A, wf%n_o)
         Ii = index_two(wf%I, i, wf%n_o)
!
         rho(Ci, BJ) = rho(Ci, BJ) + Y(KA, Ii) ! 10
!
         Ai = index_two(wf%A, i, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         IB = index_two(wf%I, wf%B, wf%n_o)
         Ji = index_two(wf%J, i, wf%n_o)
!
         rho(Ai, CK) = rho(Ai, CK) + Y(IB, Ji) ! 11
!
         Bi = index_two(wf%B, i, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         JA = index_two(wf%J, wf%A, wf%n_o)
         Ii = index_two(wf%I, i, wf%n_o)
!
         rho(Bi, CK) = rho(Bi, CK) + Y(JA, Ii) ! 12
!
         Bi = index_two(wf%B, i, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         KC = index_two(wf%K, wf%C, wf%n_o)
         Ji = index_two(wf%J, i, wf%n_o)
!
         rho(Bi, AI) = rho(Bi, AI) - two*Y(KC, Ji) ! 13
!
         Ci = index_two(wf%C, i, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         JB = index_two(wf%J, wf%B, wf%n_o)
         Ki = index_two(wf%K, i, wf%n_o)
!
         rho(Ci, AI) = rho(Ci, AI) - two*Y(JB, Ki) ! 14
!
         Ai = index_two(wf%A, i, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         KC = index_two(wf%K, wf%C, wf%n_o)
         Ii = index_two(wf%I, i, wf%n_o)
!
         rho(Ai, BJ) = rho(Ai, BJ) - two*Y(KC, Ii) ! 15
!
         Ci = index_two(wf%C, i, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         IA = index_two(wf%I, wf%A, wf%n_o)
         Ki = index_two(wf%K, i, wf%n_o)
!
         rho(Ci, BJ) = rho(Ci, BJ) - two*Y(IA, Ki) ! 16
!
         Ai = index_two(wf%A, i, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         JB = index_two(wf%J, wf%B, wf%n_o)
         Ii = index_two(wf%I, i, wf%n_o)
!
         rho(Ai, CK) = rho(Ai, CK) - two*Y(JB, Ii) ! 17
!
         Bi = index_two(wf%B, i, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         IA = index_two(wf%I, wf%A, wf%n_o)
         Ji = index_two(wf%J, i, wf%n_o)
!
         rho(Bi, CK) = rho(Bi, CK) - two*Y(IA, Ji) ! 18
!
      enddo
!
   end subroutine jacobian_sccsd_b2_sccsd
!
!
   module subroutine jacobian_sccsd_c2_sccsd(wf, rho, Z)
!!
!!    Jacobian SCCSD C2
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the C2 term,
!!
!!       t * P_IJK^ABC (2 delta_bji,AIJ Z_aBKC 
!!                      - delta_bij,AJK Z_aBIC
!!                      - delta_bji,AIK Z_aBJC),
!!
!!    and adds it to rho(ai,bj) = rho_aibj. Capital
!!    letters denote the indices of the triples amplitude 
!!    t_IJK^ABC, and t its value.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)*(wf%n_v)) :: Z ! Z(ac,ld) = Z_acld 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho ! rho(ai,bj) = rho_aibj
!
      integer(i15) :: a = 0
!
      integer(i15) :: KC = 0, KB = 0, KA = 0, JC = 0, JB = 0, JA = 0
      integer(i15) :: IC = 0, IB = 0, IA = 0, CK = 0, CJ = 0, CI = 0
      integer(i15) :: BK = 0, BJ = 0, BI = 0, aK = 0, aJ = 0, AI = 0
      integer(i15) :: aC = 0, aB = 0, aA = 0
!
!     Scale Z by the triples amplitude
!
      call dscal((wf%n_o)*(wf%n_v)**3, (wf%triples), Z, 1) 
!
!     Add the 18 contributions to rho(ai,bj) = rho_aibj 
!
      do a = 1, wf%n_v
!
         aJ = index_two(a, wf%J, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         aB = index_two(a, wf%B, wf%n_v)
         KC = index_two(wf%K, wf%C, wf%n_o)
!
         rho(aJ, AI) = rho(aJ, AI) + two*Z(aB, KC) ! 1
!
         aK = index_two(a, wf%K, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         aC = index_two(a, wf%C, wf%n_v)
         JB = index_two(wf%J, wf%B, wf%n_o)
!
         rho(aK, AI) = rho(aK, AI) + two*Z(aC, JB) ! 2
!
         aI = index_two(a, wf%I, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         aA = index_two(a, wf%A, wf%n_v)
         KC = index_two(wf%K, wf%C, wf%n_o)
!
         rho(aI, BJ) = rho(aI, BJ) + two*Z(aA, KC) ! 3
!
         aK = index_two(a, wf%K, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         aC = index_two(a, wf%C, wf%n_v)
         IA = index_two(wf%I, wf%A, wf%n_o)
!
         rho(aK, BJ) = rho(aK, BJ) + two*Z(aC, IA) ! 4
!
         aI = index_two(a, wf%I, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         aA = index_two(a, wf%A, wf%n_v)
         JB = index_two(wf%J, wf%B, wf%n_o)
!
         rho(aI, CK) = rho(aI, CK) + two*Z(aA, JB) ! 5
!
         aJ = index_two(a, wf%J, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         aB = index_two(a, wf%B, wf%n_v)
         IA = index_two(wf%I, wf%A, wf%n_o)
!
         rho(aJ, CK) = rho(aJ, CK) + two*Z(aB, IA) ! 6
!
         aJ = index_two(a, wf%J, wf%n_v)
         AK = index_two(wf%A, wf%K, wf%n_v)
         aB = index_two(a, wf%B, wf%n_v)
         IC = index_two(wf%I, wf%C, wf%n_o)
!
         rho(aJ, AK) = rho(aJ, AK) - Z(aB, IC) ! 7
!
         aK = index_two(a, wf%K, wf%n_v)
         AJ = index_two(wf%A, wf%J, wf%n_v)
         aC = index_two(a, wf%C, wf%n_v)
         IB = index_two(wf%I, wf%B, wf%n_o)
!
         rho(aK, AJ) = rho(aK, AJ) - Z(aC, IB) ! 8
!
         aI = index_two(a, wf%I, wf%n_v)
         BK = index_two(wf%B, wf%K, wf%n_v)
         aA = index_two(a, wf%A, wf%n_v)
         JC = index_two(wf%J, wf%C, wf%n_o)
!
         rho(aI, BK) = rho(aI, BK) - Z(aA, JC) ! 9
!
         aK = index_two(a, wf%K, wf%n_v)
         BI = index_two(wf%B, wf%I, wf%n_v)
         aC = index_two(a, wf%C, wf%n_v)
         JA = index_two(wf%J, wf%A, wf%n_o)
!
         rho(aK, BI) = rho(aK, BI) - Z(aC, JA) ! 10
!
         aI = index_two(a, wf%I, wf%n_v)
         CJ = index_two(wf%C, wf%J, wf%n_v)
         aA = index_two(a, wf%A, wf%n_v)
         KB = index_two(wf%K, wf%B, wf%n_o)
!
         rho(aI, CJ) = rho(aI, CJ) - Z(aA, KB) ! 11
!
         aJ = index_two(a, wf%J, wf%n_v)
         CI = index_two(wf%C, wf%I, wf%n_v)
         aB = index_two(a, wf%B, wf%n_v)
         KA = index_two(wf%K, wf%A, wf%n_o)
!
         rho(aJ, CI) = rho(aJ, CI) - Z(aB, KA) ! 12
!
         aK = index_two(a, wf%K, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         aB = index_two(a, wf%B, wf%n_v)
         JC = index_two(wf%J, wf%C, wf%n_o)
!
         rho(aK, AI) = rho(aK, AI) - Z(aB, JC) ! 13
!
         aJ = index_two(a, wf%J, wf%n_v)
         AI = index_two(wf%A, wf%I, wf%n_v)
         aC = index_two(a, wf%C, wf%n_v)
         KB = index_two(wf%K, wf%B, wf%n_o)
!
         rho(aJ, AI) = rho(aJ, AI) - Z(aC, KB) ! 14
!
         aK = index_two(a, wf%K, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         aA = index_two(a, wf%A, wf%n_v)
         IC = index_two(wf%I, wf%C, wf%n_o)
!
         rho(aK, BJ) = rho(aK, BJ) - Z(aA, IC) ! 15
!
         aI = index_two(a, wf%I, wf%n_v)
         BJ = index_two(wf%B, wf%J, wf%n_v)
         aC = index_two(a, wf%C, wf%n_v)
         KA = index_two(wf%K, wf%A, wf%n_o)
!
         rho(aI, BJ) = rho(aI, BJ) - Z(aC, KA) ! 16
!
         aJ = index_two(a, wf%J, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         aA = index_two(a, wf%A, wf%n_v)
         IB = index_two(wf%I, wf%B, wf%n_o)
!
         rho(aJ, CK) = rho(aJ, CK) - Z(aA, IB) ! 17
!
         aI = index_two(a, wf%I, wf%n_v)
         CK = index_two(wf%C, wf%K, wf%n_v)
         aB = index_two(a, wf%B, wf%n_v)
         JA = index_two(wf%J, wf%A, wf%n_o)
!
         rho(aI, CK) = rho(aI, CK) - Z(aB, JA)
!
      enddo
!
   end subroutine jacobian_sccsd_c2_sccsd
!
!
end submodule jacobian