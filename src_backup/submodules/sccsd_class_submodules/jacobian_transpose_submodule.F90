submodule (sccsd_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Contains the following family of procedures of the SCCSD class:
!!
!!    jacobian_transpose_transformation: performs the transformation by the SCCSD
!!                                       Jacobian transpose matrix A^T, placing the result in the
!!                                       incoming vector.
!!    jacobian_sccsd_x2_sccsd:           adds the contributions to the doubles transformed 
!!                                       vector arising from the triples; x = a, b, c 
!!
!
   implicit none 
!
   character(len=40) :: integral_type
!
!
contains 
!
!
   module subroutine jacobian_transpose_ccsd_transformation_sccsd(wf, b_a_i, b_aibj)
!!
!!    Jacobian transpose CCSD transformation (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
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
!!    The cluster operator is T = T_1 + T_2 + t_IJK^ABC P_IJK^ABC tau_IJK^ABC,
!!    where P_IJK^ABC performs all six permutations of the paired indices AI,
!!    BJ, and CK. The name "triples_amplitude" is given to t_IJK^ABC.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
      real(dp), dimension(wf%n_t2am, 1)   :: b_aibj 
!
      real(dp), dimension(:,:), allocatable :: sigma_a_i   ! sigma_ai   = (A^T b)_ai
      real(dp), dimension(:,:), allocatable :: sigma_aibj  ! sigma_aibj = (A^T b)_aibj packed 
!
      real(dp), dimension(:,:), allocatable :: b_ai_bj ! b_aibj unpacked 
!
      real(dp), dimension(:,:), allocatable :: b_aibj_copy
      real(dp), dimension(:,:), allocatable :: sigma_a_i_corr ! Holds the correction to sigma_ai
!
      real(dp), dimension(:,:), allocatable :: L_ia_jb ! L_iajb 
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
!
      integer(i15) :: b = 0, j = 0, a = 0, i = 0, jb = 0, ja = 0, ib = 0, ia = 0
!
      real(dp) :: begin_timer
      real(dp) :: end_timer
!
!     :: The CCSD contributions ::
!     ::::::::::::::::::::::::::::
!
!     Make a copy of the doubles contribution
!
      call wf%mem%alloc(b_aibj_copy, wf%n_t2am, 1)
      b_aibj_copy = b_aibj
!
!     Perform the CCSD Jacobian tranpose transformation (non-T3 terms)
!
      call cpu_time(begin_timer)
      call jacobian_transpose_ccsd_transformation_ccsd(wf, b_a_i, b_aibj)
      call cpu_time(end_timer)
!
      if (wf%settings%print_level == 'developer') then 
!
         write(unit_output,'(/t6,a27,f15.8)') 'Time in CC  part (seconds):', end_timer-begin_timer
!
      endif
!
      call cpu_time(begin_timer)
!
!     Copy the transformed vector (b_a_i, b_aibj) into (sigma_a_i, sigma_ai_bj)
!
      call wf%mem%alloc(sigma_a_i, wf%n_v, wf%n_o)
      sigma_a_i = zero
!
      call wf%mem%alloc(sigma_aibj, wf%n_t2am, 1)
      sigma_aibj = zero 
!
      sigma_a_i  = b_a_i
      sigma_aibj = b_aibj
!
!     Place back the copy of the doubles b 
!
      b_aibj = b_aibj_copy
      call wf%mem%dealloc(b_aibj_copy, wf%n_t2am, 1)
!
!     :: The SCCSD contributions ::
!     :::::::::::::::::::::::::::::
!
!     For the transposed transformation, we require the integrals g_iajb, for the Y and Z terms, 
!     and L_iajb, for the X terms. 
! 
!     Form g_ia_jb 
!
      call wf%mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_ia_jb)
!
!     Square up the doubles vector 
!
      call wf%mem%alloc(b_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      b_ai_bj = zero
!
      call squareup(b_aibj, b_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Scale the integral g_ia_jb by the triples amplitude 
!
      call dscal(((wf%n_o)**2)*((wf%n_v)**2), wf%triples, g_ia_jb, 1)
!
!     Add the Y and Z terms 
!
      call wf%jacobian_transpose_sccsd_b1(sigma_a_i, b_ai_bj, g_ia_jb)
      call wf%jacobian_transpose_sccsd_c1(sigma_a_i, b_ai_bj, g_ia_jb)
!
!     Form L_ia_jb = 2 * g_ia_jb - g_ia_jb(ib,ja)
!
      call wf%mem%alloc(L_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ia_jb = zero 
!
      do b = 1, wf%n_v
         do j = 1, wf%n_o
!
            jb = index_two(j, b, wf%n_o)
!
            do a = 1, wf%n_v
!
               ja = index_two(j, a, wf%n_o)
!
               do i = 1, wf%n_o
!
                  ib = index_two(i, b, wf%n_o)
                  ia = index_two(i, a, wf%n_o)
!
                  L_ia_jb(ia, jb) = two*g_ia_jb(ia, jb) - g_ia_jb(ib, ja)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the X term 
!
      call wf%jacobian_transpose_sccsd_a1(sigma_a_i, b_ai_bj, L_ia_jb)
!
      call wf%mem%dealloc(L_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(b_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Overwrite the incoming vector with the corrected one
!
      b_a_i  = sigma_a_i
      b_aibj = sigma_aibj
!
!     Final deallocations 
!
      call wf%mem%dealloc(sigma_a_i, wf%n_v, wf%n_o)
      call wf%mem%dealloc(sigma_aibj, wf%n_t2am, 1)
!
      call cpu_time(end_timer)
!
      if (wf%settings%print_level == 'developer') then 
!
         write(unit_output,'(t6,a27,f15.8/)') 'Time in SCC part (seconds):', end_timer-begin_timer
!
      endif
!
   end subroutine jacobian_transpose_ccsd_transformation_sccsd
!
!
   module subroutine jacobian_transpose_sccsd_a1_sccsd(wf, sigma, b2am, L)
!!
!!    Jacobian transpose SCCSD A1
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       t * P_IJK^ABC (b_AIBJ L_kcKC - b_AIBK L_kcJC)
!!
!!    and adds it to sigma(c,k) = sigma_ck. Capital
!!    letters denote the indices of the triples amplitude 
!!    t_IJK^ABC, and t its value.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      !  sigma(a, i) = sigma_ai 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: L    !    L(ia, jb) = t * L_iajb 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am ! b2am(ai, bj) = b_aibj 
!
!     As an index naming convention, we distinguish 'CK', 'cK', 'Ck', and 'ck' by appending a '_two'
!     to the index appearing second in the equation (reading from left to right)
!
      integer(i15) :: k = 0, c = 0, AI = 0, BJ = 0, kc = 0, KC_two = 0
      integer(i15) :: KB = 0, KA = 0, JC = 0, JB = 0, JA = 0, IC = 0, IB = 0
      integer(i15) :: IA = 0, CK = 0, CJ = 0, CI = 0, BK = 0, BI = 0, AK = 0
      integer(i15) :: AJ = 0
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            BJ = index_two(wf%B, wf%J, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            KC_two = index_two(wf%K, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + two*b2am(AI, BJ)*L(kc, KC_two) ! 1
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            CK = index_two(wf%C, wf%K, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            JB = index_two(wf%J, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + two*b2am(AI, CK)*L(kc, JB) ! 2
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            CK = index_two(wf%C, wf%K, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            IA = index_two(wf%I, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + two*b2am(BJ, CK)*L(kc, IA) ! 3
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            BK = index_two(wf%B, wf%K, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            JC = index_two(wf%J, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - b2am(AI, BK)*L(kc, JC) ! 4
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            CJ = index_two(wf%C, wf%J, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            KB = index_two(wf%K, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - b2am(AI, CJ)*L(kc, KB) ! 5
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            AK = index_two(wf%A, wf%K, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            IC = index_two(wf%I, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - b2am(BJ, AK)*L(kc, IC) ! 6
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            CI = index_two(wf%C, wf%I, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            KA = index_two(wf%K, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - b2am(BJ, CI)*L(kc, KA) ! 7
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            AJ = index_two(wf%A, wf%J, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            IB = index_two(wf%I, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - b2am(CK, AJ)*L(kc, IB) ! 8
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            BI = index_two(wf%B, wf%I, wf%n_v)
            kc = index_two(k, c, wf%n_o)
            JA = index_two(wf%J, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - b2am(CK, BI)*L(kc, JA) ! 9
!
         enddo
      enddo
!
   end subroutine jacobian_transpose_sccsd_a1_sccsd
!
!
   module subroutine jacobian_transpose_sccsd_b1_sccsd(wf, sigma, b2am, g)
!!
!!    Jacobian transpose SCCSD B1
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the B1 term,
!!
!!       t * P_IJK^ABC (b_AJcK g_IBkC + b_AIcJ g_KBkC - 2 * b_AIcK g_JBkC)
!!
!!    and adds it to sigma(c,k) = sigma_ck. Capital
!!    letters denote the indices of the triples amplitude 
!!    t_IJK^ABC, and t its value.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      !  sigma(a, i) = sigma_ai 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: g    !    g(ia, jb) = t * g_iajb 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am ! b2am(ai, bj) = b_aibj 
!
      integer(i15) :: k = 0, c = 0, AJ = 0, cK = 0, IB = 0, kC = 0
      integer(i15) :: kB = 0, kA = 0, JC = 0, JB = 0, JA = 0, IC = 0
      integer(i15) :: IA = 0, cJ = 0, cI = 0, BK = 0, BJ = 0, BI = 0
      integer(i15) :: AK = 0, AI = 0
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            AJ = index_two(wf%A, wf%J, wf%n_v)
            cK = index_two(c, wf%K, wf%n_v)
            IB = index_two(wf%I, wf%B, wf%n_o)
            kC = index_two(k, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AJ, cK)*g(IB, kC) ! 1
!
            AK = index_two(wf%A, wf%K, wf%n_v)
            cJ = index_two(c, wf%J, wf%n_v)
            IC = index_two(wf%I, wf%C, wf%n_o)
            kB = index_two(k, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AK, cJ)*g(IC, kB) ! 2
!
            BI = index_two(wf%B, wf%I, wf%n_v)
            cK = index_two(c, wf%K, wf%n_v)
            JA = index_two(wf%J, wf%A, wf%n_o)
            kC = index_two(k, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BI, cK)*g(JA, kC) ! 3
!
            BK = index_two(wf%B, wf%K, wf%n_v)
            cI = index_two(c, wf%I, wf%n_v)
            JC = index_two(wf%J, wf%C, wf%n_o)
            kA = index_two(k, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BK, cI)*g(JC, kA) ! 4
!
            CI = index_two(wf%C, wf%I, wf%n_v)
            cJ = index_two(c, wf%J, wf%n_v)
            KA = index_two(wf%K, wf%A, wf%n_o)
            kB = index_two(k, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CI, cJ)*g(KA, kB) ! 5
!
            CJ = index_two(wf%C, wf%J, wf%n_v)
            cI = index_two(c, wf%I, wf%n_v)
            KB = index_two(wf%K, wf%B, wf%n_o)
            kA = index_two(k, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CJ, cI)*g(KB, kA) ! 6
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            cJ = index_two(c, wf%J, wf%n_v)
            KB = index_two(wf%K, wf%B, wf%n_o)
            kC = index_two(k, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AI, cJ)*g(KB, kC) ! 7
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            cK = index_two(c, wf%K, wf%n_v)
            JC = index_two(wf%J, wf%C, wf%n_o)
            kB = index_two(k, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AI, cK)*g(JC, kB) ! 8
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            cI = index_two(c, wf%I, wf%n_v)
            KA = index_two(wf%K, wf%A, wf%n_o)
            kC = index_two(k, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BJ, cI)*g(KA, kC) ! 9
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            cK = index_two(c, wf%K, wf%n_v)
            IC = index_two(wf%I, wf%C, wf%n_o)
            kA = index_two(k, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BJ, cK)*g(IC, kA) ! 10
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            cI = index_two(c, wf%I, wf%n_v)
            JA = index_two(wf%J, wf%A, wf%n_o)
            kB = index_two(k, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CK, cI)*g(JA, kB) ! 11
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            cJ = index_two(c, wf%J, wf%n_v)
            IB = index_two(wf%I, wf%B, wf%n_o)
            kA = index_two(k, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CK, cJ)*g(IB, kA) ! 12
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            cK = index_two(c, wf%K, wf%n_v)
            JB = index_two(wf%J, wf%B, wf%n_o)
            kC = index_two(k, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(AI, cK)*g(JB, kC) ! 13
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            cJ = index_two(c, wf%J, wf%n_v)
            KC = index_two(wf%K, wf%C, wf%n_o)
            kB = index_two(k, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(AI, cJ)*g(KC, kB) ! 14
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            cK = index_two(c, wf%K, wf%n_v)
            IA = index_two(wf%I, wf%A, wf%n_o)
            kC = index_two(k, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(BJ, cK)*g(IA, kC) ! 15
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            cI = index_two(c, wf%I, wf%n_v)
            KC = index_two(wf%K, wf%C, wf%n_o)
            kA = index_two(k, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(BJ, cI)*g(KC, kA) ! 16
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            cJ = index_two(c, wf%J, wf%n_v)
            IA = index_two(wf%I, wf%A, wf%n_o)
            kB = index_two(k, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(CK, cJ)*g(IA, kB) ! 17
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            cI = index_two(c, wf%I, wf%n_v)
            JB = index_two(wf%J, wf%B, wf%n_o)
            kA = index_two(k, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(CK, cI)*g(JB, kA) ! 18
!
         enddo
      enddo
!
   end subroutine jacobian_transpose_sccsd_b1_sccsd
!
!
   module subroutine jacobian_transpose_sccsd_c1_sccsd(wf, sigma, b2am, g)
!!
!!    Jacobian transpose SCCSD C1
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the C1 term,
!!
!!       t * P_IJK^ABC (b_AIBk g_KcJC + b_AKBk g_JcIC - 2 * b_AIBk g_JcKC)
!!
!!    and adds it to sigma(c,k) = sigma_ck. Capital
!!    letters denote the indices of the triples amplitude 
!!    t_IJK^ABC, and t its value.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      !  sigma(a, i) = sigma_ai 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: g    !    g(ia, jb) = t * g_iajb 
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am ! b2am(ai, bj) = b_aibj 
!
      integer(i15) :: k = 0, c = 0, Kc = 0, KB = 0, KA = 0, JC = 0, JB = 0, JA = 0
      integer(i15) :: IC = 0, IB = 0, Ck = 0, CJ = 0, CI = 0, Bk = 0, BJ = 0, BI = 0
      integer(i15) :: Ak = 0, AJ = 0, AI = 0, IA = 0
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            Bk = index_two(wf%B, k, wf%n_v)
            Kc = index_two(wf%K, c, wf%n_o)
            JC = index_two(wf%J, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AI, Bk)*g(Kc, JC) ! 1
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            Ck = index_two(wf%C, k, wf%n_v)
            Jc = index_two(wf%J, c, wf%n_o)
            KB = index_two(wf%K, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AI, Ck)*g(Jc, KB) ! 2
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            Ak = index_two(wf%A, k, wf%n_v)
            Kc = index_two(wf%K, c, wf%n_o)
            IC = index_two(wf%I, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BJ, Ak)*g(Kc, IC) ! 3
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            Ck = index_two(wf%C, k, wf%n_v)
            Ic = index_two(wf%I, c, wf%n_o)
            KA = index_two(wf%K, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BJ, Ck)*g(Ic, KA) ! 4
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            Ak = index_two(wf%A, k, wf%n_v)
            Jc = index_two(wf%J, c, wf%n_o)
            IB = index_two(wf%I, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CK, Ak)*g(Jc, IB) ! 5
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            Bk = index_two(wf%B, k, wf%n_v)
            Ic = index_two(wf%I, c, wf%n_o)
            JA = index_two(wf%J, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CK, Bk)*g(Ic, JA) ! 6
!
            AK = index_two(wf%A, wf%K, wf%n_v)
            Bk = index_two(wf%B, k, wf%n_v)
            Jc = index_two(wf%J, c, wf%n_o)
            IC = index_two(wf%I, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AK, Bk)*g(Jc, IC) ! 7
!
            AJ = index_two(wf%A, wf%J, wf%n_v)
            Ck = index_two(wf%C, k, wf%n_v)
            Kc = index_two(wf%K, c, wf%n_o)
            IB = index_two(wf%I, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(AJ, Ck)*g(Kc, IB) ! 8
!
            BK = index_two(wf%B, wf%K, wf%n_v)
            Ak = index_two(wf%A, k, wf%n_v)
            Ic = index_two(wf%I, c, wf%n_o)
            JC = index_two(wf%J, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BK, Ak)*g(Ic, JC) ! 9
!
            BI = index_two(wf%B, wf%I, wf%n_v)
            Ck = index_two(wf%C, k, wf%n_v)
            Kc = index_two(wf%K, c, wf%n_o)
            JA = index_two(wf%J, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(BI, Ck)*g(Kc, JA) ! 10
!
            CJ = index_two(wf%C, wf%J, wf%n_v)
            Ak = index_two(wf%A, k, wf%n_v)
            Ic = index_two(wf%I, c, wf%n_o)
            KB = index_two(wf%K, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CJ, Ak)*g(Ic, KB) ! 11
!
            CI = index_two(wf%C, wf%I, wf%n_v)
            Bk = index_two(wf%B, k, wf%n_v)
            Jc = index_two(wf%J, c, wf%n_o)
            KA = index_two(wf%K, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) + b2am(CI, Bk)*g(Jc, KA) ! 12
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            Bk = index_two(wf%B, k, wf%n_v)
            Jc = index_two(wf%J, c, wf%n_o)
            KC = index_two(wf%K, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(AI, Bk)*g(Jc, KC) ! 13
!
            AI = index_two(wf%A, wf%I, wf%n_v)
            Ck = index_two(wf%C, k, wf%n_v)
            Kc = index_two(wf%K, c, wf%n_o)
            JB = index_two(wf%J, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(AI, Ck)*g(Kc, JB) ! 14
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            Ak = index_two(wf%A, k, wf%n_v)
            Ic = index_two(wf%I, c, wf%n_o)
            KC = index_two(wf%K, wf%C, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(BJ, Ak)*g(Ic, KC) ! 15
!
            BJ = index_two(wf%B, wf%J, wf%n_v)
            Ck = index_two(wf%C, k, wf%n_v)
            Kc = index_two(wf%K, c, wf%n_o)
            IA = index_two(wf%I, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(BJ, Ck)*g(Kc, IA) ! 16
!  
            CK = index_two(wf%C, wf%K, wf%n_v)
            Ak = index_two(wf%A, k, wf%n_v)
            Ic = index_two(wf%I, c, wf%n_o)
            JB = index_two(wf%J, wf%B, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(CK, Ak)*g(Ic, JB) ! 17
!
            CK = index_two(wf%C, wf%K, wf%n_v)
            Bk = index_two(wf%B, k, wf%n_v)
            Jc = index_two(wf%J, c, wf%n_o)
            IA = index_two(wf%I, wf%A, wf%n_o)
!
            sigma(c, k) = sigma(c, k) - two*b2am(CK, Bk)*g(Jc, IA) ! 18
!
         enddo
      enddo
!
   end subroutine jacobian_transpose_sccsd_c1_sccsd
!
!
end submodule jacobian_transpose 
