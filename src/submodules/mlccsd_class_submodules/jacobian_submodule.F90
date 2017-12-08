submodule (mlccsd_class) jacobian
!
!!
!!    Jacobian transformation submodule (MLCCSD) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2017
!!
!
   implicit none 
!
   logical :: debug   = .false.
   logical :: timings = .false.
!
   character(len=40) :: integral_type
!
contains
!
   module subroutine jacobian_mlccsd_transformation_mlccsd(wf, c_a_i, c_aibj)
!!
!!    Jacobian transformation (MLCC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
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
      class(mlccsd) :: wf 
!
!      Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_x2am, 1)   :: c_aibj ! c_aibj     
!
!      Local unpacked and reordered vectors 
!
      real(dp), dimension(:,:), allocatable :: rho_a_i         ! rho_ai   = (A c)_ai
      real(dp), dimension(:,:), allocatable :: rho_ai_bj       ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: rho_ab_ij       ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_sym   ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: c_ai_bj         ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: c_ab_ij         ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: x_ia_jb         ! rho_ai   = (A c)_aibj
!
!      Indices 
!
      integer(i15) :: a = 0, ab = 0, ai = 0, b = 0 
      integer(i15) :: bj = 0, i = 0, ij = 0, j = 0, aibj = 0
!
      integer(i15) :: offset = 0
!
!    Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: active_space = 0
!
!    Allocate and zero the transformed vector (singles part)
!
      call wf%mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!    :: CCS contributions to the singles c vector ::  
!
      call wf%read_single_amplitudes
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
!
!     :: MLCCSD contributions to transformed vector :: 
!
      call wf%jacobian_mlcc2_a1(rho_a_i, c_a_i)
!
!     Calculate number of active indices
! 
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
!     Allocate the incoming unpacked doubles vector 
!
      call wf%mem%alloc(c_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v) 
      c_ai_bj = zero
!
      call squareup(c_aibj, c_ai_bj, n_active_o*n_active_v) ! Pack out vector 
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
      do i = 1, n_active_o*n_active_v
!
         c_ai_bj(i,i) = two*c_ai_bj(i,i)
!
      enddo
!
      call wf%jacobian_mlcc2_b1(rho_a_i, c_ai_bj)
!
!     Allocate unpacked transformed vector
!
      call wf%mem%alloc(rho_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v) 
      rho_ai_bj = zero 
!
      call wf%jacobian_mlcc2_a2(rho_ai_bj, c_a_i)
      call wf%jacobian_mlccsd_b2(rho_ai_bj, c_a_i)
      call wf%jacobian_mlccsd_c2(rho_ai_bj, c_a_i)
      call wf%jacobian_mlccsd_d2(rho_ai_bj, c_a_i)
!
      call wf%jacobian_mlccsd_e2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_mlccsd_f2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_mlccsd_g2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_mlccsd_h2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_mlccsd_i2(rho_ai_bj, c_ai_bj)
!
!     Allocate temporary symmetric transformed vector 
!
      call wf%mem%alloc(rho_ai_bj_sym, n_active_o*n_active_v, n_active_o*n_active_v)
      rho_ai_bj_sym = zero
!!
      do j = 1, n_active_o
         do b = 1, n_active_v
!!
            bj = index_two(b, j, n_active_v)
!!
            do i = 1, n_active_o
               do a = 1, n_active_v
!!
                  ai = index_two(a, i, n_active_v)
!!
                  rho_ai_bj_sym(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj(bj, ai)
!!
               enddo
            enddo
         enddo
      enddo
!!
      rho_ai_bj = rho_ai_bj_sym
!
 
!     Done with temporary vector; deallocate
!  
      call wf%mem%dealloc(rho_ai_bj_sym, n_active_o*n_active_v, n_active_o*n_active_v)
! 
!     In preparation for last two terms, reorder 
!     rho_ai_bj to rho_ab_ij, and c_ai_bj to c_ab_ij
! 
      call wf%mem%alloc(rho_ab_ij, (n_active_v)**2, (n_active_o)**2)
      call wf%mem%alloc(c_ab_ij, (n_active_v)**2, (n_active_o)**2)
! 
      rho_ab_ij = zero
      c_ab_ij   = zero
! 
      do j = 1, n_active_o
         do i = 1, n_active_o
! 
            ij = index_two(i, j, n_active_o)
! 
            do b = 1, n_active_v
! 
               bj = index_two(b, j, n_active_v)
! 
               do a = 1, n_active_v
! 
                  ai = index_two(a, i, n_active_v)
                  ab = index_two(a, b, n_active_v)
! 
                  c_ab_ij(ab, ij)   = c_ai_bj(ai, bj)
                  rho_ab_ij(ab, ij) = rho_ai_bj(ai, bj)
! 
               enddo
            enddo
         enddo
      enddo
! 
      call wf%mem%dealloc(c_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
      call wf%mem%dealloc(rho_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
! 
      call wf%jacobian_mlccsd_j2(rho_ab_ij, c_ab_ij)
      call wf%jacobian_mlccsd_k2(rho_ab_ij, c_ab_ij)
! 
!     Done with reordered doubles c; deallocate 
! 
      call wf%mem%dealloc(c_ab_ij, (n_active_v)**2, (n_active_o)**2)
! 
!     Order rho_ab_ij back into rho_ai_bj & divide by 
!     the biorthonormal factor 1 + delta_ai,bj
! 
      call wf%mem%alloc(rho_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
! 
      do j = 1, n_active_o
         do b = 1, n_active_v
! 
            bj = index_two(b, j, n_active_v)
! 
            do i = 1, n_active_o
! 
               ij = index_two(i, j, n_active_o)
! 
               do a = 1, n_active_v
! 
                  ab = index_two(a, b, n_active_v)
                  ai = index_two(a, i, n_active_v)
! 
                  if ((a .eq. b) .and. (i .eq. j)) then 
! 
                     rho_ai_bj(ai, bj) = half*rho_ab_ij(ab, ij)
! 
                  else
! 
                     rho_ai_bj(ai, bj) = rho_ab_ij(ab, ij)
! 
                  endif
! 
               enddo
            enddo
         enddo
      enddo
! 
!     Done with reordered transformed vector; deallocate 
! 
      call wf%mem%dealloc(rho_ab_ij, (n_active_v)**2, (n_active_o)**2)
!
!     c_a_i -> rho_a_i
!     c_ai_bj -> rho_ai_bj
!
      c_a_i = rho_a_i
!
      c_aibj = zero
      call packin(c_aibj, rho_ai_bj, n_active_o*n_active_v)
!
!     Deallocations
!
      call wf%mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
      call wf%mem%dealloc(rho_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
!
   end subroutine jacobian_mlccsd_transformation_mlccsd
!
!
   module subroutine jacobian_mlccsd_b2_mlccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD B2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^B2 = - sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^B2,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension(:,:) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: x_c_aij ! t_ij^ac 
      real(dp), dimension(:,:), allocatable :: x_aib_k ! t_ik^ab 
! 
      real(dp), dimension(:,:), allocatable :: I_k_aij ! An intermediate
      real(dp), dimension(:,:), allocatable :: I_k_j   ! An intermediate 
!
      real(dp), dimension(:,:), allocatable :: rho_b_aij ! rho_ai_bj, unordered, term 1
      real(dp), dimension(:,:), allocatable :: rho_aib_j ! rho_ai_bj, unordered, term 2
!
      integer(i15) :: a = 0, c = 0, i = 0, j = 0, b = 0, k = 0
      integer(i15) :: ai = 0, cj = 0, aicj = 0, aij = 0, bj = 0
      integer(i15) :: bk = 0, aibk = 0, aib = 0
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: first_CCSD_o ! first active occupied index 
      integer(i15) :: first_CCSD_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
!     :: Term 1. - sum_kc F_kc t_ij^ac c_bk ::
!
!     k           - general index
!     c           - CC2 index
!     a, b, i, j  - CCSD indices
!
!     Read the amplitudes from disk 
!
      call wf%read_amplitudes
!
!     Order the amplitudes as t_c_aij = t_ij^ac 
!
      call wf%mem%alloc(x_c_aij, n_CC2_v, (n_CCSD_v)*(n_CCSD_o**2))
!
      do j = 1, n_CCSD_o
         do i = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CC2_v)
!
               aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
               do c = 1, n_CC2_v
!
                  cj = index_two(c, j, n_CC2_v)
!
                  aicj = index_packed(ai, cj)
!
                  x_c_aij(c, aij) = wf%x2am(aicj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_x2am
!
!     Form the intermediate I_k_aij = sum_c F_k_c * x_c_aij 
!
      call wf%mem%alloc(I_K_aij, wf%n_o, (n_CCSD_v)*(n_CCSD_o)**2)
!
      call dgemm('N', 'N',                   &
                  wf%n_o,                    &
                  (n_CCSD_v)*(n_CCSD_o)**2,  &
                  n_CC2_v,                   &
                  one,                       &
                  wf%fock_ia,                &
                  wf%n_o,                    &
                  x_c_aij,                   &
                  n_CC2_v,                   &   
                  zero,                      &
                  I_k_aij,                   &
                  wf%n_o)
!
      call wf%mem%dealloc(x_c_aij, n_CC2_v, (n_CCSD_v)*(n_CCSD_o**2))
!
!     Form rho_b_aij = sum_k c_b_k X_k_aij(k,aij)
!
      call wf%mem%alloc(rho_b_aij, n_CCSD_v, (n_CCSD_v)*(n_CCSD_o)**2)
!
      call dgemm('N', 'N',                   &
                  n_CCSD_v,                  &
                  (n_CCSD_v)*(n_CCSD_o)**2,  &
                  wf%n_o,                    &
                  -one,                      & 
                  c_a_i,                     & ! c_b_k
                  wf%n_v,                    &
                  I_k_aij,                   &
                  wf%n_o,                    &
                  zero,                      &
                  rho_b_aij,                 &
                  n_CCSD_v)
!
      call wf%mem%dealloc(I_k_aij, wf%n_o, (n_CCSD_v)*(n_CCSD_o)**2)
!
!     Add rho_b_aij to rho_ai_bj 
!
      do j = 1, n_CCSD_o
         do i = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CC2_v)
!
               aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
               do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CC2_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_aij(b, aij)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_b_aij, n_CCSD_v, (n_CCSD_v)*(n_CCSD_o)**2)
!
!
!     :: Term 2. - sum_kc F_kc t_ik^ab c_cj ::
!
!     c           - general index
!     k           - CC2 index
!     a, b, i, j  - CCSD indices 
!
!     Form I_k_j = sum_c F_kc c_cj = sum_c fock_ia(k,c) c_a_i(c,j)
!
      call wf%mem%alloc(I_k_j, n_CC2_o, n_CCSD_o)
!
      call dgemm('N','N',        &
                  n_CC2_o,       &
                  n_CCSD_o,      &
                  wf%n_v,        &
                  one,           &
                  wf%fock_ia,    & ! F_k_c
                  wf%n_o,        &
                  c_a_i,         & ! c_c_j
                  wf%n_v,        &
                  zero,          &
                  I_k_j,         &
                  n_CC2_o)
!
!     Order the amplitudes as x_aib_k = t_ik^ab 
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(x_aib_k, ((n_CCSD_v)**2)*(n_CCSD_o), n_CC2_o)
      x_aib_k = zero
!
      do k = 1, n_CC2_o
         do b = 1, n_CCSD_v
!
            bk = index_two(b, k, n_CC2_v)
!
            do i = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
!
                  aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
                  aibk = index_packed(ai, bk)
!
                  x_aib_k(aib, k) = wf%x2am(aibk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate doubles amplitudes    
!
      call wf%destruct_x2am
! 
!     Form rho_aib_j = - sum_k x_aib_k I_k_j
!
      call wf%mem%alloc(rho_aib_j, (n_CCSD_o)*(n_CCSD_v)**2, n_CCSD_o)
!
      call dgemm('N','N',                    &
                  (n_CCSD_o)*(n_CCSD_v)**2,  &
                  n_CCSD_o,                  &
                  n_CC2_o,                   &
                  -one,                      & 
                  x_aib_k,                   &
                  (n_CCSD_o)*(n_CCSD_v)**2,  &
                  I_k_j,                     &
                  n_CC2_o,                   &
                  zero,                      &
                  rho_aib_j,                 &
                  (n_CCSD_o)*(n_CCSD_v)**2)
!
      call wf%mem%dealloc(I_k_j, n_CC2_o, n_CCSD_o)
      call wf%mem%dealloc(x_aib_k, (n_CCSD_o)*(n_CCSD_v)**2, n_CC2_o)
!
      do a = 1, n_CCSD_v
         do i = 1, n_CCSD_o
!
            ai = index_two(a, i, n_CC2_v)
!
            do b = 1, n_CCSD_v
!
               aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
               do j = 1, n_CCSD_o
!
                  bj = index_two(b, j, n_CC2_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aib_j(aib, j)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_aib_j, (n_CCSD_o)*(n_CCSD_v)**2, n_CCSD_o)
!
   end subroutine jacobian_mlccsd_b2_mlccsd
!
!
   module subroutine jacobian_mlccsd_c2_mlccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD C2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2017
!!
!!    rho_ai_bj^C2 = sum_kcl (g_ljkc * t_ki^ac * c_bl) + (g_ljkc * t_li^bc * c_ak) 
!!                         + (g_ljkc * t_lk^ba * c_ci) 
!!                         - (L_ljkc * t_ik^ac * c_bl)- (L_ljkc * t_il^ab * c_ck)
!!                
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
      real(dp), dimension(:,:)            :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: x_ai_kc
      real(dp), dimension(:,:), allocatable :: x_kc_ai
      real(dp), dimension(:,:), allocatable :: x_ab_kl
      real(dp), dimension(:,:), allocatable :: x_aib_l
      real(dp), dimension(:,:), allocatable :: c_kc
      real(dp), dimension(:,:), allocatable :: g_kc_lj
      real(dp), dimension(:,:), allocatable :: g_lj_kc
      real(dp), dimension(:,:), allocatable :: g_lc_jk
      real(dp), dimension(:,:), allocatable :: g_kj_lc
      real(dp), dimension(:,:), allocatable :: g_kl_jc
      real(dp), dimension(:,:), allocatable :: I_lj_ai
      real(dp), dimension(:,:), allocatable :: I_kl_ji
      real(dp), dimension(:,:), allocatable :: I_bi_jk
      real(dp), dimension(:,:), allocatable :: I_ab_jc
      real(dp), dimension(:,:), allocatable :: I_ljk_i
      real(dp), dimension(:,:), allocatable :: I_lj
      real(dp), dimension(:,:), allocatable :: rho_b_jai
      real(dp), dimension(:,:), allocatable :: rho_bij_a
      real(dp), dimension(:,:), allocatable :: rho_ab_ji
      real(dp), dimension(:,:), allocatable :: rho_aib_j
!
      integer(i15) :: i = 0, k = 0, j = 0, l = 0
      integer(i15) :: a = 0, b = 0, c = 0
!
      integer(i15) :: ab = 0
      integer(i15) :: ai = 0, ak = 0, ci = 0, ck = 0, bj = 0, bl = 0, ai_CC2
      integer(i15) :: kc = 0, lc = 0, jc = 0
      integer(i15) :: jl = 0, lj = 0, kl = 0, ji = 0, kj = 0, jk = 0
      integer(i15) :: jai = 0, bij = 0, abj = 0, ljk = 0, aib = 0
      integer(i15) :: akci = 0, akbl = 0, aick = 0, aibl = 0
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: first_CCSD_o ! first active occupied index 
      integer(i15) :: first_CCSD_v ! first active virtual index
!
      integer(i15) :: last_CC2_o ! first active occupied index 
      integer(i15) :: last_CC2_v ! first active virtual index
      integer(i15) :: last_CCSD_o ! first active occupied index 
      integer(i15) :: last_CCSD_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!     Construct x_ai_kc (= t_ki^ac) for term 1
!     Will also be used as t_li^bc for term 2
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(x_ai_kc, n_CCSD_v*n_CCSD_o, n_CC2_v*n_CC2_o)
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CCSD_v)
!
            do k = 1, n_CC2_o
!   
               ak = index_two(a, k, n_CC2_v)
!
               do c = 1, n_CC2_v
!
                  kc = index_two(k, c, n_CC2_o) 
                  ci = index_two(c, i, n_CC2_v)
!
                  akci = index_packed(ak, ci) 
!
                  x_ai_kc(ai, kc) = wf%x2am(akci, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_x2am
!
!     :: Term 1 ::
!     sum_kcl (g_ljkc * t_ki^ac * c_bl)
!
!     l           - general index
!     k, c        - CC2 indices 
!     a, i, b, j  - CCSD indices
!
!     Construct g_kc_lj (=g_lj,kc)
!
      call wf%mem%alloc(g_lj_kc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_lj_kc,    &
                        1, wf%n_o,                 &
                        first_CCSD_o, last_CCSD_o, &
                        first_CC2_o, last_CC2_o,    &
                        first_CC2_v, last_CC2_v)
!
!     I_ai_lj = sum_(kc) x_ai_kc*g_kc_lj
!
      call wf%mem%alloc(I_lj_ai, (wf%n_o)*(n_CCSD_o), (n_CCSD_o)*(n_CCSD_v))
!
      call dgemm('N', 'T',                &
                  (wf%n_o)*(n_CCSD_o),    &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  n_CC2_o*n_CC2_v,        & 
                  one,                    &
                  g_lj_kc,                &
                  (wf%n_o)*(n_CCSD_o),    &
                  x_ai_kc,                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  zero,                   &
                  I_lj_ai,                &
                  (wf%n_o)*(n_CCSD_o))
!
      call wf%mem%dealloc(g_lj_kc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
!     rho_b_jai = sum_(l) I_lj_ai * c_bl
!
      call wf%mem%alloc(rho_b_jai, n_CCSD_v, (n_CCSD_o**2)*(n_CCSD_v))
!
      call dgemm('N', 'N', &
                  n_CCSD_v,                  &
                  (n_CCSD_o**2)*(n_CCSD_v),  &
                  wf%n_o,                    &
                  one,                       &
                  c_a_i,                     & ! c_b_l
                  wf%n_v,                    &
                  I_lj_ai,                   & ! I_l_jai
                  wf%n_o,                    &
                  zero,                      &
                  rho_b_jai,                 &
                  (n_CCSD_v))
!
!
      call wf%mem%dealloc(I_lj_ai, (wf%n_o)*(n_CCSD_o), (n_CCSD_o)*(n_CCSD_v))
!
!     Add terms to rho_ai_bj
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CC2_v)
!
            do j = 1, n_CCSD_o
!
               jai = index_three(j, a, i, n_CCSD_o, n_CCSD_v)
!
               do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CC2_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_jai(b, jai)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_b_jai, n_CCSD_v, (n_CCSD_o**2)*(n_CCSD_v))
!
!     :: Term 2 ::
!     sum_kcl (g_ljkc * x_li^bc * c_ak)
!
!     k           - general index
!     l, c        - CC2 index 
!     b, i, j, a  - CCSD index 
!
!     Construct g_kj_lc (= g_ljkc)
!
      call wf%mem%alloc(g_kc_lj, (wf%n_o)*(n_CC2_v), (n_CCSD_o)*(n_CC2_o))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_oo(integral_type, g_kc_lj,  &
                        1, wf%n_o,               &
                        first_CC2_v, last_CC2_v, &
                        first_CC2_o, last_CC2_o, &
                        first_CCSD_o, last_CCSD_o)
!
!     Reorder g_kc_lj to g_lc_jk
!
      call wf%mem%alloc(g_lc_jk, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(wf%n_o))
!
      do j = 1, n_CCSD_o
         do l = 1, n_CC2_o
!
            lj = index_two(l, j, n_CC2_o)
!
            do c = 1, n_CC2_v
!
               lc = index_two(l, c, n_CC2_o)
!
               do k = 1, wf%n_o
!
                  kc = index_two(k, c, wf%n_o)
                  jk = index_two(j, k, n_CCSD_o)
!
                  g_lc_jk(lc, jk) = g_kc_lj(kc, lj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_lj, (wf%n_o)*(n_CC2_v), (n_CC2_o)*(n_CCSD_o))
!
!     I_bi_jk = sum_(kc) t_bi_lc*g_lc_jk
!
      call wf%mem%alloc(I_bi_jk, (n_CCSD_o)*(n_CCSD_v), (wf%n_o)*(n_CCSD_o))
!
      call dgemm('N', 'N',                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (wf%n_o)*(n_CCSD_o),    &
                  (n_CC2_o)*(n_CC2_v),    &
                  one,                    &
                  x_ai_kc,                & ! x_bi_lc
                  (n_CCSD_o)*(n_CCSD_v),  &
                  g_lc_jk,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  zero,                   &
                  I_bi_jk,                &
                  (n_CCSD_o)*(n_CCSD_v))
!
      call wf%mem%dealloc(g_lc_jk, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(wf%n_o))
      call wf%mem%dealloc(x_ai_kc, n_CCSD_v*n_CCSD_o, n_CC2_v*n_CC2_o)
!
!     rho_bij_a = sum_(k) I_bi_jk * c_ak
!
      call wf%mem%alloc(rho_bij_a, (n_CCSD_o**2)*(n_CCSD_v), n_CCSD_v)
      call dgemm('N', 'T', &
                  (n_CCSD_o**2)*(n_CCSD_v), &
                  n_CCSD_v,                 &
                  wf%n_o,                   &
                  one,                      &
                  I_bi_jk,                  & !I_bij_k
                  (n_CCSD_o**2)*(n_CCSD_v), &
                  c_a_i,                    &
                  wf%n_v,                   &
                  zero,                     &
                  rho_bij_a,                &
                  (n_CCSD_o**2)*(n_CCSD_v))
!
      call wf%mem%dealloc(I_bi_jk, (n_CCSD_o)*(n_CCSD_v), (wf%n_o)*(n_CCSD_o))
!
!     Add terms to rho_ai_bj
!
      do i = 1, n_CCSD_o
         do b = 1, n_CCSD_v
!
            do j = 1, n_CCSD_o
!
               bj = index_two(b, j, n_CC2_v)
               bij = index_three(b, i, j, n_CCSD_v, n_CCSD_o)
!
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_bij_a(bij, a)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_bij_a, (n_CCSD_o**2)*(n_CCSD_v), n_CCSD_v) 
!
!     :: Term 3 ::
!     sum_(kcl) g_ljkc * t_lk^ba * c_ci
!
!     c           - general index 
!     l, k        - CC2 indices 
!     i, j, a, b, - CCSD indices
!
!     Construct and reorder t_kl^ab to t_ab_kl
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(x_ab_kl, n_CCSD_v**2, n_CC2_o**2)
!
      do b = 1, n_CCSD_v
         do a = 1, n_CCSD_v
!
            ab = index_two(a, b, n_CCSD_v)
!
            do k = 1, n_CC2_o
!
               ak = index_two(a, k, n_CC2_v)
!
               do l = 1, n_CC2_o
!
                  kl = index_two(k, l, n_CC2_o)
                  bl = index_two(b, l, n_CC2_v)
                  akbl = index_packed(ak, bl)
!
                  x_ab_kl(ab, kl) = wf%x2am(akbl, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_x2am
!
      call wf%mem%alloc(g_lj_kc, (n_CC2_o)*(n_CCSD_o), (n_CC2_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_lj_kc ,   &
                        first_CC2_o, last_CC2_o,   &
                        first_CCSD_o, last_CCSD_o, &
                        first_CC2_o, last_CC2_o,   &
                        1, wf%n_v)
!
!     I_ljk_i = sum_(c) g_lj_kc * c_ci
!
      call wf%mem%alloc(I_ljk_i, (n_CC2_o**2)*(n_CCSD_o), (n_CCSD_o))
!
      call dgemm('N','N',                    &
                  (n_CC2_o**2)*(n_CCSD_o),   &
                  n_CCSD_o,                  &
                  wf%n_v,                    &
                  one,                       &
                  g_lj_kc,                   & ! "g_ljk_c"
                  (n_CC2_o**2)*(n_CCSD_o),   &
                  c_a_i,                     &
                  wf%n_v,                    &
                  zero,                      &
                  I_ljk_i,                   &
                  (n_CC2_o**2)*(n_CCSD_o))
!
      call wf%mem%dealloc(g_lj_kc, (n_CC2_o)*(n_CCSD_o), (n_CC2_o)*(wf%n_v))
!
!     Reorder I_kjl_i to I_kl_ji
!
      call wf%mem%alloc(I_kl_ji, n_CC2_o**2, n_CCSD_o**2)
!
      do i = 1, n_CCSD_o
         do j = 1, n_CCSD_o
!
            ji = index_two(j, i, n_CCSD_o)
!
            do k = 1, n_CC2_o
               do l = 1, n_CC2_o
!
                  ljk = index_three(l, j, k, n_CC2_o, n_CCSD_o)
!
                  kl = index_two(k, l, n_CC2_o)
!
                  I_kl_ji(kl, ji) = I_ljk_i(ljk, i)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(I_ljk_i, (n_CC2_o**2)*(n_CCSD_o), (n_CCSD_o))
!
!     rho_ab_ji = sum_(kl) t_ab_kl * I_kl_ji
!
      call wf%mem%alloc(rho_ab_ji, n_CCSD_v**2, n_CCSD_o**2)
!
      call dgemm('N', 'N',       &
                  n_CCSD_v**2,   &
                  n_CCSD_o**2,   &
                  n_CC2_o**2,    &
                  one,           &
                  x_ab_kl,       &
                  n_CCSD_v**2,   &
                  I_kl_ji,       &
                  n_CC2_o**2,    &
                  zero,          &
                  rho_ab_ji,     &
                  n_CCSD_v**2)
!
      call wf%mem%dealloc(I_kl_ji, n_CC2_o**2, n_CCSD_o**2)
      call wf%mem%dealloc(x_ab_kl, n_CCSD_v**2, n_CC2_o**2)
!
!     Reorder into rho_ai_bj
!
      do j = 1, n_CCSD_o
         do i = 1, n_CCSD_o
!
            ji = index_two(j, i, n_CCSD_o)
!
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CC2_v)
!
               do b = 1, n_CCSD_v
!
                  ab = index_two(a, b, n_CCSD_v)
                  bj = index_two(b, j, n_CC2_v)
!  
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ab_ji(ab, ji)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_ab_ji, n_CCSD_v**2, n_CCSD_o**2) 
!
!     :: Term 4 :: 
!     - sum_(klc) L_ljkc * t_ik^ac * c_bl
!     = - 2 sum_(klc) g_ljkc * t_ik^ac * c_bl
!         + sum_(klc) g_kjlc * t_ik^ac * c_bl
!     = (4a) + (4b)
!
!     Construct t_ik^ac ordered as t_kc_ai
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(x_kc_ai, n_CC2_o*n_CC2_v, n_CCSD_o*n_CCSD_v)
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CCSD_v)
            ai_CC2 = index_two(a, i, n_CC2_v)
!
            do k = 1, n_CC2_o
               do c = 1, n_CC2_v
!
                  kc = index_two(k, c, n_CC2_o)
                  ck = index_two(c, k, n_CC2_v)
                  aick = index_packed(ai_CC2, ck)
!
                  x_kc_ai(kc, ai) = wf%x2am(aick, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_x2am
!
!     :: 4a ::
!     - 2 sum_(klc) g_ljkc * t_ik^ac * c_bl
!
!     l           - general index 
!     k, c        - CC2 indices 
!     a, i, b, j  - CCSD indices
!
!     Construct g_ljkc
!
      call wf%mem%alloc(g_lj_kc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
      integral_type ='electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_lj_kc,    &
                        1, wf%n_o,                 &
                        first_CCSD_o, last_CCSD_o, &
                        first_CC2_o, last_CC2_o,   &
                        first_CC2_v, last_CC2_v)
!
!     I_lj_ai = - 2sum_(ck) g_lj_kc *t_kc_ai
!
      call wf%mem%alloc(I_lj_ai,(wf%n_o)*(n_CCSD_o), (n_CCSD_v)*(n_CCSD_o))
!
      call dgemm('N', 'N',                &
                  (wf%n_o)*(n_CCSD_o),    &
                  (n_CCSD_v)*(n_CCSD_o),  &
                  (n_CC2_v)*(n_CC2_o),    &
                  -two,                   &
                  g_lj_kc,                &
                  (wf%n_o)*(n_CCSD_o),    &
                  x_kc_ai,                &
                  (n_CC2_v)*(n_CC2_o),    &
                  zero,                   &
                  I_lj_ai,                &
                  (wf%n_o)*(n_CCSD_o)) 
!     
      call wf%mem%dealloc(g_lj_kc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
!     :: 4b ::
!     sum_(klc) g_kjlc * t_ik^ac * c_bl
!     l           - general index 
!     k, c        - CC2 indices 
!     a, i, j, b  - CCSD indices
!
!     Construct g_kjlc
!
      call wf%mem%alloc(g_kj_lc, (n_CC2_o)*(n_CCSD_o), (wf%n_o)*(n_CC2_v))
!
      integral_type ='electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_kj_lc,    &
                        first_CC2_o, last_CC2_o,   &
                        first_CCSD_o, last_CCSD_o, &
                        1, wf%n_o,                 &
                        first_CC2_v, last_CC2_v)
!
!     Reorder g_kj_lc to g_lj_kc      
!
      call wf%mem%alloc(g_lj_kc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
      do j = 1, n_CCSD_o
         do c = 1, n_CC2_v 
            do l = 1, wf%n_o
!
               lj = index_two(l, j, wf%n_o)
!
               do k = 1, n_CC2_o
!
                  kj = index_two(k, j, n_CC2_o)
!
                  kc = index_two(k, c, n_CC2_o)  
                  lc = index_two(l, c, wf%n_o)  
!
                  g_lj_kc(lj, kc) = g_kj_lc(kj, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kj_lc, (n_CC2_o)*(n_CCSD_o), (wf%n_o)*(n_CC2_v))
!
!     I_lj_ai += sum_(kc) g_lj_kc * x_kc_ai
!
      call dgemm('N', 'N',                &
                  (wf%n_o)*(n_CCSD_o),    &
                  (n_CCSD_v)*(n_CCSD_o),  &
                  (n_CC2_v)*(n_CC2_o),    &
                  one,                    &
                  g_lj_kc,                &
                  (wf%n_o)*(n_CCSD_o),    &
                  x_kc_ai,                &
                  (n_CC2_v)*(n_CC2_o),    &
                  one,                    &
                  I_lj_ai,                &
                  (wf%n_o)*(n_CCSD_o))
!
      call wf%mem%dealloc(g_lj_kc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
      call wf%mem%dealloc(x_kc_ai, n_CC2_o*n_CC2_v, n_CCSD_o*n_CCSD_v)
!
!     rho_b_jai = sum_(l) c_bl *I_lj_ai
!
      call wf%mem%alloc(rho_b_jai, n_CCSD_v, (n_CCSD_v)*(n_CCSD_o**2))
!
      call dgemm('N', 'N',                   &
                  n_CCSD_v,                  &
                  (n_CCSD_v)*(n_CCSD_o**2),  &
                  wf%n_o,                    &
                  one,                       &
                  c_a_i,                     & ! c_bl
                  wf%n_v,                    &
                  I_lj_ai,                   & ! I_l_jai
                  wf%n_o,                    &
                  zero,                      &
                  rho_b_jai,                 &
                  (n_CCSD_v))
!
      call wf%mem%dealloc(I_lj_ai,(wf%n_o)*(n_CCSD_o), (n_CCSD_v)*(n_CCSD_o))
!
!     Add terms to rho_ai_bj
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CC2_v)
!
            do j = 1, n_CCSD_o
!
               jai = index_three(j, a, i, n_CCSD_o, n_CCSD_v)
!
               do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CC2_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_jai(b, jai)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_b_jai, n_CCSD_v, (n_CCSD_o)*(n_CCSD_v)*n_CCSD_o)
!     :: Term 5 :: 
!     - sum_(klc) L_ljkc * t_il^ab * c_ck
!     = - 2 sum_(klc) g_ljkc * t_il^ab * c_ck
!         + sum_(klc) g_kjlc * t_il^ab * c_ck
!     = (5a) + (5b)
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(x_aib_l, (n_CCSD_o)*(n_CCSD_v**2), n_CC2_o)
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CC2_v)
!
            do l = 1, n_CC2_o
               do b = 1, n_CCSD_v
!
                  bl = index_two(b, l, n_CC2_v)
!
                  aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
                  aibl = index_packed(ai, bl)
!
                  x_aib_l(aib, l) = wf%x2am(aibl, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_x2am
!
!     :: 5a ::
!     - 2 sum_(klc) g_ljkc * t_il^ab * c_ck
!
!     k, c        - general index 
!     l           - CC2 index  
!     a, i, b, j  - CCSD index  
!
!     Construct g_ljkc
!
      call wf%mem%alloc(g_lj_kc, (n_CC2_o)*(n_CCSD_o), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_lj_kc,    &
                        first_CC2_o, last_CC2_o,   &
                        first_CCSD_o, last_CCSD_o, &
                        1, wf%n_o,                 &
                        1, wf%n_v)
! 
      call wf%mem%alloc(c_kc, (wf%n_o)*(wf%n_v), 1)
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            kc = index_two(k, c, wf%n_o)
            c_kc(kc, 1) = c_a_i(c, k)
         enddo
      enddo
!
      call wf%mem%alloc(I_lj, (n_CC2_o)*(n_CCSD_o), 1)
!
      call dgemm('N', 'N',             &
                 (n_CC2_o)*(n_CCSD_o), &
                 1,                    &
                 (wf%n_o)*(wf%n_v),    &
                 -two,                 &
                 g_lj_kc,              &
                 (n_CC2_o)*(n_CCSD_o), &
                 c_kc,                 &
                 (wf%n_o)*(wf%n_v),    &
                 zero,                 &
                 I_lj,                 &
                 (n_CC2_o)*(n_CCSD_o))
!
      call wf%mem%dealloc(g_lj_kc, (n_CC2_o)*(n_CCSD_o), (wf%n_o)*(wf%n_v))
!
!     :: 5a ::
!     sum_(klc) g_kjlc * t_il^ab * c_ck
!
!     k, c        - general index 
!     l           - CC2 index  
!     a, i, b, j  - CCSD index  
!
!     Construct g_ljkc

      call wf%mem%alloc(g_kj_lc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_kj_lc,    &
                        1, wf%n_o,                 &
                        first_CCSD_o, last_CCSD_o, &
                        first_CC2_o, last_CC2_o,   &
                        1, wf%n_v)
!
!     Reorder g_kj_lc to g_lj_kc
!
      call wf%mem%alloc(g_lj_kc, (n_CC2_o)*(n_CCSD_o), (wf%n_o)*(wf%n_v))
!
      do j = 1, n_CCSD_o
         do c = 1, wf%n_v 
            do k = 1, wf%n_o
!
               kc = index_two(k, c, wf%n_o)
               kj = index_two(k, j, wf%n_o)
!
               do l = 1, n_CC2_o
!
                  lj = index_two(l, j, n_CC2_o)
                  lc = index_two(l, c, n_CC2_o)  
!
                  g_lj_kc(lj, kc) = g_kj_lc(kj, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kj_lc, (wf%n_o)*(n_CCSD_o), (n_CC2_o)*(wf%n_v))
!
      call dgemm('N', 'N',             &
                 (n_CC2_o)*(n_CCSD_o), &
                 1,                    &
                 (wf%n_o)*(wf%n_v),    &
                 one,                  &
                 g_lj_kc,              &
                 (n_CC2_o)*(n_CCSD_o), &
                 c_kc,                 &
                 (wf%n_o)*(wf%n_v),    &
                 one,                  &
                 I_lj,                 &
                 (n_CC2_o)*(n_CCSD_o))
!
      call wf%mem%dealloc(g_lj_kc, (n_CC2_o)*(n_CCSD_o), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(c_kc, (wf%n_o)*(wf%n_v), 1)
!
      call wf%mem%alloc(rho_aib_j,(n_CCSD_o)*(n_CCSD_v**2), n_CCSD_o)
!
      call dgemm('N', 'N', &
                  (n_CCSD_o)*(n_CCSD_v**2), &
                  n_CCSD_o, &
                  n_CC2_o, &
                  one, &
                  x_aib_l, &
                  (n_CCSD_o)*(n_CCSD_v**2), &
                  I_lj, &
                  n_CC2_o, &
                  zero, &
                  rho_aib_j, &
                  (n_CCSD_o)*(n_CCSD_v**2))
!
      
      call wf%mem%dealloc(I_lj, (n_CC2_o)*(n_CCSD_o), 1)
      call wf%mem%dealloc(x_aib_l, (n_CCSD_o)*(n_CCSD_v**2), n_CC2_o)
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CC2_v)
!
            do j = 1, n_CCSD_o
               do b = 1, n_CCSD_v
!                  
                  bj = index_two(b, j, n_CC2_v)
                  aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aib_j(aib, j)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_aib_j,(n_CCSD_o)*(n_CCSD_v**2), n_CCSD_o)
!
   end subroutine jacobian_mlccsd_c2_mlccsd
!
!
   module subroutine jacobian_mlccsd_d2_mlccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD D2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^D2 = - sum_kcd g_kcbd (x_ij^cd c_ak + x_kj^ad c_ci + x_ik^ca c_dj)
!!                       + sum_kcd L_kcbd (x_ik^ac c_dj + x_ij^ad c_ck)
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
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension(:,:) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: L_kc_J
      real(dp), dimension(:,:), allocatable :: L_bd_J
!
      real(dp), dimension(:,:), allocatable :: g_kc_bd
      real(dp), dimension(:,:), allocatable :: g_kc_bd_CC2
      real(dp), dimension(:,:), allocatable :: g_kd_bc
      real(dp), dimension(:,:), allocatable :: g_kb_cd
!
      real(dp), dimension(:,:), allocatable :: L_kc_bd
      real(dp), dimension(:,:), allocatable :: L_kcb_d
      real(dp), dimension(:,:), allocatable :: L_ck_bd
! 
      real(dp), dimension(:,:), allocatable :: x_cd_ij 
      real(dp), dimension(:,:), allocatable :: x_aj_kd 
      real(dp), dimension(:,:), allocatable :: x_ai_kc 
      real(dp), dimension(:,:), allocatable :: x_d_aij
!
      real(dp), dimension(:,:), allocatable :: I_bd 
      real(dp), dimension(:,:), allocatable :: I_kb_ij 
      real(dp), dimension(:,:), allocatable :: I_kcb_j 
      real(dp), dimension(:,:), allocatable :: I_kdb_i 
!
      real(dp), dimension(:,:), allocatable :: c_ck
!
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_CCSD 
      real(dp), dimension(:,:), allocatable :: rho_b_aij 
      real(dp), dimension(:,:), allocatable :: rho_a_bij 
      real(dp), dimension(:,:), allocatable :: rho_aj_bi 
!
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_b_length = 0, batch_dimension = 0
      integer(i15) :: n_batch = 0, b_first = 0, b_last = 0, b_batch = 0, b_length = 0
!
!     Indices 
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0
!
      integer(i15) :: ai = 0, aj = 0, ak = 0, ai_CCSD = 0, ai_CC2 = 0
      integer(i15) :: bi = 0, bj = 0, bd = 0, bc = 0, bj_CCSD = 0
      integer(i15) :: cd = 0, ci = 0, ck = 0
      integer(i15) :: dj = 0
      integer(i15) :: ij = 0
      integer(i15) :: kc = 0, kb = 0, kd = 0, kc_CC2
!
      integer(i15) :: aij = 0, bij = 0, kcb = 0
!
      integer(i15) :: akdj = 0, aidj = 0, akci = 0, cidj = 0, aick = 0
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: first_CCSD_o ! first active occupied index 
      integer(i15) :: first_CCSD_v ! first active virtual index
!
      integer(i15) :: last_CC2_o ! first active occupied index 
      integer(i15) :: last_CC2_v ! first active virtual index
      integer(i15) :: last_CCSD_o ! first active occupied index 
      integer(i15) :: last_CCSD_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!     Determine batch size, etc.
!     (Redo estimate once loop is done)
!
      required = max(2*(n_CCSD_v)*(wf%n_o)*(wf%n_J) +    &
                     2*(wf%n_J)*(wf%n_v)*(n_CCSD_v) +    &
                     (n_CC2_o**2)*(n_CC2_v**2),          &  ! Construction of L_bc^J 
                     (wf%n_J)*(wf%n_v)*(n_CCSD_v) +      &
                     (wf%n_o)*((wf%n_v)**2)*(n_CCSD_v) + &
                     (n_CC2_o**2)*(n_CC2_v**2))             ! Holding L_bc^J and g_aibc
!
      required = 4*required ! Words
!
      batch_dimension  = n_CCSD_v ! Batch over the virtual index b
      max_b_length = 0        ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, wf%mem%available, max_b_length, n_batch, batch_dimension)  
!
      do b_batch = 1, n_batch 
!
!        Get batching limits 
!
         call batch_limits(b_first, b_last, b_batch, max_b_length, batch_dimension)
         b_length = b_last - b_first + 1 
!
!        Construct g_kcbd for the batch
!
         call wf%mem%alloc(g_kc_bd, (wf%n_o)*(wf%n_v), (wf%n_v)*b_length)
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_vv(integral_type, g_kc_bd, &
                           1, wf%n_o,              &
                           1, wf%n_v,              &
                           b_first, b_last,        &
                           1, wf%n_v)
!
!        :: Term 1 ::
!        - sum_kcd g_kcbd x_ij^cd c_ak 
!
!        k           - general index
!        c, d        - CC2 indices
!        a, b, i, j  - CCSD indices 
!
!        Construct and order x_ij^cd as x_cd_ij 
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_cd_ij, n_CC2_v**2, n_CCSD_o**2)
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
!
               ij = index_two(i, j, n_CCSD_o)
!
               do d = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
!
                  do c = 1, n_CC2_v
!
                     ci = index_two(c, i, n_CC2_v)
                     cd = index_two(c, d, n_CC2_v)
                     cidj = index_packed(ci, dj)
!
                     x_cd_ij(cd, ij) = wf%x2am(cidj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        Reorder g_kc_bd to g_kb_cd with c and d restricted to the CC2 space 
!
         call wf%mem%alloc(g_kb_cd, (wf%n_o)*b_length, n_CC2_v**2)
!
         do c = 1, n_CC2_v
            do d = 1, n_CC2_v
!
               cd = index_two(c, d, n_CC2_v)
!
               do k = 1, wf%n_o
!
                  kc = index_two(k, c, wf%n_o)
!
                  do b = 1, b_length
!
                     bd = index_two(b, d, b_length)
                     kb = index_two(k, b, wf%n_o)
!
                     g_kb_cd(kb, cd) = g_kc_bd(kc, bd)
!
                  enddo
               enddo
            enddo
         enddo
!
!        I_kb_ij = sum(cd) g_kb_cd * x_cd_ij
!
         call wf%mem%alloc(I_kb_ij, (wf%n_o)*b_length, n_CCSD_o**2)
!
         call dgemm('N', 'N', &
                     (wf%n_o)*(b_length), &
                     n_CCSD_o**2, &
                     n_CC2_v**2, &
                     one, &
                     g_kb_cd, &
                     (wf%n_o)*(b_length), &
                     x_cd_ij, &
                     n_CC2_v**2, &
                     zero, &
                     I_kb_ij, &
                     (wf%n_o)*(b_length))
!
         call wf%mem%dealloc(g_kb_cd, (wf%n_o)*b_length, n_CC2_v**2)
         call wf%mem%dealloc(x_cd_ij, n_CC2_v**2, n_CCSD_o**2)
!
!        rho_a_bij = - sum_k c_a_k* I_kb_ij
!
         call wf%mem%alloc(rho_a_bij, n_CCSD_v, b_length*(n_CCSD_o**2))
!
         call dgemm('N', 'N',                      &
                     n_CCSD_v,                     &
                     b_length*(n_CCSD_o**2),       &
                     wf%n_o,                       &   
                     -one,                         &
                     c_a_i,                        & ! c_a_k
                     wf%n_v,                       &
                     I_kb_ij,                      & ! I_k_bij
                     wf%n_o,                       &
                     zero,                         &
                     rho_a_bij,                    &
                     n_CCSD_v)

         call wf%mem%dealloc(I_kb_ij, (wf%n_o)*b_length, n_CCSD_o**2)
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
!
                  do b = 1, b_length
!
                     bj = index_two(b + b_first - 1, j, n_CC2_v)
!
                     bij = index_three(b, i, j, b_length, n_CCSD_o)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_a_bij(a, bij)
!
                  enddo
               enddo
            enddo
         enddo 
!
         call wf%mem%dealloc(rho_a_bij, n_CCSD_v, b_length*(n_CCSD_o**2))
!
!        :: Term 2 ::
!        - sum_kcd g_kcbd x_kj^ad c_ci
!
!        c           - general index 
!        k, d        - CC2 indices 
!        a, b, i, j  - CCSD indices
!
!        Reorder g_kc_bd to g_kd_bc and constrain k and d to CC2 space
!
         call wf%mem%alloc(g_kd_bc, (n_CC2_o)*(n_CC2_v), b_length*(wf%n_v))
!
         do k = 1, n_CC2_o
            do d = 1, n_CC2_v
!
               kd = index_two(k, d, n_CC2_o)
!
               do b = 1, b_length
!
                  bd = index_two(b, d, b_length)
!
                  do c = 1, wf%n_v
!
                     bc = index_two(b, c, b_length)
                     kc = index_two(k, c, wf%n_o)
!
                     g_kd_bc(kd, bc) = g_kc_bd(kc, bd)
!
                  enddo
               enddo
            enddo
         enddo
!
!        I_kdb_i = sum_c g_kd_bc *c_ci
!
         call wf%mem%alloc(I_kdb_i, (b_length)*(n_CC2_o)*(n_CC2_v), (n_CCSD_o))
!
         call dgemm('N', 'N',                         &
                     (b_length)*(n_CC2_o)*(n_CC2_v),  &
                     (n_CCSD_o),                      &
                     wf%n_v,                          &
                     one,                             &
                     g_kd_bc,                         & ! g_kdb_c
                     (b_length)*(n_CC2_o)*(n_CC2_v),  &
                     c_a_i,                           & ! c_c_i
                     wf%n_v,                          &
                     zero,                            &  
                     I_kdb_i,                         &
                     (b_length)*(n_CC2_o)*(n_CC2_v))
!
         call wf%mem%dealloc(g_kd_bc, (n_CC2_o)*(n_CC2_v), b_length*(wf%n_v))
!
!        Construct and order  x_kj^ad to x_aj_kd
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_aj_kd, (n_CCSD_o)*(n_CCSD_v), (n_CC2_v)*(n_CC2_o))
!
         do j = 1, n_CCSD_o
            do k = 1, n_CC2_o
               do d = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
                  kd = index_two(k, d, n_CC2_o)
!
                  do a = 1, n_CCSD_v
!
                     aj = index_two(a, j, n_CCSD_v)
                     ak = index_two(a, k, n_CC2_v)
                     akdj = index_packed(ak, dj)
!
                     x_aj_kd(aj, kd) = wf%x2am(akdj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        rho_aj_bi = -sum_(kd) I_kdb_i * x_aj_kd
!
         call wf%mem%alloc(rho_aj_bi, (n_CCSD_o)*(n_CCSD_v), b_length*n_CCSD_o)
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     b_length*n_CCSD_o,      &
                     (n_CC2_v)*(n_CC2_o),    &
                     -one,                   &
                     x_aj_kd,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     I_kdb_i,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_aj_bi,              &
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(I_kdb_i, (b_length)*(n_CC2_o)*(n_CC2_v), (n_CCSD_o))
!
         call wf%mem%dealloc(x_aj_kd, (n_CCSD_o)*(n_CCSD_v), (n_CC2_v)*(n_CC2_o))
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  aj = index_two(a, j, n_CCSD_v)
!
                  do b = 1, b_length
!
                     bi = index_two(b, i, b_length)
                     bj = index_two(b + b_first - 1, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aj_bi(aj, bi) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_aj_bi, (n_CCSD_o)*(n_CCSD_v), b_length*n_CCSD_o)
!
!        :: Term 3 ::
!        - sum_kcd g_kcbd x_ik^ca c_dj
!
!        d           - general index 
!        k, c        - CC2 indices
!        a, b, i, j  - CCSD indices
!
!        Constrain k, c to CC2 space for g_kc_bd
!
         call wf%mem%alloc(g_kc_bd_CC2, (n_CC2_o)*(n_CC2_v), (b_length)*(wf%n_v))
!
         do k = 1, n_CC2_o
            do c = 1, n_CC2_v
!
               kc = index_two(k, c, wf%n_o)
               kc_CC2 = index_two(k, c, n_CC2_o)
!
               g_kc_bd_CC2(kc_CC2, :) = g_kc_bd(kc, :)
!
            enddo
         enddo
!
!        I_kcb_j = sum_(d) g_kc_bd_CC2 * c_d_j
!
         call wf%mem%alloc(I_kcb_j, (n_CC2_o)*(n_CC2_v)*b_length, n_CCSD_o)
!
         call dgemm('N', 'N',                      &
                     (n_CC2_o)*(n_CC2_v)*b_length, &
                     n_CCSD_o,                     &
                     wf%n_v,                       &
                     one,                          &
                     g_kc_bd_CC2,                  & ! g_kcb_d
                     (n_CC2_o)*(n_CC2_v)*b_length, &
                     c_a_i,                        &
                     wf%n_v,                       &
                     zero,                         &
                     I_kcb_j,                      &
                     (n_CC2_o)*(n_CC2_v)*b_length) 
!
         call wf%mem%dealloc(g_kc_bd_CC2, (n_CC2_o)*(n_CC2_v), (b_length)*(wf%n_v))
!
!        Construct and order x_ik^ca as x_ai_kc
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_ai_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_v)*(n_CC2_o))
!
         do i = 1, n_CCSD_o
            do k = 1, n_CC2_o
               do c = 1, n_CC2_v
!
                  ci = index_two(c, i, n_CC2_v)
                  kc = index_two(k, c, n_CC2_o)
!
                  do a = 1, n_CCSD_v
!
                     ai = index_two(a, i, n_CCSD_v)
                     ak = index_two(a, k, n_CC2_v)
                     akci = index_packed(ak, ci)
!
                     x_ai_kc(ai, kc) = wf%x2am(akci, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        rho_ai_bj_CCSD = -sum_(kc) x_ai_kc * I_kcb_j
!
         call wf%mem%alloc(rho_ai_bj_CCSD, n_CCSD_o*n_CCSD_v, (b_length)*(n_CCSD_o))
!
         call dgemm('N', 'N', &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (b_length)*(n_CCSD_o),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     -one,                   &
                     x_ai_kc,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     I_kcb_j,                & ! I_kc_bj
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_ai_bj_CCSD,         &
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(x_ai_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_v)*(n_CC2_o))
         call wf%mem%dealloc(I_kcb_j, (n_CC2_o)*(n_CC2_v)*b_length, n_CCSD_o)
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, b_length
!
                     bj_CCSD = index_two(b, j, b_length)
                     bj = index_two(b + b_first - 1, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD) 
!
                  enddo
               enddo
            enddo
         enddo

         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), b_length*n_CCSD_o)
!
!        Construct L_kcbd for the batch (constrain k to CC2 space)
!
         call wf%mem%alloc(L_kc_bd, (n_CC2_o)*(wf%n_v), (wf%n_v)*b_length)
!
         do c = 1, wf%n_v
            do k = 1, n_CC2_o
!
               kc = index_two(k, c, wf%n_o)
               kc_CC2 = index_two(k, c, n_CC2_o)
!
               do d = 1, wf%n_v
!
                  kd = index_two(k, d, wf%n_o)
!
                  do b = 1, b_length
!
                     bd = index_two(b, d, b_length)
                     bc = index_two(b, c, b_length)
!
                     L_kc_bd(kc_CC2, bd) = two*g_kc_bd(kc, bd) - g_kc_bd(kd, bc)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_kc_bd, (wf%n_o)*(wf%n_v), (wf%n_v)*b_length)
!
!       :: Term 4 :: 
!       sum_(kcd) L_kcbd x_ik^ac c_dj
!
!       d           - general index 
!       k, c        - CC2 indices 
!       a, b, i, j  - CCSD indices
!
         call wf%mem%alloc(L_kcb_d, n_CC2_o*n_CC2_v*b_length, wf%n_v)
!
         do c = 1, n_CC2_v
            do k = 1, n_CC2_o
!
               kc = index_two(k, c, n_CC2_o)
!
               do d = 1, wf%n_v
!
                  do b = 1, b_length
!
                     bd = index_two(b, d, b_length)
                     kcb = index_three(k, c, b, n_CC2_o, n_CC2_v)
!
                     L_kcb_d(kcb, d) = L_kc_bd(kc, bd)
!
                  enddo
               enddo
            enddo
         enddo
!
!        I_kc_bj = sum_(d) L_kc_bd *c_dj
!
         call wf%mem%alloc(I_kcb_j, (n_CC2_o)*(n_CC2_v)*(b_length), (n_CCSD_o))
!
         call dgemm('N', 'N',                         &
                     (n_CC2_o)*(n_CC2_v)*(b_length),  &
                     (n_CCSD_o),                      &
                     wf%n_v,                          &
                     one,                             &
                     L_kcb_d,                         & ! L_kcb_d
                     (n_CC2_o)*(n_CC2_v)*(b_length),  &
                     c_a_i,                           &
                     wf%n_v,                          &
                     zero,                            &
                     I_kcb_j,                         &
                     (n_CC2_o)*(n_CC2_v)*(b_length))
!
         call wf%mem%dealloc(L_kcb_d, n_CC2_o*n_CC2_v*b_length, wf%n_v)
!
!        Construct x_ik^ac ordered as x_ai_kc with
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_ai_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_v)*(n_CC2_o))
!
         do i = 1, n_CCSD_o
            do k = 1, n_CC2_o
               do c = 1, n_CC2_v
!
                  ck = index_two(c, k, n_CC2_v)
                  kc = index_two(k, c, n_CC2_o)
!
                  do a = 1, n_CCSD_v
!
                     ai = index_two(a, i, n_CCSD_v)
                     ai_CC2 = index_two(a, i, n_CC2_v)
!
                     aick = index_packed(ai_CC2, ck)
!
                     x_ai_kc(ai, kc) = wf%x2am(aick, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        rho_ai_bj_CCSD = x_ai_kc * I_kcb_j
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (b_length)*(n_CCSD_o))
!
         call dgemm('N', 'N', &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (b_length)*(n_CCSD_o),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     x_ai_kc,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     I_kcb_j,                & ! I_kc_bj
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_ai_bj_CCSD,         &
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(I_kcb_j, (n_CC2_o)*(n_CC2_v)*(b_length), (n_CCSD_o))
         call wf%mem%dealloc(x_ai_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_v)*(n_CC2_o))
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, b_length
!
                     bj_CCSD = index_two(b, j, b_length)
                     bj = index_two(b + b_first - 1, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), b_length*n_CCSD_o)
!
!        :: Term 5 ::
!        sum_kcd L_kcbd t_ij^ad c_ck
!
!        c           - general index 
!        k, d        - CC2 indices 
!        a, b, i, j  - CCSD indices
!
!
!        Reorder L_kc_bd to L_ck_bd and constrain d to CC2 space
!
         call wf%mem%alloc(L_ck_bd, (wf%n_v)*(n_CC2_o), (n_CC2_v)*(b_length))
!
         do k = 1, n_CC2_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
               kc = index_two(k, c, n_CC2_o)
!
               do d = 1, n_CC2_v
                  do b = 1, b_length
!
                     bd = index_two(b, d, b_length)
!
                     L_ck_bd(ck, bd) = L_kc_bd(kc, bd)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(L_kc_bd, (n_CC2_o)*(wf%n_v), (wf%n_v)*b_length)
!
!        I_bd = sum_(ck) L_ck_bd * c_ck
!
         call wf%mem%alloc(I_bd, 1, b_length*n_CC2_v)
!
         call dgemm('N', 'N',             &
                     1,                   &
                     b_length*n_CC2_v,    &
                     (n_CC2_o)*(wf%n_v),  &
                     one,                 &
                     c_a_i,               &
                     1,                   &
                     L_ck_bd,             &
                     (n_CC2_o)*(wf%n_v),  &
                     zero,                & 
                     I_bd,                &
                     1)
!
         call wf%mem%dealloc(L_ck_bd, (wf%n_v)*(n_CC2_o), (n_CC2_v)*(b_length))
!
!        Construct x_ij^ad ordered as x_d_aij
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_d_aij, (n_CC2_v), (n_CCSD_o**2)*(n_CCSD_v))
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do d = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
!
                  do a = 1, n_CCSD_v
!
                     ai  = index_two(a, i, n_CC2_v)
                     aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
                     aidj = index_packed(ai, dj)
!
                     x_d_aij(d, aij) = wf%x2am(aidj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        rho_b_aij = sum(d) I_bd * x_d_aij
!
         call wf%mem%alloc(rho_b_aij, b_length, (n_CCSD_v)*(n_CCSD_o**2))
!
         call dgemm('N', 'N',                   &
                     b_length,                  &
                     (n_CCSD_v)*(n_CCSD_o**2),  &
                     n_CC2_v,                   &
                     one,                       &
                     I_bd,                      & ! I_b_d
                     b_length,                  &
                     x_d_aij,                   &
                     n_CC2_v,                   &
                     zero,                      &
                     rho_b_aij,                 &
                     b_length)
!
         call wf%mem%dealloc(I_bd, 1, b_length*n_CC2_v)
         call wf%mem%dealloc(x_d_aij, (n_CC2_v), (n_CCSD_o**2)*(n_CCSD_v))
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
                  do b = 1, b_length
!
                     bj = index_two(b, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_aij(b, aij)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_b_aij, b_length, (n_CCSD_v)*(n_CCSD_o**2))
!
      enddo ! End of batches over b 
!
!
   end subroutine jacobian_mlccsd_d2_mlccsd
!
!
   module subroutine jacobian_mlccsd_e2_mlccsd(wf, rho_ai_bj, c_ai_ck)
!!
!!    Jacobian MLCCSD E2 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    rho_ai_bj^E2 = 2 sum_dlck x_bj,DL * L_KC,LD * c_ai,CK 
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(:,:) :: rho_ai_bj
      real(dp), dimension(:,:) :: c_ai_ck
!
      real(dp), dimension(:,:), allocatable :: x_dl_bj
      real(dp), dimension(:,:), allocatable :: g_kc_ld
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: I_ck_bj
      real(dp), dimension(:,:), allocatable :: c_ai_ck_CCSD
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_CCSD
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0, i = 0, j = 0, k = 0, l = 0
      integer(i15) :: ck = 0, dl = 0, kc = 0, kd = 0, lc = 0, ld = 0
      integer(i15) :: ai = 0, bj = 0, ai_CCSD = 0, bj_CCSD = 0
      integer(i15) :: dlbj = 0
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: first_CCSD_o ! first active occupied index 
      integer(i15) :: first_CCSD_v ! first active virtual index
!
      integer(i15) :: last_CC2_o ! first active occupied index 
      integer(i15) :: last_CC2_v ! first active virtual index
      integer(i15) :: last_CCSD_o ! first active occupied index 
      integer(i15) :: last_CCSD_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!     Read X2 amplitudes from disk
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(x_dl_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
      do l = 1, n_CC2_o
         do d = 1, n_CC2_v
!
            dl = index_two(d, l, n_CC2_v)
!
            do j = 1, n_CCSD_o
               do b = 1, n_CCSD_v
!
                  bj       = index_two(b, j, n_CC2_v)
                  bj_CCSD  = index_two(b, j, n_CCSD_v)
                  dlbj = index_packed(dl, bj)
!
                  x_dl_bj(dl, bj_CCSD) = wf%x2am(dlbj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_x2am
!
      call wf%mem%alloc(g_kc_ld, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld,  &
                        first_CC2_o, last_CC2_o, &
                        first_CC2_v, last_CC2_v, &
                        first_CC2_o, last_CC2_o, &
                        first_CC2_v, last_CC2_v)
!
!     Construct L_kc,ld ordered as L_ck_dl
!
      call wf%mem%alloc(L_ck_dl, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
      L_ck_dl = zero
!
      do l = 1, n_CC2_o
         do d = 1, n_CC2_v
!
            dl = index_two(d, l, n_CC2_v)
            ld = index_two(l, d, n_CC2_o)
!
            do k = 1, n_CC2_o
!
               kd = index_two(k, d, n_CC2_o)
!
               do c = 1, n_CC2_v
!
                  ck = index_two(c, k, n_CC2_v)
                  kc = index_two(k, c, n_CC2_o)
                  lc = index_two(l, c, n_CC2_o)
!
                  L_ck_dl(ck, dl) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     Intermediate I_ck_bj = sum_dl L_ck_dl * x_dl_bj
!
      call wf%mem%alloc(I_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
      call dgemm('N', 'N',                &
                  (n_CC2_o)*(n_CC2_v),    &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  one,                    &
                  L_ck_dl,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  x_dl_bj,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  zero,                   &
                  I_ck_bj,                &        
                  (n_CC2_o)*(n_CC2_v))
!
      call wf%mem%dealloc(x_dl_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
      call wf%mem%dealloc(L_ck_dl, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     rho_ai_bj = 2 * sum_ck c_ai_CK * I_CK_bj
!
      call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
      call wf%mem%alloc(c_ai_ck_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v)) 
!
      do a = 1, n_CCSD_v
         do i = 1, n_CCSD_o
            ai_CCSD = index_two(a, i, n_CCSD_v)
            ai = index_two(a, i, n_CC2_v)
            c_ai_ck_CCSD(ai_CCSD, :) = c_ai_ck(ai, :)
         enddo
      enddo
!
      call dgemm('N', 'N',                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  two,                    &
                  c_ai_ck_CCSD,           &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  I_ck_bj,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  zero,                   &
                  rho_ai_bj_CCSD,         &        
                  (n_CCSD_o)*(n_CCSD_v))
!
      call wf%mem%dealloc(I_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
      call wf%mem%dealloc(c_ai_ck_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v)) 
!
      do i = 1, n_CCSD_o
         do j = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CC2_v)
               ai_CCSD = index_two(a, i, n_CCSD_v)
!
               do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CC2_v)
                  bj_CCSD = index_two(b, j, n_CCSD_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
   end subroutine jacobian_mlccsd_e2_mlccsd
!
!
      module subroutine jacobian_mlccsd_f2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD F2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^F2 =   - sum_(CKDL) x_ai,CK * L_KC,LD * c_bL,Dj 
!!                        - sum_(CKDL) x_ai,Dj * L_KC,LD * c_bL,CK
!!                        - sum_(CKDL) x_ai_bL * L_KC,LD * c_CK,Dj
!!
!!
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!
         real(dp), dimension(:,:), allocatable :: L_KC_J
         real(dp), dimension(:,:), allocatable :: g_KC_LD
!
         real(dp), dimension(:,:), allocatable :: L_CK_DL
         real(dp), dimension(:,:), allocatable :: L_D_LCK
         real(dp), dimension(:,:), allocatable :: L_L_CKD
!
         real(dp), dimension(:,:), allocatable :: c_DL_bj
!
         real(dp), dimension(:,:), allocatable :: x_ai_CK
         real(dp), dimension(:,:), allocatable :: x_aij_D
         real(dp), dimension(:,:), allocatable :: x_aib_L
 !               
         real(dp), dimension(:,:), allocatable :: I_D_b
         real(dp), dimension(:,:), allocatable :: I_L_j
         real(dp), dimension(:,:), allocatable :: I_CK_bj
!
         real(dp), dimension(:,:), allocatable :: rho_aij_b
         real(dp), dimension(:,:), allocatable :: rho_aib_j
         real(dp), dimension(:,:), allocatable :: rho_ai_bj_CCSD
!
         integer(i15) :: a = 0, b = 0, c = 0, d = 0
         integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
         integer(i15) :: ai = 0, bj = 0, bk = 0, bl= 0, ck = 0, cl = 0, dj = 0, dl = 0
         integer(i15) :: ai_CC2 = 0, ai_CCSD = 0, bj_CCSD = 0
         integer(i15) :: kc = 0, kd = 0, ld = 0, lc = 0 
!
         integer(i15) :: aij = 0, aib = 0, lck = 0, ckd = 0
!
         integer(i15) :: bldj = 0, aidj = 0, bkcl = 0, aibl = 0, aick = 0
!
!        Active space variables
!
         integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
         integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
         integer(i15) :: first_CC2_o ! first active occupied index 
         integer(i15) :: first_CC2_v ! first active virtual index
!
         integer(i15) :: last_CC2_o ! first active occupied index 
         integer(i15) :: last_CC2_v ! first active virtual index
!
!        Calculate first/last indeces
! 
         call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
!
         call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
         call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
         last_CC2_o = first_CC2_o + n_CC2_o - 1
         last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
!        :: Construct L_kc,ld ::
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!   
         call wf%mem%alloc(L_CK_DL, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
         L_CK_DL = zero
!
!        Construct L_KC,LD ordered as L_CK_DL
!             
         do C = 1, n_CC2_v   
            do K = 1, n_CC2_o
!
               KC = index_two(k, c, n_CC2_o)
               CK = index_two(c, k, n_CC2_v)
!
               do D = 1, n_CC2_v
                  do L = 1, n_CC2_o
!
                     LC = index_two(L, C, n_CC2_o)
                     LD = index_two(L, D, n_CC2_o)
                     DL = index_two(D, L, n_CC2_v)
                     KD = index_two(K, D, n_CC2_o)
!
                     L_CK_DL(CK, DL) = two*g_KC_LD(KC, LD) - g_KC_LD(KD, LC)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        :: Term 1: - sum_CKLD t_ai,CK * L_KC,LD * c_bL,Dj ::
!
!        Reorder c_bL_Dj as c_DL_bj, 
!
         call wf%mem%alloc(c_DL_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         c_DL_bj = zero
!
         do L = 1, n_CC2_o
            do j = 1, n_CCSD_o
               do D = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
                  dl = index_two(d, l, n_CC2_v)
!
                  do b = 1, n_CCSD_v   
!
                     bj = index_two(b, j, n_CCSD_v)    
                     bl = index_two(b, l, n_CC2_v)
!  
                     c_dl_bj(dl, bj) = c_ai_bj(bl, dj)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(I_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!        I_ck_bj = sum_dl L_CK_DL*c_DL_bj = sum_DL L_KC,LD*c_bL,Dj
!
         call dgemm('N', 'N',                &  
                     (n_CC2_o)*(n_CC2_v),    &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     L_ck_dl,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     c_dl_bj,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     I_ck_bj,                &        
                     (n_CC2_o)*(n_CC2_v))
!
         call wf%mem%dealloc(c_DL_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(L_CK_DL, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_ai_ck, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         do i = 1, n_CCSD_o   
            do a = 1, n_CCSD_v
!
               AI_CC2   = index_two(a, i, n_CC2_v)
               ai       = index_two(a, i, n_CCSD_v)
!
               do K = 1, n_CC2_o
                  do C = 1, n_CC2_v
! 
                     CK = index_two(C, K, n_CC2_v)
                     AICK = index_packed(AI_CC2, CK)
!
                     x_ai_ck(ai, ck) = wf%x2am(AICK, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        rho_ai_bj = sum_ck t_ai_ck*I_ck_bj
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N', 'N',                &
                    (n_CCSD_o)*(n_CCSD_v),   &
                    (n_CCSD_o)*(n_CCSD_v),   &
                    (n_CC2_o)*(n_CC2_v),     &
                    -one,                    &
                    x_ai_ck,                 &
                    (n_CCSD_o)*(n_CCSD_v),   &
                    I_ck_bj,                 &
                    (n_CC2_o)*(n_CC2_v),     &
                    zero,                    &
                    rho_ai_bj_CCSD,          &        
                    (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(I_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(x_ai_ck, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
!        :: Term 2: - sum_ckdl x_ai,dj * L_kc,ld * c_bl,ck
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!
!        Construct L_ck,dl reordered as L_d_clk
!
         call wf%mem%alloc(L_D_LCK, n_CC2_v, (n_CC2_v)*((n_CC2_o)**2))
!
         do K = 1, n_CC2_o
            do L = 1, n_CC2_o
               do C = 1, n_CC2_v
!
                  LCK = index_three(L, C, K, n_CC2_o, n_CC2_v)
!
                  KC = index_two(K, C, n_CC2_o)
                  LC = index_two(L, C, n_CC2_o)
!
                  do D = 1, n_CC2_v
!
                     LD = index_two(L, D, n_CC2_o)
                     KD = index_two(K, D, n_CC2_o)
!
                     L_D_LCK(D, LCK) = two*g_KC_LD(KC, LD) - g_KC_LD(KD, LC)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        I_d_b = sum_clk L_d_LCK * c_b_LCK
!        Here dgemm is tricked to believe that c_bL_CK is c_b_LCK
! 
         call wf%mem%alloc(I_D_b, n_CC2_v, n_CCSD_v)
!
         call dgemm('N', 'T',                   &
                     n_CC2_v,                   &
                     n_CCSD_v,                  &
                     ((n_CC2_o)**2)*(n_CC2_v),  &
                     one,                       &
                     L_d_lck,                   &
                     n_CC2_v,                   &
                     c_ai_bj,                   & ! c_b_lck
                     n_CC2_v,                   &  
                     zero,                      &
                     I_d_b,                     &
                     n_CC2_v)
!
         call wf%mem%dealloc(L_D_LCK, n_CC2_v, (n_CC2_v)*((n_CC2_o)**2))
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_aij_d, (n_CCSD_v)*((n_CCSD_o)**2), n_CC2_v)
!
!        Reorder X2 amplitudes
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
               do D = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
!
                  do a = 1, n_CCSD_v
!
                     ai = index_two(a, i, n_CC2_v)
!
                     aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
                     aidj = index_packed(ai, dj)
!
                     x_aij_D(aij, D) = wf%x2am(aidj,1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
         call wf%mem%alloc(rho_aij_b, (n_CCSD_v)*((n_CCSD_o)**2), n_CCSD_v)
!
!        rho_aij_b = sum_d t_aij_d*Y_d_b
!
         call dgemm('N','N',                       &
                     (n_CCSD_v)*((n_CCSD_o)**2),   &
                     n_CCSD_v,                     &
                     n_CC2_v,                      &
                     -one,                         &
                     x_aij_d,                      &
                     (n_CCSD_v)*((n_CCSD_o)**2),   &
                     I_d_b,                        &
                     n_CC2_v,                      &
                     zero,                         &
                     rho_aij_b,                    &
                     (n_CCSD_v)*((n_CCSD_o)**2))
!
         call wf%mem%dealloc(x_aij_d, (n_CCSD_v)*((n_CCSD_o)**2), n_CC2_v)
         call wf%mem%dealloc(I_D_b, n_CC2_v, n_CCSD_v)
!
!        Adding term 2 to rho_ai_bj
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
                  ai = index_two(a, i, n_CC2_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_aij_b, (n_CCSD_v)*((n_CCSD_o)**2), n_CCSD_v)
!
!        :: Term 3: - sum_(CKDL) x_ai,bL * L_KC,LD * c_CK,Dj ::
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!   
         call wf%mem%alloc(L_L_CKD, (n_CC2_o), (n_CC2_o)*((n_CC2_v)**2))
!
!        Construct L_kc,dl ordered as L_l_ckd
!             
         do C = 1, n_CC2_v
            do K = 1, n_CC2_o
!
               KC = index_two(K, C, n_CC2_o)
!
               do D = 1, n_CC2_v
!
                  CKD = index_three(C, K, D, n_CC2_v, n_CC2_o)
!
                  do L = 1, n_CC2_o
!
                     LC = index_two(L, C, n_CC2_o)
                     LD = index_two(L, D, n_CC2_o)
                     KD = index_two(K, D, n_CC2_o)
!
                     L_L_CKD(L, CKD) = two*g_KC_LD(KC, LD) - g_KC_LD(KD, LC)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         call wf%mem%alloc(I_L_j, n_CC2_o, n_CCSD_o)
!
!        I_L_j = sum_(CKD) L_L_CKD * c_CKD_j 
!  
         call dgemm('N', 'N',                   &
                     n_CC2_o,                   &
                     n_CCSD_o,                  &
                     ((n_CC2_v)**2)*(n_CC2_o),  &
                     one,                       &
                     L_L_CKD,                   &
                     n_CC2_o,                   &
                     c_ai_bj,                   & ! c_ai_bj(ck,dl)= c_ckd_l
                     ((n_CC2_v)**2)*(n_CC2_o),  &
                     zero,                      &
                     I_L_j,                     &
                     n_CC2_o)
!
         call wf%mem%dealloc(L_L_CKD, (n_CC2_o), (n_CC2_o)*((n_CC2_v)**2))
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_aib_L, (n_CCSD_o)*((n_CCSD_v)**2), n_CC2_o)
!
         do L = 1, n_CC2_o
            do i = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
!
                  do b = 1, n_CCSD_v
!
                     bl = index_two(b, l, n_CC2_v)
!
                     aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
                     aibl = index_packed(ai, bl)
!
                     x_aib_L(aib, L) = wf%x2am(aibl, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        rho_ai_bj_CCSD = sum_L x_aib_L * I_L_j
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N','N',                    &
                     (n_CCSD_o)*(n_CCSD_v**2),  &
                     n_CCSD_o,                  &
                     n_CC2_o,                   &
                     -one,                      &
                     x_aib_L,                   &
                     (n_CCSD_o)*(n_CCSD_v**2),  &
                     I_l_j,                     &
                     n_CC2_o,                   &
                     zero,                      &
                     rho_ai_bj_CCSD,            &
                     (n_CCSD_o)*(n_CCSD_v**2))
!
         call wf%mem%dealloc(x_aib_L, (n_CCSD_o)*((n_CCSD_v)**2), n_CC2_o)
         call wf%mem%dealloc(I_L_j, n_CC2_o, n_CCSD_o)
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
      end subroutine jacobian_mlccsd_f2_mlccsd
!
!
      module subroutine jacobian_mlccsd_g2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD G2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^G2 =  - sum_ckdl x_bL,Dj * L_KC,LD * c_ai,CK 
!!                       - sum_ckdl x_CK_bL * L_KC,LD * c_ai,Dj 
!!                       - sum_ckld x_CK,Dj * L_KC,LD * c_ai,bL 
!!
!!
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!         
         real(dp), dimension(:,:), allocatable :: g_KC_LD
!
         real(dp), dimension(:,:), allocatable :: L_CK_DL
         real(dp), dimension(:,:), allocatable :: L_D_CLK
         real(dp), dimension(:,:), allocatable :: L_L_CKD
!
         real(dp), dimension(:,:), allocatable :: c_ai_ck
         real(dp), dimension(:,:), allocatable :: c_aib_l
         real(dp), dimension(:,:), allocatable :: c_aij_d
!
         real(dp), dimension(:,:), allocatable :: x_dl_bj
         real(dp), dimension(:,:), allocatable :: x_clk_b
         real(dp), dimension(:,:), allocatable :: x_ckd_j
 !               
         real(dp), dimension(:,:), allocatable :: I_ck_bj
         real(dp), dimension(:,:), allocatable :: I_d_b
         real(dp), dimension(:,:), allocatable :: I_l_j
!
         real(dp), dimension(:,:), allocatable :: rho_aij_b
         real(dp), dimension(:,:), allocatable :: rho_ai_bj_CCSD
!
         integer(i15) :: a = 0, b = 0, c = 0, d = 0
         integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
         integer(i15) :: ai = 0, bj = 0, bk = 0, bl = 0, ck = 0, cl = 0, dj = 0, dl = 0
         integer(i15) :: kc = 0, lc = 0, kd = 0, ld = 0
         integer(i15) :: ai_CCSD = 0, ai_CC2 = 0, bj_CCSD = 0
!
         integer(i15) :: aib = 0, aij = 0, ckd = 0, clk = 0
!
         integer(i15) :: ckbl = 0, ckdj = 0, bldj = 0
!
!        Active space variables
!
         integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
         integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
         integer(i15) :: first_CC2_o ! first active occupied index 
         integer(i15) :: first_CC2_v ! first active virtual index
!
         integer(i15) :: last_CC2_o ! first active occupied index 
         integer(i15) :: last_CC2_v ! first active virtual index
!
!        Calculate first/last indeces
! 
         call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
!
         call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
         call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
         last_CC2_o = first_CC2_o + n_CC2_o - 1
         last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
!        :: Term 1: - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck  ::
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!   
         call wf%mem%alloc(L_CK_DL, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        Construct L_kc_ld ordered as L_ck_dl
!             
         do C = 1, n_CC2_v
            do K = 1, n_CC2_o
!
               KC = index_two(K, C, n_CC2_o)
               CK = index_two(C, K, n_CC2_v)
!
               do D = 1, n_CC2_v
                  do L = 1, n_CC2_o
!
                     LC = index_two(L, C, n_CC2_o)
                     LD = index_two(L, D, n_CC2_o)
                     DL = index_two(D, L, n_CC2_v)
                     KD = index_two(K, D, n_CC2_o)
!
                     L_CK_DL(CK, DL) = two*g_KC_LD(KC, LD) - g_KC_LD(KD, LC)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        Reorder t_bl_dj as t_dl_bj
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_DL_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
         do L = 1, n_CC2_o
            do j = 1, n_CCSD_o
               do D = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
                  dl = index_two(d, l, n_CC2_v)
!
                  do b = 1, n_CCSD_v    
!
                     bj = index_two(b, j, n_CCSD_v)    
                     bl = index_two(b, l, n_CC2_v)
!  
                     bldj = index_packed(bl, dj)
!
                     x_dl_bj(dl, bj) = wf%x2am(bldj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
         call wf%mem%alloc(I_CK_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!        I_CK_bj = sum_dl x_bL,Dj * L_KC,LD = sum_dl L_CK_DL t_DL_bj 
!
         call dgemm('N', 'N',                &
                     (n_CC2_o)*(n_CC2_v),    &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     L_CK_DL,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     x_DL_bj,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     I_CK_bj,                &         
                     (n_CC2_o)*(n_CC2_v))
!
         call wf%mem%dealloc(x_DL_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(L_CK_DL, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        rho_ai_bj =+ - sum_CK c_ai,CK X_CK_bj
!
         call wf%mem%alloc(c_ai_CK, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         do a = 1, n_CCSD_v
            do i = 1, n_CCSD_o
!
               ai     = index_two(a, i, n_CCSD_v)
               ai_CC2 = index_two(a, i, n_CC2_v)
!
               c_ai_CK(ai,:) = c_ai_bj(ai_CC2, :)
!
            enddo
         enddo
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N', 'N',                &
                    (n_CCSD_o)*(n_CCSD_v),   &
                    (n_CCSD_o)*(n_CCSD_v),   &
                    (n_CC2_o)*(n_CC2_v),     &
                    -one,                    &
                    c_ai_ck,                 &
                    (n_CCSD_o)*(n_CCSD_v),   &
                    I_ck_bj,                 &
                    (n_CC2_o)*(n_CC2_v),     &
                    zero,                    &
                    rho_ai_bj_CCSD,          &        
                    (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(c_ai_CK, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
         call wf%mem%dealloc(I_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!        Add to rho_ai_bj
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
!        :: Term 2: - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!
!        Reorder L_ck_dl to L_d_clk
!
         call wf%mem%alloc(L_D_CLK, n_CC2_v, (n_CC2_v)*((n_CC2_o)**2))
!
         do K = 1, n_CC2_o
            do L = 1, n_CC2_o
               do C = 1, n_CC2_v
!
                  CLK = index_three(C, L, K, n_CC2_v, n_CC2_o)
!
                  KC = index_two(K, C, n_CC2_o)
                  LC = index_two(L, C, n_CC2_o)
!
                  do D = 1, n_CC2_v
!
                     LD = index_two(L, D, n_CC2_o)
                     KD = index_two(K, D, n_CC2_o)
!
                     L_D_CLK(D, CLK) = two*g_KC_LD(KC, LD) - g_KC_LD(KD, LC)
!
                  enddo
               enddo
            enddo
         enddo
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        Reorder x_ck,bl as x_clk_b
!        
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_CLK_b, (n_CC2_v)*((n_CC2_o)**2), n_CCSD_v)
!
         do K = 1, n_CC2_o
            do L = 1, n_CC2_o
               do C = 1, n_CC2_v
!  
                  CK = index_two(C, K, n_CC2_v)
!
                  CLK = index_three(C, L, K, n_CC2_v, n_CC2_o)
!
                  do b = 1, n_CCSD_v
!
                     BL = index_two(b, l, n_CC2_v)
!
                     CKBL = index_packed(CK, BL)
!
                     x_CLK_b(CLK, b) = wf%x2am(ckbl, 1)
!
                 enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
!        I_D_b = sum_CLK L_D_CLK * c_CLK_b 
!
         call wf%mem%alloc(I_D_b, n_CC2_v, n_CCSD_v)
!
         call dgemm('N', 'N',                &
                     n_CC2_v,                &
                     n_CCSD_v,               &
                     ((n_CC2_o)**2)*n_CC2_v, &
                     one,                    &
                     L_D_CLK,                &
                     n_CC2_v,                &
                     x_CLK_b,                &
                     ((n_CC2_o)**2)*n_CC2_v, &
                     zero,                   &
                     I_D_b,                  &
                     n_CC2_v)
!
         call wf%mem%dealloc(x_CLK_b, (n_CC2_v)*((n_CC2_o)**2), n_CCSD_v)
         call wf%mem%dealloc(L_D_CLK, n_CC2_v, (n_CC2_v)*((n_CC2_o)**2)) 
!
         call wf%mem%alloc(c_aij_D, (n_CCSD_v)*((n_CCSD_o)**2), n_CC2_v)
!
!        Reorder c_ai_dj to c_aij_d
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
               do d = 1, n_CC2_v
!
                  dj = index_two(d, j, n_CC2_v)
!
                  do a = 1, n_CCSD_v
!
                     ai = index_two(a, i, n_CC2_v)
!
                     aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
                     c_aij_d(aij, d) = c_ai_bj(ai, dj)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(rho_aij_b, (n_CCSD_v)*((n_CCSD_o)**2), n_CCSD_v)
!
!        rho_aij_b = sum_d c_aij_d * I_d_b
!
         call dgemm('N','N',                       &
                     (n_CCSD_v)*((n_CCSD_o)**2),   &
                     n_CCSD_v,                     &
                     n_CC2_v,                      &
                     -one,                         &
                     c_aij_D,                      &
                     (n_CCSD_v)*((n_CCSD_o)**2),   &
                     I_d_b,                        &
                     n_CC2_v,                      &
                     zero,                         &
                     rho_aij_b,                    &
                     (n_CCSD_v)*((n_CCSD_o)**2))
!
         call wf%mem%dealloc(c_aij_D, (n_CCSD_v)*((n_CCSD_o)**2), n_CC2_v)
         call wf%mem%dealloc(I_D_b, n_CC2_v, n_CCSD_v)
!
!        Adding term 2 to rho_ai_bj
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  aij = index_three(a, i, j, n_CCSD_v, n_CCSD_o)
!
                  ai = index_two(a, i, n_CC2_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_aij_b, (n_CCSD_v)*((n_CCSD_o)**2), n_CCSD_v)
!
!        :: Term 3: - sum_ckld x_ck,dj * L_kc,ld * c_ai,bl ::
!

!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!   
         call wf%mem%alloc(L_L_CKD, n_CC2_o, (n_CC2_o)*((n_CC2_v)**2))
!
!        Construct L_kc_ld ordered as  L_l_ckd
!                 
         do C = 1, n_CC2_v
            do K = 1, n_CC2_o
!
               KC = index_two(K, C, n_CC2_o)
!
               do D = 1, n_CC2_v
!
                  CKD = index_three(C, K, D, n_CC2_v, n_CC2_o)
!
                  do L = 1, n_CC2_o
!
                     LC = index_two(L, C, n_CC2_o)
                     LD = index_two(L, D, n_CC2_o)
                     KD = index_two(K, D, n_CC2_o)
!
                     L_L_CKD(L, CKD) = two*g_KC_LD(KC, LD) - g_KC_LD(KD, LC)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        Reorder x_ck,dj to x_ckd_j 
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_ckd_j, (n_CC2_v**2)*(n_CC2_o), n_CCSD_o)
!
         do K = 1, n_CC2_o
            do C = 1, n_CC2_v
!               
               CK = index_two(C, K, n_CC2_v)
!             
               do D = 1, n_CC2_v
! 
                  CKD = index_three(C, K, D, n_CC2_v, n_CC2_o)
!
                  do j = 1, n_CCSD_o
!
                     DJ = index_two(D, j, n_CC2_v)
!
                     CKDJ = index_packed(CK, DJ)
!
                     x_CKD_j(CKD, j) = wf%x2am(CKDJ, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
         call wf%mem%alloc(I_L_j, n_CC2_o, n_CCSD_o)
!
!        I_L_j = sum_ckd L_L_CKD*x_CKD_j
!
         call dgemm('N', 'N',                &
                     n_CC2_o,                &
                     n_CCSD_o,               &
                     ((n_CC2_v)**2)*n_CC2_o, &
                     one,                    &
                     L_l_ckd,                &
                     n_CC2_o,                &
                     x_ckd_j,                &
                     ((n_CC2_v)**2)*n_CC2_o, &
                     zero,                   &
                     I_l_j,                  &
                     n_CC2_o)
!
         call wf%mem%dealloc(L_L_CKD, n_CC2_o, (n_CC2_o)*((n_CC2_v)**2))
         call wf%mem%dealloc(x_ckd_j, (n_CC2_v**2)*(n_CC2_o), n_CCSD_o)
!
!        rho_aib_j = sum_L c_aib_L*I_L_j
!
         call wf%mem%alloc(c_aib_L, n_CCSD_o*(n_CCSD_v**2), n_CC2_o)
!
         do L = 1, n_CC2_o
            do b = 1, n_CCSD_v
!
               BL = index_two(B, L, n_CC2_v)
!
               do i = 1, n_CCSD_o
                  do a = 1, n_CCSD_v
!
                     AI = index_two(A, I, n_CC2_v)
!
                     aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
                     c_aib_L(aib, L) = c_ai_bj(AI, BL)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N','N',                 &
                     n_CCSD_o*(n_CCSD_v**2), &
                     n_CCSD_o,               &
                     n_CC2_o,                &
                     -one,                   &
                     c_aib_L,                &
                     n_CCSD_o*(n_CCSD_v**2), &
                     I_L_j,                  &
                     n_CC2_o,                &
                     zero,                   &
                     rho_ai_bj_CCSD,         &
                     n_CCSD_o*(n_CCSD_v**2))
!
         call wf%mem%dealloc(c_aib_L, n_CCSD_o*(n_CCSD_v**2), n_CC2_o)
         call wf%mem%dealloc(I_L_j, wf%n_o, wf%n_o)
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
! 
      end subroutine jacobian_mlccsd_g2_mlccsd
!
!
      module subroutine jacobian_mlccsd_h2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD H2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^H2 =  sum_CKDL x_Ci,aK * g_KC,LD * c_bL,Dj 
!!                     + sum_CKDL x_Cj,aL * g_KC,LD * c_bK,Di
!!                
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!         
         real(dp), dimension(:,:), allocatable :: g_kc_ld
         real(dp), dimension(:,:), allocatable :: g_lc_kd
!
         real(dp), dimension(:,:), allocatable :: x_ai_kc
         real(dp), dimension(:,:), allocatable :: x_aj_lc
!
         real(dp), dimension(:,:), allocatable :: c_ld_bj
         real(dp), dimension(:,:), allocatable :: c_kd_bi
!
         real(dp), dimension(:,:), allocatable :: I_ai_ld
         real(dp), dimension(:,:), allocatable :: I_aj_kd
!
         real(dp), dimension(:,:), allocatable :: rho_aj_bi
         real(dp), dimension(:,:), allocatable :: rho_ai_bj_CCSD
!
         integer(i15) :: a = 0, b = 0, c = 0, d = 0
         integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
         integer(i15) :: ak = 0, ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, bl = 0, ci = 0, cj = 0, di = 0, dj = 0
         integer(i15) :: kc = 0, kd = 0, ld = 0, lc = 0
         integer(i15) :: ai_CCSD = 0, bj_CCSD = 0
         integer(i15) :: akci = 0, alcj = 0
!
!        Active space variables
!
         integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
         integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
         integer(i15) :: first_CC2_o ! first active occupied index 
         integer(i15) :: first_CC2_v ! first active virtual index
!
         integer(i15) :: last_CC2_o ! first active occupied index 
         integer(i15) :: last_CC2_v ! first active virtual index
!
!        Calculate first/last indeces
! 
         call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
!
         call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
         call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
         last_CC2_o = first_CC2_o + n_CC2_o - 1
         last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
!        :: Term 1: sum_ckld t_ci,ak * g_kc,ld * c_bl,dj ::
!
!        Construct g_kc_ld
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!
!        x_ak,ci ordered as x_ai_kc
!  
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_ai_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         do i = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CCSD_v)
!
               do c = 1, n_CC2_v
                  do k = 1, n_CC2_o
!
                     CI = index_two(c, i, n_CC2_v)
                     AK = index_two(a, k, n_CC2_v)
                     KC = index_two(k, c, n_CC2_o)
!
                     AKCI = index_packed(AK, CI)
!
                     x_ai_KC(ai, KC) = wf%x2am(AKCI, 1)
! 
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!  
         call wf%mem%alloc(I_ai_ld, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!        I_ai_ld = sum_ck t_ai_kc*g_kc_ld
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     x_ai_KC,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     g_kc_LD,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     I_ai_LD,                &   
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(x_ai_KC, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         call wf%mem%alloc(c_LD_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!        Reorder c_bl,dj as c_ld_bj
!
         do L = 1, n_CC2_o
            do b = 1, n_CCSD_v
!
               BL = index_two(b, l, n_CC2_v)
!
               do j = 1, n_CCSD_o
!
                  bj = index_two(b, j, n_CCSD_v)
!
                  do D = 1, n_CC2_v
!
                     DJ = index_two(D, J, n_CC2_v)
                     LD = index_two(l, d, n_CC2_o)
!
                     c_LD_bj(LD, bj) = c_ai_bj(BL, DJ)
!
                  enddo
               enddo
            enddo
         enddo
!
!        rho_ai_bj += sum_ld X_ai_ld*c_ld_bj
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     I_ai_ld,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     c_ld_bj,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_ai_bj_CCSD,         &     
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(c_LD_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(I_ai_ld, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
!        :: Term 2: sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!
!        Construct g_kc_ld
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!
!        Reorder g_kc_ld to g_lc_kd 
!
         call wf%mem%alloc(g_LC_KD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         do C = 1, n_CC2_v
            do D = 1, n_CC2_v
               do K = 1, n_CC2_o
!
                  KC = index_two(k, c, n_CC2_o)
                  KD = index_two(k, d, n_CC2_o)
!
                  do L = 1, n_CC2_o
!
                     LC = index_two(L, C, n_CC2_o)
                     LD = index_two(L, D, n_CC2_o)
!
                     g_LC_KD(LC, KD) = g_KC_LD(KC, LD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!       x_al,cj ordered as x_aj_lc
!  
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_aj_lc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         do j = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               aj = index_two(a, j, n_CCSD_v)
!
               do C = 1, n_CC2_v
!
                  CJ = index_two(C, J, n_CC2_v)
!
                  do L = 1, n_CC2_o
!
                     AL = index_two(A, L, n_CC2_v)
                     LC = index_two(L, C, n_CC2_o)
!
                     ALCJ = index_packed(AL, CJ)
!
                     x_aj_LC(aj, LC) = wf%x2am(ALCJ, 1)
! 
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_x2am
!
         call wf%mem%alloc(I_aj_kd, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!        I_aj_kd = sum_lc x_aj_lc * g_lc_kd
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     x_aj_lc,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     g_lc_kd,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &  
                     I_aj_kd,                &
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(g_LC_KD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
         call wf%mem%dealloc(x_aj_lc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!  
!        Reorder c_bk,di as c_kd_bi
!
         call wf%mem%alloc(c_kd_bi, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
         do i = 1, n_CCSD_o
            do k = 1, n_CC2_o
               do d = 1, n_CC2_v
!
                  kd = index_two(k, d, n_CC2_o)
                  di = index_two(d, i, n_CC2_v)
!
                  do b = 1, n_CCSD_v
!
                     bk = index_two(b, k, n_CC2_v)
                     bi = index_two(b, i, n_CCSD_v)
!
                     c_kd_bi(kd, bi) = c_ai_bj(bk, di)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(rho_aj_bi, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
!        rho_aj_bi = sum_kd  Y_aj_kd * c_kd_bi
!
         call dgemm('N', 'N',                & 
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     I_aj_kd,                &  
                     (n_CCSD_o)*(n_CCSD_v),  &
                     c_kd_bi,                &
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_aj_bi,              &   
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(c_kd_bi, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(I_aj_kd, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!        Reorder into rho_ai_bj
!
         do i = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CC2_v)
!
               do j = 1, n_CCSD_o
!
                  aj = index_two(a, j, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bi = index_two(b, i, n_CCSD_v)
                     bj = index_two(b, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aj_bi(aj, bi) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_aj_bi, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
      end subroutine jacobian_mlccsd_h2_mlccsd
!
!
      module subroutine jacobian_mlccsd_i2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD I2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^I2 =  sum_C F_bC * c_ai,Cj - sum_K F_jK * c_ai,bK
!!                     + sum_ck L_bj,KC * c_ai,CK 
!!                     - sum_ck ( g_KC,bj * c_aK,Ci + g_Ki,bC * c_aK,Cj ) 
!!                
!!       Batch over c to construct  g_ki_bC
!! 
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!
         real(dp), dimension(:,:), allocatable :: c_aij_c 
         real(dp), dimension(:,:), allocatable :: c_aib_k
         real(dp), dimension(:,:), allocatable :: c_ai_ck
         real(dp), dimension(:,:), allocatable :: c_aj_ck
!
         real(dp), dimension(:,:), allocatable :: rho_aij_b
         real(dp), dimension(:,:), allocatable :: rho_aib_j
         real(dp), dimension(:,:), allocatable :: rho_aj_bi
         real(dp), dimension(:,:), allocatable :: rho_ai_bj_CCSD
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
         integer(i15) :: ai_CCSD = 0, ai_CC2 = 0, bj_CCSD = 0
!
         integer(i15) :: aij = 0, aib = 0
!
         integer(i15) :: offset = 0
!
         integer(i15) :: required = 0, available = 0, max_batch_length = 0, n_batch = 0, batch_dimension = 0
         integer(i15) :: c_batch = 0, c_first = 0, c_last = 0, c_length = 0
!
!        Active space variables
!
         integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
         integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
         integer(i15) :: first_CC2_o ! first active occupied index 
         integer(i15) :: first_CC2_v ! first active virtual index
         integer(i15) :: first_CCSD_o ! first active occupied index 
         integer(i15) :: first_CCSD_v ! first active virtual index
!
         integer(i15) :: last_CC2_o ! first active occupied index 
         integer(i15) :: last_CC2_v ! first active virtual index
         integer(i15) :: last_CCSD_o ! first active occupied index 
         integer(i15) :: last_CCSD_v ! first active virtual index
!
!        Calculate first/last indeces
! 
         call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
         call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
         call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
         call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
         last_CC2_o = first_CC2_o + n_CC2_o - 1
         last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
         last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
         last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!        :: sum_c F_bC * c_ai,Cj ::
!
!        Reorder c_ai,cj to c_aij_c
!
         call wf%mem%alloc(c_aij_c, (n_CC2_v)*((n_CC2_o)**2), n_CC2_v)
!
         do j = 1, n_CC2_o
            do i = 1, n_CC2_o
               do a = 1, n_CC2_v
!
                  ai = index_two(a, i, n_CC2_v)
!
                  aij = index_three(a, i, j, n_CC2_v, n_CC2_o)
!
                  do c = 1, n_CC2_v
!
                     cj = index_two(c, j, n_CC2_v)
!
                     c_aij_c(aij, c) = c_ai_bj(ai, cj)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(rho_aij_b, (n_CC2_v)*((n_CC2_o)**2), n_CC2_v)
!
!        rho_ai_bj += sum_c F_bc * c_ai,cj = sum_c c_aij_c(aij,c) F_ab(b,c) = sum_c c_aij_c(aij,c) F_ab^T(c,b)
!
         call dgemm('N','T',                       & 
                     (n_CC2_v)*((n_CC2_o)**2),     &
                     n_CC2_v,                      & 
                     n_CC2_v,                      & 
                     one,                          & 
                     c_aij_c,                      & 
                     (n_CC2_v)*((n_CC2_o)**2),     &
                     wf%fock_ab,                   & 
                     wf%n_v,                       & 
                     zero,                         & 
                     rho_aij_b,                    & 
                     (n_CC2_v)*((n_CC2_o)**2))
!
         call wf%mem%dealloc(c_aij_c, (n_CC2_v)*((n_CC2_o)**2), n_CC2_v)
!
!        Reorder rho_aij_b into rho_ai_bj
!
         do i = 1, n_CC2_o
            do j = 1, n_CC2_o
               do a = 1, n_CC2_v
!  
                  ai  = index_two(a, i, n_CC2_v)
!
                  aij = index_three(a, i, j, n_CC2_v, n_CC2_o)
!
                  do b = 1, n_CC2_v
!
                     bj = index_two(b, j, n_CC2_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b)
!
                  enddo
               enddo
            enddo
         enddo
!
!
         call wf%mem%dealloc(rho_aij_b, (n_CCSD_v)*((n_CCSD_o)**2), n_CCSD_v)
!
!       ::  - sum_k F_jk * c_ai,bk  ::
!
!        rho_ai_bj += - sum_k F_jk * c_ai,bk = - sum_k c_aib_k(aib,k) F_ij(k,j)^T 
!
         call wf%mem%alloc(c_aib_k, (n_CC2_o)*((n_CC2_v)**2), n_CC2_o)
!
         do b = 1, n_CC2_v
            do i = 1, n_CC2_o
               do a = 1, n_CC2_v
!
                  ai = index_two(a, i, n_CC2_v)
!
                  aib = index_three(a, i, b, n_CC2_v, n_CC2_o)
!
                  do k = 1, n_CC2_o
!
                     bk = index_two(b, k, n_CC2_v)
!
                     c_aib_k(aib, k) = c_ai_bj(ai, bk)
!
                  enddo
               enddo
            enddo
         enddo
!         
         call dgemm('N', 'N',                      & 
                     (n_CC2_o)*((n_CC2_v)**2),     &
                     n_CC2_o,                      & 
                     n_CC2_o,                      & 
                     -one,                         & 
                     c_aib_k,                      & 
                     (n_CC2_o)*((n_CC2_v)**2),     &
                     wf%fock_ij,                   & 
                     wf%n_o,                       & 
                     one,                          & 
                     rho_ai_bj,                    & 
                     (n_CC2_o)*((n_CC2_v)**2))
!
         call wf%mem%dealloc(c_aib_k, (n_CC2_o)*((n_CC2_v)**2), n_CC2_o)
!
!        ::   sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj ) ::
!            
!         sum_ck ( g_bj,kc*(2*c_ai,ck - c_ak,ci) - g_bc,kj*c_ai,ck - g_ki,bc*c_ak,cj ) 
!
!        Construct g_bj,kc 
!
         call wf%mem%alloc(g_bj_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_vo_ov(integral_type, g_bj_kc, &
                           first_CCSD_v, last_CCSD_v, &
                           first_CCSD_o, last_CCSD_o, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!
!        Reordering g_bj_kc to g_ck_bj
!
         call wf%mem%alloc(g_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
         do j = 1, n_CCSD_o
            do b = 1, n_CCSD_v
!
               bj = index_two(b, j, n_CCSD_v) 
!
               do k = 1, n_CC2_o
                  do c = 1, n_CC2_v
!
                     ck = index_two(c, k, n_CC2_v)
                     kc = index_two(k, c, n_CC2_o)
!
                     g_ck_bj(ck, bj) = g_bj_kc(bj, kc)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_bj_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!       rho_ai_bj += sum_ck 2*c_ai_ck * g_ck_bj
!
         call wf%mem%alloc(c_ai_ck, (n_CCSD_v)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
         do a = 1, n_CCSD_v
            do i = 1, n_CCSD_o
!
               ai     = index_two(a, i, n_CCSD_v)
               ai_CC2 = index_two(a, i, n_CC2_v)
!
               c_ai_ck(ai, :) = c_ai_bj(ai_CC2, :)
!
            enddo
         enddo
!
         call wf%mem%alloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     two,                    &
                     c_ai_ck,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     g_ck_bj,                &   
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_ai_bj_CCSD,    &
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(c_ai_ck, (n_CCSD_v)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
!        Reorder c_ak,ci to c_ai_ck
!
         call wf%mem%alloc(c_ai_ck, (n_CCSD_v)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
!
         do i = 1, n_CCSD_o
            do k = 1, n_CC2_o
               do a = 1, n_CCSD_v
!
                  ak = index_two(a, k, n_CC2_v)
                  ai = index_two(a, i, n_CCSD_v)
! 
                  do c = 1, n_CC2_v
!
                     ck = index_two(c, k, n_CC2_v)
                     ci = index_two(c, i, n_CC2_v) 
!
                     c_ai_ck(ai, ck) = c_ai_bj(ak, ci)
!
                  enddo
               enddo
            enddo
         enddo
!
!        rho_ai_bj_CCSD += - sum_ck g_ck_bj*c_ai_ck
!
         call dgemm('N', 'N',                &
                     (n_CCSD_v)*(n_CCSD_o),  &
                     (n_CCSD_v)*(n_CCSD_o),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     -one,                   &
                     c_ai_ck,                &
                     (n_CCSD_v)*(n_CCSD_o),  &
                     g_ck_bj,                &   
                     (n_CC2_o)*(n_CC2_v),    &
                     one,                    &
                     rho_ai_bj_CCSD,         &
                     (n_CCSD_v)*(n_CCSD_o))
!
         call wf%mem%dealloc(c_ai_ck, (n_CCSD_v)*(n_CCSD_o), (n_CC2_o)*(n_CC2_v))
         call wf%mem%dealloc(g_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!        Add to rho
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%alloc(g_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!       Start batching over c ! MUST BE FIXED
! 
        required = 2*(wf%n_J)*((wf%n_v)**2) &
                 + 4*(wf%n_J)*(wf%n_v)*(wf%n_o) &
                 + 2*(wf%n_J)*((wf%n_o)**2)
!    
        required = 4*required         ! In words
!
        batch_dimension  = n_CC2_v ! Batch over the virtual index a
        max_batch_length = 0      ! Initilization of unset variables 
        n_batch          = 0
!
        call num_batch(required, wf%mem%available, max_batch_length, n_batch, batch_dimension)           
!
!       Loop over the number of a batches 
!
        do c_batch = 1, n_batch
!
!          For each batch, get the limits for the a index 
!
           call batch_limits(c_first, c_last, c_batch, max_batch_length, batch_dimension)
           c_length = c_last - c_first + 1

!
!           Construct g_bc_kj       
!
            call wf%mem%alloc(g_bc_kj, c_length*(n_CCSD_v), (n_CC2_o)*(n_CCSD_o))
!           
            integral_type = 'electronic_repulsion'
            call wf%get_vv_oo(integral_type, g_bc_kj,    &
                              first_CCSD_v, last_CCSD_v, &
                              c_first, c_last,           &
                              first_CC2_o, last_CC2_o,   &
                              first_CCSD_o, last_CCSD_o)
!
!
!           Reorder g_bc_kj
!
            do c = 1, c_length
               do b = 1, n_CCSD_v
 !
                  bc = index_two(b, c, n_CCSD_v)
 !
                  do j = 1, n_CCSD_o
 !
                     bj = index_two(b, j, n_CCSD_v)
 !
                     do k = 1, n_CC2_o
 !
                        kj = index_two(k, j, n_CC2_o)
                        ck = index_two(c + c_first - 1, k, n_CC2_v)
 !    
                        g_ck_bj(ck, bj) = g_bc_kj(bc, kj)
 !
                     enddo
                  enddo
               enddo
            enddo
!
            call wf%mem%dealloc(g_bc_kj, c_length*(n_CCSD_v), (n_CC2_o)*(n_CCSD_o))
        enddo
!
!        rho_ai_bj += - sum_ck c_ai_ck * g_ck_bj       
!
         call wf%mem%alloc(c_ai_ck, n_CCSD_o*n_CCSD_v, n_CC2_o*n_CC2_v)
!
         do a = 1, n_CCSD_v
            do i = 1, n_CCSD_o
               ai       = index_two(a, i, n_CCSD_v)
               ai_CC2   = index_two(a, i, n_CC2_v)
               c_ai_ck(ai,:) = c_ai_bj(ai_CC2, :)
            enddo
         enddo
!
         call wf%mem%alloc(rho_ai_bj_CCSD, n_CCSD_o*n_CCSD_v, n_CCSD_o*n_CCSD_v)
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     -one,                   &
                     c_ai_ck,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     g_ck_bj,                &   
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_ai_bj_CCSD,         &
                     (n_CCSD_o)*(n_CCSD_v))

      call wf%mem%dealloc(c_ai_ck, n_CCSD_o*n_CCSD_v, n_CC2_o*n_CC2_v)
!
!        Add to rho
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CC2_v)
                  ai_CCSD = index_two(a, i, n_CCSD_v)
!
                  do b = 1, n_CCSD_v
!
                     bj = index_two(b, j, n_CC2_v)
                     bj_CCSD = index_two(b, j, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj_CCSD(ai_CCSD, bj_CCSD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_bj_CCSD, (n_CCSD_v)*(n_CCSD_o), (n_CCSD_v)*(n_CCSD_o))
!
!
!        Reorder  c_ak,cj to c_aj_ck
!
         call wf%mem%alloc(c_aj_ck, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
         do j = 1, n_CCSD_o
            do k = 1, n_CC2_o
               do a = 1, n_CCSD_v
!
                  ak = index_two(a, k, n_CC2_v)
                  aj = index_two(a, j, n_CCSD_v)
! 
                  do c = 1, n_CC2_v
!
                     ck = index_two(c, k, n_CC2_v)
                     cj = index_two(c, j, n_CC2_v) 
!
                     c_aj_ck(aj, ck) = c_ai_bj(ak, cj)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(rho_aj_bi, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
         call dgemm('N', 'N',                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     (n_CC2_o)*(n_CC2_v),    &
                     -one,                   &
                     c_aj_ck,                &
                     (n_CCSD_o)*(n_CCSD_v),  &
                     g_ck_bj,                &  !g_ck_bi(ck,bi) = g_bc,ki 
                     (n_CC2_o)*(n_CC2_v),    &
                     zero,                   &
                     rho_aj_bi,              &
                     (n_CCSD_o)*(n_CCSD_v))
!
         call wf%mem%dealloc(g_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
         call wf%mem%dealloc(c_aj_ck, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!        Reorder rho_aj_bi into rho_ai_bj
!
         do i = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!  
               ai  = index_two(a, i, n_CC2_v)
!                 
               do j = 1, n_CCSD_o
!                
                  aj = index_two(a, j, n_CCSD_v)
!
                  do b = 1, n_CCSD_v          
!
                     bj = index_two(b, j, n_CC2_v)
                     bi = index_two(b, i, n_CCSD_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aj_bi(aj, bi)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_aj_bi, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!  
!
      end subroutine jacobian_mlccsd_i2_mlccsd
!
!
      module subroutine jacobian_mlccsd_j2_mlccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD J2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl 
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!             
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ab_ij
         real(dp), dimension(:,:) :: c_ab_ij
!
         real(dp), dimension(:,:), allocatable :: g_kc_ld
         real(dp), dimension(:,:), allocatable :: g_kl_cd
!
         real(dp), dimension(:,:), allocatable :: x_cd_ij
         real(dp), dimension(:,:), allocatable :: x_ab_kl
!
         real(dp), dimension(:,:), allocatable :: c_ab_kl
         real(dp), dimension(:,:), allocatable :: c_cd_ij
!
         real(dp), dimension(:,:), allocatable :: I_kl_ij
!
         real(dp), dimension(:,:), allocatable :: rho_ab_ij_CCSD
!
         integer(i15) :: a = 0, b = 0, c = 0, d = 0
         integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
         integer(i15) :: ab = 0, cd = 0
         integer(i15) :: ai = 0, bj = 0, bl = 0, ak = 0, ci = 0, dj = 0
         integer(i15) :: kl = 0, ij = 0
         integer(i15) :: kc = 0, ld = 0
         integer(i15) :: ab_CC2 = 0, ij_CC2 = 0
!
         integer(i15) :: aibj = 0, akbl = 0, cidj = 0
!
!        Active space variables
!
         integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
         integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
         integer(i15) :: first_CC2_o ! first active occupied index 
         integer(i15) :: first_CC2_v ! first active virtual index
         integer(i15) :: first_CCSD_o ! first active occupied index 
         integer(i15) :: first_CCSD_v ! first active virtual index
!
         integer(i15) :: last_CC2_o ! first active occupied index 
         integer(i15) :: last_CC2_v ! first active virtual index
         integer(i15) :: last_CCSD_o ! first active occupied index 
         integer(i15) :: last_CCSD_v ! first active virtual index
!
!        Calculate first/last indeces
! 
         call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
         call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
         call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
         call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
         last_CC2_o = first_CC2_o + n_CC2_o - 1
         last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
         last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
         last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!        Constructing g_kc_ld
!
         call wf%mem%alloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_KC_LD,  &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v, &
                           first_CC2_o, last_CC2_o, &
                           first_CC2_v, last_CC2_v)
!
         call wf%mem%alloc(g_KL_CD, (n_CC2_o)**2, (n_CC2_v)**2)
!
!        Reorder g_kc_ld to g_kl_cd
!
         do C = 1, n_CC2_v
            do D = 1, n_CC2_v
!
               CD = index_two(C, D, n_CC2_v)
!
               do K = 1, n_CC2_o
!
                  KC = index_two(K, C, n_CC2_o)
!
                  do L = 1, n_CC2_o
! 
                     KL = index_two(K, L, n_CC2_o)
                     LD = index_two(L, D, n_CC2_o)
!
                     g_KL_CD(KL, CD) = g_KC_LD(KC, LD)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        Reordered X2 amplitudes  
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_cd_ij, (n_CC2_v)**2, (n_CCSD_o)**2)
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
!
               ij = index_two(i, j, n_CCSD_o)
!
               do c = 1, n_CC2_v
!
                  ci = index_two(c, i, n_CC2_v)
!
                  do d = 1, n_CC2_v
!
                     dj = index_two(d, j, n_CC2_v)
                     cd = index_two(c, d, n_CC2_v)
!  
                     cidj = index_packed(ci, dj)
!
                     x_cd_ij(cd, ij) = wf%x2am(cidj, 1)
!
                  enddo
               enddo
            enddo
         enddo
         

         call wf%destruct_x2am
!
!        I_kl_ij = g_kl_cd * t_cd_ij
!
         call wf%mem%alloc(I_kl_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
!
         call dgemm('N', 'N',       &
                     (n_CC2_o)**2,  &
                     (n_CCSD_o)**2, &
                     (n_CC2_v)**2,  &
                     one,           &  
                     g_kl_cd,       &
                     (n_CC2_o)**2,  &
                     x_cd_ij,       &
                     (n_CC2_v)**2,  &
                     zero,          &
                     I_kl_ij,       &
                     (n_CC2_o)**2)
!
         call wf%mem%dealloc(x_cd_ij, (n_CC2_v)**2, (n_CCSD_o)**2)
         call wf%mem%alloc(c_ab_kl, n_CCSD_v**2, n_CC2_o**2)
!
         do a = 1, n_CCSD_v
            do b = 1, n_CCSD_v
!
               ab     = index_two(a, b, n_CCSD_v)
               ab_CC2 = index_two(a, b, n_CC2_v)
!
               c_ab_kl(ab, :) = c_ab_ij(ab_CC2, :)
!   
            enddo
         enddo
!
!        rho_ab_ij += c_ab_kl * X_kl_ij
!
         call wf%mem%alloc(rho_ab_ij_CCSD, n_CCSD_v**2, n_CCSD_o**2)
!

         call dgemm('N', 'N',       &
                     (n_CCSD_v)**2, &
                     (n_CCSD_o)**2, &
                     (n_CC2_o)**2,  &
                     one,           &  
                     c_ab_kl,       &
                     (n_CCSD_v)**2, &
                     I_kl_ij,       &
                     (n_CC2_o)**2,  &
                     zero,          &
                     rho_ab_ij_CCSD,&
                     (n_CCSD_v)**2)
!
         call wf%mem%dealloc(c_ab_kl, n_CCSD_v**2, n_CC2_o**2)
         call wf%mem%dealloc(I_kl_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
!
!        I_kl_ij = g_kl_cd * c_cd_ij
!
         call wf%mem%alloc(c_cd_ij, n_CC2_v**2, n_CCSD_o**2)
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
!
               ij     = index_two(i, j, n_CCSD_o)
               ij_CC2 = index_two(i, j, n_CC2_o)
!
               c_cd_ij(:, ij) = c_ab_ij(:, ij_CC2)
!  
            enddo
         enddo
!
         call wf%mem%alloc(I_kl_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
!
         call dgemm('N', 'N',       &
                     (n_CC2_o)**2,  &
                     (n_CCSD_o)**2, &
                     (n_CC2_v)**2,  &
                     one,           &  
                     g_kl_cd,       &
                     (n_CC2_o)**2,  &
                     c_cd_ij,       &
                     (n_CC2_v)**2,  &
                     zero,          &
                     I_kl_ij,       &
                     (n_CC2_o)**2)
!
         call wf%mem%dealloc(c_cd_ij, n_CC2_v**2, n_CCSD_o**2)    
         call wf%mem%dealloc(g_kl_cd, (n_CC2_o)**2, (n_CC2_v)**2)
!
!        Reordered X2 amplitudes        
!
         call wf%read_amplitudes
!
         call wf%mem%alloc(x_ab_kl, (n_CCSD_v)**2, (n_CC2_o)**2)
!
         do k = 1, n_CC2_o
            do l = 1, n_CC2_o
!
               kl = index_two(k, l, n_CC2_o)
!
               do a = 1, n_CCSD_v
!
                  ak = index_two(a, k, n_CC2_v)
!
                  do b = 1, n_CCSD_v
!
                     bl = index_two(b, l, n_CC2_v)
                     ab = index_two(a, b, n_CCSD_v)
!  
                     akbl = index_packed(ak, bl)
!
                     x_ab_kl(ab, kl) = wf%x2am(akbl, 1)
!
                  enddo
               enddo
            enddo
         enddo

         call wf%destruct_x2am
!
!        rho_ab_ij += t_ab_kl * X_kl_ij
!
         call dgemm('N', 'N',       &
                     (n_CCSD_v)**2, &
                     (n_CCSD_o)**2, &
                     (n_CC2_o)**2,  &
                     one,           &  
                     x_ab_kl,       &
                     (n_CCSD_v)**2, &
                     I_kl_ij,       &
                     (n_CC2_o)**2,  &
                     one,           &
                     rho_ab_ij_CCSD,&
                     (n_CCSD_v)**2)
!
         call wf%mem%dealloc(I_kl_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
         call wf%mem%dealloc(x_ab_kl, (n_CCSD_v)**2, (n_CC2_o)**2)
!
         do a = 1, n_CCSD_v
            do b = 1, n_CCSD_v
!
               ab = index_two(a, b, n_CCSD_v)
               ab_CC2 = index_two(a, b, n_CC2_v)
!
               do i = 1, n_CCSD_o
                  do j = 1, n_CCSD_o
!
                     ij = index_two(i, j, n_CCSD_o)
                     ij_CC2 = index_two(i, j, n_CC2_o)
!
                     rho_ab_ij(ab_CC2, ij_CC2) = rho_ab_ij(ab_CC2, ij_CC2) +  rho_ab_ij_CCSD(ab, ij)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ab_ij_CCSD, n_CCSD_v**2, n_CCSD_o**2)
!
      end subroutine jacobian_mlccsd_j2_mlccsd
!
!
      module subroutine jacobian_mlccsd_k2_mlccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian MLCCSD K2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^K2 =    sum_kl g_Ki,Lj * c_aK,bL 
!!                       + sum_cd g_aC,bD * c_Ci,Dj
!! 
!!       For the last term we batch over a and b and 
!!       add each batch to rho_ai_bj 
!!               
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ab_ij
         real(dp), dimension(:,:) :: c_ab_ij
!
         real(dp), dimension(:,:), allocatable :: L_Ki_J
         real(dp), dimension(:,:), allocatable :: L_ca_J
         real(dp), dimension(:,:), allocatable :: L_ac_J
         real(dp), dimension(:,:), allocatable :: L_db_J
         real(dp), dimension(:,:), allocatable :: L_bd_J
!
         real(dp), dimension(:,:), allocatable :: g_ki_lj
         real(dp), dimension(:,:), allocatable :: g_kl_ij
         real(dp), dimension(:,:), allocatable :: g_ac_bd
         real(dp), dimension(:,:), allocatable :: g_ab_cd
!
         real(dp), dimension(:,:), allocatable :: c_ab_kl
         real(dp), dimension(:,:), allocatable :: c_cd_ij
!
         real(dp), dimension(:,:), allocatable :: rho_batch_ab_ij
         real(dp), dimension(:,:), allocatable :: rho_ab_ij_CCSD
!
         integer(i15) :: i = 0, j = 0, k = 0, l = 0  
         integer(i15) :: a = 0, b = 0, c = 0, d = 0  
!
         integer(i15) :: ab = 0, db = 0, ca = 0, cd = 0, full_ab = 0, ac = 0, bd = 0
         integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0, full_ij = 0
         integer(i15) :: ab_CC2 = 0, ij_CC2 = 0
!
!        Batching and memory handling variables
!
         integer(i15) :: a_n_batch = 0, a_first = 0, a_last = 0, a_length = 0, a_max_length = 0, a_batch = 0
         integer(i15) :: b_n_batch = 0, b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0
!
         integer(i15) :: required = 0, available = 0
!
!        Active space variables
!
         integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
         integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
         integer(i15) :: first_CC2_o ! first active occupied index 
         integer(i15) :: first_CC2_v ! first active virtual index
         integer(i15) :: first_CCSD_o ! first active occupied index 
         integer(i15) :: first_CCSD_v ! first active virtual index
!
         integer(i15) :: last_CC2_o ! first active occupied index 
         integer(i15) :: last_CC2_v ! first active virtual index
         integer(i15) :: last_CCSD_o ! first active occupied index 
         integer(i15) :: last_CCSD_v ! first active virtual index
!
!        Calculate first/last indeces
! 
         call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
         call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
         call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
         call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
         last_CC2_o = first_CC2_o + n_CC2_o - 1
         last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
         last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
         last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
         call wf%mem%alloc(g_Ki_Lj, (n_CCSD_o)*(n_CC2_o), (n_CCSD_o)*(n_CC2_o))
!
         integral_type = 'electronic_repulsion'
         call wf%get_oo_oo(integral_type, g_Ki_Lj,    &
                           first_CC2_o, last_CC2_o,   &
                           first_CCSD_o, last_CCSD_o, &
                           first_CC2_o, last_CC2_o,   &
                           first_CCSD_o, last_CCSD_o)
!
!
!        Reorder g_ki_lj to g_kl_ij
!
         call wf%mem%alloc(g_kl_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
!
         do j = 1, n_CCSD_o
            do i = 1, n_CCSD_o
!
               ij = index_two(i, j, n_CCSD_o)
!
               do k = 1, n_CC2_o
!
                  ki = index_two(k, i, n_CC2_o)
!
                  do l = 1, n_CC2_o
! 
                     kl = index_two(k, l, n_CC2_o)
                     lj = index_two(l, j, n_CC2_o)
!
                     g_KL_ij(kl, ij) = g_Ki_Lj(ki, lj) 
!                     
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_Ki_Lj, (n_CCSD_o)*(n_CC2_o), (n_CCSD_o)*(n_CC2_o))
!
!        rho_ab_ij += sum_kl g_ki,lj * c_ak,bl = sum_kl c_ab_ij(ab,kl) g_kl_ij(kl,ij)  
!
         call wf%mem%alloc(c_ab_kl, n_CCSD_v**2, n_CC2_o**2)
!
         do a = 1, n_CCSD_v
            do b = 1, n_CCSD_v
!
               ab       = index_two(a, b, n_CCSD_v)
               ab_CC2   = index_two(a, b, n_CC2_v)
!
               c_ab_kl(ab, :) = c_ab_ij(ab_CC2, :)
! 
            enddo
         enddo
!
         call wf%mem%alloc(rho_ab_ij_CCSD, n_CCSD_v**2, n_CCSD_o**2)
!
         call dgemm('N', 'N',          & 
                     (n_CCSD_v)**2,    &
                     (n_CCSD_o)**2,    &
                     (n_CC2_o)**2,     &
                     one,              &
                     c_ab_kl,          &
                     (n_CCSD_v)**2,    &
                     g_kl_ij,          & 
                     (n_CC2_o)**2,     &
                     zero,             &
                     rho_ab_ij_CCSD,   &
                     (n_CCSD_v)**2)
!
         call wf%mem%dealloc(g_kl_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
         call wf%mem%dealloc(c_ab_kl, n_CCSD_v**2, n_CC2_o**2)
!
         do a = 1, n_CCSD_v
            do b = 1, n_CCSD_v
!
               ab = index_two(a, b, n_CCSD_v)
               ab_CC2 = index_two(a, b, n_CC2_v)
!
               do i = 1, n_CCSD_o
                  do j = 1, n_CCSD_o
!
                     ij = index_two(i, j, n_CCSD_o)
                     ij_CC2 = index_two(i, j, n_CC2_o)
!
                     rho_ab_ij(ab_CC2, ij_CC2) = rho_ab_ij(ab_CC2, ij_CC2) +  rho_ab_ij_CCSD(ab, ij)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ab_ij_CCSD, n_CCSD_v**2, n_CCSD_o**2)
!
!        Prepare for batching over a and b
!
!        ::  sum_cd g_aC,bD * c_Ci,Dj ::
!
         required = max(3*(wf%n_v)**2*(wf%n_J) + 2*(wf%n_v)*(wf%n_o)*(wf%n_J),      & ! Needed to get L_db_J
                     (wf%n_v)**4 + 2*(wf%n_v)**2*(wf%n_J))                            ! Needed to get g_ac_bd
!
         required = required*4  ! Words
!
         a_max_length = 0
         call num_two_batch(required, wf%mem%available, a_max_length, a_n_batch, n_CCSD_v)
!
!        Initialize some variables for batching
!
         a_first  = 0
         a_last   = 0
         a_length = 0
!
!        Start looping over a-batches
!
         do a_batch = 1, a_n_batch
!   
            call batch_limits(a_first, a_last, a_batch, a_max_length, n_CCSD_v)
            a_length = a_last - a_first + 1     
!
!           Start looping over batches of b
!
            b_first  = 0
            b_last   = 0
            b_length = 0
!
            b_max_length = a_max_length
!
            do b_batch = 1, a_n_batch
!
               call batch_limits(b_first ,b_last ,b_batch, b_max_length, n_CCSD_v)
               b_length = b_last - b_first + 1 
!
!              Allocate g_ac_bd = g_acbd
!
               call wf%mem%alloc(g_ac_bd, (n_CC2_v)*a_length, (n_CC2_v)*b_length)
!
               integral_type = 'electronic_repulsion'
               call wf%get_vv_vv(integral_type, g_ac_bd, &
                                 a_first, a_last,        &
                                 first_CC2_v, last_CC2_v,&
                                 b_first, b_last,        &
                                 first_CC2_v, last_CC2_v)
!
!              sum_cd g_ac,bd * c_ci,dj = sum_cd g_ac,bd c_cd,ij = sum_cd g_ab_cd c_cd_ij  
!
!              Reorder g_ca_db into g_ab_cd 
!              (Here, g_ab_cd = g_acbd = g_ca_db.)
!
               call wf%mem%alloc(g_ab_cd, a_length*b_length, (n_CC2_v)**2) 
!
               do b = 1, b_length
                  do a = 1, a_length
!
                     ab = index_two(a, b, a_length)
!
                     do d = 1, n_CC2_v
!
                        bd = index_two(b, d, b_length)
!
                        do c = 1, n_CC2_v
!
                           ac = index_two(a, c, a_length)
                           cd = index_two(c, d, n_CC2_v)
!
                           g_ab_cd(ab, cd) = g_ac_bd(ac, bd) ! = g_acbd 
!
                        enddo
                     enddo
                  enddo
               enddo
!
               call wf%mem%dealloc(g_ac_bd, (n_CC2_v)*a_length, (n_CC2_v)*b_length) 
!
               call wf%mem%alloc(rho_batch_ab_ij, a_length*b_length, (n_CC2_o)**2)
!
!              rho_ab_ij += sum_cd g_ac,bd * c_ci,dj = sum_cd g_ab_cd(ab, cd) c_ab_ij(cd, ij) 
!
               call wf%mem%alloc(c_cd_ij, n_CC2_v**2, n_CCSD_o**2)
!
               do i = 1, n_CCSD_o
                  do j = 1, n_CCSD_o
!
                     ij       = index_two(i, j, n_CCSD_o)
                     ij_CC2   = index_two(i, j, n_CC2_o)
!
                     c_cd_ij(:, ij) = c_ab_ij(:, ij_CC2)
!
                  enddo
               enddo
!
               call dgemm('N', 'N',            &  
                            a_length*b_length, &
                            (n_CCSD_o)**2,     &  
                            (n_CC2_v)**2,      &  
                            one,               &  
                            g_ab_cd,           &
                            a_length*b_length, &
                            c_cd_ij,           & 
                            (n_CC2_v)**2,      &  
                            zero,              &
                            rho_batch_ab_ij,   &
                            a_length*b_length)
!               
               call wf%mem%dealloc(g_ab_cd, a_length*b_length, (n_CC2_v)**2)
               call wf%mem%dealloc(c_cd_ij, n_CC2_v**2, n_CCSD_o**2)
!
!              Reorder into rho_ab_ij
!
               do b = 1, b_length
                  do a = 1, a_length
!
                     ab = index_two(a, b, a_length)
!
                     full_ab = index_two(a + a_first - 1, b + b_first - 1, n_CC2_v)
!
                     do i = 1, n_CCSD_o
                        do j = 1, n_CCSD_o
!
                           ij = index_two(i, j, n_CCSD_o)
                           full_ij = index_two(i, j, n_CC2_o)
!
                           rho_ab_ij(full_ab, full_ij) = rho_ab_ij(full_ab, full_ij) + rho_batch_ab_ij(ab, ij)
!
                        enddo
                     enddo
                  enddo
               enddo
!
               call wf%mem%dealloc(rho_batch_ab_ij,  a_length*b_length, (n_CCSD_o)**2) 
!
            enddo ! End batches of b 
         enddo ! End batches of a 
!         
      end subroutine jacobian_mlccsd_k2_mlccsd
!
!
end submodule jacobian