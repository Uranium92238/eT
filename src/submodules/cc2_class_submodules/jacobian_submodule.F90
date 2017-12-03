submodule (cc2_class) jacobian
!
!!
!!    jacobian transformation submodule (cc2)) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, September 2017
!!
!!    contains the following family of procedures of the cc2 class:
!!
!!    jacobian_transformation: Directs the transformation by a.
!!
!!    cc2 contributions to jacobi transformation
!!
!!    jacobian_cc2_a1
!!    jacobian_cc2_b1
!!    jacobian_cc2_a2
!!    jacobian_cc2_b2
!!
!!    Upper case indices are general indices, lower case indices are restricted
!!    to the cc2 orbital space.
!! 
!
   implicit none 
!
   logical :: debug   = .false.
   logical :: timings = .false.
!
   character(len=40) :: integral_type
!
!
contains
!
!
   module subroutine jacobian_cc2_transformation_cc2(wf, c_a_i, c_aibj)
!!
!!    Jacobian transformation (cc2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Directs the transformation by the ccSD Jacobi matrix,
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
      class(cc2) :: wf 
!
!     Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_s2am, 1)   :: c_aibj ! c_aibj     
!
!     Local unpacked and reordered vectors 
!
      real(dp), dimension(:,:), allocatable :: rho_a_i         ! rho_ai   = (A c)_ai
      real(dp), dimension(:,:), allocatable :: rho_ai_bj       ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_sym   ! rho_ai   = (A c)_aibj
      real(dp), dimension(:,:), allocatable :: c_ai_bj         ! rho_ai   = (A c)_aibj
!
!     Indices 
!
      integer(i15) :: a = 0, ab = 0, ai = 0, b = 0 
      integer(i15) :: bj = 0, i = 0, ij = 0, j = 0, aibj = 0
!
!     Allocate and zero the transformed vector (singles part)
!
      call wf%mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!     :: ccS contributions to the singles c vector ::  
! 
      call wf%initialize_amplitudes
      call wf%read_amplitudes
!
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
!
!     :: cc2 contributions to transformed vector :: 
!
      call wf%jacobian_cc2_a1(rho_a_i, c_a_i)
!
!     Allocate the incoming unpacked doubles vector 
!
      call wf%mem%alloc(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
      c_ai_bj = zero
!
      call squareup(c_aibj, c_ai_bj, (wf%n_o)*(wf%n_v)) ! Pack out vector 
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
      do i = 1, (wf%n_o)*(wf%n_v)
!
         c_ai_bj(i,i) = two*c_ai_bj(i,i)
!
      enddo
!
!     - B1 term -
!
      call wf%jacobian_cc2_b1(rho_a_i, c_ai_bj)
!
!
!     Allocate unpacked transformed vector
!
      call wf%mem%alloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
      rho_ai_bj = zero 
!
!     - A2 term -
!   
      call wf%jacobian_cc2_a2(rho_ai_bj, c_a_i)
!
!     Last term is already symmetric (B2). Perform the symmetrization 
!     rho_ai_bj = P_ij^ab rho_ai_bj now, for convenience 
!
!     Allocate temporary symmetric transformed vector 
!
      call wf%mem%alloc(rho_ai_bj_sym, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
      rho_ai_bj_sym = zero
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
                  rho_ai_bj_sym(ai, bj) = rho_ai_bj(ai, bj) + rho_ai_bj(bj, ai)
!
               enddo
            enddo
         enddo
      enddo
!
      rho_ai_bj = rho_ai_bj_sym

!
!     Done with temporary vector; deallocate
! 
      call wf%mem%dealloc(rho_ai_bj_sym, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     - B2 term -
!
      call wf%jacobian_cc2_b2(rho_ai_bj, c_ai_bj)
!
!     Divide rho_ai_bj by biorthonormal, and save to c_aibj     
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = index_two(a, i, wf%n_v)
!
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  if (ai == bj) rho_ai_bj(ai, bj) = half*rho_ai_bj(ai, bj)
!
                   c_aibj(aibj, 1) = rho_ai_bj(ai, bj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      c_a_i = rho_a_i
!
      call wf%mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
   end subroutine jacobian_cc2_transformation_cc2
!
!
   module subroutine jacobian_cc2_a1_cc2(wf, rho_a_i, c_a_i)
!!
!!    jacobian tem a1 (cc2))
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, june 2017
!!
!!    calculates the a1 contribution to the jacobi transformation,
!!
!!      rho_ai^a1 = sum_ckbj  (L_kc,jb u_ki^ca c_bj  - g_kc,jb u_ki^cb c_aj - g_kc,jb u_kj^ab c_ci)
!!
!!     with, 
!!
!!    u_ik^ac = 2*s_ik^ac - 2*s_ik^ca,
!!
!!    which is constructed while batching over c
!!
      implicit none
!
      class(cc2) :: wf 
!
!     incoming vectors c and rho 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i 
!
!     allocatables
!
      real(dp), dimension(:,:), allocatable :: u_a_jbk
      real(dp), dimension(:,:), allocatable :: u_bkc_i
      real(dp), dimension(:,:), allocatable :: u_ai_ck
      real(dp), dimension(:,:), allocatable :: L_ck_bj
      real(dp), dimension(:,:), allocatable :: g_jb_kc
      real(dp), dimension(:,:), allocatable :: X_a_c
      real(dp), dimension(:,:), allocatable :: X_j_i
      real(dp), dimension(:,:), allocatable :: X_ai_bj
!
!     indices
!
      integer(i15) :: k = 0, c = 0, j = 0, a = 0, i = 0, b = 0
      integer(i15) :: kc = 0, ck = 0, ai = 0, ak = 0, ic = 0, jb = 0, bj = 0, jc = 0, aj = 0
      integer(i15) :: bi = 0, kb = 0, ci = 0, bk = 0
      integer(i15) :: jkc = 0, bkc = 0, jbk = 0
      integer(i15) :: ajbk = 0, akbj = 0, ciak = 0, ckai = 0, ckbi = 0, cibk = 0
!
      call wf%read_amplitudes
!           
!     :: Term a1.3 ::
!     - sum_(bjck) g_kc,jb u_kj^ab c_ci
!
      call wf%mem%alloc(u_a_jbk, (wf%n_v), (wf%n_o**2)*(wf%n_v)) ! u_ak_bj
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do a = 1, wf%n_v
!
              aj = index_two(a, j, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ak = index_two(a, k, wf%n_v)
                  bk = index_two(b, k, wf%n_v)
!
                  akbj = index_packed(ak, bj)
                  ajbk = index_packed(aj, bk)
                  jbk = index_three(j, b, k, wf%n_o, wf%n_v)
!
                  u_a_jbk(a, jbk) = two*wf%s2am(akbj, 1) - wf%s2am(ajbk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_s2am
!
!     g_jb_kc = g_jb,kc 
!
      call wf%mem%alloc(g_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_jb_kc)
!
!     X_a_c = - sum_(bjk)u_a_jbk * g_jbk_c
!
      call wf%mem%alloc(X_a_c, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o**2), &
                  -one,                 &
                  u_a_jbk,              &
                  wf%n_v,               &
                  g_jb_kc,              &
                  (wf%n_v)*(wf%n_o**2), &
                  zero,                 &
                  X_a_c,                &
                  wf%n_v)
!
      call wf%mem%dealloc(u_a_jbk, (wf%n_v), (wf%n_o**2)*(wf%n_v))  
!
!     rho_a_i =+ sum_(c) X_a_c * c_c_i 
!
      call dgemm('N', 'N',  &
                  wf%n_v,   &
                  wf%n_o,   &
                  wf%n_v,   &
                  one,      &
                  X_a_c,    &
                  wf%n_v,   &
                  c_a_i,    & ! c_c_i
                  wf%n_v,   &
                  one,      &
                  rho_a_i,  &
                  wf%n_v)
!
      call wf%mem%dealloc(X_a_c, wf%n_v, wf%n_v)
!           
!     :: Term a.2 ::
!     - sum_(bjck) g_jb,kc u_ki^cb c_aj
!
!     construct u_ck,bi ordered as u_bkc_i
!
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(u_bkc_i, (wf%n_v**2)*(wf%n_o), (wf%n_o)) ! u_ak_bj
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do b = 1, wf%n_v
!
              bk = index_two(b, k, wf%n_v)
!
               do i = 1, wf%n_o
!
                  ci = index_two(c, i, wf%n_v)
                  bi = index_two(b, i, wf%n_v)
!
                  ckbi = index_packed(ck, bi)
                  cibk = index_packed(ci, bk)
                  bkc = index_three(b, k, c, wf%n_v, wf%n_o)
!
                  u_bkc_i(bkc, i) = two*wf%s2am(ckbi, 1) - wf%s2am(cibk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_s2am
!
!     X_j_i = -sum_(bkc) g_j_bkc * u_bkc_i
!
      call wf%mem%alloc(X_j_i, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',            &
                  wf%n_o,             &
                  wf%n_o,             &
                  wf%n_o*(wf%n_v**2), &
                  -one,               &
                  g_jb_kc,            &
                  wf%n_o,             &
                  u_bkc_i,            &
                  wf%n_o*(wf%n_v**2), &
                  zero,               &
                  X_j_i,              &
                  wf%n_o)
!
      call wf%mem%alloc(u_bkc_i, (wf%n_v**2)*(wf%n_o), (wf%n_o)) ! u_ak_bj
!
!     rho_a_i =+ sum_(c) X_j_i * c_a_j 
!
      call dgemm('N', 'N',  &
                  wf%n_v,   &
                  wf%n_o,   &
                  wf%n_o,   &
                  one,      &
                  c_a_i,    & ! c_a_j
                  wf%n_v,   &
                  X_j_i,    & 
                  wf%n_o,   &
                  one,      &
                  rho_a_i,  &
                  wf%n_v)
!
      call wf%mem%dealloc(X_j_i, wf%n_o, wf%n_o)
!          
!     :: Term a1.1 ::
!     sum_(ckbj) L_jb,kc u_ki^ca c_bj
!
!     construct L_kc,bj ordered as L_kc_bj
!
      call wf%mem%alloc(L_ck_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do k = 1, wf%n_o
        do c = 1, wf%n_v
!
          kc = index_two(k, c, wf%n_o)
          ck = index_two(c, k, wf%n_v)
!
          do j = 1, wf%n_o
!
            jc = index_two(j, c, wf%n_o)
!
            do b = 1, wf%n_v
!
              kb = index_two(k, b, wf%n_o)
              jb = index_two(j, b, wf%n_o)
              bj = index_two(b, j, wf%n_v)
!
              L_ck_bj(ck, bj) = two*g_jb_kc(kc, jb) - g_jb_kc(jc, kb)
!
            enddo
          enddo
        enddo
      enddo
!
!     construct u_ai,ck ordered as u_ai_ck
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(u_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) 
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do a = 1, wf%n_v
!
              ak = index_two(a, k, wf%n_v)
!
               do i = 1, wf%n_o
!
                  ci = index_two(c, i, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
!
                  ckai = index_packed(ck, ai)
                  ciak = index_packed(ci, ak)
!
                  u_ai_ck(ai, ck) = two*wf%s2am(ckai, 1) - wf%s2am(ciak, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_s2am
!
!     X_ai_bj = sum_(ck) u_ai_ck * L_ck_bj
!
      call wf%mem%alloc(X_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'N',    &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  L_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ai_bj,           &
                  (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(u_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%mem%dealloc(L_ck_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) 
!
      call dgemm('N', 'N',            &
                  (wf%n_o)*(wf%n_v),  &
                  1,                  &
                  (wf%n_o)*(wf%n_v),  &
                  one,                &
                  X_ai_bj,            &
                  (wf%n_o)*(wf%n_v),  &
                  c_a_i,              & ! c_bj
                  (wf%n_o)*(wf%n_v),  &
                  one,                &
                  rho_a_i,            &
                  (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(X_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_cc2_a1_cc2
!
!
   module subroutine jacobian_cc2_b1_cc2(wf, rho_a_i, c_ai_bj)
!!
!!    jacobian tem b1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, june 2017
!!
!!    calculates the b1 contribution to the jacobi transformation,
!!
!!       b1:   sum_ck F_kc*(2c_ai,ck - c_ak,ci) 
!!           - sum_ckj L_jikc * c_aj,ck + sum_cbk L_abkc * c_bi,ck
!!
!!
!!    L_abkc is constructed while batching over a.
!!
      implicit none
!  
      class(cc2) :: wf
!
      real(dp), dimension(:,:)            :: c_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
!     batching variables 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, offset = 0
      integer(i15) :: a_n_batch = 0, a_first = 0, a_last = 0, a_batch = 0, a_length = 0
!
!     allocatables
! 
      real(dp), dimension(:,:), allocatable :: d_ai_kc 
      real(dp), dimension(:,:), allocatable :: c_bkc_i
      real(dp), dimension(:,:), allocatable :: g_ji_kc
      real(dp), dimension(:,:), allocatable :: g_ab_kc
      real(dp), dimension(:,:), allocatable :: L_jck_i
      real(dp), dimension(:,:), allocatable :: L_ab_kc
!
      integer(i15) :: k = 0, c = 0, a = 0, i = 0, j = 0, b = 0
      integer(i15) :: ak = 0, ck = 0, ai = 0, ci = 0, ji = 0, ki = 0, kc = 0, jc = 0, ac = 0, ab = 0, bi = 0, kb = 0
      integer(i15) :: jck = 0, bkc
!
!     :: Term 1 ::
!     sum_ck F_kc*(2c_ai,ck - c_ak,ci) 
!
      call wf%mem%alloc(d_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     construct 2c_ai,ck - c_ak,ci = d_ai_kc
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do k = 1, wf%n_o
!
               ak = index_two(a, k, wf%n_v)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  kc = index_two(k, c, wf%n_o)
                  ci = index_two(c, i, wf%n_v)
!
                  d_ai_kc(ai, kc) = two*c_ai_bj(ai, ck) - c_ai_bj(ak, ci)
!
               enddo
            enddo
         enddo
      enddo
!
!
      call dgemm('N', 'N',            &
                  (wf%n_o)*(wf%n_v),  &
                  1,                  &
                  (wf%n_o)*(wf%n_v),  &
                  one,                &
                  d_ai_kc,            &
                  (wf%n_o)*(wf%n_v),  &
                  wf%fock_ia,         & !f_kc
                  (wf%n_o)*(wf%n_v),  &
                  one,                &
                  rho_a_i,            &
                  (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(d_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 2 ::
!     - sum_ckj L_ji,kc * c_aj,ck = - sum_ckj (2*g_jikc - g_kijc) * c_aj,ck (L_ji,kc ordered as L_jck_i)
!
      call wf%mem%alloc(g_ji_kc, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_ji_kc)
!
!     construct L_jikc ( = two*g_ji_kc - g_ki_jc )
!
      call wf%mem%alloc(L_jck_i, (wf%n_o**2)*(wf%n_v), wf%n_o)
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do c = 1, wf%n_v
               do k = 1, wf%n_o
                  kc = index_two(k, c, wf%n_o)
                  ki = index_two(k, i, wf%n_o)
                  jc = index_two(j, c, wf%n_o)
                  ji = index_two(j, i, wf%n_o)
                  jck = index_three(j, c, k, wf%n_o, wf%n_v)
!
                  L_jck_i(jck, i) = two*g_ji_kc(ji, kc) - g_ji_kc(ki, jc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ji_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     add term to rho_a_i
!
      call dgemm('N', 'N',            &
                  wf%n_v,             &
                  wf%n_o,             &
                  wf%n_v*wf%n_o**2,   &
                  -one,               &
                  c_ai_bj,            &
                  wf%n_v,             &
                  L_jck_i,            &
                  wf%n_v*wf%n_o**2,   &
                  one,                &
                  rho_a_i,            &
                  wf%n_v) 
!
      call wf%mem%dealloc(L_jck_i, (wf%n_o**2)*wf%n_v, wf%n_o)
!
!    :: Term 3 ::
!    sum_cbk L_abkc * c_bi,ck
!
!    Prepare for batching over a
!
     required = max(2*(wf%n_v**3)*(wf%n_o), &
                    (wf%n_v**3)*(wf%n_o) + wf%n_o*wf%n_v*(wf%n_j) + (wf%n_v**2)*(wf%n_j))
!
     required = required*4  ! Words
!
     max_length = 0
     call num_batch(required, wf%mem%available, max_length, a_n_batch, wf%n_v)
!
!    initialize some variables for batching
!
     a_first  = 0
     a_last   = 0
     a_length = 0
!
!    Start looping over a-batches
!
     do a_batch = 1, a_n_batch
!  
        call batch_limits(a_first ,a_last ,a_batch, max_length, wf%n_v)
!
        a_length = a_last - a_first + 1
!
        call wf%mem%alloc(g_ab_kc, wf%n_v*a_length, (wf%n_v)*(wf%n_o))
!
        integral_type = 'electronic_repulsion'
        call wf%get_vv_ov(integral_type, g_ab_kc, a_first, a_last, 1, wf%n_v, 1, wf%n_o, 1, wf%n_v)
!
!       construct L_ab_kc = 2*g_ab_kc - g_ac_kb 
!
        call wf%mem%alloc(L_ab_kc, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
        do a = 1, a_length
           do b = 1, wf%n_v
              do c = 1, wf%n_v
                 do k = 1, wf%n_o
!
                    ab = index_two(a, b, a_length)
                    kc = index_two(k, c, wf%n_o)
                    kb = index_two(k, b, wf%n_o)
                    ac = index_two(a, c, a_length)
!
                    L_ab_kc(ab, kc) = two*g_ab_kc(ab, kc) - g_ab_kc(ac, kb)
!
                 enddo
              enddo
           enddo
        enddo
!
        call wf%mem%dealloc(g_ab_kc, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
!
!       Reorder c_bi_ck to c_bkc_i
!
        call wf%mem%alloc(c_bkc_i, (wf%n_v**2)*(wf%n_o), wf%n_o)
        do i = 1, wf%n_o
           do b = 1, wf%n_v
              do c = 1, wf%n_v
                 do k = 1, wf%n_o
!
                    bi = index_two(b, i, wf%n_v)
                    ck = index_two(c, k, wf%n_v)
                    bkc = index_three(b, k, c, wf%n_v, wf%n_o)
!
                    c_bkc_i(bkc, i) = c_ai_bj(bi, ck)
!
                 enddo
              enddo
           enddo
        enddo
!
!       add contribution for current batch to rho
!

        call dgemm('N', 'N',                 &
                    a_length,                &
                    wf%n_o,                  &
                    (wf%n_v**2)*wf%n_o,      &
                    one,                     &
                    L_ab_kc,                 &
                    a_length,                &
                    c_bkc_i,                 &
                    (wf%n_v**2)*wf%n_o,      &
                    one,                     &
                    rho_a_i(a_first, 1),     &
                    wf%n_v) 
!
        call wf%mem%dealloc(L_ab_kc, (wf%n_v)*a_length, (wf%n_o)*(wf%n_v))
        call wf%mem%dealloc(c_bkc_i, (wf%n_v**2)*wf%n_o, wf%n_o)
!
      enddo ! batching over a
!
!
   end subroutine jacobian_cc2_b1_cc2
!
!
   module subroutine jacobian_cc2_a2_cc2(wf, rho_ai_bj, c_a_i)
!!
!!    jacobian tem a2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, june 2017
!!
!!    calculates the a2 contribution to the jacobi transformation,
!!
!!       a2:   sum_c g_ai,bc * c_cj - sum_K g_ai,kj * c_bk.
!!
!!    g_ai,bc is constructed in batches of c.
!!
      implicit none
!  
      class(cc2) :: wf
!
      real(dp), dimension(: , :)          :: rho_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
!
!     batching variables 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: c_n_batch = 0, c_first = 0, c_last = 0, c_batch = 0, c_length = 0
!
!     allocatables
!
      real(dp), dimension(:,:), allocatable :: g_ai_bc
      real(dp), dimension(:,:), allocatable :: g_ai_Kj
      real(dp), dimension(:,:), allocatable :: g_aij_K
      real(dp), dimension(:,:), allocatable :: rho_aij_b
      real(dp), dimension(:,:), allocatable :: rho_aib_j
!
!     indices 
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, K = 0
      integer(i15) :: Kj = 0, ai = 0, aij = 0, bj = 0
!
!     :: Term 1 ::
!     sum_c g_aib_c * c_c_j(= sum_c g_ai,bc * c_cj)
!
      required = max(2*wf%n_v*(wf%n_v)*(wf%n_J) + 2*(wf%n_v)*(wf%n_o)*(wf%n_J), &
                     wf%n_v*(wf%n_v)*(wf%n_J) + (wf%n_v**2)*(wf%n_v)*wf%n_o)
!
!
      required = required*4  ! Words
!
      max_length = 0
      call num_batch(required, wf%mem%available, max_length, c_n_batch, wf%n_v)
!
!     Initialize some variables for batching
!
      c_first  = 0
      c_last   = 0
      c_length = 0
!
!     Start looping over batches of c
!
      do c_batch = 1, c_n_batch
!
         call batch_limits(c_first, c_last, c_batch, max_length, wf%n_v)
!
!        construct g_ai,bc ordered as g_ai_bc batching over c
!
         c_length = c_last - c_first + 1
!
         call wf%mem%alloc(g_ai_bc, wf%n_v*wf%n_o, wf%n_v*c_length)
!
         integral_type = 'electronic_repulsion'
         call wf%get_vo_vv(integral_type, g_ai_bc, 1, wf%n_v, 1, wf%n_o, 1, wf%n_v, c_first, c_last)
!
!        Add contribution tho rho
!        
          call dgemm('N', 'N', &
                      (wf%n_v**2)*wf%n_o, &
                      wf%n_o, &
                      c_length, &
                      one, &
                      g_ai_bc, &
                      ((wf%n_v)**2)*wf%n_o, &
                      c_a_i(c_first, 1), &
                      wf%n_v, &
                      one, &
                      rho_ai_bj, &
                      (wf%n_v**2)*wf%n_o)
!
         call wf%mem%dealloc(g_ai_bc, wf%n_v*wf%n_o, wf%n_v*c_length)

!
      enddo ! batching over c
! 
!     :: Term 2 ::
!     - sum_K g_aij_K * c_K_b(= - sum_K g_ai,Kj * c_bK)
! 
!     construct g_ai,Kj ordered as g_aij_K
! 
      call wf%mem%alloc(g_ai_kj, wf%n_v*wf%n_o, wf%n_o**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_oo(integral_type, g_ai_kj)
! 
!     Reorder g_ai_Kj to g_aij_K
! 
      call wf%mem%alloc(g_aij_K, wf%n_v*wf%n_o**2, wf%n_o)
      g_aij_K = zero
      do a = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
               ai = index_two(a, i, wf%n_v)
               do K = 1, wf%n_o
                  Kj = index_two(K, j, wf%n_o)
                  g_aij_K(aij, K) = g_ai_Kj(ai, Kj)
               enddo
            enddo
         enddo
      enddo
      call wf%mem%dealloc(g_ai_Kj, wf%n_v*wf%n_o, wf%n_o*(wf%n_o))
! 
      call wf%mem%alloc(rho_aij_b, wf%n_v*wf%n_o**2, wf%n_v)
! 
!      Add term to rho_ai_bj
! 
      call dgemm('N', 'T',            &
                  wf%n_v*wf%n_o**2,   &
                  wf%n_v,             &
                  wf%n_o,             &
                  -one,               &
                  g_aij_K,            &
                  wf%n_v*wf%n_o**2,   &
                  c_a_i,              &
                  wf%n_v,             &
                  zero,               &
                  rho_aij_b,          &
                  wf%n_v*wf%n_o**2)
! 
      call wf%mem%dealloc(g_aij_K, wf%n_v*wf%n_o**2, wf%n_o)
! 
! 
!     Reorder into rho_ai_bj
! 
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do j = 1, wf%n_o
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b)
               enddo
            enddo
         enddo
      enddo
      call wf%mem%dealloc(rho_aij_b, wf%n_v*wf%n_o**2, wf%n_v)
!
!
   end subroutine jacobian_cc2_a2_cc2
!
!
   module subroutine jacobian_cc2_b2_cc2(wf, rho_ai_bj, c_ai_bj)
!!
!!    jacobian tem b2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, june 2017
!!
!!    calculates the b2 contribution to the jacobi transformation,
!!
!!       b2:   ε_ij^ab*c_ai,bj.
!!
!!
      implicit none
!  
      class(cc2) :: wf
!
      real(dp), dimension(:,:)   :: c_ai_bj 
      real(dp), dimension(:,:)   :: rho_ai_bj
!
!     Local routine variables
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0
!
       do i = 1, wf%n_o
!
          do a = 1, wf%n_v
!
             ai = index_two(a, i, wf%n_v)
!
             do j = 1, wf%n_o
!
                do b = 1, wf%n_v
!
                   bj = index_two(b, j, wf%n_v)
!
                   rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) +  c_ai_bj(ai, bj)*(wf%fock_diagonal(wf%n_o + a, 1) &
                                                  - wf%fock_diagonal(i, 1) &
                                                  + wf%fock_diagonal(wf%n_o + b, 1) &
                                                  - wf%fock_diagonal(j, 1))
!
                enddo
             enddo
          enddo
       enddo
!
   end subroutine jacobian_cc2_b2_cc2
!
!
    module subroutine cvs_rho_aibj_projection_cc2(wf, vec_aibj)
!!
!!    Rho projection for cVS (cc2),
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
      implicit none
!
      class(cc2) :: wf
      real(dp), dimension(:, :) :: vec_aibj
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital
!
      do i = 1, wf%n_o
       do j = 1, wf%n_o
!
          core_orbital = .false.
          do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
             if ((i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) .or. &
                (j .eq. wf%core_excited_state_specifications%index_core_mo(core, 1))) core_orbital = .true.
!
          enddo
!
          if (.not. core_orbital) then
             do a = 1, wf%n_v
                do b = 1, wf%n_v
                   ai = index_two(a, i, wf%n_v)
                   bj = index_two(b, j, wf%n_v)
                   aibj = index_packed(ai, bj)

                   vec_aibj(aibj, 1) = zero
!
                enddo
             enddo
          endif
       enddo
    enddo
!
  end subroutine cvs_rho_aibj_projection_cc2
!
!
!
end submodule jacobian