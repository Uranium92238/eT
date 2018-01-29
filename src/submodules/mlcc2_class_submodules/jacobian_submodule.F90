submodule (mlcc2_class) jacobian
!
!!
!!    Jacobian transformation submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!    jacobian_transformation: Directs the transformation by A.
!!
!!    MLCC2 contributions to jacobi transformation
!!
!!    jacobian_mlcc2_a1
!!    jacobian_mlcc2_b1
!!    jacobian_mlcc2_a2
!!    jacobian_mlcc2_b2
!!
!!    Upper case indices are general indices, lower case indices are restricted
!!    to the CC2 orbital space.
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

   module subroutine jacobian_mlcc2_transformation_mlcc2(wf, c_a_i, c_aibj)
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
      class(mlcc2) :: wf 
!
!     Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_x2am, 1)   :: c_aibj ! c_aibj     
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
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
!
!     Allocate and zero the transformed vector (singles part)
!
      call wf%mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!     :: CCS contributions to the singles c vector ::  
!
      call wf%read_single_amplitudes
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
!
!     :: MLCC2 contributions to transformed vector :: 
!
      call wf%jacobian_mlcc2_a1(rho_a_i, c_a_i)
!
!     Calculatenumber of active indices
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
!     - B1 term -
!
      call wf%jacobian_mlcc2_b1(rho_a_i, c_ai_bj)
!
!
!     Allocate unpacked transformed vector
!
      call wf%mem%alloc(rho_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v) 
      rho_ai_bj = zero 
!
!     - A2 term -
!  
      call wf%jacobian_mlcc2_a2(rho_ai_bj, c_a_i)
!
!     Last term is already symmetric (B2). Perform the symmetrization 
!     rho_ai_bj = P_ij^ab rho_ai_bj now, for convenience 
!
!     Allocate temporary symmetric transformed vector 
!
      call wf%mem%alloc(rho_ai_bj_sym, n_active_o*n_active_v, n_active_o*n_active_v) 
      rho_ai_bj_sym = zero
!
      do j = 1, n_active_o
         do b = 1, n_active_v
!
            bj = index_two(b, j, n_active_v)
!
            do i = 1, n_active_o
               do a = 1, n_active_v
!
                  ai = index_two(a, i, n_active_v)
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
      call wf%jacobian_mlcc2_b2(rho_ai_bj, c_ai_bj)
!
!     Divide rho_ai_bj by biorthonormal, and save to c_aibj     
!
      do a = 1, n_active_v
         do i = 1, n_active_o
!
            ai = index_two(a, i, n_active_v)
!
            do b = 1, n_active_v
               do j = 1, n_active_o
!
                  bj = index_two(b, j, n_active_v)
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
      call wf%mem%dealloc(rho_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
      call wf%mem%dealloc(c_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
!
      c_a_i = rho_a_i
!
      call wf%mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
   end subroutine jacobian_mlcc2_transformation_mlcc2
!  
!
   module subroutine jacobian_mlcc2_a1_mlcc2(wf, rho_a_i, c_a_i)
!!
!!    Jacobian tem A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Calculates the A1 contribution to the jacobi transformation,
!!
!!       A1: 2*sum_BJck u_ik^ac*g_kc,JB*c_BJ - sum_Bjck u_kj^ca*g_kc,jB*c_BI
!!            - sum_BJck u_ik^ac*g_kB,Jc*c_BJ - sum_bJck u_ki^cb*g_kc,Jb*c_AJ, 
!!
!!     with, 
!!
!!    u_ik^ac = 2*s_ik^ac - 2*s_ik^ca,
!!
!!    which is constructed while batching over c
!!
      implicit none
!
      class(mlcc2) :: wf 
!
!     Incoming vectors c and rho 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i 
!
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: L_ai_J  
      real(dp), dimension(:,:), allocatable :: L_kc_J
      real(dp), dimension(:,:), allocatable :: L_BJ_J
      real(dp), dimension(:,:), allocatable :: L_jB_J
      real(dp), dimension(:,:), allocatable :: L_Jc_J
      real(dp), dimension(:,:), allocatable :: L_kB_J
      real(dp), dimension(:,:), allocatable :: u_ai_kc
      real(dp), dimension(:,:), allocatable :: g_ai_kc
      real(dp), dimension(:,:), allocatable :: u_bkc_i
      real(dp), dimension(:,:), allocatable :: g_kc_BJ
      real(dp), dimension(:,:), allocatable :: g_kc_JB
      real(dp), dimension(:,:), allocatable :: g_jkc_B
      real(dp), dimension(:,:), allocatable :: g_jB_kc
      real(dp), dimension(:,:), allocatable :: g_kB_Jc
      real(dp), dimension(:,:), allocatable :: c_JB
      real(dp), dimension(:,:), allocatable :: X_kc
      real(dp), dimension(:,:), allocatable :: X_a_B
      real(dp), dimension(:,:), allocatable :: X_J_i
      real(dp), dimension(:,:), allocatable :: rho_ai_active_space
!
!     Indices
!
      integer(i15) :: k = 0, c = 0, j = 0, a = 0, i = 0, b = 0
      integer(i15) :: kc = 0, ck = 0, ai = 0, ak = 0, ic = 0, jb = 0, BJ = 0, Jc = 0
      integer(i15) :: bi = 0, kB = 0, ci = 0
      integer(i15) :: jkc = 0, bkc = 0
      integer(i15) :: akci = 0, aick = 0
!
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
      call wf%read_amplitudes
!
      call wf%mem%alloc(u_ai_kc, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      do i = 1, n_active_o
         do a = 1, n_active_v
!
            ai = index_two(a, i, n_active_v)
!
            do c = 1, n_active_v
!
               do k = 1, n_active_o
!
                  kc = index_two(k, c, n_active_o)
                  ck = index_two(c, k, n_active_v)
                  ak = index_two(a, k, n_active_v)
                  ci = index_two(c, i, n_active_v)
!
                  aick = index_packed(ai, ck)
                  akci = index_packed(ak, ci)
!
                  u_ai_kc(ai,kc) = two*wf%x2am(aick, 1) - wf%x2am(akci, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_cc2_double_amplitudes
!           
!     :: Term 1 ::
!     2* sum_BJck u_ai_kc * g_kc_BJ * c_BJ ( = sum_BJck u_ik^ac * g_kc,JB * c_BJ )
!
!     Construct g_kc,JB ordered as g_kc_BJ
!     We will use read_cholesky_ai to get reordered cholesky vectors of ia-type.
!     This is possible because L^J_ia is unchanged on T1-transformation
!
      call wf%mem%alloc(g_kc_JB, n_active_o*n_active_v, (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_JB,        &
                        first_active_o, last_active_o, &
                        first_active_v, last_active_v, &
                        1, wf%n_o,                     &
                        1, wf%n_v)
!
      call wf%mem%alloc(c_JB, (wf%n_o)*(wf%n_v), 1)
      do B = 1, wf%n_v
         do J = 1, wf%n_o
!
            JB = index_two(J, B, wf%n_o)
!
            c_JB(JB, 1) = c_a_i(B, J)
!
         enddo
      enddo
!
!      Add contribution to rho_a_i
!
!      sum_BJ g_kc_BJ*c_BJ = X_kc
!
       call wf%mem%alloc(X_kc, n_active_o*n_active_v, 1)
       call dgemm('N', 'T',                    &
                   (n_active_o)*n_active_v,    &
                   1,                          &
                   (wf%n_o)*(wf%n_v),          &
                   one,                        &
                   g_kc_JB,                    &
                   (n_active_o)*n_active_v,    &
                   c_JB,                       &
                   1,                          &
                   zero,                       &
                   X_kc,                       &
                   (n_active_o)*n_active_v)
!
       call wf%mem%dealloc(g_kc_JB, n_active_o*n_active_v, (wf%n_o)*(wf%n_v))
       call wf%mem%dealloc(c_JB, (wf%n_o)*(wf%n_v), 1)
!  
       call wf%mem%alloc(rho_ai_active_space, n_active_v, n_active_o)
!
!      rho_a_i += 2 * sum_ck u_ai_kc * X_kc
!
       call dgemm('N', 'N',                    &
                   (n_active_o)*(n_active_v),  &
                   1,                          &
                   n_active_o*n_active_v,      &
                   two,                        &
                   u_ai_kc,                    &
                   (n_active_o)*(n_active_v),  &
                   X_kc,                       &
                   (n_active_o)*n_active_v,    &
                   zero,                       &
                   rho_ai_active_space,        &
                   (n_active_o)*(n_active_v))
!
       call wf%mem%dealloc(X_kc, n_active_o*n_active_v, 1)
!
       do i = 1, n_active_o
          do a = 1, n_active_v
             rho_a_i(a + first_active_v - 1, i + first_active_o - 1) = rho_a_i(a + first_active_v - 1, i + first_active_o - 1)&
                                                                       + rho_ai_active_space(a,i)
          enddo
       enddo
!

       call wf%mem%dealloc(rho_ai_active_space, n_active_v, n_active_o)
!
!       :: Term 2 ::
!        - sum_Bjck u_a_jkc * g_jkc_B * c_BI (= - sum_Bjck u_kj^ca * g_kc,jB * c_BI)
!
!       Construct g_kc_jB
!
        call wf%mem%alloc(g_kc_jB, n_active_o*n_active_v, n_active_o*wf%n_v)
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_kc_jB,     &
                        first_active_o, last_active_o, &
                        first_active_v, last_active_v, &
                        first_active_o, last_active_o, &
                        1, wf%n_v)
!
!       Reorder g_kc_jB to g_jkc_B
!
        call wf%mem%alloc(g_jkc_B, n_active_v*n_active_o**2, wf%n_v)
        do c = 1, n_active_v
           do k = 1, n_active_o
              kc = index_two(k, c, n_active_o)
              do j = 1, n_active_o
                 jkc = index_three(j, k, c, n_active_o, n_active_o)
                 do B = 1, wf%n_v
                    jB = index_two(j, B, n_active_o)
                    g_jkc_B(jkc, B) = g_kc_jB(kc, jB)
                 enddo
              enddo
           enddo
        enddo
!
        call wf%mem%dealloc(g_kc_jB, n_active_o*n_active_v, n_active_o*wf%n_v)
!
!       Add contribution to rho_a_i
!
!       sum u_a_jkc*g_jkc_B = X_a_B
!
        call wf%mem%alloc(X_a_B, n_active_v, wf%n_v)
        call dgemm('N', 'N',                          &
                    (n_active_v),                     &
                    wf%n_v,                           &
                    n_active_o*n_active_v*n_active_o, &
                    one,                              &
                    u_ai_kc,                          &
                    (n_active_v),                     &
                    g_jkc_B,                          &
                    n_active_o*n_active_v*n_active_o, &
                    zero,                             &
                    X_a_B,                            &
                    (n_active_v))
!
        call wf%mem%dealloc(g_jkc_B, n_active_v*n_active_o**2, wf%n_v)
!
!       rho_a_i += - sum_ck  X_a_B*c_B_I
!
        call dgemm('N', 'N',                    &
                    (n_active_v),               &
                    wf%n_o,                     &     
                    wf%n_v,                     &
                    -one,                       &
                    X_a_B,                      &
                    (n_active_v),               &
                    c_a_i,                      &
                    wf%n_v,                     &
                    one,                        &
                    rho_a_i(first_active_v, 1), &
                    wf%n_v)
!
        call wf%mem%dealloc(X_a_B, n_active_v, wf%n_v)
!
!        :: Term 3 ::
!        - sum_BJck u_ai_kc * g_kc_BJ * c_BJ ( = - sum_BJck u_ik^ac * g_kB,Jc * c_BJ )
!
!        construct g_kB_Jc
!
         call wf%mem%alloc(g_kB_Jc, n_active_o*wf%n_v, n_active_v*wf%n_o)
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_kB_Jc, &
                           first_active_o, last_active_o, &
                           1, wf%n_v,                     &
                           1, wf%n_o,                     &
                           first_active_v, last_active_v)
!
!        Reorder g_kB_Jc to g_kc_BJ
!           
         call wf%mem%alloc(g_kc_BJ, n_active_o*n_active_v, (wf%n_o)*(wf%n_v))
         g_kc_BJ = zero
!
         do c = 1, n_active_v
            do k = 1, n_active_o
               kc = index_two(k, c, n_active_o)
               do B = 1, wf%n_v
                  kB = index_two(k, B, n_active_o)
                  do J = 1, wf%n_o
                     BJ = index_two(B, J, wf%n_v)
                     Jc = index_two(J, c, wf%n_o)
                     g_kc_BJ(kc, BJ) = g_kB_Jc(kB, Jc) 
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_kB_Jc, n_active_o*wf%n_v, n_active_v*wf%n_o)
!
!        Add term to rho_a_i
!
!        sum_BJ g_kc_BJ*c_BJ = X_kc
!
         call wf%mem%alloc(X_kc, n_active_o*n_active_v, 1)
         call dgemm('N', 'N',                    &
                     (n_active_o)*n_active_v,    &
                     1,                          &
                     (wf%n_o)*(wf%n_v),          &
                     one,                        &
                     g_kc_BJ,                    &
                     (n_active_o)*n_active_v,    &
                     c_a_i,                      &
                     (wf%n_o)*(wf%n_v),          &
                     zero,                       &
                     X_kc,                       &
                     (n_active_o)*n_active_v)
!
         call wf%mem%dealloc(g_kc_BJ, n_active_o*n_active_v, (wf%n_o)*(wf%n_v))
!  
         call wf%mem%alloc(rho_ai_active_space, n_active_v, n_active_o)
!
!        rho_a_i += 2 * sum_ck u_ai_kc * X_kc
!
         call dgemm('N', 'N',                    &
                     (n_active_o)*(n_active_v),  &
                     1,                          &
                     n_active_o*n_active_v,      &
                     -one,                       &
                     u_ai_kc,                    &
                     (n_active_o)*(n_active_v),  &
                     X_kc,                       &
                     (n_active_o)*n_active_v,    &
                     zero,                       &
                     rho_ai_active_space,        &
                     (n_active_o)*(n_active_v))
!
         call wf%mem%dealloc(X_kc, n_active_o*n_active_v, 1)
!
         do i = 1, n_active_o
            do a = 1, n_active_v
               rho_a_i(a + first_active_v - 1, i + first_active_o - 1) = rho_a_i(a + first_active_v - 1, i + first_active_o - 1)&
                                                                         + rho_ai_active_space(a,i)
            enddo
         enddo
!
         call wf%mem%dealloc(rho_ai_active_space, n_active_v, n_active_o)
!
!        :: Term 4 ::
!        - sum_bJck c_A_J * g_J_bkc u_bkc_i*( =  - sum_bJck u_ki^cb * g_Jb,kc * c_AJ)
!
!        Reorder u_ai_kc(bi, kc) to u_bkc_i(bkc, i)
!
         call wf%mem%alloc(u_bkc_i, n_active_o*n_active_v*n_active_v, n_active_o)
!
         do i = 1, n_active_o
            do b = 1, n_active_v
               bi = index_two(b, i, n_active_v)
               do c = 1, n_active_v
                  do k = 1, n_active_o
                     kc = index_two(k, c, n_active_o)
                     bkc = index_three(b, k, c, n_active_v, n_active_o)
                     u_bkc_i(bkc, i) = u_ai_kc(bi, kc)
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(u_ai_kc, (n_active_o)*(n_active_v), (n_active_o)*n_active_v)
!
!        Construct g_Jb_kc
!
         call wf%mem%alloc(g_Jb_kc, n_active_v*wf%n_o, n_active_v*n_active_o)
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_Jb_kc, &
                           1, wf%n_o,                     &
                           first_active_v, last_active_v, &
                           first_active_o, last_active_o, &
                           first_active_v, last_active_v)
!
!        Add contribution to rho_a_i
!
!        sum_bck g_Jb_kc*u_bkc_i = X_J_i
!
         call wf%mem%alloc(X_J_i, wf%n_o, n_active_o)
!
         call dgemm('N', 'N',                      &
                     wf%n_o,                       &
                     n_active_o,                   &
                     (n_active_v**2)*(n_active_o), &
                     one,                          &
                     g_Jb_kc,                      &
                     wf%n_o,                       &
                     u_bkc_i,                      &
                     (n_active_v**2)*(n_active_o), &
                     zero,                         &
                     X_J_i,                        &
                     wf%n_o)
!
         call wf%mem%dealloc(g_Jb_kc, n_active_v*wf%n_o, n_active_v*n_active_o)
         call wf%mem%dealloc(u_bkc_i, n_active_o*n_active_v*n_active_v, n_active_o)
!
!        sum_J C_A_J * X_J_i
!

         call dgemm('N', 'N',                      &
                     wf%n_v,                       &
                     n_active_o,                   &
                     wf%n_o,                       &
                     -one,                         &
                     c_a_i,                        &
                     wf%n_v,                       &
                     X_J_i,                        &
                     wf%n_o,                       &
                     one,                          &
                     rho_a_i(1, first_active_o),   &
                     wf%n_v)
!
         call wf%mem%dealloc(X_J_i, wf%n_o, n_active_o)
!
   end subroutine jacobian_mlcc2_a1_mlcc2
!
!
   module subroutine jacobian_mlcc2_b1_mlcc2(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian tem B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Calculates the B1 contribution to the jacobi transformation,
!!
!!       B1:   sum_ck F_kc*(2c_ai,ck - c_ak,ci) 
!!           - sum_ckj L_jIkc * c_aj,ck + sum_cbk L_Abkc * c_bi,ck
!!
!!
!!    L_Abkc is constructed while batching over A.
!!
      implicit none
!  
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:)            :: c_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
!     Variables for batching
!
      integer(i15) :: required = 0 
      integer(i15) :: current_A_batch = 0
!
      type(batching_index) :: batch_A 
!
!     Allocatables
! 
      real(dp), dimension(:,:), allocatable :: c_ai_ck ! = c_ak,ci
      real(dp), dimension(:,:), allocatable :: c_bkc_i
      real(dp), dimension(:,:), allocatable :: rho_ai_active_space
      real(dp), dimension(:,:), allocatable :: fock_ck_active_space
      real(dp), dimension(:,:), allocatable :: L_jI_J
      real(dp), dimension(:,:), allocatable :: L_kc_J
      real(dp), dimension(:,:), allocatable :: L_Ab_J
      real(dp), dimension(:,:), allocatable :: g_jI_kc
      real(dp), dimension(:,:), allocatable :: g_Ab_kc
      real(dp), dimension(:,:), allocatable :: L_jck_I
      real(dp), dimension(:,:), allocatable :: L_Ab_kc
!
      integer(i15) :: k = 0, c = 0, a = 0, i = 0, j = 0, b = 0
      integer(i15) :: ak = 0, ck = 0, ai = 0, ci = 0, jI = 0, kI = 0, kc = 0, jc = 0, Ac = 0, Ab = 0, bi = 0, kb = 0
      integer(i15) :: jck = 0, bkc
!
      logical :: reorder
!
!     Active space variables
!  
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
!     :: Term 1 ::
!     sum_ck F_kc*(2c_ai,ck - c_ak,ci) 
!
!     Restrict Fock matrix ia-block to active space ordered as F_ck
!
      call wf%mem%alloc(fock_ck_active_space, n_active_v, n_active_o) 
!
      do  k = 1, n_active_o
         do c = 1, n_active_v
!
            fock_ck_active_space(c, k) = wf%fock_ia(k + first_active_o - 1, c + first_active_v - 1)
!
         enddo
      enddo
!
      call wf%mem%alloc(c_ai_ck, n_active_o*n_active_v, n_active_o*n_active_v)
!
!     Construct 2c_ai,ck - c_ak,ci
!
      do i = 1, n_active_o
         do a = 1, n_active_v
!
            ai = index_two(a, i, n_active_v)
!
            do k = 1, n_active_o
!
               ak = index_two(a, k, n_active_v)
!
               do c = 1, n_active_v
!
                  ck = index_two(c, k, n_active_v)
                  ci = index_two(c, i, n_active_v)
!
                  c_ai_ck(ai, ck) = two*c_ai_bj(ai, ck) - c_ai_bj(ak, ci)
!
               enddo
            enddo
         enddo
      enddo
!
!     Adding term to rho     
!
      call wf%mem%alloc(rho_ai_active_space, n_active_v, n_active_o)
!
      call dgemm('N', 'N', &
                  n_active_o*n_active_v, &
                  1,                     &
                  n_active_o*n_active_v, &
                  one,                  &
                  c_ai_ck,               &
                  n_active_o*n_active_v, &
                  fock_ck_active_space,  &
                  n_active_o*n_active_v, &
                  zero,                   &
                  rho_ai_active_space,   &
                  n_active_o*n_active_v)
!
      call wf%mem%dealloc(c_ai_ck, n_active_o*n_active_v, n_active_o*n_active_v)
!
      call wf%mem%dealloc(fock_ck_active_space, n_active_o, n_active_v)
! 
      do a = 1, n_active_v
         do i = 1, n_active_o
            rho_a_i(a + first_active_v - 1,i + first_active_o - 1) = rho_a_i(a + first_active_v - 1,i + first_active_o - 1) &
                                                                   +  rho_ai_active_space(a,i)
         enddo
      enddo
!
      call wf%mem%dealloc(rho_ai_active_space, n_active_v, n_active_o)
!
!     :: Term 2 ::
!     - sum_ckj L_jI,kc * c_aj,ck = - sum_ckj (2*g_jIkc - g_kIjc) * c_aj,ck (L_jI,kc ordered as L_jck_I)
!
      call wf%mem%alloc(g_jI_kc, n_active_o*wf%n_o, n_active_v*n_active_o)
!  
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_jI_kc,        &
                        first_active_o, last_active_o, &
                        1, wf%n_o,                     &
                        first_active_o, last_active_o, &
                        first_active_v, last_active_v)
!
!     Construct L_jIkc ( = two*g_jI_kc - g_kI_jc )
!
      call wf%mem%alloc(L_jck_I, n_active_o*n_active_v*n_active_o, wf%n_o)
!
      do j = 1, n_active_o
         do i = 1, wf%n_o
            do c = 1, n_active_v
               do k = 1, n_active_o
                  kc = index_two(k, c, n_active_o)
                  kI = index_two(k, I, n_active_o)
                  jc = index_two(j, c, n_active_o)
                  jI = index_two(j, I, n_active_o)
                  jck = index_three(j, c, k, n_active_o, n_active_v)
!
                  L_jck_I(jck, I) = two*g_jI_kc(jI, kc) - g_jI_kc(kI, jc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_jI_kc, n_active_o*wf%n_o, n_active_v*n_active_o)
!
!     Add term to rho_a_i
!
      call dgemm('N', 'N',                    &
                  n_active_v,                 &
                  wf%n_o,                     &
                  n_active_v*n_active_o**2,   &
                  -one,                       &
                  c_ai_bj,                    &
                  n_active_v,                 &
                  L_jck_I,                    &
                  n_active_v*n_active_o**2,   &
                  one,                        &
                  rho_a_i(first_active_v,1),  &
                  wf%n_v)
!
      call wf%mem%dealloc(L_jck_I, n_active_o*n_active_v*n_active_o, wf%n_o)
!
!     :: Term 3 ::
!     sum_cbk L_Abkc * c_bi,ck
!!
!     Prepare for batching over A
!
      required = max(wf%get_vvov_required_mem(), 2*(wf%n_v**3)*(wf%n_o))
!
!     Initialize batching variable 
!
      call batch_A%init(wf%n_v)
      call wf%mem%num_batch(batch_A, required)         
!
!     Loop over the number of a batches 
!
      do current_A_batch = 1, batch_A%num_batches
!
!        Determine the limits for the current a-batch 
!
         call batch_a%determine_limits(current_A_batch)
!
         flush(unit_output)
        call wf%mem%alloc(g_Ab_kc, n_active_v*(batch_A%length), n_active_o*n_active_v)
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_ov(integral_type, g_Ab_kc, &
                           batch_A%first, batch_A%last, &
                           first_active_v, last_active_v, &
                           first_active_o, last_active_o, &
                           first_active_v, last_active_v)
!
!        Construct L_Ab_kc = 2*g_Ab_kc - g_Ac_kb 
!
         call wf%mem%alloc(L_Ab_kc, n_active_v*(batch_A%length), n_active_o*n_active_v)
         do A = 1, batch_A%length
            do b = 1, n_active_v
               do c = 1, n_active_v
                  do k = 1, n_active_o
!
                     Ab = index_two(A, b, batch_A%length)
                     kc = index_two(k, c, n_active_o)
                     kb = index_two(k, b, n_active_o)
                     Ac = index_two(A, c, batch_A%length)
!
                     L_Ab_kc(Ab, kc) = two*g_Ab_kc(Ab, kc) - g_Ab_kc(Ac, kb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_Ab_kc, n_active_v*batch_A%length, n_active_o*n_active_v)
!
!        Reorder c_bi_ck to c_bkc_i
!
         call wf%mem%alloc(c_bkc_i, (n_active_v**2)*n_active_o, n_active_o)
         do i = 1, n_active_o
            do b = 1, n_active_v
               do c = 1, n_active_v
                  do k = 1, n_active_o
!
                     bi = index_two(b, i, n_active_v)
                     ck = index_two(c, k, n_active_v)
                     bkc = index_three(b, k, c, n_active_v, n_active_o)
!
                     c_bkc_i(bkc, i) = c_ai_bj(bi, ck)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Add contribution for current batch to rho
!

         call dgemm('N', 'N',                                &
                     batch_A%length,                         &
                     n_active_o,                             &
                     (n_active_v**2)*n_active_o,             &
                     one,                                    &
                     L_Ab_kc,                                &
                     batch_A%length,                         &
                     c_bkc_i,                                &
                     (n_active_v**2)*n_active_o,             &
                     one,                                    &
                     rho_a_i(batch_A%first, first_active_o), &
                     wf%n_v) 
!
         call wf%mem%dealloc(L_Ab_kc, n_active_v*(batch_A%length), n_active_o*n_active_v)
         call wf%mem%dealloc(c_bkc_i, (n_active_v**2)*n_active_o, n_active_o)
!
      enddo ! Batching over A
!
!
   end subroutine jacobian_mlcc2_b1_mlcc2
!
!
   module subroutine jacobian_mlcc2_a2_mlcc2(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian tem A2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Calculates the A2 contribution to the jacobi transformation,
!!
!!       A2:   sum_C g_ai,bC * c_Cj - sum_K g_ai,Kj * C_bK.
!!
!!    g_ai,bC is constructed in batches of C.
!!
      implicit none
!  
      class(mlcc2) :: wf
!
      real(dp), dimension(: , :)          :: rho_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
!
!     Batching variables 
!
      integer(i15) :: required = 0, offset = 0
      integer(i15) :: current_c_batch = 0
!
      type(batching_index) :: batch_c 

!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: g_ai_bC
      real(dp), dimension(:,:), allocatable :: g_ai_Kj
      real(dp), dimension(:,:), allocatable :: g_aij_K
      real(dp), dimension(:,:), allocatable :: rho_aij_b
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: L_bC_J
      real(dp), dimension(:,:), allocatable :: L_Kj_J
      real(dp), dimension(:,:), allocatable :: rho_aib_j
!
!     Indices 
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, K = 0
      integer(i15) :: Kj = 0, ai = 0, aij = 0, bj = 0
!
!     Active space variables
!  
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
!     :: Term 1 ::
!     sum_C g_aib_C * c_C_j(= sum_C g_ai,bC * c_Cj)
!
      required = wf%get_vvvo_required_mem()
!
!     Initialize batching variable 
!
      call batch_c%init(wf%n_v)
      call wf%mem%num_batch(batch_c, required)         
!
!     Loop over the number of a batches 
!
      do current_c_batch = 1, batch_c%num_batches
!
!        Determine the limits for the current a-batch 
!
         call batch_c%determine_limits(current_c_batch)
!
         call wf%mem%alloc(g_ai_bC, n_active_v*n_active_o, n_active_v*batch_c%length)
!
         integral_type = 'electronic_repulsion'
         call wf%get_vo_vv(integral_type, g_ai_bC, &
                           first_active_v, last_active_v, &
                           first_active_o, last_active_o, &
                           first_active_v, last_active_v, &
                           batch_c%first, batch_c%last)
!
!        Add contribution tho rho
!        
          call dgemm('N', 'N', &
                      (n_active_v**2)*n_active_o, &
                      n_active_o, &
                      batch_c%length, &
                      one, &
                      g_ai_bC, &
                      ((n_active_v)**2)*n_active_o, &
                      c_a_i(batch_c%first, first_active_o), &
                      wf%n_v, &
                      one, &
                      rho_ai_bj, &
                      (n_active_v**2)*n_active_o)
!
         call wf%mem%dealloc(g_ai_bC, n_active_v*n_active_o, n_active_v*batch_c%length)

!
      enddo ! batching over c
      write(unit_output,*)'1'
      flush(unit_output)

! 
!     :: Term 2 ::
!     - sum_K g_aij_K * C_K_b(= - sum_K g_ai,Kj * C_bK)
! 
!     Construct g_ai,Kj ordered as g_aij_K
! 
      call wf%mem%alloc(g_ai_Kj, n_active_v*n_active_o, n_active_o*(wf%n_o))
      write(unit_output,*)'2.1'
      flush(unit_output)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_oo(integral_type, g_ai_Kj, &
                        first_active_v, last_active_v, &
                        first_active_o, last_active_o, &
                        1, wf%n_o,                     &
                        first_active_o, last_active_o)
      write(unit_output,*)'2.2'
      flush(unit_output)
! 
!     Reorder g_ai_Kj to g_aij_K
! 
      call wf%mem%alloc(g_aij_K, n_active_v*n_active_o**2, wf%n_o)
      g_aij_K = zero
      write(unit_output,*)'2.3'
      flush(unit_output)
      do a = 1, n_active_v
         do j = 1, n_active_o
            do i = 1, n_active_o
               aij = index_three(a, i, j, n_active_v, n_active_o)
               ai = index_two(a, i, n_active_v)
               do K = 1, wf%n_o
                  Kj = index_two(K, j, wf%n_o)
                  g_aij_K(aij, K) = g_ai_Kj(ai, Kj)
               enddo
            enddo
         enddo
      enddo
      write(unit_output,*)'2'
      flush(unit_output)
      call wf%mem%dealloc(g_ai_Kj, n_active_v*n_active_o, n_active_o*(wf%n_o))
! 
      call wf%mem%alloc(rho_aij_b, n_active_v*n_active_o**2, n_active_v)
! 
!      Add term to rho_ai_bj
! 
      call dgemm('N', 'T',                    &
                  n_active_v*n_active_o**2,   &
                  n_active_v,                 &
                  wf%n_o,                     &
                  -one,                       &
                  g_aij_K,                    &
                  n_active_v*n_active_o**2,   &
                  c_a_i(first_active_v,1),    &
                  wf%n_v,                     &
                  zero,                       &
                  rho_aij_b,                  &
                  n_active_v*n_active_o**2)
! 
      call wf%mem%dealloc(g_aij_K, n_active_v*n_active_o**2, wf%n_o)
! 
! 
!     Reorder into rho_ai_bj
! 
    write(unit_output,*)'3'
      flush(unit_output)
      do a = 1, n_active_v
         do b = 1, n_active_v
            do i = 1, n_active_o
               do j = 1, n_active_o
                  ai = index_two(a, i, n_active_v)
                  bj = index_two(b, j, n_active_v)
                  aij = index_three(a, i, j, n_active_v, n_active_o)
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b)
               enddo
            enddo
         enddo
      enddo
      call wf%mem%dealloc(rho_aij_b, n_active_v*n_active_o**2, n_active_v)
      write(unit_output,*)'4'
      flush(unit_output)
!
   end subroutine jacobian_mlcc2_a2_mlcc2
!
!
   module subroutine jacobian_mlcc2_b2_mlcc2(wf, rho_ai_bj, c_ai_bj)
!!
!!    Jacobian tem B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Calculates the B2 contribution to the jacobi transformation,
!!
!!       B2:   ε_ij^ab*c_ai,bj.
!!
!!
      implicit none
!  
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:)   :: c_ai_bj 
      real(dp), dimension(:,:)   :: rho_ai_bj
!
!     Local routine variables
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0
!
!     Active space variables
!  
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
       do i = 1, n_active_o
!
          do a = 1, n_active_v
!
             ai = index_two(a, i, n_active_v)
!
             do j = 1, n_active_o
!
                do b = 1, n_active_v
!
                   bj = index_two(b, j, n_active_v)
!
                   rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) +  c_ai_bj(ai, bj)*(wf%fock_diagonal(wf%n_o + a + first_active_v - 1, 1) &
                                                  - wf%fock_diagonal(i + first_active_o - 1, 1) &
                                                  + wf%fock_diagonal(wf%n_o + b + first_active_v - 1, 1) &
                                                  - wf%fock_diagonal(j + first_active_o - 1, 1))
!
                enddo
             enddo
          enddo
       enddo
!
   end subroutine jacobian_mlcc2_b2_mlcc2
!
!
    module subroutine cvs_rho_aibj_projection_mlcc2(wf, vec_aibj)
!!
!!    Rho projection for CVS (MLCC2),
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(:, :) :: vec_aibj
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital
!
!     Active space variables
! 
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
! 
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
! 
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1         
!
      do i = 1, n_active_o
       do j = 1, n_active_o
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
             do a = 1, n_active_v
                do b = 1, n_active_v
                   ai = index_two(a, i, n_active_v)
                   bj = index_two(b, j, n_active_v)
                   aibj = index_packed(ai, bj)

                   vec_aibj(ai, bj) = zero
                enddo
             enddo
          endif
       enddo
    enddo
!
  end subroutine cvs_rho_aibj_projection_mlcc2
!
end submodule jacobian