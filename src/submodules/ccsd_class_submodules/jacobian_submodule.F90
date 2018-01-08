submodule (ccsd_class) jacobian
!
!!
!!    Jacobian submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    jacobian_transformation: performs the transformation by the CCSD
!!                             Jacobian matrix A, placing the result in the
!!                             incoming vector. 
!!
!!    jacobian_ccsd_x1:        adds the X1 term to the transformed singles vector; x = a, b, c, d
!!    jacobian_ccsd_x2:        adds the X2 term to the transformed doubles vector; x = a, b, ..., k
!!
!
   implicit none 
!
   character(len=40) :: integral_type
!
contains
!
!
   module subroutine jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
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
!     Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: c_aibj ! c_aibj     
!
!     Local unpacked and reordered vectors 
!
      real(dp), dimension(:,:), allocatable :: rho_a_i   ! rho_ai   = (A c)_ai
      real(dp), dimension(:,:), allocatable :: rho_ai_bj ! rho_ai_bj = (A c)_aibj
!
      real(dp), dimension(:,:), allocatable :: c_ai_bj ! Unpacked c_aibj
      real(dp), dimension(:,:), allocatable :: c_ab_ij ! c_ai_bj, reordered
!
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_sym ! Symmetrized rho_ai_bj, temporary
      real(dp), dimension(:,:), allocatable :: rho_ab_ij     ! rho_ai_bj, reordered
!
!     Indices 
!
      integer(i15) :: a = 0, ab = 0, ai = 0, b = 0 
      integer(i15) :: bj = 0, i = 0, ij = 0, j = 0, aibj = 0
!
!     Timings 
!
      real(dp) :: begin_timer, end_timer
!
      real(dp) :: ccsd_a1_time
      real(dp) :: ccsd_b1_time 
      real(dp) :: ccsd_c1_time
      real(dp) :: ccsd_d1_time
      real(dp) :: ccs_a1_time 
      real(dp) :: ccs_b1_time 
!
      real(dp) :: ccsd_a2_time
      real(dp) :: ccsd_b2_time 
      real(dp) :: ccsd_c2_time 
      real(dp) :: ccsd_d2_time 
      real(dp) :: ccsd_e2_time  
      real(dp) :: ccsd_f2_time
      real(dp) :: ccsd_g2_time 
      real(dp) :: ccsd_h2_time 
      real(dp) :: ccsd_i2_time 
      real(dp) :: ccsd_j2_time 
      real(dp) :: ccsd_k2_time 
!
!     Allocate and zero the transformed vector (singles part)
!
      call wf%mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!     :: CCS contributions to the singles c vector ::  
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
      call cpu_time(end_timer)
      ccs_a1_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
      call cpu_time(end_timer)
      ccs_b1_time = end_timer - begin_timer
!
!     :: CCSD contributions to the transformed singles vector :: 
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_a1(rho_a_i, c_a_i)
      call cpu_time(end_timer)
      ccsd_a1_time = end_timer - begin_timer
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
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_b1(rho_a_i, c_ai_bj) 
      call cpu_time(end_timer)
      ccsd_b1_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_c1(rho_a_i, c_ai_bj)
      call cpu_time(end_timer)
      ccsd_c1_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_d1(rho_a_i, c_ai_bj)
      call cpu_time(end_timer)
      ccsd_d1_time = end_timer - begin_timer
!
!     :: CCSD contributions to the transformed doubles vector ::  
!
!     Allocate unpacked transformed vector 
!
      call wf%mem%alloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      rho_ai_bj = zero 
!
!     Contributions from singles vector c 
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_a2(rho_ai_bj, c_a_i)
      call cpu_time(end_timer)
      ccsd_a2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_b2(rho_ai_bj, c_a_i)
      call cpu_time(end_timer)
      ccsd_b2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_c2(rho_ai_bj, c_a_i)
      call cpu_time(end_timer)
      ccsd_c2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_d2(rho_ai_bj, c_a_i)
      call cpu_time(end_timer)
      ccsd_d2_time = end_timer - begin_timer
!
!     Done with singles vector c; overwrite it with 
!     transformed vector for exit
!
      call dcopy((wf%n_o)*(wf%n_v), rho_a_i, 1, c_a_i, 1)   
!
!     Contributions from doubles vector c
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_e2(rho_ai_bj, c_ai_bj)
      call cpu_time(end_timer)
      ccsd_e2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_f2(rho_ai_bj, c_ai_bj)
      call cpu_time(end_timer)
      ccsd_f2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_g2(rho_ai_bj, c_ai_bj) 
      call cpu_time(end_timer)
      ccsd_g2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_h2(rho_ai_bj, c_ai_bj)
      call cpu_time(end_timer)
      ccsd_h2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_i2(rho_ai_bj, c_ai_bj)
      call cpu_time(end_timer)
      ccsd_i2_time = end_timer - begin_timer
!
!     Last two terms are already symmetric (J2 and K2). Perform the symmetrization 
!     rho_ai_bj = P_ij^ab rho_ai_bj now, for convenience 
!
!     Allocate temporary symmetric transformed vector 
!
      call wf%mem%alloc(rho_ai_bj_sym, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
!     In preparation for last two terms, reorder 
!     rho_ai_bj to rho_ab_ij, and c_ai_bj to c_ab_ij
! 
      call wf%mem%alloc(rho_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
      call wf%mem%alloc(c_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
! 
      rho_ab_ij = zero
      c_ab_ij   = zero
! 
      do j = 1, wf%n_o
         do i = 1, wf%n_o
! 
            ij = index_two(i, j, wf%n_o)
! 
            do b = 1, wf%n_v
! 
               bj = index_two(b, j, wf%n_v)
! 
               do a = 1, wf%n_v
! 
                  ai = index_two(a, i, wf%n_v)
                  ab = index_two(a, b, wf%n_v)
! 
                  c_ab_ij(ab, ij)   = c_ai_bj(ai, bj)
                  rho_ab_ij(ab, ij) = rho_ai_bj(ai, bj)
! 
               enddo
            enddo
         enddo
      enddo
! 
      call wf%mem%dealloc(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_j2(rho_ab_ij, c_ab_ij)
      call cpu_time(end_timer)
      ccsd_j2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%jacobian_ccsd_k2(rho_ab_ij, c_ab_ij)
      call cpu_time(end_timer)
      ccsd_k2_time = end_timer - begin_timer
! 
!     Done with reordered doubles c; deallocate 
! 
      call wf%mem%dealloc(c_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
! 
!     Order rho_ab_ij back into rho_ai_bj & divide by 
!     the biorthonormal factor 1 + delta_ai,bj
! 
      call wf%mem%alloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      rho_ai_bj = zero 
! 
      do j = 1, wf%n_o
         do b = 1, wf%n_v
! 
            bj = index_two(b, j, wf%n_v)
! 
            do i = 1, wf%n_o
! 
               ij = index_two(i, j, wf%n_o)
! 
               do a = 1, wf%n_v
! 
                  ab = index_two(a, b, wf%n_v)
                  ai = index_two(a, i, wf%n_v)
! 
                  if (a .eq. b .and. i .eq. j) then 
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
      call wf%mem%dealloc(rho_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     Overwrite the incoming doubles c vector & pack in
!
      c_aibj = zero
      call packin(c_aibj, rho_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Remaining deallocations 
!
      call wf%mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
      call wf%mem%dealloc(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Print timings
!
      if (wf%settings%print_level == 'developer') then 
!
         write(unit_output,'(t3,a/)') 'Breakdown of CCSD Jacobian timings:'
!
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCS  A1 (seconds):', ccs_a1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCS  B1 (seconds):', ccs_b1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD A1 (seconds):', ccsd_a1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD B1 (seconds):', ccsd_b1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD C1 (seconds):', ccsd_c1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD D1 (seconds):', ccsd_d1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD A2 (seconds):', ccsd_a2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD B2 (seconds):', ccsd_b2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD C2 (seconds):', ccsd_c2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD D2 (seconds):', ccsd_d2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD E2 (seconds):', ccsd_e2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD F2 (seconds):', ccsd_f2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD G2 (seconds):', ccsd_g2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD H2 (seconds):', ccsd_h2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD I2 (seconds):', ccsd_i2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD J2 (seconds):', ccsd_j2_time
         write(unit_output,'(t6,a26,f14.8/)') 'Time in CCSD K2 (seconds):', ccsd_k2_time
!
         flush(unit_output)
!
      endif
!
   end subroutine jacobian_ccsd_transformation_ccsd
!
!
  module subroutine jacobian_ccsd_a1_ccsd(wf, rho_a_i, c_a_i)
!!
!!    Jacobian CCSD A1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!!
!!    rho_ai^A1 = sum_ckdl L_lckd (u_li^ca c_dk  - t_li^cd c_ak - t_lk^ad c_ci)
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!    where c_a_i(a,i) = c_ai above. 
!!
      implicit none 
!
      class(ccsd) :: wf
!
!     Vectors sent to the routine 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i   ! c_ai 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_lc_kd ! g_lckd 
      real(dp), dimension(:,:), allocatable :: L_lc_dk ! L_lckd reordered
      real(dp), dimension(:,:), allocatable :: L_lkd_c ! L_lckd
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: X_lc
      real(dp), dimension(:,:), allocatable :: X_i_k
      real(dp), dimension(:,:), allocatable :: X_a_c
!
!     Amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_i_lcd ! t_li^cd 
      real(dp), dimension(:,:), allocatable :: t_a_lkd ! t_lk^ad 
      real(dp), dimension(:,:), allocatable :: u_ai_lc ! u_li^ca reordered 
!
!     Indices 
!
      integer(i15) :: a = 0, c = 0, d = 0, i = 0, k = 0, l = 0
      integer(i15) :: ai = 0, al = 0, dk = 0, ci = 0, cl = 0, di = 0
      integer(i15) :: kc = 0, kd = 0, ld = 0, lc = 0
      integer(i15) :: lkd = 0, lcd = 0, lad
      integer(i15) :: cial = 0, clai = 0, cldi = 0
!
!
!     :: Term 1: sum_ckdl L_lckd u_li^ca c_dk ::
!
!     g_lc_kd = g_lckd 
!
      call wf%mem%alloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_lc_kd)
!
      call wf%mem%alloc(L_lc_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_lc_dk = zero 
!
!     L_lc_dk = L_lckd = 2*g_lckd - g_ldkc = 2*g_lc_kd(lc,kd) - g_lc_kd(ld,kc)
!
      do k = 1, wf%n_o
         do d = 1, wf%n_v
!
            dk = index_two(d, k, wf%n_v)
            kd = index_two(k, d, wf%n_o)
!
            do c = 1, wf%n_v
!
               kc = index_two(k, c, wf%n_o)
!
               do l = 1, wf%n_o
!
                  lc = index_two(l, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
!
                  L_lc_dk(lc, dk) = two*g_lc_kd(lc, kd) - g_lc_kd(ld, kc) ! L_lckd 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_lc_kd 
!
      call wf%mem%dealloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     X_lc = sum_kd L_lckd c_dk = sum_kd L_lc_dk c_dk 
!
      call wf%mem%alloc(X_lc, (wf%n_o)*(wf%n_v), 1)
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
!     Form u_ai_lc = u_li^ca = 2 * t_li^ca - t_il^ca = 2 * t2am(clai,1) - t2am(cial,1)
!
!     Allocate & read double amplitudes 
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(u_ai_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      u_ai_lc = zero
!
      do c = 1, wf%n_v
         do l = 1, wf%n_o
!
            lc = index_two(l, c, wf%n_o)
            cl = index_two(c, l, wf%n_v)
!
            do i = 1, wf%n_o
!
               ci = index_two(c, i, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  al = index_two(a, l, wf%n_v)
!
                  clai = index_packed(cl, ai)
                  cial = index_packed(ci, al)
!
                  u_ai_lc(ai, lc) = two*(wf%t2am(clai, 1)) - wf%t2am(cial, 1)
!
               enddo
            enddo
         enddo
      enddo
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
      call wf%mem%dealloc(u_ai_lc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%mem%dealloc(X_lc, (wf%n_v)*(wf%n_o), 1)
!
!
!     :: Term 2. - sum_ckdl L_lckd t_li^cd c_ak ::
!
!     Reorder amplitudes to t_lcd_i = t_li^cd 
!
      call wf%mem%alloc(t_i_lcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      t_i_lcd = zero
!
      do i = 1, wf%n_o
         do d = 1, wf%n_v
!
            di = index_two(d, i, wf%n_v)
!
            do c = 1, wf%n_v
               do l = 1, wf%n_o
!
                  cl = index_two(c, l, wf%n_v)
!
                  cldi = index_packed(cl, di)
!
                  lcd = index_three(l, c, d, wf%n_o, wf%n_v)
!
                  t_i_lcd(i, lcd) = wf%t2am(cldi, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Done with doubles amplitudes: deallocate 
!
      call wf%destruct_double_amplitudes
!
!     Calculate X_i_k = sum_cdl L_lcd_k t_i_lcd
!
      call wf%mem%alloc(X_i_k, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_i_lcd,              &
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
      call wf%mem%dealloc(X_i_k, wf%n_o, wf%n_o)
!
!     :: Term 3: - sum_ckdl L_lckd t_lk^ad c_ci ::
!
!     Reorder to L_lkd_c = L_lckd = L_k_lcd
!
      call wf%mem%alloc(L_lkd_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      L_lkd_c = zero
!
      do c = 1, wf%n_v
         do d = 1, wf%n_v
            do k = 1, wf%n_o
               do l = 1, wf%n_o
!
                  lkd = index_three(l, k, d, wf%n_o, wf%n_o)
                  lc = index_two(l, c, wf%n_o)
                  dk = index_two(d, k, wf%n_v)
!
                  L_lkd_c(lkd, c) = L_lc_dk(lc, dk) ! L_lckd 
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(L_lc_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder amplitudes to t_a_lkd = t_lk^ad 
!
      call wf%mem%alloc(t_a_lkd, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      t_a_lkd = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
!
            do l = 1, wf%n_o
!
               lkd = index_three(l, k, d, wf%n_o, wf%n_o)
!
               do a = 1, wf%n_v
!
                  lad = index_three(l, a, d, wf%n_o, wf%n_v)
!
                  t_a_lkd(a, lkd) = t_i_lcd(k, lad)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(t_i_lcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Calculate X_a_c = sum_kdl t_a_lkd L_lkd_c 
!
      call wf%mem%alloc(X_a_c, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  t_a_lkd,              &
                  wf%n_v,               &
                  L_lkd_c,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_a_c,                &
                  wf%n_v)
!
      call wf%mem%dealloc(L_lkd_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      call wf%mem%dealloc(t_a_lkd, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
      call wf%mem%dealloc(X_a_c, wf%n_v, wf%n_v)
!
   end subroutine jacobian_ccsd_a1_ccsd
!
!
module subroutine jacobian_ccsd_b1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!    rho_ai^B1 = sum_bj F_jb (2*c_ai_bj  -  c_aj_bi) 
!!              = sum_bj F_jb v_ai_jb
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!    where c_a_i(a,i) = c_ai above. 
!!

      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))   :: c_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      real(dp), dimension(:,:), allocatable :: v_ai_jb 
      real(dp), dimension(:,:), allocatable :: F_bj
!
      integer(i15) :: a = 0, b = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ai = 0, aj = 0, bi = 0, jb = 0, bj = 0
!
!     Construct v_ai_jb = 2*c_aibj - c_ajbi
!
      call wf%mem%alloc(v_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      v_ai_jb = zero
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
            jb = index_two(j, b, wf%n_o)
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
                  v_ai_jb(ai, jb) = two*c_ai_bj(ai, bj) - c_ai_bj(aj, bi)        
!
               enddo
            enddo
         enddo
      enddo
!
!     sum_bj F_jb*v_ai_jb
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
      call wf%mem%dealloc(v_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_b1_ccsd
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD C1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!    rho_ai^C1 = - sum_bjk L_jikb c_aj_bk 
!!              = - sum_bjk (2*g_jikb - g_kijb) c_aj_bk 
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!    where c_ai_bj(ai,bj) = c_aibj above. 
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))   :: c_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      real(dp), dimension(:,:), allocatable :: g_ji_kb
      real(dp), dimension(:,:), allocatable :: L_jbk_i
      real(dp), dimension(:,:), allocatable :: c_a_jbk
!
      integer(i15) :: i = 0, j = 0, k = 0
      integer(i15) :: a = 0, b = 0 
!
      integer(i15) :: ji = 0, ik = 0, ki = 0
      integer(i15) :: jb = 0, kb = 0
      integer(i15) :: aj = 0, bk = 0
!
      integer(i15) :: jbk = 0
!
      integer(i15) :: ajbk = 0
!
!     Construct the integral g_ji_kb = sum_J L_ji_J * L_kb_J
!
      call wf%mem%alloc(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_ji_kb)
!
!     Constructing L_jikb = 2*g_jikb - g_jbki
!
      call wf%mem%alloc(L_jbk_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)  
      L_jbk_i = zero
!   
      do b = 1, wf%n_v
         do k = 1, wf%n_o
!
            kb = index_two(k, b, wf%n_o)
!
            do j = 1, wf%n_o
!
               jb = index_two(j, b, wf%n_o)
               jbk = index_three(j, b, k, wf%n_o, wf%n_v)
!
               do i = 1, wf%n_o
!
                  ji = index_two(j, i, wf%n_o)
                  ki = index_two(k, i, wf%n_o)
!
                  L_jbk_i(jbk, i) = two*g_ji_kb(ji, kb) - g_ji_kb(ki,jb)

!
               enddo
            enddo
         enddo
      enddo 
!        
      call wf%mem%dealloc(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     - sum_bjk L_jbk_i * c_a_jbk
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  c_ai_bj,                & ! c_aj_bk = c_a_jbk
                  wf%n_v,                 &
                  L_jbk_i,                &
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  rho_a_i,                &
                  wf%n_v)
!
         call wf%mem%dealloc(L_jbk_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
!
   module subroutine jacobian_ccsd_d1_ccsd(wf, rho_a_i, c_bi_cj)
!!
!!    Jacobian CCSD D1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!    rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!    where c_bi_cj(bi,cj) = c_bicj above. 
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)) :: c_bi_cj 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i 
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
      real(dp), dimension(:,:), allocatable :: L_a_cjb
!
      integer(i15) :: a = 0, b = 0, c = 0
      integer(i15) :: i = 0, j = 0 
!
      integer(i15) :: jb = 0, jc = 0, bi = 0, cj = 0 
      integer(i15) :: ca = 0, ab = 0, ac = 0
!
      integer(i15) :: cjb = 0
!
      integer(i15) :: bicj = 0 
!
!     Prepare for batching over index a
! 
      required = wf%get_vvov_required_mem() + (wf%n_v**3)*(wf%n_o)*dp
!
!     Initialize batching variable 
!
      call batch_a%init(wf%n_v)
      call wf%mem%num_batch(batch_a, required)         
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
         call wf%mem%alloc(g_ab_jc, (batch_a%length)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_ov(integral_type, &
                           g_ab_jc,       &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Construct L_abjc ordered as L_a_cjb
!
         call wf%mem%alloc(L_a_cjb, (batch_a%length), ((wf%n_v)**2)*(wf%n_o))
         L_a_cjb = zero
!       
         do b = 1, wf%n_v
            do j = 1, wf%n_o
               do c = 1, wf%n_v
!
                  cjb = index_three(c, j, b, wf%n_v, wf%n_o)
!
                  jb = index_two(j, b, wf%n_o)
                  jc = index_two(j, c, wf%n_o)
!
                  do a = 1, batch_a%length 
!
                     ac = index_two(a, c, batch_a%length )
                     ab = index_two(a, b, batch_a%length )
!
                     L_a_cjb(a, cjb) = two*g_ab_jc(ab, jc) - g_ab_jc(ac, jb)
!
                  enddo
!
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_ab_jc, (batch_a%length)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
         call dgemm('N', 'N',                   & 
                     batch_a%length,            &
                     wf%n_o,                    &
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     L_a_cjb,                   &
                     batch_a%length,            &
                     c_bi_cj,                   & ! c_cjb_i
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     rho_a_i(batch_a%first, 1), &
                     wf%n_v)
!
         call wf%mem%dealloc(L_a_cjb, batch_a%length, (wf%n_o)*(wf%n_v)**2)
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
!!    where c_a_i(a,i) = c_ai above. 
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_ai_kj ! g_aikj 
      real(dp), dimension(:,:), allocatable :: g_k_aij ! g_aikj reordered 
      real(dp), dimension(:,:), allocatable :: g_ai_bc ! g_aibc, batching over b
      real(dp), dimension(:,:), allocatable :: g_bc_ai ! g_aibc = g_bcai, batching over b 
!
      real(dp), dimension(:,:), allocatable :: rho_b_aij ! rho_ai_bj, term 1 (see below)
      real(dp), dimension(:,:), allocatable :: rho_aib_j ! rho_ai_bj, term 2 (batching over b)
!
      integer(i15) :: a = 0, ai = 0, aij = 0, b = 0, bj = 0, i = 0, j = 0
      integer(i15) :: k = 0, kj = 0, c = 0, cb = 0, aib = 0, bc = 0
!
      logical :: reorder 
!
      real(dp) :: begin_reorder
      real(dp) :: end_reorder
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
      call wf%mem%alloc(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_oo(integral_type, g_ai_kj)
!
!     Reorder to g_k_aij = g_aikj = g_ai_kj
!
      call wf%mem%alloc(g_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
      g_k_aij = zero
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               do k = 1, wf%n_o
!
                  kj = index_two(k, j, wf%n_o)
!  
                  g_k_aij(k, aij) = g_ai_kj(ai, kj) ! g_aikj
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
!     Calculate rho_b_aij = - sum_k c_bk g_aikj = - sum_k c_a_i(b, k) g_k_aij(k, aij)
!
      call wf%mem%alloc(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                &
                  wf%n_v,               &
                  g_k_aij,              &
                  wf%n_o,               &
                  zero,                 &
                  rho_b_aij,            &
                  wf%n_v)
!
      call wf%mem%dealloc(g_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
!     Add rho_b_aij to rho_ai_bj 
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_aij(b, aij)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     
!     :: Term 2. rho_ai_bj =+ sum_c g_aibc c_cj ::
!
!     We do the matrix multiplication as g_aib_c c_cj,
!     batching over the b index.
!
      call wf%mem%alloc(rho_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o) ! rho_ai_bj formed in batch 
      rho_aib_j = zero
!
      required = wf%get_vvvo_required_mem()
!
!     Initialize batching variable  
!
      call batch_b%init(wf%n_v)
      call wf%mem%num_batch(batch_b, required)
!
      do current_b_batch = 1, batch_b%num_batches 
!
!        Determine limits for the current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
         call wf%mem%alloc(g_ai_bc, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_b%length))
!
         integral_type = 'electronic_repulsion'
         call wf%get_vo_vv(integral_type, &
                           g_ai_bc,       &
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
                     rho_aib_j(aib_offset,1),            &
                     (wf%n_o)*(wf%n_v)**2)
!
          call wf%mem%dealloc(g_ai_bc, (wf%n_v)*(wf%n_o), wf%n_v*(batch_b%length))
!
      enddo ! End of batches over b
!
!     Now that rho_aib_j has been made in batches, add it to rho_ai_bj 
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
                  aib = index_three(a, i, b, wf%n_v, wf%n_o)
!  
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aib_j(aib, j)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call wf%mem%dealloc(rho_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
   end subroutine jacobian_ccsd_a2_ccsd
!
!
   module subroutine jacobian_ccsd_b2_ccsd(wf, rho_ai_bj, c_a_i)
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
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: t_c_aij ! t_ij^ac 
      real(dp), dimension(:,:), allocatable :: t_aib_k ! t_ik^ab 
! 
      real(dp), dimension(:,:), allocatable :: X_k_aij ! An intermediate
      real(dp), dimension(:,:), allocatable :: X_k_j   ! An intermediate 
!
      real(dp), dimension(:,:), allocatable :: rho_b_aij ! rho_ai_bj, unordered, term 1
      real(dp), dimension(:,:), allocatable :: rho_aib_j ! rho_ai_bj, unordered, term 2
!
      integer(i15) :: a = 0, c = 0, i = 0, j = 0, b = 0, k = 0
      integer(i15) :: ai = 0, cj = 0, aicj = 0, aij = 0, bj = 0
      integer(i15) :: bk = 0, aibk = 0, aib = 0
!
!     :: Term 1. - sum_kc F_kc t_ij^ac c_bk ::
!
!     Allocate & read the amplitudes from disk 
!
      call wf%read_double_amplitudes
!
!     Order the amplitudes as t_c_aij = t_ij^ac 
!
      call wf%mem%alloc(t_c_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               do c = 1, wf%n_v
!
                  cj = index_two(c, j, wf%n_v)
!
                  aicj = index_packed(ai, cj)
!
                  t_c_aij(c, aij) = wf%t2am(aicj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Form the intermediate X_k_aij = sum_c F_k_c t_c_aij 
!
      call wf%mem%alloc(X_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  wf%fock_ia,           &
                  wf%n_o,               &
                  t_c_aij,              &
                  wf%n_v,               &
                  zero,                 &
                  X_k_aij,              &
                  wf%n_o)
!
      call wf%mem%dealloc(t_c_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     Form rho_b_aij = sum_k c_a_i(b,k) X_k_aij(k,aij)
!
      call wf%mem%alloc(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
                  rho_b_aij,            &
                  wf%n_v)
!
      call wf%mem%dealloc(X_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
!     Add rho_b_aij to rho_ai_bj 
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_aij(b, aij)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 2. - sum_kc F_kc t_ik^ab c_cj :: 
!
!     Form X_k_j = sum_c F_kc c_cj = sum_c fock_ia(k,c) c_a_i(c,j)
!
      call wf%mem%alloc(X_k_j, wf%n_o, wf%n_o)
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
!     Order the amplitudes as t_aib_k = t_ik^ab 
!
      call wf%mem%alloc(t_aib_k, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
      t_aib_k = zero
!
      do k = 1, wf%n_o
         do b = 1, wf%n_v
!
            bk = index_two(b, k, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                  aibk = index_packed(ai, bk)
!
                  t_aib_k(aib, k) = wf%t2am(aibk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate doubles amplitudes    
!
      call wf%destruct_double_amplitudes
! 
!     Form rho_aib_j = - sum_k t_aib_k X_k_j
!     (Interpret rho_ai_bj as rho_aib_j)
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 & 
                  t_aib_k,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  X_k_j,                &
                  wf%n_o,               &
                  one,                  &
                  rho_ai_bj,            &
                  (wf%n_o)*(wf%n_v)**2)
!
      call wf%mem%dealloc(X_k_j, wf%n_o, wf%n_o)
      call wf%mem%dealloc(t_aib_k, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
   end subroutine jacobian_ccsd_b2_ccsd
!
!
 module subroutine jacobian_ccsd_c2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD C2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^C2 = sum_kcl g_ljkc (t_ki^ac c_bl + t_li^bc c_ak + t_lk^ba c_ci)
!!                 - sum_kcl L_ljkc (t_il^ab c_ck + t_ik^ac c_bl)
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_lj_kc ! g_ljkc
      real(dp), dimension(:,:), allocatable :: g_kj_lc ! g_ljkc 
      real(dp), dimension(:,:), allocatable :: L_lj_ck ! L_ljkc 
!
      real(dp), dimension(:,:), allocatable :: X_lj_ai ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_kj_bi ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_kjl_i ! An intermediate, term 3
      real(dp), dimension(:,:), allocatable :: X_kl_ij ! X_kjl_i reordered
      real(dp), dimension(:,:), allocatable :: X_lj    ! An intermediate, term 4
      real(dp), dimension(:,:), allocatable :: Y_lj_ai ! An intermediate, term 5
!
      real(dp), dimension(:,:), allocatable :: t_kc_ai ! t_ki^ac
      real(dp), dimension(:,:), allocatable :: t_lc_bi ! t_li^bc
      real(dp), dimension(:,:), allocatable :: t_ba_kl ! t_lk^ba
      real(dp), dimension(:,:), allocatable :: t_aib_l ! t_il^ab 
      real(dp), dimension(:,:), allocatable :: t_ck_ai ! t_ik^ac
!
      real(dp), dimension(:,:), allocatable :: rho_b_jai ! rho_ai_bj, term 1 & 5
      real(dp), dimension(:,:), allocatable :: rho_a_jbi ! rho_ai_bj, term 2
      real(dp), dimension(:,:), allocatable :: rho_ba_ij ! rho_ai_bj, term 3 
!
      integer(i15) :: ci = 0, i = 0, k = 0, kc = 0, c = 0, ak = 0, ai = 0, akci = 0
      integer(i15) :: a = 0, j = 0, b = 0, bj = 0, jai = 0, kj = 0, jbi = 0, l = 0, lc = 0
      integer(i15) :: lj = 0, blci = 0, bl = 0, bi = 0, ij = 0, kjl = 0, kl = 0, ba = 0
      integer(i15) :: blak = 0, jk = 0, ck = 0, aib = 0, aibl = 0
!
!     :: Term 1. sum_kcl g_ljkc t_ki^ac c_bl :: 
!
!     Form g_lj_kc = g_ljkc
!
      call wf%mem%alloc(g_lj_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_ov(integral_type, g_lj_kc)
!
!     Read the amplitudes from disk
!
      call wf%read_double_amplitudes
!
!     Order as t_kc_ai = t_ki^ac 
!
      call wf%mem%alloc(t_kc_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_kc_ai = zero 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do c = 1, wf%n_v
!
               ci = index_two(c, i, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ak = index_two(a, k, wf%n_v)
                  kc = index_two(k, c, wf%n_o)
!
                  akci = index_packed(ak, ci)
!
                  t_kc_ai(kc, ai) = wf%t2am(akci, 1) ! t_ki^ac 
!
               enddo
            enddo
         enddo
      enddo
!
!     Form X_lj_ai = sum_ck g_lj_kc t_kc_ai
!
      call wf%mem%alloc(X_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
!
!     Calculate rho_b_jai = sum_l c_bl X_lj_ai 
!     (Interpret the X array as an X_l_jai object in the matrix multiplication)
!
      call wf%mem%alloc(rho_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_a_i,                &
                  wf%n_v,               &
                  X_lj_ai,              & ! "X_l_jai"
                  wf%n_o,               & 
                  zero,                 &
                  rho_b_jai,            &
                  wf%n_v)
!
      call wf%mem%dealloc(X_lj_ai, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Add rho_b_jai to rho_ai_bj 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do j = 1, wf%n_o
!
               jai = index_three(j, a, i, wf%n_o, wf%n_v)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_jai(b, jai)
!
               enddo
            enddo
         enddo
      enddo
!
!     Final deallocations for term 1 (keep g_lj_kc for later use)
!
      call wf%mem%dealloc(rho_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 2. sum_kcl g_ljkc t_li^bc c_ak ::
!
!     Reorder to g_kj_lc = g_lj_kc = g_ljkc     
!
      call wf%mem%alloc(g_kj_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      g_kj_lc = zero 
!
      do c = 1, wf%n_v
         do l = 1, wf%n_o 
!
            lc = index_two(l, c, wf%n_o)
!
            do j = 1, wf%n_o
!
               lj = index_two(l, j, wf%n_o)
!
               do k = 1, wf%n_o
!
                  kc = index_two(k, c, wf%n_o)
                  kj = index_two(k, j, wf%n_o)
!
                  g_kj_lc(kj, lc) = g_lj_kc(lj, kc) ! g_ljkc
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_lj_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_kj_bi = sum_lc g_kj_lc t_lc_bi
!
      call wf%mem%alloc(X_kj_bi, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(t_kc_ai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
!     Calculate rho_a_jbi = sum_k c_ak X_kj_bi 
!
      call wf%mem%alloc(rho_a_jbi, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_a_i,                &
                  wf%n_v,               &
                  X_kj_bi,              & ! "X_k_jbi"
                  wf%n_o,               &
                  zero,                 &
                  rho_a_jbi,            &
                  wf%n_v)
!
      call wf%mem%dealloc(X_kj_bi, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Add rho_a_jbi to rho_ai_bj 
!
      do i = 1, wf%n_o
         do b = 1, wf%n_v
            do j = 1, wf%n_o
!
               jbi = index_three(j, b, i, wf%n_o, wf%n_v)
!
               bj = index_two(b, j, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_a_jbi(a, jbi)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocations for term 2 (keep g_kj_lc = g_ljkc)
!
      call wf%mem%dealloc(rho_a_jbi, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 3. sum_kcl g_ljkc t_lk^ba c_ci :: 
!
!     Form the intermediate X_kjl_i = sum_c g_ljkc c_ci = sum_c g_kj_lc c_c_i
!
!     Note: interpret g_kj_lc as g_kjl_c in matrix multiplication.
!
      call wf%mem%alloc(X_kjl_i, (wf%n_o)**3, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**3, &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  g_kj_lc,     & ! "g_kjl_c"
                  (wf%n_o)**3, &
                  c_a_i,       &
                  wf%n_v,      &
                  zero,        &
                  X_kjl_i,     &
                  (wf%n_o)**3)
!
!     Reorder to X_kl_ij = X_kjl_i
!
      call wf%mem%alloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
      X_kl_ij = zero
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do l = 1, wf%n_o
               do k = 1, wf%n_o
!
                  kl = index_two(k, l, wf%n_o)
!
                  kjl = index_three(k, j, l, wf%n_o, wf%n_o)
!
                  X_kl_ij(kl, ij) = X_kjl_i(kjl, i) 
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(X_kjl_i, (wf%n_o)**3, wf%n_o)
!
!     Order amplitudes as t_kl_ba = t_lk^ba 
!
      call wf%mem%alloc(t_ba_kl, (wf%n_v)**2, (wf%n_o)**2)
      t_ba_kl = zero 
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            ba = index_two(b, a, wf%n_v)
!
            do l = 1, wf%n_o
!
               bl = index_two(b, l, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ak = index_two(a, k, wf%n_v)
                  kl = index_two(k, l, wf%n_o) 
!
                  blak = index_packed(bl, ak)
!
                  t_ba_kl(ba, kl) = wf%t2am(blak, 1) ! t_lk^ba 
!
               enddo
            enddo
         enddo
      enddo
!
!     Calculate rho_ba_ij = sum_kcl g_ljkc t_lk^ba c_ci 
!                         = sum kl ( sum_c g_ljkc c_ci ) t_lk^ba 
!                         = sum_kl t_ba_kl X_kl_ij
!
      call wf%mem%alloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
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
      call wf%mem%dealloc(t_ba_kl, (wf%n_v)**2, (wf%n_o)**2)
      call wf%mem%dealloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Add rho_ba_ij into rho_ai_bj 
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               do b = 1, wf%n_v
!
                  ba = index_two(b, a, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!  
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_ba_ij(ba, ij)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocations for term 3 (keep g_kj_lc = g_ljkc)
!
      call wf%mem%dealloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!
!     :: Term 4. - sum_kcl L_ljkc t_il^ab c_ck ::
!
!     Form L_lj_ck = L_ljkc = 2 * g_ljkc - g_lckj = 2 * g_ljkc - g_kjlc = g_kj_lc(kj,lc) - g_kj_lc(lj,kc)
!
      call wf%mem%alloc(L_lj_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
            kc = index_two(k, c, wf%n_o)
!
            do j = 1, wf%n_o
!
               kj = index_two(k, j, wf%n_o)
!
               do l = 1, wf%n_o
!
                  lj = index_two(l, j, wf%n_o)
                  lc = index_two(l, c, wf%n_o)
!
                  L_lj_ck(lj, ck) = two*g_kj_lc(kj, lc) - g_kj_lc(lj, kc) ! L_ljkc
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kj_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Calculate the intermediate X_lj = sum_ck L_lj_ck c_ck 
!
      call wf%mem%alloc(X_lj, (wf%n_o)**2, 1)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)**2,       &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lj_ck,           &
                  (wf%n_o)**2,       &
                  c_a_i,             & ! Interpret as "c_ck"
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lj,              &
                  (wf%n_o)**2)
!
!     Order the amplitudes as t_aib_l = t_il^ab 
!
      call wf%mem%alloc(t_aib_l, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      t_aib_l = zero 
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
                  aibl = index_packed(ai, bl)
!
                  t_aib_l(aib, l) = wf%t2am(aibl, 1) ! t_il^ab 
!
               enddo
            enddo
         enddo
      enddo
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
                  t_aib_l,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  X_lj,                 & ! "X_l_j"
                  wf%n_o,               &
                  one,                  &
                  rho_ai_bj,            & ! "rho_aib_j"
                  (wf%n_o)*(wf%n_v)**2)
!
      call wf%mem%dealloc(t_aib_l, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call wf%mem%dealloc(X_lj, (wf%n_o)**2, 1)
!
!     :: Term 5. - sum_kcl L_ljkc t_ik^ac c_bl ::
!
!     Square up the amplitudes, t_ck_ai(ck, ai) = t_ki^ca = t_ik^ac 
!
      call wf%mem%alloc(t_ck_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ck_ai = zero
!
      call squareup(wf%t2am, t_ck_ai, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate Y_lj_ai = sum_kc L_ljkc t_ik^ac = sum_kc L_lj_ck t_ck_ai 
!
      call wf%mem%alloc(Y_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lj_ck,           &
                  (wf%n_o)**2,       &
                  t_ck_ai,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_lj_ai,           &
                  (wf%n_o)**2)
!
      call wf%mem%dealloc(t_ck_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate rho_b_jai =+ - sum_l c_bl Y_lj_ai
!
!     Note: interpret Y_lj_ai as Y_l_jai in the matrix multiplication
!
      call wf%mem%alloc(rho_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                &
                  wf%n_v,               &
                  Y_lj_ai,              & ! "Y_l_jai" 
                  wf%n_o,               &
                  zero,                 &
                  rho_b_jai,            &
                  wf%n_v)
!
!     Add rho_b_jai to rho_ai_bj 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do j = 1, wf%n_o
!
               jai = index_three(j, a, i, wf%n_o, wf%n_v)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_b_jai(b, jai)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocations and cleanup
!
      call wf%mem%dealloc(Y_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(L_lj_ck, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(rho_b_jai, (wf%n_v), (wf%n_v)*(wf%n_o)**2)
!
      call wf%destruct_double_amplitudes 
!
   end subroutine jacobian_ccsd_c2_ccsd
!
!
  module subroutine jacobian_ccsd_d2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD D2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
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
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_bd_kc ! g_kcbd 
      real(dp), dimension(:,:), allocatable :: g_cd_kb ! g_kcbd reordered
      real(dp), dimension(:,:), allocatable :: g_ckb_d ! g_kcbd reordered 
      real(dp), dimension(:,:), allocatable :: L_ckb_d ! L_kcbd = 2 g_kcbd - g-kdbc
!
      real(dp), dimension(:,:), allocatable :: t_ij_cd ! t_ij^cd
      real(dp), dimension(:,:), allocatable :: t_dk_aj ! t_kj^ad 
      real(dp), dimension(:,:), allocatable :: t_ai_ck ! t_ik^ca 
      real(dp), dimension(:,:), allocatable :: t_d_aij ! t_ij^ad 
!
      real(dp), dimension(:,:), allocatable :: X_ij_kb ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_k_ijb ! The above intermediate, reordered 
      real(dp), dimension(:,:), allocatable :: X_i_dkb ! An intermediate, term 2
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
!  ! Move into batching loop
      call wf%read_double_amplitudes
!
!     Determine batch size, etc.
!     (Redo estimate once loop is done)
!
      required = wf%get_vvov_required_mem() + &
                  dp*(max((wf%n_o)*(wf%n_v**3), &
                     (wf%n_o**2)*(wf%n_v**2) + (wf%n_o**3)*(wf%n_v), &
                      2*(wf%n_o**3)*(wf%n_v), &
                      3*(wf%n_o**2)*(wf%n_v**2)))
!
!     Initialize batching variable 
!
      call batch_b%init(wf%n_v)
      call wf%mem%num_batch(batch_b, required)
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
         call wf%mem%alloc(g_bd_kc, (wf%n_v)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_ov(integral_type, &
                           g_bd_kc,       &
                           batch_b%first, &
                           batch_b%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Reorder to g_cd_kb = g_kc_db = g_kcbd 
!
         call wf%mem%alloc(g_cd_kb, (wf%n_v)**2, (wf%n_o)*(batch_b%length))
!
         do b = 1, batch_b%length
            do k = 1, wf%n_o
!
               kb = index_two(k, b, wf%n_o)
!
               do d = 1, wf%n_v
!
                  bd = index_two(b, d, batch_b%length)
!
                  do c = 1, wf%n_v
!
                     kc = index_two(k, c, wf%n_o)
                     cd = index_two(c, d, wf%n_v)
!
                     g_cd_kb(cd, kb) = g_bd_kc(bd, kc) ! g_kcbd 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_bd_kc, (wf%n_v)*(batch_b%length), (wf%n_o)*(wf%n_v))
!
!        Order amplitudes as t_ij_cd = t_ij^cd 
!
         call wf%mem%alloc(t_ij_cd, (wf%n_o)**2, (wf%n_v)**2)
         t_ij_cd = zero
!
         do d = 1, wf%n_v
            do c = 1, wf%n_v
!
               cd = index_two(c, d, wf%n_v)
!
               do j = 1, wf%n_o
!
                  dj = index_two(d, j, wf%n_v)
!
                  do i = 1, wf%n_o
!
                     ci = index_two(c, i, wf%n_v)
                     ij = index_two(i, j, wf%n_o)
!
                     cidj = index_packed(ci, dj)
!
                     t_ij_cd(ij, cd) = wf%t2am(cidj, 1) ! t_ij^cd
!
                  enddo
               enddo
            enddo
         enddo
!
!        Form intermediate X_ij_kb = sum_cd g_kcdb t_ij^cd 
!                                  = sum_cd t_ij_cd g_cd_kb
!
         call wf%mem%alloc(X_ij_kb, (wf%n_o)**2, (wf%n_o)*(batch_b%length))
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
         call wf%mem%dealloc(t_ij_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!        sum_kcd g_kcbd t_ij^cd c_ak = sum_k X_ij_kb c_ak 
!        Reorder to X_k_ijb = X_ij_kb 
!
         call wf%mem%alloc(X_k_ijb, (wf%n_o), (batch_b%length)*(wf%n_o)**2)
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
         call wf%mem%dealloc(X_ij_kb, (wf%n_o)**2, (wf%n_o)*(batch_b%length))
!
!        Form rho_a_ijb = - sum_k c_ak X_k_ijb = - sum_k c_a_i(a,k) X_k_ijb(k, ijb)
!
         call wf%mem%alloc(rho_a_ijb, wf%n_v, (batch_b%length)*(wf%n_o)**2)
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
         call wf%mem%dealloc(X_k_ijb, wf%n_o, (batch_b%length)*(wf%n_o)**2)
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
         call wf%mem%dealloc(rho_a_ijb, wf%n_v, (batch_b%length)*(wf%n_o)**2)
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
         call wf%mem%alloc(X_i_dkb, wf%n_o, (wf%n_v)*(wf%n_o)*(batch_b%length))
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
                     X_i_dkb,                            &
                     wf%n_o)
!
!        sum_kcd g_kcbd t_kj^ad c_ci = sum_kd (sum_c c_ci g_kcbd) t_kj^ad
!                                    = sum_kd X_i_dkb t_kj^ad 
!                                    = sum_kd X_ib_dk t_dk_aj
!
!        Reorder to X_ib_dk = X_i_dkb 
!
         call wf%mem%alloc(X_ib_dk, (wf%n_o)*(batch_b%length), (wf%n_v)*(wf%n_o))
!
         do k = 1, wf%n_o
            do d = 1, wf%n_v
!
               dk = index_two(d, k, wf%n_v)
!
               do b = 1, batch_b%length
!
                  dkb = index_three(d, k, b, wf%n_v, wf%n_o)
!
                  do i = 1, wf%n_o
!
                     ib = index_two(i, b, wf%n_o)
!
                     X_ib_dk(ib, dk) = X_i_dkb(i, dkb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(X_i_dkb, wf%n_o, (wf%n_o)*(wf%n_v)*(batch_b%length))
!
!        Order the amplitudes as t_dk_aj = t_kj^ad 
!
         call wf%mem%alloc(t_dk_aj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_dk_aj = zero
!
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               aj = index_two(a, j, wf%n_v)
!
               do k = 1, wf%n_o
!
                  ak = index_two(a, k, wf%n_v)
!
                 do d = 1, wf%n_v
!
                     dj = index_two(d, j, wf%n_v)
                     dk = index_two(d, k, wf%n_v)
!
                     akdj = index_packed(ak, dj)
!
                     t_dk_aj(dk, aj) = wf%t2am(akdj, 1) ! t_kj^ad 
!
                  enddo
               enddo
            enddo
         enddo
!
!        Calculate rho_ib_aj = - sum_kcd g_kcbd t_kj^ad c_ci
!                            = - sum_dk X_ib_dk t_dk_aj
!
         call wf%mem%alloc(rho_ib_aj, (wf%n_o)*(batch_b%length), (wf%n_o)*(wf%n_v))
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
         call wf%mem%dealloc(X_ib_dk, (wf%n_o)*(batch_b%length), (wf%n_o)*(wf%n_v))
         call wf%mem%dealloc(t_dk_aj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call wf%mem%dealloc(rho_ib_aj, (wf%n_o)*(batch_b%length), (wf%n_v)*(wf%n_o))
!
!
!        :: Term 3. - sum_kcd g_kcbd t_ik^ca c_dj ::
!
!        sum_d g_kcbd c_dj = sum_d g_cd_kb c_dj
!
!        Reorder integrals to g_cd_kb to g_ckb_d
!
         call wf%mem%alloc(g_ckb_d, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_v)
         g_ckb_d = zero
!
         do d = 1, wf%n_v
            do b = 1, batch_b%length
               do k = 1, wf%n_o
!
                  kb = index_two(k, b, wf%n_o)
!
                  do c = 1, wf%n_v
!
                     cd = index_two(c, d, wf%n_v)
!
                     ckb = index_three(c, k, b, wf%n_v, wf%n_o)
!
                     g_ckb_d(ckb, d) = g_cd_kb(cd, kb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_cd_kb, (wf%n_v)**2, (wf%n_o)*(batch_b%length))
!
!        Form the intermediate X_ckb_j = sum_d g_kcbd c_dj = sum_d g_ckb_d c_d_j 
!
         call wf%mem%alloc(X_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_ckb_d,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_a_i,                              &
                     wf%n_v,                             &
                     zero,                               &
                     X_ckb_j,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ca
!
         call wf%mem%alloc(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
         t_ai_ck = zero 
!
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               do i = 1, wf%n_o
!
                  ci = index_two(c, i, wf%n_v)
!
                  do a = 1, wf%n_v
!
                     ak = index_two(a, k, wf%n_v)
                     ai = index_two(a, i, wf%n_v)
!
                     ciak = index_packed(ci, ak)
!
                     t_ai_ck(ai, ck) = wf%t2am(ciak, 1) ! t_ik^ca 
!
                  enddo
               enddo
            enddo
         enddo
!
!        Form rho_aib_j = -sum_kcd g_kcbd t_ik^ca c_dj = sum_ck t_ai_ck X_ckb_j 
!
!        Note: X_ckb_j is interpreted as X_ck_bj in the matrix multiplication.
!        Note: rho_aib_j is interpreted as rho_ai_bj in the matrix multiplication.
!
         call wf%mem%alloc(rho_aib_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
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
         call wf%mem%dealloc(X_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
         call wf%mem%dealloc(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
         call wf%mem%dealloc(rho_aib_j, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_o)
!
!
!        :: Term 4.  sum_kcd L_kcbd t_ik^ac c_dj :: 
!
!        sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj 
!
!        Form L_ckb_d = L_kcbd = 2 * g_kcbd - g_kdbc = 2 * g_ckb_d(ckb, d) - g_ckb_d(dkb, c)
!
         call wf%mem%alloc(L_ckb_d, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_v)
         L_ckb_d = zero
!
         do d = 1, wf%n_v
            do b = 1, batch_b%length
               do k = 1, wf%n_o
!
                  dkb = index_three(d, k, b, wf%n_v, wf%n_o)
!
                  do c = 1, wf%n_v
!
                     ckb = index_three(c, k, b, wf%n_v, wf%n_o)
!
                     L_ckb_d(ckb, d) = two*g_ckb_d(ckb, d) - g_ckb_d(dkb, c)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_ckb_d, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_v)
!
!        Form the intermediate Y_ckb_j = sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj 
!
         call wf%mem%alloc(Y_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), & 
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     L_ckb_d,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_a_i,                              &
                     wf%n_v,                             &
                     zero,                               &
                     Y_ckb_j,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ac = t2am(aick, 1)
!
         call wf%mem%alloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_ai_ck = zero 
!
         call squareup(wf%t2am, t_ai_ck, (wf%n_o)*(wf%n_v))
!
!        Form rho_aib_j =  sum_ck t_ai_ck Y_ckb_j 
!
!        Note: we interpret Y_ckb_j as Y_ck_bj in the matrix multiplication
!        Note: we interpret rho_aib_j as rho_ai_bj in the matrix multiplication
!
         call wf%mem%alloc(rho_aib_j, (wf%n_o)*(wf%n_v)*(batch_b%length), wf%n_o)
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
         call wf%mem%dealloc(Y_ckb_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
         call wf%mem%dealloc(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
         call wf%mem%dealloc(rho_aib_j, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_o)
!
!
!        :: Term 5.  sum_kcd L_kcbd t_ij^ad c_ck ::
!
!        Form the intermediate X_1,bd = sum_ck c_ck L_kcbd = sum_ck c_1,ck L_ckb_d
!
!        Note: L_ckb_d is interpreted as L_ck_bd in the matrix multiplication 
!        Note: c_a_i is interpreted as c_1,ai in the matrix multiplication 
!
         call wf%mem%alloc(X_bd, 1, (batch_b%length)*(wf%n_v))
!  
         call dgemm('N','N',                    &
                     1,                         &
                     (batch_b%length)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v),         &
                     one,                       &
                     c_a_i,                     & ! "c_1,ai"
                     1,                         &
                     L_ckb_d,                   & ! "L_ck_bd"
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     X_bd,                      &
                     1)
!
!        Order amplitudes as t_d_aij = t_ij^ad 
!
         call wf%mem%alloc(t_d_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
         t_d_aij = zero
!
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  do d = 1, wf%n_v
!
                     dj = index_two(d, j, wf%n_v)
!
                     aidj = index_packed(ai, dj)
!
                     t_d_aij(d, aij) = wf%t2am(aidj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Form rho_b_aij =  sum_kcd L_kcbd t_ij^ad c_ck
!                       =  sum_d X_bd t_d_aij
!
!        Note: X_bd is interpreted as X_b_d in the matrix multiplication
!
         call wf%mem%alloc(rho_b_aij, batch_b%length, (wf%n_v)*(wf%n_o)**2)
!
         call dgemm('N','N',               &
                     batch_b%length,       &
                     (wf%n_v)*(wf%n_o)**2, &
                     wf%n_v,               &
                     one,                  & 
                     X_bd,                 & ! "X_b_d"
                     batch_b%length,       &
                     t_d_aij,              &
                     wf%n_v,               &
                     zero,                 &
                     rho_b_aij,            &
                     batch_b%length)
!
         call wf%mem%dealloc(X_bd, 1, (batch_b%length)*(wf%n_v))
         call wf%mem%dealloc(t_d_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
         call wf%mem%dealloc(rho_b_aij, batch_b%length, (wf%n_v)*(wf%n_o)**2)
         call wf%mem%dealloc(L_ckb_d, (wf%n_v)*(wf%n_o)*(batch_b%length), wf%n_v)
!
     enddo ! End of batches over b 
!
!    Destroy amplitudes from memory
!
     call wf%destruct_double_amplitudes
!
   end subroutine jacobian_ccsd_d2_ccsd
!
!
    module subroutine jacobian_ccsd_e2_ccsd(wf, rho_ai_bj, c_ai_ck)
!!
!!    Jacobian CCSD E2 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    rho_ai_bj^E2 = 2 sum_dlck t_bj,dl * L_kc,ld * c_ai,ck 
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_ck
!
      real(dp), dimension(:,:), allocatable :: t_dl_bj
      real(dp), dimension(:,:), allocatable :: g_kc_ld
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: X_ck_bj
!
      integer(i15) :: c = 0, d = 0, k = 0, l = 0
      integer(i15) :: ck = 0, dl = 0, kc = 0, kd = 0, lc = 0, ld = 0
!
!     Read T2 amplitudes from disk
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_dl_bj = zero
!
      call squareup(wf%t2am, t_dl_bj, (wf%n_o)*(wf%n_v))
!
      call wf%destruct_double_amplitudes
!
!     Construct g_kcld 
!
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!
!     Construct L_kc,ld ordered as L_ck_dl
!
      call wf%mem%alloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ck_dl = zero
!
      do l = 1, wf%n_o
         do d = 1, wf%n_v
!
            dl = index_two(d, l, wf%n_v)
            ld = index_two(l, d, wf%n_o)
!
            do k = 1, wf%n_o
!
               kd = index_two(k, d, wf%n_o)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  kc = index_two(k, c, wf%n_o)
                  lc = index_two(l, c, wf%n_o)
!
                  L_ck_dl(ck, dl) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Intermediate X_ck_bj = sum_dl L_ck_dl * t_dl_bj
!
      call wf%mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD F2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
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
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_kc_ld
!
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: L_d_lck
      real(dp), dimension(:,:), allocatable :: L_l_ckd
!
      real(dp), dimension(:,:), allocatable :: c_dl_bj
      real(dp), dimension(:,:), allocatable :: c_clk_b
      real(dp), dimension(:,:), allocatable :: c_ckd_j
!
      real(dp), dimension(:,:), allocatable :: t_ai_ck
      real(dp), dimension(:,:), allocatable :: t_aij_d
      real(dp), dimension(:,:), allocatable :: t_aib_l
 !               
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_l_j
      real(dp), dimension(:,:), allocatable :: X_ck_bj
!
      real(dp), dimension(:,:), allocatable :: rho_aij_b
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

      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!   
      call wf%mem%alloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ck_dl = zero
!
!        Construct L_kc,ld ordered as L_ck_dl
!             
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            kc = index_two(k, c, wf%n_o)
            ck = index_two(c, k, wf%n_v)
!
            do d = 1, wf%n_v
               do l = 1, wf%n_o
!
                  lc = index_two(l, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
                  dl = index_two(d, l, wf%n_v)
                  kd = index_two(k, d, wf%n_o)
!
                  L_ck_dl(ck, dl) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 1: - sum_ckld t_ai,ck * L_kc,ld * c_bl,dj ::
!
!     Reorder c_bl_dj as c_dl_bj
!
      call wf%mem%alloc(c_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      c_dl_bj = zero
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
                  c_dl_bj(dl,bj) = c_ai_bj(bl,dj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(c_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ai_ck = zero
!
      call squareup(wf%t2am, t_ai_ck, (wf%n_o)*(wf%n_v))
      call wf%destruct_double_amplitudes
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
      call wf%mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 2: - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!
!     Construct L_ck,dl reordered as L_d_clk
!
      call wf%mem%alloc(L_d_lck, wf%n_v, (wf%n_v)*((wf%n_o)**2))
      L_d_lck = zero
!
      do k = 1, wf%n_o
         do l = 1, wf%n_o
            do c = 1, wf%n_v
!
               lck = index_three(l, c, k, wf%n_o, wf%n_v)
!
               kc = index_two(k, c, wf%n_o)
               lc = index_two(l, c, wf%n_o)
!
               do d = 1, wf%n_v
!
                  ld = index_two(l, d, wf%n_o)
                  kd = index_two(k, d, wf%n_o)
!
                  L_d_lck(d, lck) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Y_d_b = sum_clk L_d_lck * c_b_lck
!     Here dgemm is tricked to believe that c_bl_ck is c_b_lck 
!
      call wf%mem%alloc(Y_d_b, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  wf%n_v,               &
                  ((wf%n_o)**2)*wf%n_v, &
                  one,                  &
                  L_d_lck,              &
                  wf%n_v,               &
                  c_ai_bj,              &
                  wf%n_v,               &
                  zero,                 &
                  Y_d_b,                &
                  wf%n_v)
!
      call wf%mem%dealloc(L_d_lck, wf%n_v, (wf%n_v)*((wf%n_o)**2)) 
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      t_aij_d = zero
!
!     Reorder T2 amplitudes
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do d = 1, wf%n_v
!
               dj = index_two(d, j, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  aidj = index_packed(ai,dj)
!
                  t_aij_d(aij, d) = wf%t2am(aidj,1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_double_amplitudes
!
      call wf%mem%alloc(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!     rho_aij_b = sum_d t_aij_d*Y_d_b
!
      call dgemm('N','N',                 &
                  ((wf%n_o)**2)*(wf%n_v), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  t_aij_d,                &
                  ((wf%n_o)**2)*(wf%n_v), &
                  Y_d_b,                  &
                  wf%n_v,                 &
                  zero,                   &
                  rho_aij_b,              &
                  ((wf%n_o)**2)*(wf%n_v))
!
      call wf%mem%dealloc(t_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      call wf%mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     Adding term 2 to rho_ai_bj
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               ai = index_two(a, i, wf%n_v)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b) 
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!     :: Term 3: - sum_ckdl t_ai,bl * L_kc,ld * c_ck,dj ::
!
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!   
      call wf%mem%alloc(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2))
      L_l_ckd = zero
!
!     Construct L_kc,dl ordered as L_l_ckd
!             
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            kc = index_two(k, c, wf%n_o)
!
            do d = 1, wf%n_v
!
               ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
               do l = 1, wf%n_o
!
                  lc = index_two(l, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
                  kd = index_two(k, d, wf%n_o)
!
                  L_l_ckd(l, ckd) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(Z_l_j, wf%n_o, wf%n_o)
!
!     Z_l_j = sum_ckd L_l_ckd * c_ckd_l 
!  
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  one,                  &
                  L_l_ckd,              &
                  wf%n_o,               &
                  c_ai_bj,              & ! c_ai_bj(ck,dl)= c_ckd_l
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_l_j,                &
                  wf%n_o)
!
      call wf%mem%dealloc(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2))
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_aib_l, (wf%n_o)*((wf%n_v)), wf%n_o*(wf%n_v))
      t_aib_l = zero
      call squareup(wf%t2am, t_aib_l, wf%n_o*(wf%n_v))

!
      call wf%destruct_double_amplitudes
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
      call wf%mem%dealloc(t_aib_l, (wf%n_o)*((wf%n_v)), wf%n_o*(wf%n_v))
      call wf%mem%dealloc(Z_l_j, wf%n_o, wf%n_o)
! 
   end subroutine jacobian_ccsd_f2_ccsd
!
!
   module subroutine jacobian_ccsd_g2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!    Jacobian CCSD G2 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
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
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!         
      real(dp), dimension(:,:), allocatable :: g_kc_ld
!
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: L_d_clk
      real(dp), dimension(:,:), allocatable :: L_l_ckd
!
      real(dp), dimension(:,:), allocatable :: c_ai_ck
      real(dp), dimension(:,:), allocatable :: c_aib_l
      real(dp), dimension(:,:), allocatable :: c_aij_d
!
      real(dp), dimension(:,:), allocatable :: t_dl_bj
      real(dp), dimension(:,:), allocatable :: t_clk_b
      real(dp), dimension(:,:), allocatable :: t_ckd_j
 !               
      real(dp), dimension(:,:), allocatable :: X_ck_bj
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_l_j
!
      real(dp), dimension(:,:), allocatable :: rho_aij_b
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
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!   
      call wf%mem%alloc(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ck_dl = zero
!
!     Construct L_kc_ld ordered as L_ck_dl
!             
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            kc = index_two(k, c, wf%n_o)
            ck = index_two(c, k, wf%n_v)
!
            do d = 1, wf%n_v
               do l = 1, wf%n_o
!
                  lc = index_two(l, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
                  dl = index_two(d, l, wf%n_v)
                  kd = index_two(k, d, wf%n_o)
!
                  L_ck_dl(ck, dl) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_bl_dj as t_dl_bj
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
                  t_dl_bj(dl, bj) = wf%t2am(bldj, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_double_amplitudes
!
      call wf%mem%alloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Term 2: - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!
!     Reorder L_ck_dl to L_d_clk
!
      call wf%mem%alloc(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2))
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
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ck,bl as t_clk_b
!        
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
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
                  t_clk_b(clk, b) = wf%t2am(ckbl, 1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_double_amplitudes
!
!     Y_d_b = sum_clk L_d_clk * c_clk_b 
!
      call wf%mem%alloc(Y_d_b, wf%n_v, wf%n_v)
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
      call wf%mem%dealloc(t_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      call wf%mem%dealloc(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2)) 
!
      call wf%mem%alloc(c_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      c_aij_d = zero
!
!     Reorder c_ai_dj to c_aij_d
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do d = 1, wf%n_v
!
               dj = index_two(d, j, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                  c_aij_d(aij, d) = c_ai_bj(ai, dj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_double_amplitudes
!
      call wf%mem%alloc(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!     rho_aij_b = sum_d c_aij_d * Y_d_b
!
      call dgemm('N','N',                 &
                  ((wf%n_o)**2)*(wf%n_v), &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  c_aij_d,                &
                  ((wf%n_o)**2)*(wf%n_v), &
                  Y_d_b,                  &
                  wf%n_v,                 &
                  zero,                   &
                  rho_aij_b,              &
                  ((wf%n_o)**2)*(wf%n_v))
!
      call wf%mem%dealloc(c_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      call wf%mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     Adding term 2 to rho_ai_bj
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               ai = index_two(a, i, wf%n_v)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_aij_b(aij, b) 
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!     :: Term 3: - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl ::
!
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!   
      call wf%mem%alloc(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2))
      L_l_ckd = zero
!
!     Construct L_kc_ld ordered as  L_l_ckd
!             
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            kc = index_two(k, c, wf%n_o)
!
            do d = 1, wf%n_v
               do l = 1, wf%n_o
!
                  lc = index_two(l, c, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
                  kd = index_two(k, d, wf%n_o)
!
                  ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
                  L_l_ckd(l, ckd) = two*g_kc_ld(kc, ld) - g_kc_ld(kd, lc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ck,dj to t_ckd_j 
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_ckd_j, ((wf%n_v))*(wf%n_o), wf%n_o*(wf%n_v))
      t_ckd_j = zero
      call squareup(wf%t2am, t_ckd_j, wf%n_o*(wf%n_v))
!
      call wf%destruct_double_amplitudes
!
      call wf%mem%alloc(Z_l_j, wf%n_o, wf%n_o)
!
!     Z_l_j = sum_ckd L_l_ckd*t_ckd_j
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  one,                  &
                  L_l_ckd,              &
                  wf%n_o,               &
                  t_ckd_j,              &
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_l_j,                &
                  wf%n_o)
!
      call wf%mem%dealloc(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2)) 
      call wf%mem%dealloc(t_ckd_j, ((wf%n_v))*(wf%n_o), wf%n_o*(wf%n_v))
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
      call wf%mem%dealloc(Z_l_j, wf%n_o, wf%n_o)
! 
   end subroutine jacobian_ccsd_g2_ccsd
!
      module subroutine jacobian_ccsd_h2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD H2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^H2 =  sum_ckdl t_ci,ak * g_kc,ld * c_bl,dj 
!!                     + sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!!                
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
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
         call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_kc_ld)
!
!        t_ak,ci ordered as t_ai_kc
!  

         call wf%read_double_amplitudes
!
         call wf%mem%alloc(t_ai_kc, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         t_ai_kc = zero
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               do c = 1, wf%n_v
                  do k = 1, wf%n_o
!
                     ci = index_two(c, i, wf%n_v)
                     ak = index_two(a, k, wf%n_v)
                     kc = index_two(k, c, wf%n_o)
!
                     akci = index_packed(ak, ci)
!
                     t_ai_kc(ai, kc) = wf%t2am(akci, 1)
! 
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_double_amplitudes
!  
         call wf%mem%alloc(X_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call wf%mem%dealloc(t_ai_kc, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%mem%alloc(c_ld_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         c_ld_bj = zero
!
!        Reorder c_bl,dj as c_ld_bj
!
         do l = 1, wf%n_o
            do b = 1, wf%n_v
!
               bl = index_two(b, l, wf%n_v)
!
               do j = 1, wf%n_o
!
                  bj = index_two(b, j, wf%n_v)
!
                  do d = 1, wf%n_v
!
                     dj = index_two(d, j, wf%n_v)
                     ld = index_two(l, d, wf%n_o)
!
                     c_ld_bj(ld, bj) = c_ai_bj(bl, dj)
!
                  enddo
               enddo
            enddo
         enddo
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
         call wf%mem%dealloc(c_ld_bj, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         call wf%mem%dealloc(X_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        :: Term 2: sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!
!        Construct g_kc_ld
!
         call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_ov(integral_type, g_kc_ld)
!
!        Reorder g_kc_ld to g_lc_kd 
!
         call wf%mem%alloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         g_lc_kd = zero
!
         do c = 1, wf%n_v
            do d = 1, wf%n_v
               do k = 1, wf%n_o
!
                  kc = index_two(k, c, wf%n_o)
                  kd = index_two(k, d, wf%n_o)
!
                  do l = 1, wf%n_o
!
                     lc = index_two(l, c, wf%n_o)
                     ld = index_two(l, d, wf%n_o)
!
                     g_lc_kd(lc, kd) = g_kc_ld(kc, ld)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        t_al,cj ordered as t_aj_lc
!  

         call wf%read_double_amplitudes
!
         call wf%mem%alloc(t_aj_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_aj_lc = zero
!
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               aj = index_two(a, j, wf%n_v)
!
               do c = 1, wf%n_v
!
                  cj = index_two(c, j, wf%n_v)
!
                  do l = 1, wf%n_o
!
                     al = index_two(a, l, wf%n_v)
                     lc = index_two(l, c, wf%n_o)
!
                     alcj = index_packed(al, cj)
!
                     t_aj_lc(aj, lc) = wf%t2am(alcj, 1)
! 
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_double_amplitudes
!
         call wf%mem%alloc(Y_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call wf%mem%dealloc(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call wf%mem%dealloc(t_aj_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!  
!        Reorder c_bk,di as c_kd_bi
!
         call wf%mem%alloc(c_kd_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         c_kd_bi = zero
!
         do i = 1, wf%n_o
            do k = 1, wf%n_o
               do d = 1, wf%n_v
!
                  kd = index_two(k, d, wf%n_o)
                  di = index_two(d, i, wf%n_v)
!
                  do b = 1, wf%n_v
!
                     bk = index_two(b, k, wf%n_v)
                     bi = index_two(b, i, wf%n_v)
!
                     c_kd_bi(kd, bi) = c_ai_bj(bk, di)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%alloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call wf%mem%dealloc(c_kd_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call wf%mem%dealloc(Y_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder into rho_ai_bj
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               do j = 1, wf%n_o
!
                  aj = index_two(a, j, wf%n_v)
!
                  do b = 1, wf%n_v
!
                     bi = index_two(b, i, wf%n_v)
                     bj = index_two(b, j, wf%n_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aj_bi(aj, bi) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      end subroutine jacobian_ccsd_h2_ccsd
!
!
   module subroutine jacobian_ccsd_i2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!     Jacobian CCSD I2 
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
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
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!
      real(dp), dimension(:,:), allocatable :: c_aij_c 
      real(dp), dimension(:,:), allocatable :: c_aib_k
      real(dp), dimension(:,:), allocatable :: c_ai_ck
      real(dp), dimension(:,:), allocatable :: c_aj_ck
!
      real(dp), dimension(:,:), allocatable :: rho_aij_b
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
!     Reorder c_ai,cj to c_aij_c
!
      call wf%mem%alloc(c_aij_c, (wf%n_v)*((wf%n_o)**2), wf%n_v)
      c_aij_c = zero
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               do c = 1, wf%n_v
!
                  cj = index_two(c, j, wf%n_v)
!
                  c_aij_c(aij, c) = c_ai_bj(ai, cj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%alloc(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!     rho_ai_bj += sum_c F_bc * c_ai,cj = sum_c c_aij_c(aij,c) F_ab(b,c) = sum_c c_aij_c(aij,c) F_ab^T(c,b)
!
      call dgemm('N','T',                 & 
                  (wf%n_v)*((wf%n_o)**2), &
                  wf%n_v,                 & 
                  wf%n_v,                 & 
                  one,                    & 
                  c_aij_c,                & 
                  (wf%n_v)*((wf%n_o)**2), &
                  wf%fock_ab,             & 
                  wf%n_v,                 & 
                  zero,                   & 
                  rho_aij_b,              & 
                  (wf%n_v)*((wf%n_o)**2))
!
      call wf%mem%dealloc(c_aij_c, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!     Reorder rho_aij_b into rho_ai_bj
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!  
               ai  = index_two(a, i, wf%n_v)
!
               aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aij_b(aij, b)
!
               enddo
            enddo
         enddo
      enddo
!
!
      call wf%mem%dealloc(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
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
      call wf%mem%alloc(g_bj_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_ov(integral_type, g_bj_kc)
!
!     Reordering g_bj_kc to g_ck_bj
!
      call wf%mem%alloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ck_bj = zero
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v) 
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  kc = index_two(k, c, wf%n_o)
!
                  g_ck_bj(ck, bj) = g_bj_kc(bj, kc)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_bj_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%alloc(c_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      c_ai_ck = zero
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
            do a = 1, wf%n_v
!
               ak = index_two(a, k, wf%n_v)
               ai = index_two(a, i, wf%n_v)
! 
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  ci = index_two(c, i, wf%n_v) 
!
                  c_ai_ck(ai, ck) = c_ai_bj(ak, ci)
!
               enddo
            enddo
         enddo
      enddo
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
      call wf%mem%dealloc(c_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ck_bj = zero
!
!     Start batching over c
! 
      required = wf%get_vvoo_required_mem()
!    
!     Initialize batching variable 
!
      call batch_c%init(wf%n_v)
      call wf%mem%num_batch(batch_c, required)         
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
         call wf%mem%alloc(g_bc_kj, (wf%n_v)*(batch_c%length), (wf%n_o)**2)
         g_bc_kj = zero
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_oo(integral_type, &
                           g_bc_kj,       &
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
         call wf%mem%dealloc(g_bc_kj, (wf%n_v)*(batch_c%length), (wf%n_o)**2)
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
      call wf%mem%alloc(c_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      c_aj_ck = zero
!
      do j = 1, wf%n_o
         do k = 1, wf%n_o
            do a = 1, wf%n_v
!
               ak = index_two(a, k, wf%n_v)
               aj = index_two(a, j, wf%n_v)
! 
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
                  cj = index_two(c, j, wf%n_v) 
!
                  c_aj_ck(aj, ck) = c_ai_bj(ak, cj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%alloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%mem%dealloc(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(c_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder rho_aj_bi into rho_ai_bj
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!  
            ai  = index_two(a, i, wf%n_v)
!                 
            do j = 1, wf%n_o
!                
               aj = index_two(a, j, wf%n_v)
!
               do b = 1, wf%n_v          
!
                  bj = index_two(b, j, wf%n_v)
                  bi = index_two(b, i, wf%n_v)
!
                  rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aj_bi(aj, bi)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_i2_ccsd
!
!
   module subroutine jacobian_ccsd_j2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!     Jacobian CCSD J2 
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl 
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!             
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: c_ab_ij
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
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!
      call wf%mem%alloc(g_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
      g_kl_cd = zero
!
!     Reorder g_kc_ld to g_kl_cd
!
      do c = 1, wf%n_v
         do d = 1, wf%n_v
!
            cd = index_two(c, d, wf%n_v)
!
            do k = 1, wf%n_o
!
               kc = index_two(k, c, wf%n_o)
!
               do l = 1, wf%n_o
! 
                  kl = index_two(k, l, wf%n_o)
                  ld = index_two(l, d, wf%n_o)
!
                  g_kl_cd(kl, cd) = g_kc_ld(kc, ld)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reordered T2 amplitudes        
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(t_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
      t_ab_ij = zero
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do b = 1, wf%n_v
!
               bj = index_two(b, j, wf%n_v)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  ab = index_two(a, b, wf%n_v)
!  
                  aibj = index_packed(ai, bj)
!
                  t_ab_ij(ab, ij) = wf%t2am(aibj, 1)
!
               enddo
            enddo
         enddo
      enddo

      call wf%destruct_double_amplitudes
!
!     X_kl_ij = g_kl_cd * t_cd_ij
!
      call wf%mem%alloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
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
!        rho_ab_ij += c_ab_kl * X_kl_ij
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
      call wf%mem%dealloc(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
      call wf%mem%dealloc(g_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
      call wf%mem%dealloc(t_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
   end subroutine jacobian_ccsd_j2_ccsd
!
!
   module subroutine jacobian_ccsd_k2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!    Jacobian CCSD K2 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
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
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: c_ab_ij
!
      real(dp), dimension(:,:), allocatable :: g_ki_lj
      real(dp), dimension(:,:), allocatable :: g_kl_ij
      real(dp), dimension(:,:), allocatable :: g_ac_bd
      real(dp), dimension(:,:), allocatable :: g_ab_cd
!
      real(dp), dimension(:,:), allocatable :: rho_batch_ab_ij
!
      integer(i15) :: i = 0, j = 0, k = 0, l = 0  
      integer(i15) :: a = 0, b = 0, c = 0, d = 0  
!
      integer(i15) :: ab = 0, bd = 0, cd = 0, full_ab = 0, ac = 0
      integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0
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
      call wf%mem%alloc(g_ki_lj, (wf%n_o)**2, (wf%n_o)**2)
! 
      integral_type = 'electronic_repulsion'
      call wf%get_oo_oo(integral_type, g_ki_lj)
!
!     Reorder g_ki_lj to g_kl_ij
!
      call wf%mem%alloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
      g_kl_ij = zero
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do k = 1, wf%n_o
!
               ki = index_two(k, i, wf%n_o)
!
               do l = 1, wf%n_o
! 
                  kl = index_two(k, l, wf%n_o)
                  lj = index_two(l, j, wf%n_o)
!
                  g_kl_ij(kl, ij) = g_ki_lj(ki, lj) 
!                     
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ki_lj, (wf%n_o)**2, (wf%n_o)**2)
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
      call wf%mem%dealloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Prepare for batching over a and b
!
!     ::  sum_cd g_ac,bd * c_ci,dj ::
!
      required = wf%get_vvvv_required_mem() + dp*(wf%n_v**4)
!
!     Initialize batching variables 
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call wf%mem%num_two_batch(batch_a, batch_b, required)
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
            call wf%mem%alloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!           g_ca_db = sum_J L_ca_J*L_db_J
!     
            integral_type = 'electronic_repulsion'
            call wf%get_vv_vv(integral_type, &
                                 g_ac_bd,       &
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
!           Reorder g_ca_db into g_ab_cd 
!           (Here, g_ab_cd = g_acbd = g_ca_db.)
!
            call wf%mem%alloc(g_ab_cd, (batch_a%length)*(batch_b%length), (wf%n_v)**2) 
            g_ab_cd = zero
!
            do b = 1, batch_b%length
               do a = 1, batch_a%length
!
                  ab = index_two(a, b, batch_a%length)
!
                  do d = 1, wf%n_v
!
                     bd = index_two(b, d, batch_b%length)
!
                     do c = 1, wf%n_v
!
                        ac = index_two(a, c, batch_a%length)
                        cd = index_two(c, d, wf%n_v)
!
                        g_ab_cd(ab, cd) = g_ac_bd(ac, bd) ! = g_acbd 
!
                     enddo
                  enddo
               enddo
            enddo
!
            call wf%mem%dealloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length)) 
!
            call wf%mem%alloc(rho_batch_ab_ij, (batch_a%length)*(batch_b%length), (wf%n_o)**2)
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
            call wf%mem%dealloc(g_ab_cd, (batch_a%length)*(batch_b%length), (wf%n_v)**2)
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
            call wf%mem%dealloc(rho_batch_ab_ij,  (batch_a%length)*(batch_b%length), (wf%n_o)**2) 
!
         enddo ! End batches of b 
      enddo ! End batches of a 
!         
   end subroutine jacobian_ccsd_k2_ccsd
!
!
end submodule jacobian