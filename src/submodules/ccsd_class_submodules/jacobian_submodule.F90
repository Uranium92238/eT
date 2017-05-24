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
!
   implicit none 
!
!
contains
!
!
   module subroutine jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
!!
!!    Jacobian transformation (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
      implicit none
!
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: c_aibj ! c_aibj     
!
      real(dp), dimension(:,:), allocatable :: rho_a_i  ! rho_ai   = (A c)_ai
      real(dp), dimension(:,:), allocatable :: rho_aibj ! rho_aibj = (A c)_aibj
!
      call allocator(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!     CCS contributions from the singles c vector 
!
      call wf%jacobian_ccs_a1(c_a_i, rho_a_i)
      call wf%jacobian_ccs_b1(c_a_i, rho_a_i)
!
!     CCSD contributions from the singles c vector 
!
      call wf%jacobian_ccsd_a1(c_a_i, rho_a_i)
!
      call deallocator(rho_a_i, wf%n_v, wf%n_o)
!
   end subroutine jacobian_ccsd_transformation_ccsd
!
!
   module subroutine jacobian_ccsd_a1_ccsd(wf, c_a_i, rho_a_i)
!!
!!    Jacobian CCSD A1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!!
!!    rho_ai^A1 = sum_ckdl L_lckd (u_li^ca c_dk  - t_li^cd c_ak - t_lk^ad c_ci)
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i   ! c_ai 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      real(dp), dimension(:,:), allocatable :: L_lc_J  ! L_lc^J
      real(dp), dimension(:,:), allocatable :: g_lc_kd ! g_lckd 
      real(dp), dimension(:,:), allocatable :: L_lc_dk ! L_lckd reordered
!
      real(dp), dimension(:,:), allocatable :: u_ai_lc ! u_li^ca reordered 
!
      real(dp), dimension(:,:), allocatable :: X_lc ! An intermediate
!
      integer(i15) :: a = 0, ai = 0, al = 0, c = 0, ci = 0, cial = 0, clai = 0
      integer(i15) :: d = 0, dk = 0, i = 0, k = 0, kc = 0, kd = 0, l = 0, lc = 0
      integer(i15) :: ld = 0, cl = 0
!
!
!     :: Term 1: sum_ckdl L_lckd u_li^ca c_dk ::
!
!
!     Construct L_lc_dk = L_lckd = 2 * g_lckd - g_ldkc 
!
      call allocator(L_lc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_lc_J)
!
      call allocator(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     g_lc_kd = g_lckd 
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_lc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_lc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_lc_kd,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_lc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call allocator(L_lc_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_lc_dk = zero 
!
!     L_lc_dk = L_lckd
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
      call deallocator(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     X_lc = sum_kd L_lckd c_dk = sum_kd L_lc_dk c_dk 
!
      call allocator(X_lc, (wf%n_o)*(wf%n_v), 1)
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
      call wf%initialize_amplitudes ! Allocate t amplitudes, then set them to zero 
      call wf%read_amplitudes       ! Read the converged amplitudes from disk 
!
      call allocator(u_ai_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      do c = 1, wf%n_v
         do l = 1, wf%n_o
!
            lc = index_two(l, c, wf%n_o)
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
      call dgemm('N','N',            &
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
   end subroutine jacobian_ccsd_a1_ccsd
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, c_aibj, rho_a_i)
!!
!!       Jacobian CCSD B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!       rho_ai^B1 = sum_bj F_jb (2*c_aibj  -  c_ajbi) = sum_bj F_jb v_ai_bj
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj   ! c_aibj
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
         real(dp), dimension(:,:), allocatable :: v_ai_bj
         real(dp), dimension(:,:), allocatable :: F_bj
         real(dp), dimension(:,:), allocatable :: rho_ai
!
         integer(i15) :: a = 0, b = 0
         integer(i15) :: i = 0, j = 0
!
         integer(i15) :: ai = 0, aj = 0, bi = 0, bj = 0
!
         integer(i15) :: aibj = 0, ajbi = 0

!
!        Construct v_ai_bj = 2*c_aibj - c_ajbi
!        Reorder F_j_b to F_bj 
!
         call allocator(v_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call allocator(F_bj, (wf%n_o)*(wf%n_v), 1)
         call allocator(rho_ai, (wf%n_o)*(wf%n_v), 1)
!
         do j = 1, wf%n_o
            do b = 1, wf%n_v
!
               bj = index_two(b, j, wf%n_v)
               F_bj = wf%fock_ia(j, b)
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
                     aibj = index_packed(ai, bj)
                     ajbi = index_packed(aj, bi)
!
                     v_ai_bj(ai, bj) =  two*c_aibj(aibj,1) - c_aibj(ajbi,1)        
!
                  enddo
               enddo
            enddo
         enddo
!
!        sum_bj F_jb v_ai_bj 
! 
         call dgemm('N', 'N',           &
                     (wf%n_o)*(wf%n_v), &
                     1,                 &
                     (wf%n_o)*(wf%n_v), &
                     one,               &
                     v_ai_bj,           &
                     (wf%n_o)*(wf%n_v), &
                     F_bj,              &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &        
                     rho_ai,            &
                     (wf%n_o)*(wf%n_v)) 
!
         call deallocator(v_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(F_bj, (wf%n_o)*(wf%n_v), 1)
!
!        Reorder into rho_a_i
! 
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               rho_a_i(a,i) = rho_ai(ai, 1)
!
            enddo
         enddo

         call deallocator(rho_ai, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, c_aibj, rho_a_i)
!!
!!       Jacobian CCSD C1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!       rho_ai^C1 = - sum_bjk L_jikb c_aj_bk = - sum_bjk (2*g_jikb - g_jbki) c_aj_bk
!!                 = - sum_bjk (2*g_jikb - g_kijb) c_aj_bk 
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj   ! c_aibj
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
         real(dp), dimension(:,:), allocatable :: L_ji_J
         real(dp), dimension(:,:), allocatable :: L_kb_J
         real(dp), dimension(:,:), allocatable :: g_ji_kb
         real(dp), dimension(:,:), allocatable :: L_jkb_i
         real(dp), dimension(:,:), allocatable :: c_a_jkb
!
         integer(i15) :: i = 0, j = 0, k = 0
         integer(i15) :: a = 0, b = 0 
!
         integer(i15) :: ji = 0, ik = 0
         integer(i15) :: jb = 0, kb = 0
         integer(i15) :: aj = 0, bk = 0
!
         integer(i15) :: jkb = 0
!
         integer(i15) :: ajbk = 0
!
!        Construct integral g_ji_kb = sum_J L_ji_J * L_kb_J
!
         call allocator(L_ji_J, (wf%n_o)**2, wf%n_J)
         call wf%get_cholesky_ij(L_ji_J)
!
         call allocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
         call wf%get_cholesky_ia(L_kb_J)
!
         call allocator(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
         call dgemm('N', 'T',           &
                     (wf%n_o)**2,       &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_J),          &
                     one,               &
                     L_ji_J,            &
                     (wf%n_o)**2,       & 
                     L_kb_J,            &      
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     g_ji_kb,           &
                     (wf%n_o)**2)
!
         call deallocator(L_ji_J, (wf%n_o)**2, wf%n_J)
         call deallocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Construct L_jikb ordered as L_jkb_i
!
         call allocator(L_jkb_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)
!
         
            
         do b = 1, wf%n_v
            do k = 1, wf%n_o
!
               kb =index_two(k, b, wf%n_o)
!
               do j = 1, wf%n_o
!
                  jb = index_two(j, b, wf%n_o)
                  jkb = index_three(j, k, b, wf%n_o, wf%n_o)
!
                  do i = 1, wf%n_o
!
                     ji = index_two(j, i, wf%n_o)
                     ik = index_two(i, k, wf%n_o)
!
                     L_jkb_i(jkb, i) = two*g_ji_kb(ji, kb) - g_ji_kb(ik,jb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!        Reorder c_ajbk as c_a_jkb
!
         call allocator(c_a_jkb, wf%n_v, (wf%n_v)*((wf%n_o)**2))
!
         do b = 1, wf%n_v
            do k = 1, wf%n_o
!
               bk = index_two(b, k, wf%n_v)
!
               do j = 1, wf%n_o
!
                  jkb = index_three(j, k, b, wf%n_o, wf%n_o)
!
                  do a = 1, wf%n_v
!
                     aj = index_two(a, j, wf%n_v)
!
                     ajbk = index_packed(aj, bk)
!
                     c_a_jkb(a, jkb) = c_aibj(ajbk,1)
!
                  enddo
               enddo
            enddo
         enddo
!
!        - sum_bjk L_jkb_i * c_a_jkb
!
         call dgemm('N', 'N',                &
                     wf%n_v,                 &
                     wf%n_o,                 &
                     (wf%n_v)*((wf%n_o)**2), &
                     -one,                   &
                     c_a_jkb,                &
                     wf%n_v,                 &
                     L_jkb_i,                &
                     (wf%n_v)*((wf%n_o)**2), &
                     zero,                   &
                     rho_a_i,                &
                     wf%n_v)
!
         call deallocator(c_a_jkb, wf%n_v, (wf%n_v)*((wf%n_o)**2))
         call deallocator(L_jkb_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)
         
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
    module subroutine jacobian_ccsd_d1_ccsd(wf, c_bicj, rho_a_i)
!!
!!    Jacobian CCSD D1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!    rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_t2am, 1)   :: c_bicj  ! c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
!     Variables for batching
!
      integer(i15) :: required = 0, available = 0
      integer(i15) :: batch_dimension = 0, max_batch_length = 0, n_batch = 0
      integer(i15) :: a_batch = 0, a_first = 0, a_last = 0, a_length = 0
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_ba_J ! Reordered L_ab_J
      real(dp), dimension(:,:), allocatable :: L_jc_J
      real(dp), dimension(:,:), allocatable :: g_ba_jc ! Reordered g_ab_jc
      real(dp), dimension(:,:), allocatable :: L_a_bjc
      real(dp), dimension(:,:), allocatable :: c_bjc_i
!
      integer(i15) :: a = 0, b = 0, c = 0
      integer(i15) :: i = 0, j = 0 
!
      integer(i15) :: jb = 0, jc = 0, bi = 0, cj = 0 
      integer(i15) :: ba = 0, ca = 0
!
      integer(i15) :: bjc = 0
!
      integer(i15) :: bicj = 0 
!
!     Prepare for batching over index a
! 
      rho_a_i = zero
! 
      required = max(2*(wf%n_J)*((wf%n_v)**2) + 2*(wf%n_J)*(wf%n_v)*(wf%n_o), &
                       (wf%n_J)*((wf%n_v)**2) + (wf%n_J)*(wf%n_v)*(wf%n_o) + ((wf%n_v)**3)*(wf%n_o), &
                        2*((wf%n_v)**3)*(wf%n_o))
!     
      required = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
         call allocator(L_ba_J, (wf%n_v)*a_length, wf%n_J)
         call wf%get_cholesky_ab(L_ba_J, a_first, a_last, a_length*(wf%n_v), .true.)
!
         call allocator(L_jc_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call wf%get_cholesky_ia(L_jc_J)
!
!        g_abjc = sum_J L_ab_J * L_jc_J ordered as g_ba_jc
!
         call allocator(g_ba_jc, a_length*(wf%n_v), (wf%n_v)*(wf%n_o))
!  
         call dgemm('N', 'T', &
                     (wf%n_v)*a_length, &
                     (wf%n_v)*(wf%n_o), &
                     wf%n_J,            &
                     one,               &
                     L_ba_J,            &
                     (wf%n_v)*a_length, &
                     L_jc_J,            &
                     (wf%n_v)*(wf%n_o), &
                     zero,              &
                     g_ba_jc,           &
                     (wf%n_v)*a_length)
!
         call deallocator(L_jc_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call deallocator(L_ba_J, (wf%n_v)*a_length, wf%n_J)
!
!        Construct L_abjc ordered as L_a_bjc
!        Reorder c_bicj to c_bjc_i
!
         call allocator(L_a_bjc, a_length, ((wf%n_v)**2)*(wf%n_o))
         call allocator(c_bjc_i, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
!
         
         do c = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bjc = index_three(b, j, c, wf%n_v, wf%n_o)
                  jb = index_two(j, b, wf%n_o)
!
                  do i = 1, wf%n_o
!
                     bi = index_two(b, i, wf%n_v)
                     cj = index_two(c, j, wf%n_v)
                     jc = index_two(j, c, wf%n_o)
!
                     bicj = index_packed(bi, cj)
!
                     c_bjc_i(bjc, i) = c_bicj(bicj, 1)
!
                  enddo
                  do a = 1, a_length
!
                     ca = index_two(c, a, wf%n_v)
                     ba = index_two(b, a, wf%n_v)
!
                     L_a_bjc(a, bjc) = two*g_ba_jc(ba, jc) - g_ba_jc(ca, jb)
!
                  enddo
               enddo
            enddo
         enddo
         call deallocator(g_ba_jc, a_length*(wf%n_v), (wf%n_v)*(wf%n_o))
!
         call dgemm('N', 'N',                & 
                     a_length,               &
                     wf%n_o,                 &
                     (wf%n_o)*((wf%n_v)**2), &
                     one,                    &
                     L_a_bjc,                &
                     a_length,               &
                     c_bjc_i,                &
                     (wf%n_o)*((wf%n_v)**2), &
                     one,                    &
                     rho_a_i(a_first, 1),    &
                     wf%n_v)
!
         call deallocator(c_bjc_i, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
         call deallocator(L_a_bjc, a_length, ((wf%n_v)**2)*(wf%n_o))
!
      enddo
!
   end subroutine jacobian_ccsd_d1_ccsd
!
!
end submodule jacobian