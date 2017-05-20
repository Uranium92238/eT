submodule (cc3_class) omega
!
!!
!!    Omega submodule (CC3)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Contains the following subroutines:
!!
!!    construct_omega: constructs the projection vector by adding the CC3 
!!                     and CCSD contributions to omega1 and omega2. 
!!    omega_integrals: constructs the integrals needed for the CC3 contributions,
!!                     and saves them to disk.
!!    omega_e1:        adds the E1 term to omega1 (for a given i, j, and k).
!!    omega_f2:        adds the F2 term to omega2 (for a given i, j, and k).
!!    omega_g2:        adds the G2 term to omega2 (for a given i, j, and k).
!!
!!    The submodule was based on the MLCC3 routines implemented by Rolf H. Myhre
!!    and Henrik Koch in the Dalton quantum chemistry program. 
!!
!!    Note: the construction of omega has not yet been optimized, and is rather
!!    slow. Among other things, the factor of 1/6 (i >= j >= k) that can be saved
!!    in the construction of t_ijk^abc has not been implemented. Todo. 
!!
!
   implicit none 
!
contains
!
   module subroutine construct_omega_cc3(wf)
!!
!!    Construct Omega (CC3)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Directs the calculation of the projection vector (omega1, omega2)
!!    for the CC3 level of theory.
!!
      implicit none 
!
      class(cc3) :: wf 
!
      real(dp), dimension(:,:), allocatable :: t_abc ! Triples, t_ijk^abc, fixed ijk
!
      real(dp), dimension(:,:), allocatable :: omega_ai_bj ! Unpacked holder of CC3 contributions
                                                           ! to the omega2 vector
!
      integer(i15) :: i = 0, j = 0, k = 0 
      integer(i15) :: a = 0, b = 0, ai = 0, bj = 0, aibj = 0
!
!     Set the omega vector to zero 
!
      wf%omega1 = zero
      wf%omega2 = zero    
!
!     :: CC3 integrals :: 
!
!     Calculate and save to disk the integrals needed 
!     to form the triples amplitudes & CC3 omega contributions     
!
      call wf%omega_integrals
!
!     :: CC3 contributions to omega ::
!
!     We construct the triples amplitudes t_abc (= t_ijk^abc) 
!     for a given set of occupied indices, i, j, k, adding
!     the appropriate CC3 omega contributions. 
!     
      call allocator(t_abc, (wf%n_v)**3, 1)
      t_abc = zero
!
      call allocator(omega_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      omega_ai_bj = zero 
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
            do k = 1, wf%n_o 
!
!              Calculate the triples amplitudes:
!
!              Calculate W_abc = P_ijk^abc ( sum_d t_ij^ad g_bdck - sum_l t_il^ab g_ljck )
!              and divide by orbital energy difference, t_abc = - W_abc / e_abc
!
               t_abc = zero
               call wf%calc_triples(t_abc,i,j,k)
!
!              Add the CC3 omega terms, using the calculated triples amplitudes:
!
               call wf%omega_e1(t_abc,i,j,k)
!
               call wf%omega_f2(omega_ai_bj,t_abc,i,j,k)
               call wf%omega_g2(omega_ai_bj,t_abc,i,j,k)
!
            enddo
         enddo
      enddo
!
!     Pack in squared omega into the wavefunction's packed omega vector 
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
                  aibj = index_packed(ai, bj)
!
                  if (ai .ge. bj) then 
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega_ai_bj(ai,bj) + &
                                                               omega_ai_bj(bj,ai)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
      call deallocator(t_abc, (wf%n_v)**3, 1)
!
!     :: CCSD contributions to omega :: 
!
!     Construct singles contributions (CCSD)
!
      call wf%omega_a1
      call wf%omega_b1
      call wf%omega_c1
      call wf%omega_ccs_a1
!
!     Construct doubles contributions (CCSD)
!
      call wf%omega_a2
      call wf%omega_b2
      call wf%omega_c2
      call wf%omega_d2
      call wf%omega_e2 
!
   end subroutine construct_omega_cc3
!
!
   subroutine omega_integrals_cc3(wf)
!!
!!     Omega Integrals (CC3)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!     Calculates g_bdck and saves to disk in the order bcd, k (one record per d and k)
!!     Calculates g_ljck and saves to disk in the order lc, jk (one record per j and k)
!!     Calculates g_jbkc and saves to disk in the order bc, kj (one record per k and j)
!!     Calculates g_ilkc and saves to disk in the order cl, ik (one record per i and k)
!!     Calculates g_dbkc and saves to disk in the order bcd, k (one record per d and k)
!!
      implicit none 
!
      class(cc3) :: wf
!
      integer(i15) :: unit_bdck = 0, unit_ljck = 0
      integer(i15) :: unit_jbkc = 0, unit_ilkc = 0, unit_dbkc = 0
!
      integer(i15) :: rec_number = 0, ioerror = 0 
!
      real(dp), dimension(:,:), allocatable :: L_bd_J 
      real(dp), dimension(:,:), allocatable :: L_ck_J
      real(dp), dimension(:,:), allocatable :: L_lj_J
      real(dp), dimension(:,:), allocatable :: L_kc_J
!
      real(dp), dimension(:,:), allocatable :: g_bd_ck ! g_bdck 
      real(dp), dimension(:,:), allocatable :: g_bcd_k ! g_bdck or g_dbkc 
!
      real(dp), dimension(:,:), allocatable :: g_lj_ck ! g_ljck
      real(dp), dimension(:,:), allocatable :: g_lc_jk ! g_ljck 
!
      real(dp), dimension(:,:), allocatable :: g_jb_kc ! g_jbkc 
      real(dp), dimension(:,:), allocatable :: g_bc_kj ! g_jbkc 
!
      real(dp), dimension(:,:), allocatable :: g_il_kc ! g_ilkc 
      real(dp), dimension(:,:), allocatable :: g_cl_ik ! g_ilkc
!
      real(dp), dimension(:,:), allocatable :: g_bd_kc ! g_dbkc
!
      integer(i15) :: d = 0, c = 0, b = 0, k = 0, i = 0
      integer(i15) :: bcd = 0, ck = 0, j = 0, jk = 0, l = 0
      integer(i15) :: lc = 0, lj = 0, bc = 0, jb = 0, kc = 0
      integer(i15) :: kj = 0, cl = 0, il = 0, ik = 0
!
      logical :: reorder = .false.
!
!     :: g_bdck integrals ::
!
      call generate_unit_identifier(unit_bdck)
      open(unit=unit_bdck, file='bdck', action='write', status='replace', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)**2, iostat=ioerror)
!
      do d = 1, wf%n_v 
!
!        Read Cholesky vectors L_bd_J and L_ck_J 
!
         call allocator(L_bd_J, (wf%n_v), wf%n_J)
         L_bd_J = zero
!
         reorder = .false.
         call wf%get_cholesky_ab(L_bd_J, d, d, wf%n_v, reorder)
!
         call allocator(L_ck_J, (wf%n_v)*(wf%n_o), wf%n_J)
         L_ck_J = zero
!
         call wf%get_cholesky_ai(L_ck_J)
!
!        Form the integral g_bd_ck = g_bdck
!
         call allocator(g_bd_ck, wf%n_v, (wf%n_v)*(wf%n_o))
         g_bd_ck = zero
!
         call dgemm('N','T',            &
                     wf%n_v,            &
                     (wf%n_v)*(wf%n_o), &
                     wf%n_J,            &
                     one,               &
                     L_bd_J,            &
                     wf%n_v,            &
                     L_ck_J,            &
                     (wf%n_v)*(wf%n_o), &
                     zero,              &
                     g_bd_ck,           &
                     wf%n_v)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ck_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call deallocator(L_bd_J, (wf%n_v), wf%n_J)
!
!        Reorder the integral as g_bcd_k
!
         call allocator(g_bcd_k, (wf%n_v)**2, wf%n_o)
         g_bcd_k = zero
!
         do c = 1, wf%n_v
            do b = 1, wf%n_v
               do k = 1, wf%n_o 
!
                  ck  = index_two(c, k, wf%n_v)
                  bcd = index_two(b, c, wf%n_v) ! d is fixed 
!
                  g_bcd_k(bcd,k) = g_bd_ck(b, ck)
!
               enddo
            enddo
         enddo
!
!        Deallocate g_bd_ck 
!
         call deallocator(g_bd_ck, wf%n_v, (wf%n_v)*(wf%n_o))
!
!        Write the integral g_bcd_k to disk 
!
         do k = 1, wf%n_o
!
            rec_number = index_two(d, k, wf%n_v)
            write(unit_bdck, rec=rec_number, iostat=ioerror) (g_bcd_k(i,k), i = 1, (wf%n_v)**2)
!
         enddo
!
         call deallocator(g_bcd_k, (wf%n_v)**2, wf%n_o)
!
      enddo ! End of do's over d 
!
!     Close bdck file 
!
      close(unit_bdck)
!
!     :: g_dbkc integrals :: 
!
      call generate_unit_identifier(unit_dbkc)
      open(unit=unit_dbkc, file='dbkc', action='write', status='replace', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)**2, iostat=ioerror)
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call wf%get_cholesky_ia(L_kc_J)
!
      do d = 1, wf%n_v
!
!        Read Cholesky vector L_bd_J(bd,J) = L_db^J 
!
         call allocator(L_bd_J, (wf%n_v), wf%n_J)
         L_bd_J = zero
!
         reorder = .true. ! bd is actually db => reorder
         call wf%get_cholesky_ab(L_bd_J, d, d, wf%n_v, reorder)
!
!        Form the integral g_bd_kc = g_dbkc
!
         call allocator(g_bd_kc, wf%n_v, (wf%n_v)*(wf%n_o))
         g_bd_kc = zero
!
         call dgemm('N','T',            &
                     wf%n_v,            &
                     (wf%n_v)*(wf%n_o), &
                     wf%n_J,            &
                     one,               &
                     L_bd_J,            &
                     wf%n_v,            &
                     L_kc_J,            &
                     (wf%n_v)*(wf%n_o), &
                     zero,              &
                     g_bd_kc,           &
                     wf%n_v)
!
!        Deallocate Cholesky vector
!
         call deallocator(L_bd_J, (wf%n_v), wf%n_J)
!
!        Reorder the integral as g_bcd_k
!
         call allocator(g_bcd_k, (wf%n_v)**2, wf%n_o)
         g_bcd_k = zero
!
         do c = 1, wf%n_v
            do b = 1, wf%n_v
               do k = 1, wf%n_o 
!
                  kc  = index_two(k, c, wf%n_o)
                  bcd = index_two(b, c, wf%n_v) ! d is fixed 
!
                  g_bcd_k(bcd,k) = g_bd_kc(b, kc) ! g_dbkc
!
               enddo
            enddo
         enddo
!
!        Deallocate g_bd_ck 
!
         call deallocator(g_bd_kc, wf%n_v, (wf%n_v)*(wf%n_o))
!
!        Write the integral g_bcd_k to disk (= g_dbkc)
!
         do k = 1, wf%n_o
!
            rec_number = index_two(d, k, wf%n_v)
            write(unit_dbkc, rec=rec_number, iostat=ioerror) &
                        (g_bcd_k(i,k), i = 1, (wf%n_v)**2)
!
         enddo
!
         call deallocator(g_bcd_k, (wf%n_v)**2, wf%n_o)
!
      enddo ! End of do's over d 
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      close(unit_dbkc)
!
!     :: g_ljck integrals ::
!
      call generate_unit_identifier(unit_ljck)
      open(unit=unit_ljck, file='ljck', action='write', status='replace', &
           access='direct', form='unformatted', recl=dp*(wf%n_o)*(wf%n_v), iostat=ioerror)
!
      call allocator(L_lj_J, (wf%n_o)**2, wf%n_J)
      call allocator(L_ck_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
      call wf%get_cholesky_ij(L_lj_J)
      call wf%get_cholesky_ai(L_ck_J)
!
!     Form g_lj_ck 
!
      call allocator(g_lj_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
      call dgemm('N','T',            &
                  (wf%n_o)**2,       &
                  (wf%n_v)*(wf%n_o), &
                  wf%n_J,            &
                  one,               &
                  L_lj_J,            &
                  (wf%n_o)**2,       &
                  L_ck_J,            &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  g_lj_ck,           &
                  (wf%n_o)**2)
!
!     Deallocate Cholesky vectors 
!
      call deallocator(L_lj_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_ck_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder to g_lc_jk
!
      call allocator(g_lc_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      g_lc_jk = zero
!
      do c = 1, wf%n_v
         do l = 1, wf%n_o 
!
            lc = index_two(l, c, wf%n_o)
!
            do k = 1, wf%n_o
!
               ck = index_two(c, k, wf%n_v)
!
               do j = 1, wf%n_o
!
                  lj = index_two(l, j, wf%n_o)
                  jk = index_two(j, k, wf%n_o)
!
                  g_lc_jk(lc, jk) = g_lj_ck(lj, ck)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integral 
!
      call deallocator(g_lj_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
!     Write the integrals to the ljck file,
!     ordered as lc, jk  
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
!
            jk = index_two(j, k, wf%n_o)
            rec_number = jk
            write(unit_ljck, rec=rec_number, iostat=ioerror) (g_lc_jk(i,jk), i = 1, (wf%n_o)*(wf%n_v))
!
         enddo
      enddo
!
      call deallocator(g_lc_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
!     Close ljck file 
!
      close(unit_ljck)
!
!     :: g_jbkc integrals ::
!
      call generate_unit_identifier(unit_jbkc)
      open(unit=unit_jbkc, file='jbkc', action='write', status='replace', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)**2, iostat=ioerror)
!
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_kc_J)
!
      call allocator(g_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_jb_kc = zero
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
                  g_jb_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate Cholesky vector L_kc_J
!
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder g_jb_kc to g_bc_kj
!
      call allocator(g_bc_kj, (wf%n_v)**2, (wf%n_o)**2)
      g_bc_kj = zero
!
      do c = 1, wf%n_v
         do b = 1, wf%n_v
!
            bc = index_two(b, c, wf%n_v)
!
            do j = 1, wf%n_o
!
               jb = index_two(j, b, wf%n_o)
!
               do k = 1, wf%n_o
!
                  kc = index_two(k, c, wf%n_o)
                  kj = index_two(k, j, wf%n_o)
!
                  g_bc_kj(bc,kj) = g_jb_kc(jb,kc) ! g_jbkc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate unordered integral g_jb_kc
!
      call deallocator(g_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Write the integrals to the jbkc file,
!     ordered as bc, kj  
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
!
            kj = index_two(k, j, wf%n_o)
            rec_number = kj
            write(unit_jbkc, rec=rec_number, iostat=ioerror) (g_bc_kj(i,kj), i = 1, (wf%n_v)**2)
!
         enddo
      enddo
!
      call deallocator(g_bc_kj, (wf%n_v)**2, (wf%n_o)**2)
!
!     Close jbkc file 
!
      close(unit_jbkc)
!
!     :: g_ilkc integrals ::
!
      call generate_unit_identifier(unit_ilkc)
      open(unit=unit_ilkc, file='ilkc', action='write', status='replace', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
      call allocator(L_lj_J, (wf%n_o)**2, wf%n_J)
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ij(L_lj_J)
      call wf%get_cholesky_ia(L_kc_J)
!
      call allocator(g_il_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      g_il_kc = zero
!
      call dgemm('N','T',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_lj_J,            &
                  (wf%n_o)**2,       &
                  L_kc_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_il_kc,           &
                  (wf%n_o)**2)
!
      call deallocator(L_lj_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Reorder g_il_kc to g_cl_ik 
!
      call allocator(g_cl_ik, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
      g_cl_ik = zero 
!
      do k = 1, wf%n_o
         do i = 1, wf%n_o
!
            ik = index_two(i, k, wf%n_o)
!
            do l = 1, wf%n_o
!
               il = index_two(i, l, wf%n_o)
!
               do c = 1, wf%n_v
!
                  kc = index_two(k, c, wf%n_o)
                  cl = index_two(c, l, wf%n_v)
!
                  g_cl_ik(cl, ik) = g_il_kc(il, kc)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_il_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Write the g_ilkc integrals to the ilkc file,
!     ordered as cl, ik
!
      do i = 1, wf%n_o
         do k = 1, wf%n_o
!
            ik = index_two(i, k, wf%n_o)
            rec_number = ik
            write(unit_ilkc, rec=rec_number, iostat=ioerror) (g_cl_ik(j,ik), j = 1, (wf%n_v)*(wf%n_o))
!
         enddo
      enddo
!
!     Close the ilkc file 
! 
      close(unit_ilkc)
!
      call deallocator(g_cl_ik, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
   end subroutine omega_integrals_cc3
!
!
   subroutine calc_triples_cc3(wf,w_abc,i,j,k)
!!
!!    Calculate Triples (CC3)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculate W_abc = P_ijk^abc ( sum_d t_ij^ad g_bdck - sum_l t_il^ab g_ljck )
!!    and divide by orbital energy difference, t_abc = - W_abc / e_abc. On exit,
!!    w_abc = t_abc, the CC3 triples amplitude t_ijk^abc for a given set of occupied
!!    indices i, j, and k.
!!
      implicit none 
!
      class(cc3) :: wf 
!
      real(dp), dimension((wf%n_v)**3, 1) :: w_abc
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(:,:), allocatable :: t_a_b ! Holds t_pq^ab, where p,q is i,j, or k 
!
      real(dp), dimension(:,:), allocatable :: t_ab_l ! Holds t_pl^ab, where p is i, j, or k 
!
      real(dp), dimension(:,:), allocatable :: g_bc_d ! g_ckbd for the given k 
                                                      ! g_bc_d(bc,d) = g_ckbd (holds for j & i also)
!
      real(dp), dimension(:,:), allocatable :: g_l_c ! g_ljck for given j and k 
!
      real(dp), dimension(:,:), allocatable :: w ! Holds differently ordered contributions to w_abc
!
      real(dp) :: e_abc ! Orbital energy difference e_ijk^abc 
      real(dp) :: e_ijk ! e_i + e_j + e_k 
!
      integer(i15) :: a = 0, d = 0, ai = 0, dj = 0, aidj = 0, l = 0
      integer(i15) :: abc = 0, acb = 0, b = 0, c = 0, aidk = 0, dk = 0
      integer(i15) :: bj = 0, bjdk = 0, bca = 0, bac = 0, jk = 0, ab = 0
      integer(i15) :: aibl = 0, bl = 0, kj = 0, al = 0, albj = 0, ik = 0
      integer(i15) :: ki = 0, cba = 0, albk = 0, bk = 0, ij = 0, ji = 0
      integer(i15) :: ck = 0, djck = 0
!
      integer(i15) :: unit_bdck = -1, unit_ljck = -1
      integer(i15) :: ioerror = 0, iostat = 0, rec_number = 0
!
!
!     :: The g_ckbd terms, i.e., P_ijk^abc sum_d t_ij^ad g_bdck :: 
!
!
!     Term 1. sum_d t_ij^ad g_bdck = sum_d t_a_b(a,d) g_bc_d(bc,d)
!
      call allocator(t_a_b, wf%n_v, wf%n_v)
      t_a_b = zero
!
      do a = 1, wf%n_v
         do d = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
            dj = index_two(d, j, wf%n_v)
!
            aidj = index_packed(ai, dj)
!
            t_a_b(a, d) = wf%t2am(aidj, 1) ! t_a_b(a,d) = t_ij^ad
!
         enddo
      enddo
!
!     Read g_bdck, ordered as g_bc_d for k fixed 
!
      call generate_unit_identifier(unit_bdck)
      open(unit=unit_bdck, file='bdck', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)**2, iostat=ioerror)
!
      call allocator(g_bc_d, (wf%n_v)**2, wf%n_v)
      g_bc_d = zero
!
      do d = 1, wf%n_v
!
         rec_number = index_two(d, k, wf%n_v)
         read(unit_bdck, rec=rec_number, iostat=ioerror) (g_bc_d(l,d), l = 1, (wf%n_v)**2)
!
      enddo
!
      call dgemm('N','T',      &
                  wf%n_v,      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  one,         &
                  t_a_b,       &
                  wf%n_v,      &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  one,         &
                  w_abc,       &
                  wf%n_v)
!
!     Term 2. sum_d t_ij^db g_adck = sum_d g_bc_d(ac,d) t_a_d(d,b) -> ordered (ac,b)
!
!     Note: t_a_d(a,d) = t_ij^ad,
!        => t_a_d(d,b) = t_ij^db 
!
!           g_bc_d(bc,d) = g_bdck
!        => g_bc_d(ac,d) = g_adck 
!
      call allocator(w, (wf%n_v)**3, 1)
      w = zero
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_v,      &
                  one,         &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  t_a_b,       &
                  wf%n_v,      &
                  zero,        &
                  w,           &
                  (wf%n_v)**2) ! Note: wait to add to w_abc
!
!     Term 3. sum_d t_ik^ad g_bjcd
! 
      t_a_b = zero 
!
      do a = 1, wf%n_v
         do d = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
            dk = index_two(d, k, wf%n_v)
!
            aidk = index_packed(ai, dk)
!
            t_a_b(a, d) = wf%t2am(aidk, 1) ! t_a_b(a,d) = t_ik^ad
!
         enddo
      enddo
!
!     g_bc_d(bc,d) = g_bdck 
!  => g_bc_d(cb,d) = g_cdbk = g_bkcd -> switch k with j 
!
      g_bc_d = zero
!
      do d = 1, wf%n_v
!
         rec_number = index_two(d, j, wf%n_v)
         read(unit_bdck, rec=rec_number, iostat=ioerror) (g_bc_d(l,d), l = 1, (wf%n_v)**2)
!
      enddo
!
!     Now, sum_d t_ik^ad g_bjcd = sum_d t_a_b(a,d) g_bc_d(cb,d) -> ordered as a,cb 
!
      call dgemm('N','T',      &
                  wf%n_v,      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  one,         &
                  t_a_b,       &
                  wf%n_v,      &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  one,         & ! Accumulate from previous term
                  w,           &
                  (wf%n_v))
!
!     Add w(acb) to w_abc(abc)
!
      do b = 1, wf%n_v
         do c = 1, wf%n_v
            do a = 1, wf%n_v
!
               acb = index_three(a, c, b, wf%n_v, wf%n_v)
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
!
               w_abc(abc,1) = w_abc(abc,1) + w(acb,1)
!
            enddo
         enddo
      enddo
!
!     Term 4. sum_d t_ik^dc g_bjad = sum_d g_bc_d(ab,d) t_a_b(d,c) -> ordered as ab, c
!
!     Note: g_bc_d(cb,d) = g_bjcd
!        => g_bc_d(ab,d) = g_bjad
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_v,      &
                  one,         &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  t_a_b,       &
                  wf%n_v,      &
                  one,         &
                  w_abc,       &
                  (wf%n_v)**2)
!
!     Term 5. sum_d t_jk^bd g_cdai
!
      t_a_b = zero
!
      do b = 1, wf%n_v
         do d = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
            dk = index_two(d, k, wf%n_v)
!
            bjdk = index_packed(bj, dk)
!
            t_a_b(b, d) = wf%t2am(bjdk, 1) ! t_a_b(b,d) = t_jk^bd
!
         enddo
      enddo
!
!     g_bc_d(bc,d) = g_bdck 
!  => g_bc_d(ca,d) = g_cdak -> switch k with i 
!
      g_bc_d = zero
!
      do d = 1, wf%n_v
!
         rec_number = index_two(d, i, wf%n_v)
         read(unit_bdck, rec=rec_number, iostat=ioerror) (g_bc_d(l,d), l = 1, (wf%n_v)**2)
!
      enddo
!
!     sum_d t_jk^bd g_cdai = sum_d t_a_b(b,d) g_bc_d(ca,d) -> ordered as b,ca 
!
      w = zero
!
      call dgemm('N','T',      &
                  wf%n_v,      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  one,         &
                  t_a_b,       &
                  wf%n_v,      &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  zero,        &
                  w,           &
                  (wf%n_v))
!
!     Add w(bca) to w_abc(abc)
!
      do a = 1, wf%n_v
         do c = 1, wf%n_v
            do b = 1, wf%n_v
!
               bca = index_three(b, c, a, wf%n_v, wf%n_v)
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
!
               w_abc(abc,1) = w_abc(abc,1) + w(bca,1)
!
            enddo
         enddo
      enddo
!
!     Term 6. sum_d t_jk^dc g_bdai = sum_d g_bc_d(ba,d) t_a_b(d,c) -> ordered as ba,c 
!
!     Note: g_bc_d(ca,d) = g_cdai
!        => g_bc_d(ba,d) = g_bdai
!
      t_a_b = zero
!
      do c = 1, wf%n_v
         do d = 1, wf%n_v
!
            dj = index_two(d, j, wf%n_v)
            ck = index_two(c, k, wf%n_v)
!
            djck = index_packed(dj, ck)
!
            t_a_b(d, c) = wf%t2am(djck, 1) ! t_a_b(d,c) = t_jk^dc 
!
         enddo
      enddo
!
      w = zero
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_v,      &
                  one,         &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  t_a_b,       &
                  wf%n_v,      &
                  zero,        &
                  w,           &
                  (wf%n_v)**2)
!
!     Add w(bac) to w_abc(abc)
!
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               bac = index_three(b, a, c, wf%n_v, wf%n_v)
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
!
               w_abc(abc,1) = w_abc(abc,1) + w(bac,1)
!
            enddo
         enddo
      enddo
!
!     Close integral file
!
      close(unit_bdck)
!
!     Deallocations 
!
      call deallocator(g_bc_d, (wf%n_v)**2, wf%n_v)
      call deallocator(t_a_b, wf%n_v, wf%n_v)
!
!     :: The g_ljck terms, i.e., - P_ijk^abc sum_l t_il^ab g_ljck :: 
!
!     Read the integrals g_ljck from file 
!
      call generate_unit_identifier(unit_ljck)
      open(unit=unit_ljck, file='ljck', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_o)*(wf%n_v), iostat=ioerror)
!
!     Read into g_l_c, ordered as lc, jk, such that g_l_c(l,c) = g_cklj 
!
      call allocator(g_l_c, wf%n_o, wf%n_v)
      g_l_c = zero
!
      jk = index_two(j, k, wf%n_o)
      rec_number = jk
      read(unit_ljck, rec=rec_number, iostat=ioerror) g_l_c
!
!     Set the amplitudes
!
      call allocator(t_ab_l, (wf%n_v)**2, wf%n_o)
      t_ab_l = zero
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            ab = index_two(a, b, wf%n_v)
!
            do l = 1, wf%n_o
!
               ai = index_two(a, i, wf%n_v)
               bl = index_two(b, l, wf%n_v)
!
               aibl = index_packed(ai, bl)
!  
               t_ab_l(ab,l) = wf%t2am(aibl, 1) ! t_ab_l(ab,l) = t_il^ab 
!
            enddo
         enddo
      enddo
!
!     Term 1. - sum_l t_il^ab g_cklj = - sum_l t_ab_l(ab,l) g_l_c(l,c) -> ordered as ab, c
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_o,      &
                  -one,        &
                  t_ab_l,      &
                  (wf%n_v)**2, &
                  g_l_c,       &
                  wf%n_o,      &
                  one,         &
                  w_abc,       &
                  (wf%n_v)**2)
!
!     Term 2. sum_l t_il^ac g_bjlk = sum_l t_ab_l(ac,l) * ? 
!
!     Note: g_l_c(l,c) = g_cklj
!        => g_l_c(l,b) = g_bklj -> switch k and j 
!
      g_l_c = zero
!
      kj = index_two(k, j, wf%n_o)
      rec_number = kj
      read(unit_ljck, rec=rec_number, iostat=ioerror) g_l_c
!
!     Now g_l_c(l,b) = g_bjlk, giving 
!
!        sum_l t_il^ac g_bjlk = sum_l t_ab_l(ac,l) g_l_c(l,b) -> ordered as ac,b
!
      w = zero
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_o,      &
                  -one,        &
                  t_ab_l,      &
                  (wf%n_v)**2, &
                  g_l_c,       &
                  wf%n_o,      &
                  zero,        &
                  w,           &
                  (wf%n_v)**2)
!
!     Add w(acb) to w_abc(abc)
!
      do b = 1, wf%n_v
         do c = 1, wf%n_v
            do a = 1, wf%n_v
!
               acb = index_three(a, c, b, wf%n_v, wf%n_v)
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
!
               w_abc(abc,1) = w_abc(abc,1) + w(acb,1)
!
            enddo
         enddo
      enddo
!
!     Term 3. - sum_l t_lj^ab g_ckli
!
      t_ab_l = zero 
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            ab = index_two(a, b, wf%n_v)
!
            do l = 1, wf%n_o
!
               al = index_two(a, l, wf%n_v)
               bj = index_two(b, j, wf%n_v)
!
               albj = index_packed(al, bj)
!  
               t_ab_l(ab,l) = wf%t2am(albj, 1) ! t_ab_l(ab,l) = t_lj^ab 
!
            enddo
         enddo
      enddo
!
!     Note: g_l_c(l,b) = g_bjlk
!        => g_l_c(l,c) = g_cjlk, switch kj -> ik 
!
      g_l_c = zero
!
      ik = index_two(i, k, wf%n_o)
      rec_number = ik
      read(unit_ljck, rec=rec_number, iostat=ioerror) g_l_c ! g_l_c(l,c) = g_ckli
!
!     - sum_l t_lj^ab g_ckli = - sum_l t_ab_l(ab,l) g_l_c(l,c) -> ordered as ab, c
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_o,      &
                  -one,        &
                  t_ab_l,      &
                  (wf%n_v)**2, &
                  g_l_c,       &
                  wf%n_o,      &
                  one,         &
                  w_abc,       &
                  (wf%n_v)**2)
!
!     Term 4. - sum_l t_lj^cb g_ailk
!
!     Note: g_l_c(l,c) = g_ckli 
!        => g_l_c(l,a) = g_akli, switch i and k 
!
      g_l_c = zero
!
      ki = index_two(k, i, wf%n_o)
      rec_number = ki
      read(unit_ljck, rec=rec_number, iostat=ioerror) g_l_c ! g_l_c(l,a) = g_ailk 
!
!     - sum_l t_lj^cb g_ailk = - sum_l t_ab_l(cb,l) g_l_c(l,a) -> ordered as cb, a
!
      w = zero
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_o,      &
                  -one,        &
                  t_ab_l,      &
                  (wf%n_v)**2, &
                  g_l_c,       &
                  wf%n_o,      &
                  zero,        &
                  w,           &
                  (wf%n_v)**2)
!
!     Add w(cba) to w_abc(abc)
!
      do b = 1, wf%n_v
         do c = 1, wf%n_v
            do a = 1, wf%n_v
!
               cba = index_three(c, b, a, wf%n_v, wf%n_v)
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
!
               w_abc(abc,1) = w_abc(abc,1) + w(cba,1)
!
            enddo
         enddo
      enddo
!
!     Term 5. - sum_l t_lk^ac g_bjli 
!
      t_ab_l = zero 
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            ab = index_two(a, b, wf%n_v)
!
            do l = 1, wf%n_o
!
               al = index_two(a, l, wf%n_v)
               bk = index_two(b, k, wf%n_v)
!
               albk = index_packed(al, bk)
!  
               t_ab_l(ab,l) = wf%t2am(albk, 1) ! t_ab_l(ab,l) = t_lk^ab
!
            enddo
         enddo
      enddo
!
!     Note: g_l_c(l,c) = g_cilk 
!        => g_l_c(l,b) = g_bilk, switch ki with ij 
!
      g_l_c = zero
!
      ij = index_two(i, j, wf%n_o)
      rec_number = ij
      read(unit_ljck, rec=rec_number, iostat=ioerror) g_l_c ! g_l_c(l,b) = g_bjli 
!
!     - sum_l t_lk^ac g_bjli = - sum_l t_ab_l(ac,l) g_l_c(l,b) -> ordered as ac, b
!
      w = zero
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_v,      &
                  wf%n_o,      &
                  -one,        &
                  t_ab_l,      &
                  (wf%n_v)**2, &
                  g_l_c,       &
                  wf%n_o,      &
                  zero,        &
                  w,           &
                  (wf%n_v)**2)
!
!     Add w(acb) to w_abc(abc)
!
      do b = 1, wf%n_v
         do c = 1, wf%n_v
            do a = 1, wf%n_v
!
               acb = index_three(a, c, b, wf%n_v, wf%n_v)
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
!
               w_abc(abc,1) = w_abc(abc,1) + w(acb,1)
!
            enddo
         enddo
      enddo
!
!     Term 6. - sum_l t_lk^bc g_ailj 
!
!     Note: g_l_c(l,b) = g_bjli
!        => g_l_c(l,a) = g_ajli, switch i and j 
!
      g_l_c = zero
!
      ji = index_two(j, i, wf%n_o)
      rec_number = ji
      read(unit_ljck, rec=rec_number, iostat=ioerror) g_l_c ! g_l_c(l,a) = g_ailj 
!
!     - sum_l t_lk^bc g_ailj = - sum_l g_l_c(l,a) t_ab_l(bc,l) -> transposed, comes out a, bc 
!
      call dgemm('T','T',      &
                  wf%n_v,      &
                  (wf%n_v)**2, &
                  wf%n_o,      &
                  -one,        &
                  g_l_c,       &
                  wf%n_o,      &
                  t_ab_l,      &
                  (wf%n_v)**2, &
                  one,         &
                  w_abc,       &
                  wf%n_v)
!
!     Close integral file ljck 
!
      close(unit_ljck)
!
!     Divide by orbital energy difference (w_abc -> t_abc)
!
      e_ijk = wf%fock_diagonal(i, 1) + wf%fock_diagonal(j, 1) + wf%fock_diagonal(k, 1)
!
      if (i .eq. j .and. j .eq. k) then
!
         w_abc = zero 
!
      else
!
         do c = 1, wf%n_v
            do b = 1, wf%n_v
               do a = 1, wf%n_v 
!  
                  abc = index_three(a, b, c, wf%n_v, wf%n_v)
!  
                  if (a .eq. b .and. b .eq. c) then 
!  
                     w_abc(abc, 1) = zero
!  
                  else
!  
                     e_abc = -one/(wf%fock_diagonal(wf%n_o + a, 1) + &
                                wf%fock_diagonal(wf%n_o + b, 1) + &
                                wf%fock_diagonal(wf%n_o + c, 1) - e_ijk)
!  
                     w_abc(abc, 1) = e_abc*w_abc(abc, 1)
!  
                  endif
!  
               enddo
            enddo
         enddo
      endif
!
!     Deallocations 
!
      call deallocator(t_ab_l, (wf%n_v)**2, wf%n_o)
      call deallocator(g_l_c, wf%n_o, wf%n_v)
      call deallocator(w, (wf%n_v)**3, 1)
!
   end subroutine calc_triples_cc3
!
!
   subroutine omega_e1_cc3(wf,t_abc,i,j,k)
!!
!!    Omega E1 (CC3)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the E1 term,
!! 
!!       sum_bc (t_ijk^abc - t_ijk^cba) L_jbkc,
!!
!!    for a given set of i, j, and k, and adds the contribution to 
!!    the singles projection vector (omega1).
!!
      implicit none 
!
      class(cc3) :: wf 
!
      real(dp), dimension((wf%n_v)**3, 1) :: t_abc 
!
      real(dp), dimension(:,:), allocatable :: g_b_c ! g_jbkc 
      real(dp), dimension(:,:), allocatable :: L_bc  ! L_jbkc = 2 g_jbkc - g_jckb
!
      real(dp), dimension(:,:), allocatable :: u_a_bc ! t_ijk^abc - t_ijk^cba 
!
      integer(i15) :: unit_jbkc = 0
      integer(i15) :: rec_number = 0, ioerror = 0
!
      integer(i15), intent(in) :: i, j, k
!
      integer(i15) :: b = 0, c = 0, a = 0, abc = 0, cba = 0, bc = 0
!
!     Read the integrals g_jbkc, ordered as bc for a given j and k 
!
      call allocator(g_b_c, wf%n_v, wf%n_v)
      g_b_c = zero 
!
      call generate_unit_identifier(unit_jbkc)
      open(unit=unit_jbkc, file='jbkc', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)**2, iostat=ioerror)
!
      rec_number = index_two(k, j, wf%n_o)
      read(unit_jbkc, rec=rec_number, iostat=ioerror) g_b_c
!
!     Calculate L_bc(bc, 1) = L_jbkc = 2 g_b_c(b,c) - g_b_c(c,b) 
!
      call allocator(L_bc, (wf%n_v)**2, 1)
      L_bc = zero 
!
      do c = 1, wf%n_v
         do b = 1, wf%n_v
!
            bc = index_two(b, c, wf%n_v)
            L_bc(bc, 1) = two*g_b_c(b,c) - g_b_c(c,b)
!
         enddo
      enddo
!
!     Deallocate g_b_c 
!
      call deallocator(g_b_c, wf%n_v, wf%n_v)
!
!     Close jbkc file 
!
      close(unit_jbkc)
!
!     Form u_a_bc(a,bc) = t_abc(abc, 1) - t_abc(cba, 1)
!
      call allocator(u_a_bc, wf%n_v, (wf%n_v)**2)
      u_a_bc = zero 
!
      do c = 1, wf%n_v
         do b = 1, wf%n_v
!
            bc = index_two(b, c, wf%n_v)
!
            do a = 1, wf%n_v
!
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
               cba = index_three(c, b, a, wf%n_v, wf%n_v)
!
               u_a_bc(a, bc) = t_abc(abc, 1) - t_abc(cba, 1)
!
            enddo
         enddo
      enddo
!
!     Calculate the E1 contribution,
!
!        sum_bc (t_ijk^abc - t_ijk^cba) L_jbkc 
!                          = sum_bc u_a_bc(a,bc) L_bc(bc,1)
!
      call dgemm('N','N',         &
                  wf%n_v,         &
                  1,              &
                  (wf%n_v)**2,    &
                  one,            & 
                  u_a_bc,         &
                  wf%n_v,         &
                  L_bc,           &
                  (wf%n_v)**2,    &
                  one,            &
                  wf%omega1(1,i), &
                  wf%n_v)
!
!     Deallocations
!
      call deallocator(u_a_bc, wf%n_v, (wf%n_v)**2)
      call deallocator(L_bc, (wf%n_v)**2, 1)
!
   end subroutine omega_e1_cc3
!
!
   subroutine omega_f2_cc3(wf,omega_ai_bj,t_abc,i,j,k)
!!
!!    Omega F2 (CC3)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the F2 term,
!! 
!!       sum_c (t_ijk^abc - t_ijk^cba) F_kc,
!!
!!    for a given set of i, j, and k, and adds the contribution to 
!!    the doubles projection vector (omega2), element ai_bj.
!!
      implicit none 
!
      class(cc3) :: wf 
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: omega_ai_bj 
      real(dp), dimension((wf%n_v)**3, 1) :: t_abc
!
      integer(i15), intent(in) :: i, j, k 
!
      real(dp), dimension(:,:), allocatable :: u_ab_c ! t_ijk^abc - t_ijk^cba 
!
      real(dp), dimension(:,:), allocatable :: F_c ! F_kc 
!
      real(dp), dimension(:,:), allocatable :: omega_ab ! sum_c u_ab_c * F_c 
!
      integer(i15) :: a = 0, b = 0, c = 0, ab = 0, abc = 0, cba = 0
      integer(i15) :: ai = 0, bj = 0, ab = 0
!
!     Set u_ab_c(ab,c) = t_ijk^abc - t_ijk^cba 
!
      call allocator(u_ab_c, (wf%n_v)**2, wf%n_v)
      u_ab_c = zero 
!
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               ab = index_two(a, b, wf%n_v)
!
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
               cba = index_three(c, b, a, wf%n_v, wf%n_v)
!
               u_ab_c(ab,c) = t_abc(abc, 1) - t_abc(cba, 1)
!
            enddo
         enddo
      enddo
!
!     Set F_c(c, 1) = F_kc 
!
      call allocator(F_c, wf%n_v, 1)
      F_c = zero 
!
      do c = 1, wf%n_v
!
         F_c(c, 1) = wf%fock_ia(k, c)
!
      enddo
!
!     Allocate and compute the F2 contribution,
!
!        sum_c (t_ijk^abc - t_ijk^cba) F_kc = sum_c u_ab_c(ab,c) F_c(c, 1) 
!
      call allocator(omega_ab, (wf%n_v)**2, 1)
      omega_ab = zero 
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  1,           &
                  wf%n_v,      &
                  one,         &
                  u_ab_c,      &
                  (wf%n_v)**2, &
                  F_c,         &
                  wf%n_v,      &
                  zero,        &
                  omega_ab,    &
                  (wf%n_v)**2)
!
      call deallocator(F_c, wf%n_v, 1)
      call deallocator(u_ab_c, (wf%n_v)**2, wf%n_v)
!
!     Add the contribution into the unpacked 
!     projection vector 
!
      do a = 1, wf%n_v
!
         ai = index_two(a, i, wf%n_v)
!
         do b = 1, wf%n_v
!
            ab = index_two(a, b, wf%n_v)
            bj = index_two(b, j, wf%n_v)
!
            omega_ai_bj(ai, bj) = omega_ai_bj(ai, bj) + omega_ab(ab, 1)
!
         enddo
      enddo
!
      call deallocator(omega_ab, (wf%n_v)**2, 1)
!
   end subroutine omega_f2_cc3
!
!
   subroutine omega_g2_cc3(wf,omega_ai_bj,t_abc,i,j,k)
!!
!!    Omega G2 (CC3)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the G2 term,
!! 
!!       omega(al,bj) = - sum_c  (2 t_ijk^abc - t_ijk^cba - t_ijk^acb) g_ilkc
!!       omega(ai,dj) = + sum_bc (2 t_ijk^abc - t_ijk^cba - t_ijk^acb) g_dbkc
!!
!!    for a given set of i, j, and k, and adds the contribution to 
!!    the doubles projection vector (omega2).
!!
      implicit none 
!
      class(cc3) :: wf 
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: omega_ai_bj 
      real(dp), dimension((wf%n_v)**3, 1) :: t_abc 
!
      integer(i15), intent(in) :: i, j, k 
!
      integer(i15) :: unit_ilkc = 0, unit_dbkc = 0
      integer(i15) :: rec_number = 0, ioerror = 0
!
      real(dp), dimension(:,:), allocatable :: g_c_l ! g_ilkc 
!
      real(dp), dimension(:,:), allocatable :: g_bc_d ! g_dbkc 
!
      real(dp), dimension(:,:), allocatable :: v_abc ! (2 t_ijk^abc - t_ijk^cba - t_ijk^acb)
!
      real(dp), dimension(:,:), allocatable :: omega_ab_l ! Holds omega(al,bj)
      real(dp), dimension(:,:), allocatable :: omega_a_d  ! Holds omega(ai,dj)
!
      integer(i15) :: abc = 0, cba = 0, acb = 0, a = 0, b = 0, c = 0
      integer(i15) :: ab = 0, al = 0, bj = 0, l = 0, d = 0, ai = 0, dj = 0
!
!     Form v_abc(abc, 1) = 2 t_ijk^abc - t_ijk^cba - t_ijk^acb
!
      call allocator(v_abc, (wf%n_v)**3, 1)
      v_abc = zero 
!
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               abc = index_three(a, b, c, wf%n_v, wf%n_v)
               cba = index_three(c, b, a, wf%n_v, wf%n_v)
               acb = index_three(a, c, b, wf%n_v, wf%n_v)
!
               v_abc(abc, 1) = two*t_abc(abc, 1) - t_abc(cba, 1) &
                                                 - t_abc(acb, 1)
!
            enddo
         enddo
      enddo
!
!     Read in g_ilkc integrals, g_c_l(c,l) = g_ilkc 
!
      call generate_unit_identifier(unit_ilkc)
      open(unit=unit_ilkc, file='ilkc', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
      call allocator(g_c_l, wf%n_v, wf%n_o)
      g_c_l = zero
!
      rec_number = index_two(i, k, wf%n_o)
      read(unit_ilkc, rec=rec_number, iostat=ioerror) g_c_l
!
      close(unit_ilkc)
!
!     Construct the first contribution 
!
      call allocator(omega_ab_l, (wf%n_v)**2, wf%n_o)
      omega_ab_l = zero 
!
!     - sum_c  (2 t_ijk^abc - t_ijk^cba - t_ijk^acb) g_ilkc
!     = - sum_c v_abc g_c_l -> comes out ab, l
! 
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  wf%n_o,      &
                  wf%n_v,      &
                  -one,        &
                  v_abc,       &
                  (wf%n_v)**2, &
                  g_c_l,       &
                  wf%n_v,      &
                  zero,        &
                  omega_ab_l,  &
                  (wf%n_v)**2)
!
!     Add the contribution to the projection vector 
!
      do l = 1, wf%n_o
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               ab = index_two(a, b, wf%n_v)
!
               al = index_two(a, l, wf%n_v)
               bj = index_two(b, j, wf%n_v)
!
               omega_ai_bj(al,bj) = omega_ai_bj(al,bj) + omega_ab_l(ab,l)
!
            enddo
         enddo
      enddo
!
      call deallocator(g_c_l, wf%n_v, wf%n_o)
!
!     Read the g_dbkc integrals, ordered as bcd, k,
!     into g_bc_d, such that g_bc_d(bc, d) = g_dbkc 
!
      call generate_unit_identifier(unit_dbkc)
      open(unit=unit_dbkc, file='dbkc', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)**2, iostat=ioerror)
!
      call allocator(g_bc_d, (wf%n_v)**2, wf%n_v)
      g_bc_d = zero
!
      do d = 1, wf%n_v
!
         rec_number = index_two(d, k, wf%n_v)
         read(unit_dbkc, rec=rec_number, iostat=ioerror) (g_bc_d(l,d), l = 1, (wf%n_v)**2)
!
      enddo
!
      close(unit_dbkc)
!
!     Calculate the second contribution, 
!     sum_bc (2 t_ijk^abc - t_ijk^cba - t_ijk^acb) g_dbkc = sum_bc v_abc g_bc_d -> a,d 
!
      call allocator(omega_a_d, wf%n_v, wf%n_v)
      omega_a_d = zero 
!
      call dgemm('N','N',      & 
                  wf%n_v,      &
                  wf%n_v,      &
                  (wf%n_v)**2, &
                  one,         &
                  v_abc,       &
                  wf%n_v,      &
                  g_bc_d,      &
                  (wf%n_v)**2, &
                  zero,        &
                  omega_a_d,   &
                  wf%n_v)
!
!     Add the second contribution to the omega vector 
!
      do d = 1, wf%n_v
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
            dj = index_two(d, j, wf%n_v)
!
            omega_ai_bj(ai,dj) = omega_ai_bj(ai,dj) + omega_a_d(a,d)
!
         enddo
      enddo
!
      call deallocator(g_bc_d, (wf%n_v)**2, wf%n_v)
      call deallocator(omega_ab_l, (wf%n_v)**2, wf%n_o)
      call deallocator(omega_a_d, wf%n_v, wf%n_v)
      call deallocator(v_abc, (wf%n_v)**3, 1)
!
   end subroutine omega_g2_cc3
!
!
end submodule omega