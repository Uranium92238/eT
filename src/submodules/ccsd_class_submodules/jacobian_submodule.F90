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
   module subroutine calculate_orbital_differences_ccsd(wf,orbital_diff)
!!
!!       Calculate Orbital Differences (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
         integer(i15) :: a = 0, i = 0, b = 0, j = 0
         integer(i15) :: ai = 0, bj = 0
         integer(i15) :: aibj = 0
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               orbital_diff(ai, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1)
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     aibj = index_packed(ai, bj)
!
                     orbital_diff((wf%n_o)*(wf%n_v)+aibj, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1) &
                                                               + wf%fock_diagonal(b + wf%n_o, 1) - wf%fock_diagonal(j, 1)

!
                  enddo
               enddo
            enddo
         enddo
         call vec_print(orbital_diff,wf%n_parameters,1)
!
   end subroutine calculate_orbital_differences_ccsd
!
!
   module subroutine transform_trial_vectors_ccsd(wf, first_trial, last_trial)
!!
!!    Transformation Trial Vectors (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Each trial vector in first_trial to last_trial is read from file and
!!    transformed before the transformed vector is written to file.
!!
      implicit none
!
      class(ccsd) :: wf
!
      integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      real(dp), dimension(:,:), allocatable :: c_a_i
      real(dp), dimension(:,:), allocatable :: c_aibj
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
      integer(i15) :: trial = 0 
!
      write(unit_output,*)'In vector transformation'
      flush(unit_output)
!
!     Allocate c_a_i
!
      call allocator(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero 
!
      call allocator(c_aibj, wf%n_t2am, 1)
      c_aibj = zero 
!
!     Open trial vector and transformed vector files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*wf%n_parameters, iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*wf%n_parameters, iostat=ioerror)
!
!     For each trial vector: read, transform and write
!   
      write(unit_output,*)'Jacobian transformation'     
      flush(unit_output)    
!  
      do trial = first_trial, last_trial
!
         read(unit_trial_vecs, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
         call wf%jacobian_ccsd_transformation(c_a_i, c_aibj)
!
!        Write transformed vector to file
!
         write(unit_rho, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
      enddo
!
      close(unit_trial_vecs) 
      close(unit_rho)                                
!
!     Deallocate c_a_i
!
      call deallocator(c_a_i, wf%n_v, wf%n_o)
      call deallocator(c_aibj, wf%n_t2am, 1)
!
   end subroutine transform_trial_vectors_ccsd
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
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i       ! c_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: c_aibj      ! c_aibj     
!
      real(dp), dimension(:,:), allocatable :: rho_a_i   ! rho_ai   = (A c)_ai
      real(dp), dimension(:,:), allocatable :: rho_ai_bj ! rho_ai_bj = (A c)_aibj
!
      real(dp), dimension(:,:), allocatable :: c_ai_bj       ! Unpacked c_aibj
      real(dp), dimension(:,:), allocatable :: rho_ai_bj_sym ! Symmetrized rho_ai_bj, temporary
      real(dp), dimension(:,:), allocatable :: rho_ab_ij     ! rho_ai_bj, reordered
      real(dp), dimension(:,:), allocatable :: c_ab_ij       ! c_ai_bj, reordered
!
      integer(i15) :: a = 0, ab = 0, ai = 0, b = 0, bj = 0, i = 0, ij = 0, j = 0, aibj = 0
!
!     Allocate and zero the transformed vector (singles part)
!
      call allocator(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!     :: CCS contributions to the singles c vector ::  
!
      call wf%jacobian_ccs_a1(c_a_i, rho_a_i)
      call wf%jacobian_ccs_b1(c_a_i, rho_a_i)
!
!     :: CCSD contributions to the transformed singles vector :: 
!
      call wf%jacobian_ccsd_a1(c_a_i, rho_a_i)
!
!     Allocate the incoming unpacked doubles vector 
!
      call allocator(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call wf%jacobian_ccsd_b1(c_ai_bj, rho_a_i) 
      call wf%jacobian_ccsd_c1(c_ai_bj, rho_a_i)
      call wf%jacobian_ccsd_d1(c_ai_bj, rho_a_i)
!
!     :: CCSD contributions to the transformed doubles vector ::  
!
!     Allocate unpacked transformed vector 
!
      call allocator(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      rho_ai_bj = zero 
!
!     Contributions from singles vector c 
!
      call wf%jacobian_ccsd_a2(rho_ai_bj, c_a_i)
      call wf%jacobian_ccsd_b2(rho_ai_bj, c_a_i)
      call wf%jacobian_ccsd_c2(rho_ai_bj, c_a_i)
      call wf%jacobian_ccsd_d2(rho_ai_bj, c_a_i)
!
!     Done with singles vector c; overwrite it with 
!     transformed vector for exit
!
      call dcopy((wf%n_o)*(wf%n_v), rho_a_i, 1, c_a_i, 1)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_f2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_g2(rho_ai_bj, c_ai_bj) 
      call wf%jacobian_ccsd_h2(rho_ai_bj, c_ai_bj)
      call wf%jacobian_ccsd_i2(rho_ai_bj, c_ai_bj)
!
!     Last two terms are already symmetric (j2 and k2). Perform the symmetrization 
!     rho_ai_bj = P_ij^ab rho_ai_bj now, for convenience 
!
!     Allocate temporary symmetric transformed vector 
!
      call allocator(rho_ai_bj_sym, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
      call deallocator(rho_ai_bj_sym, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     In preparation for last two terms, reorder 
!     rho_ai_bj to rho_ab_ij, and c_ai_bj to c_ab_ij
!
      call allocator(rho_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
      call allocator(c_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
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
      call deallocator(c_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%jacobian_ccsd_j2(rho_ab_ij, c_ab_ij)
      call wf%jacobian_ccsd_k2(rho_ab_ij, c_ab_ij)
!
!     Done with reordered doubles c; deallocate 
!
      call deallocator(c_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     Order rho_ab_ij back into rho_ai_bj & divide by 
!     the biorthonormal factor 1 + delta_ai,bj
!
      call allocator(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call deallocator(rho_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     Overwrite the incoming doubles c vector & pack in
!
      c_aibj = zero
      call packin(c_aibj, rho_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Remaining deallocations 
!
      call deallocator(rho_a_i, wf%n_v, wf%n_o)
      call deallocator(rho_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      real(dp), dimension(:,:), allocatable :: X_k_i   ! An intermediate 
      real(dp), dimension(:,:), allocatable :: L_k_lcd ! L_lckd = L_lc_dk 
!
      real(dp), dimension(:,:), allocatable :: t_lcd_i ! t_li^cd 
!
      real(dp), dimension(:,:), allocatable :: t_a_lkd ! t_lk^ad 
      real(dp), dimension(:,:), allocatable :: L_lkd_c ! L_lckd
!
      real(dp), dimension(:,:), allocatable :: X_a_c ! An intermediate 
!
      integer(i15) :: a = 0, ai = 0, al = 0, c = 0, ci = 0, cial = 0, clai = 0
      integer(i15) :: d = 0, dk = 0, i = 0, k = 0, kc = 0, kd = 0, l = 0, lc = 0
      integer(i15) :: ld = 0, cl = 0, lcd = 0, cldi = 0, di = 0, aldk = 0, lkd = 0
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
!     L_lc_dk = L_lckd = 2*g_lckd - g_ldkc = 2*g_lc_kd(lc,kd)-g_lc_kd(ld,kc)
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
      call wf%initialize_amplitudes        ! Allocate t amplitudes, then set them to zero 
      call wf%read_double_amplitudes       ! Read the converged amplitudes from disk 
!
      call allocator(u_ai_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
!     Deallocations (keep L_lc_dk = L_lckd)
!
      call deallocator(u_ai_lc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call deallocator(X_lc, (wf%n_v)*(wf%n_o), 1)
!
!
!     :: Term 2. - sum_ckdl L_lckd t_li^cd c_ak ::
!
!
!     Reorder to L_k_lcd = L_lckd = L_lc_dk 
!
      call allocator(L_k_lcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
      L_k_lcd = zero
!
      do d = 1, wf%n_v
         do c = 1, wf%n_v
            do l = 1, wf%n_o
!
               lc = index_two(l, c, wf%n_o)
!
               lcd = index_three(l, c, d, wf%n_o, wf%n_v)
!
               do k = 1, wf%n_o 
!
                  dk = index_two(d, k, wf%n_v)
!
                  L_k_lcd(k, lcd) = L_lc_dk(lc, dk) ! L_lckd 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(L_lc_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder amplitudes to t_lcd_i = t_li^cd 
!
      call allocator(t_lcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      t_lcd_i = zero
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
                  t_lcd_i(lcd, i) = wf%t2am(cldi, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Calculate X_k_i = sum_cdl L_k_lcd t_lcd_i
!
      call allocator(X_k_i, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  L_k_lcd,              &
                  wf%n_o,               &
                  t_lcd_i,              &
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_k_i,                &
                  wf%n_o)
!
      call deallocator(t_lcd_i, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Calculate rho_a_i =+ - sum_k c_a_i(a,k) X_k_i(k,i)
!
      call dgemm('N','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  c_a_i,   &
                  wf%n_v,  &
                  X_k_i,   &
                  wf%n_o,  &
                  one,     &
                  rho_a_i, &
                  wf%n_v)
!
!     Deallocations (keep L_k_lcd = L_lckd)
! 
      call deallocator(X_k_i, wf%n_o, wf%n_o)
!
!
!     :: Term 3: - sum_ckdl L_lckd t_lk^ad c_ci ::
!
!     Reorder to L_lkd_c = L_lckd = L_k_lcd
!
      call allocator(L_lkd_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      L_lkd_c = zero
!
      do c = 1, wf%n_v
         do d = 1, wf%n_v
            do k = 1, wf%n_o
               do l = 1, wf%n_o
!
                  lkd = index_three(l, k, d, wf%n_o, wf%n_o)
                  lcd = index_three(l, c, d, wf%n_o, wf%n_v)
!
                  L_lkd_c(lkd, c) = L_k_lcd(k, lcd) ! L_lckd 
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(L_k_lcd, wf%n_o, (wf%n_o)*(wf%n_v)**2)
!
!     Reorder amplitudes to t_a_lkd = t_lk^ad 
!
      call allocator(t_a_lkd, wf%n_v, (wf%n_v)*(wf%n_o)**2)
      t_a_lkd = zero
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
!
            dk = index_two(d, k, wf%n_v)
!
            do l = 1, wf%n_o
!
               lkd = index_three(l, k, d, wf%n_o, wf%n_o)
!
               do a = 1, wf%n_v
!
                  al = index_two(a, l, wf%n_v)
!
                  aldk = index_packed(al, dk)
!
                  t_a_lkd(a, lkd) = wf%t2am(aldk, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Done with doubles amplitudes: deallocate 
!
      call wf%destruct_amplitudes
!
!     Calculate X_a_c = sum_kdl t_a_lkd L_lkd_c 
!
      call allocator(X_a_c, wf%n_v, wf%n_v)
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
!     Calculate rho_a_i =+ - sum_c X_a_c(a,c) c_a_i(c,i)
!
      call dgemm('N','N',  &
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
!     Deallocations 
!
      call deallocator(L_lkd_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      call deallocator(X_a_c, wf%n_v, wf%n_v)
      call deallocator(t_a_lkd, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
   end subroutine jacobian_ccsd_a1_ccsd
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, c_ai_bj, rho_a_i)
!!
!!       Jacobian CCSD B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017 
!!
!!       rho_ai^B1 = sum_bj F_jb (2*c_aibj  -  c_aj_bi) = sum_bj F_jb v_ai_bj
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))   :: c_ai_bj 
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
!        Construct v_ai_bj = 2*c_aibj - c_ajbi
!        Reorder F_j_b to F_bj 
!
         call allocator(v_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         v_ai_bj = zero
!
         call allocator(F_bj, (wf%n_o)*(wf%n_v), 1)
         F_bj = zero
!
         call allocator(rho_ai, (wf%n_o)*(wf%n_v), 1)
!
         do j = 1, wf%n_o
            do b = 1, wf%n_v
!
               bj = index_two(b, j, wf%n_v)
!
               F_bj(bj,1) = wf%fock_ia(j, b)
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
                     v_ai_bj(ai, bj) = two*c_ai_bj(ai, bj) - c_ai_bj(aj, bi)        
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
               rho_a_i(a,i) = rho_a_i(a,i) + rho_ai(ai, 1)
!
            enddo
         enddo

         call deallocator(rho_ai, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, c_ai_bj, rho_a_i)
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
         real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))   :: c_ai_bj 
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
         L_ji_J = zero
         call wf%get_cholesky_ij(L_ji_J)
!
         call allocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_kb_J = zero
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
         L_jkb_i = zero
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
         c_a_jkb = zero
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
                     c_a_jkb(a, jkb) = c_ai_bj(aj,bk)
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
                     one,                    &
                     rho_a_i,                &
                     wf%n_v)
!
         call deallocator(c_a_jkb, wf%n_v, (wf%n_v)*((wf%n_o)**2))
         call deallocator(L_jkb_i, (wf%n_v)*((wf%n_o)**2), wf%n_o)
         
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
    module subroutine jacobian_ccsd_d1_ccsd(wf, c_bi_cj, rho_a_i)
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
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))   :: c_bi_cj 
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
         L_ba_J = zero
         call wf%get_cholesky_ab(L_ba_J, a_first, a_last, a_length*(wf%n_v), .true.)
!
         call allocator(L_jc_J, (wf%n_v)*(wf%n_o), wf%n_J)
         L_jc_J = zero
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
         L_a_bjc = zero
         c_bjc_i = zero
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
                     c_bjc_i(bjc, i) = c_bi_cj(bi, cj)
!
                  enddo
!
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
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
!
      real(dp), dimension(:,:), allocatable :: L_ai_J ! L_ai^J
      real(dp), dimension(:,:), allocatable :: L_kj_J ! L_kj^J 
      real(dp), dimension(:,:), allocatable :: L_cb_J ! L_bc^J, batching over b
!
      real(dp), dimension(:,:), allocatable :: g_ai_kj ! g_aikj 
      real(dp), dimension(:,:), allocatable :: g_k_aij ! g_aikj reordered 
      real(dp), dimension(:,:), allocatable :: g_ai_cb ! g_aibc, batching over b
      real(dp), dimension(:,:), allocatable :: g_aib_c ! g_aibc, reordered, batching over b
!
      real(dp), dimension(:,:), allocatable :: rho_b_aij ! rho_ai_bj, term 1 (see below)
      real(dp), dimension(:,:), allocatable :: rho_aib_j ! rho_ai_bj, term 2 (batching over b)
!
      integer(i15) :: a = 0, ai = 0, aij = 0, b = 0, bj = 0, i = 0, j = 0
      integer(i15) :: k = 0, kj = 0, c = 0, cb = 0, aib = 0
!
      logical :: reorder 
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, batch_dimension = 0
      integer(i15) :: n_batch = 0, b_begin = 0, b_end = 0, b_batch = 0, batch_length = 0
!
      integer(i15) :: aib_offset = 0
!
!
!     :: Term 1. - sum_k g_aikj c_bk ::
!
!     Calculate g_ai_kj 
!
      call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call allocator(L_kj_J, (wf%n_o)**2, wf%n_J)
!
      call wf%get_cholesky_ai(L_ai_J)
      call wf%get_cholesky_ij(L_kj_J)
!
      call allocator(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)**2,       &
                  wf%n_J,            &
                  one,               &
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_kj_J,            &
                  (wf%n_o)**2,       &
                  zero,              &
                  g_ai_kj,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate Cholesky's (keep L_ai_J for later)
!
      call deallocator(L_kj_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Reorder to g_k_aij = g_aikj = g_ai_kj
!
      call allocator(g_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
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
      call deallocator(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate rho_b_aij =+ - sum_k c_bk g_aikj = - sum_k c_a_i(b, k) g_k_aij(k, aij)
!
      call allocator(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N','N',               &
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
      call deallocator(g_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
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
      call deallocator(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     
!     :: Term 2. rho_ai_bj =+ sum_c g_aibc c_cj ::
!
!     We do the matrix multiplication as g_aib_c c_cj,
!     batching over the b index.
!
      call allocator(rho_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o) ! rho_ai_bj formed in batch 
      rho_aib_j = zero
!
!     We hold L_bc^J and g_aibc in two orderings. The construction of L_bc^J
!     requires 2*n_v*n_o*n_J + n_J*n_v**2 extra memory.
!
      required = max(2*(wf%n_v)*(wf%n_o)*(wf%n_J) + &
                     2*(wf%n_J)*(wf%n_v)**2,        &  ! Constr of L_bc^J 
                     (wf%n_J)*(wf%n_v)**2 +         &
                     (wf%n_o)*(wf%n_v)**3)             ! Holding L_bc^J and g_aibc
!
      required = 4*required ! Words
!
      available = get_available()
! 
      batch_dimension  = wf%n_v ! Batch over the virtual index b
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)  
!
      do b_batch = 1, n_batch 
!
!        Get batching limits 
!
         call batch_limits(b_begin, b_end, b_batch, max_batch_length, batch_dimension)
         batch_length = b_end - b_begin + 1 
!
!        Get L_cb_J = L_bc^J for the current batch of b's
!
         call allocator(L_cb_J, (wf%n_v)*batch_length, wf%n_J)
!
         reorder = .true.
         call wf%get_cholesky_ab(L_cb_J, b_begin, b_end, (wf%n_v)*batch_length, reorder)
!
!        Calculate g_ai_cb = g_aibc
!
         call allocator(g_ai_cb, (wf%n_v)*(wf%n_o), (wf%n_v)*batch_length)
!
         call dgemm('N','T',                &
                     (wf%n_v)*(wf%n_o),     &
                     (wf%n_v)*batch_length, &
                     wf%n_J,                &
                     one,                   &
                     L_ai_J,                &
                     (wf%n_v)*(wf%n_o),     &
                     L_cb_J,                &
                     (wf%n_v)*batch_length, &
                     zero,                  &
                     g_ai_cb,               &
                     (wf%n_v)*(wf%n_o))
!
!        Reorder to g_aib_c 
!
         call allocator(g_aib_c, (wf%n_v)*(wf%n_o)*batch_length, wf%n_v)
         g_aib_c = zero 
!
         do c = 1, wf%n_v
            do b = 1, batch_length
!
               cb = index_two(c, b, wf%n_v)
!
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     ai = index_two(a, i, wf%n_v)
!
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     g_aib_c(aib, c) = g_ai_cb(ai, cb)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ai_cb, (wf%n_v)*(wf%n_o), (wf%n_v)*batch_length)
!
!        Calculate the contribution to rho_aib_j = sum_c g_aib_c c_cj
!
         aib_offset = index_three(1, 1, b_begin, wf%n_v, wf%n_o)
!
         call dgemm('N','N',                         &
                     (wf%n_v)*(wf%n_o)*batch_length, &
                     wf%n_o,                         &
                     wf%n_v,                         &
                     one,                            &
                     g_aib_c,                        &
                     (wf%n_v)*(wf%n_o)*batch_length, &
                     c_a_i,                          &
                     wf%n_v,                         &
                     one,                            &
                     rho_aib_j(aib_offset,1),        &
                     (wf%n_o)*(wf%n_v)**2)
!
         call deallocator(g_aib_c, (wf%n_v)*(wf%n_o)*batch_length, wf%n_v)
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
      call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
      call deallocator(rho_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
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
!!    Bug fix 28 May: changed the sign of the entire term, which was
!!    in disagreement with equations.
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
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
!
!     :: Term 1. sum_kc F_kc t_ij^ac c_bk ::
!
!     Read the amplitudes from disk 
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
!     Order the amplitudes as t_c_aij = t_ij^ac 
!
      call allocator(t_c_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
      call allocator(X_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N','N',               &
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
      call deallocator(t_c_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     Form rho_b_aij = sum_k c_a_i(b,k) X_k_aij(k,aij)
!
      call allocator(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N','N',               &
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
      call deallocator(X_k_aij, wf%n_o, (wf%n_v)*(wf%n_o)**2)
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
      call deallocator(rho_b_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!
!     :: Term 2. sum_kc F_kc t_ik^ab c_cj :: 
!
!     Form X_k_j = sum_c F_kc c_cj = sum_c fock_ia(k,c) c_a_i(c,j)
!
      call allocator(X_k_j, wf%n_o, wf%n_o)
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
      call allocator(t_aib_k, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
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
!     Form rho_aib_j = sum_k t_aib_k X_k_j 
!
      call allocator(rho_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
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
                  zero,                 &
                  rho_aib_j,            &
                  (wf%n_o)*(wf%n_v)**2)
!
      call deallocator(X_k_j, wf%n_o, wf%n_o)
      call deallocator(t_aib_k, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Add rho_aib_j to rho_ai_bj 
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
      call deallocator(rho_aib_j, (wf%n_o)*(wf%n_v)**2, wf%n_o)
!
!     Deallocate amplitudes    
!
      call wf%destruct_amplitudes
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
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
!
      real(dp), dimension(:,:), allocatable :: L_lj_J ! L_lj^J 
      real(dp), dimension(:,:), allocatable :: L_kc_J ! L_kc^J 
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
!   
!     :: Term 1. sum_kcl g_ljkc t_ki^ac c_bl :: 
!
!     Form g_lj_kc = g_ljkc
!
      call allocator(L_lj_J, (wf%n_o)**2, wf%n_J)
      call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ij(L_lj_J)
      call wf%get_cholesky_ia(L_kc_J)
!
      call allocator(g_lj_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
                  g_lj_kc,           &
                  (wf%n_o)**2)
!
      call deallocator(L_lj_J, (wf%n_o)**2, wf%n_J)
      call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Read the amplitudes from disk
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
!     Order as t_kc_ai = t_ki^ac 
!
      call allocator(t_kc_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call allocator(X_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
      call deallocator(t_kc_ai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
!     Calculate rho_b_jai = sum_l c_bl X_lj_ai 
!     (Interpret the X array as an X_l_jai object in the matrix multiplication)
!
      call allocator(rho_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_a_i,                &
                  wf%n_v,               &
                  X_lj_ai,              &
                  wf%n_o,               & ! "X_l_jai"
                  zero,                 &
                  rho_b_jai,            &
                  wf%n_v)
!
      call deallocator(X_lj_ai, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
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
      call deallocator(rho_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!     :: Term 2. sum_kcl g_ljkc t_li^bc c_ak ::
!
!     Reorder to g_kj_lc = g_lj_kc = g_ljkc     
!
      call allocator(g_kj_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
      call deallocator(g_lj_kc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Reorder amplitudes as t_lc_bi = t_li^bc 
!
      call allocator(t_lc_bi, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
      t_lc_bi = zero 
!
      do i = 1, wf%n_o
         do b = 1, wf%n_v
!
            bi = index_two(b, i, wf%n_v)
!
            do c = 1, wf%n_v
!
               ci = index_two(c, i, wf%n_v)
!
               do l = 1, wf%n_o
!
                  bl = index_two(b, l, wf%n_v)
                  lc = index_two(l, c, wf%n_o)
!
                  blci = index_packed(bl, ci)
!
                  t_lc_bi(lc, bi) = wf%t2am(blci, 1) ! t_li^bc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Form the intermediate X_kj_bi = sum_lc g_kj_lc t_lc_bi
!
      call allocator(X_kj_bi, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_kj_lc,           &
                  (wf%n_o)**2,       &
                  t_lc_bi,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_kj_bi,           &
                  (wf%n_o)**2)
!
      call deallocator(t_lc_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate rho_a_jbi = sum_k c_ak X_kj_bi 
!
      call allocator(rho_a_jbi, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
      call dgemm('N','N',               &
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
      call deallocator(X_kj_bi, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
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
      call deallocator(rho_a_jbi, wf%n_v, (wf%n_v)*(wf%n_o)**2)
!
!
!     :: Term 3. sum_kcl g_ljkc t_lk^ba c_ci :: 
!
!
!     Form the intermediate X_kjl_i = sum_c g_ljkc c_ci = sum_c g_kj_lc c_c_i
!
!     Note: interpret g_kj_lc as g_kjl_c in matrix multiplication.
!
      call allocator(X_kjl_i, (wf%n_o)**3, wf%n_o)
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
      call allocator(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
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
      call deallocator(X_kjl_i, (wf%n_o)**3, wf%n_o)
!
!     Order amplitudes as t_kl_ba = t_lk^ba 
!
      call allocator(t_ba_kl, (wf%n_v)**2, (wf%n_o)**2)
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
      call allocator(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
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
      call deallocator(t_ba_kl, (wf%n_v)**2, (wf%n_o)**2)
      call deallocator(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
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
      call deallocator(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!
!     :: Term 4. - sum_kcl L_ljkc t_il^ab c_ck ::
!
!     Form L_lj_ck = L_ljkc = 2 * g_ljck - g_lkcj = 2 * g_kj_lc(kj,lc) - g_kj_lc(jk,lc)
!
      call allocator(L_lj_ck, (wf%n_o)**2, (wf%n_v)*(wf%n_o))
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do j = 1, wf%n_o
!
               kj = index_two(k, j, wf%n_o)
               jk = index_two(j, k, wf%n_o)
!
               do l = 1, wf%n_o
!
                  lj = index_two(l, j, wf%n_o)
                  lc = index_two(l, c, wf%n_o)
!
                  L_lj_ck(lj, ck) = two*g_kj_lc(kj,lc) - g_kj_lc(jk,lc) ! L_ljkc
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_kj_lc, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
!     Calculate the intermediate X_lj = sum_ck L_lj_ck c_ck 
!
      call allocator(X_lj, (wf%n_o)**2, 1)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_lj_ck,           &
                  (wf%n_o)**2,       &
                  c_a_i,             & ! Interpret as "c_ai"
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_lj,              &
                  (wf%n_o)**2)
!
!     Order the amplitudes as t_aib_l = t_il^ab 
!
      call allocator(t_aib_l, (wf%n_o)*(wf%n_v)**2, wf%n_o)
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
      call deallocator(t_aib_l, (wf%n_o)*(wf%n_v)**2, wf%n_o)
      call deallocator(X_lj, (wf%n_o)**2, 1)
!
!     :: Term 5. - sum_kcl L_ljkc t_ik^ac c_bl ::
!
!     Square up the amplitudes, t_ck_ai(ck, ai) = t_ki^ca = t_ik^ac 
!
      call allocator(t_ck_ai, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_ck_ai = zero
!
      call squareup(wf%t2am, t_ck_ai, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate Y_lj_ai = sum_kc L_ljkc t_ik^ac = sum_kc L_lj_ck t_ck_ai 
!
      call allocator(Y_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
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
!     Calculate rho_b_jai =+ - sum_l c_bl Y_lj_ai
!
!     Note: interpret Y_lj_ai as Y_l_jai in the matrix multiplication
!
      call allocator(rho_b_jai, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
      call deallocator(Y_lj_ai, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      call deallocator(L_lj_ck, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
      call deallocator(rho_b_jai, (wf%n_v), (wf%n_v)*(wf%n_o)**2)
!
      call wf%destruct_amplitudes 
!
   end subroutine jacobian_ccsd_c2_ccsd
!
!
   module subroutine jacobian_ccsd_d2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD D2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^D2 = - ( sum_kcd g_kcbd (t_ij^cd c_ak + t_kj^ad c_ci + t_ik^ca c_dj)
!!                       - sum_kcd L_kcbd (t_ik^ac c_dj + t_ij^ad c_ck))
!!
!!    Note: the code is structured so that we batch over the index b,
!!          where the integrals are made as g_kc_db = g_kcbd and held
!!          in some ordering or other throughout a given batch (i.e.,
!!          all five terms are constructed gradually in the batches).
!!
!!    Bug fix 28 May: changed the sign of the entire term, which was
!!    in disagreement with equations.
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
!
      real(dp), dimension(:,:), allocatable :: L_kc_J ! L_kc^J
      real(dp), dimension(:,:), allocatable :: L_db_J ! L_bd^J
!
      real(dp), dimension(:,:), allocatable :: g_kc_db ! g_kcbd 
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
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, batch_dimension = 0
      integer(i15) :: n_batch = 0, b_begin = 0, b_end = 0, b_batch = 0, batch_length = 0
!
!     Indices 
!
      integer(i15) :: b = 0, c = 0, cd = 0, ci = 0, dj = 0, cidj = 0, d = 0, db = 0
      integer(i15) :: k = 0, j = 0, kb = 0, kc = 0, i = 0, ij = 0, ijb = 0, kb = 0
      integer(i15) :: a = 0, ai = 0, bj = 0, ib = 0, dkb = 0, dk = 0, akdj = 0, ak = 0
      integer(i15) :: aj = 0, ck = 0, ckb = 0, ciak = 0, aib = 0, dkb = 0, aidj = 0, aij = 0
!
!     Read amplitudes from disk
! 
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
!     Determine batch size, etc.
!     (Redo estimate once loop is done)
!
      required = max(2*(wf%n_v)*(wf%n_o)*(wf%n_J) + &
                     2*(wf%n_J)*(wf%n_v)**2,        &  ! Constr of L_bc^J 
                     (wf%n_J)*(wf%n_v)**2 +         &
                     (wf%n_o)*(wf%n_v)**3)             ! Holding L_bc^J and g_aibc
!
      required = 4*required ! Words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index b
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)  
!
      do b_batch = 1, n_batch 
!
!        Get batching limits 
!
         call batch_limits(b_begin, b_end, b_batch, max_batch_length, batch_dimension)
         batch_length = b_end - b_begin + 1 
!
!
!        :: Term 1. sum_kcd g_kcbd t_ij^cd c_ak ::
!
!        Form g_kc_db = g_kcbd
!
         call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)     ! Consider moving outside batching loop 
         call allocator(L_db_J, (wf%n_v)*batch_length, wf%n_J) ! L_bd^J 
!
         call wf%get_cholesky_ia(L_kc_J)
!
         reorder = .true.
         call wf%get_cholesky_ab(L_db_J, b_begin, b_end, (wf%n_v)*batch_length, reorder)
!
         call allocator(g_kc_db, (wf%n_o)*(wf%n_v), (wf%n_v)*batch_length)
!
         call dgemm('N','T',                &
                     (wf%n_o)*(wf%n_v),     &
                     (wf%n_v)*batch_length, &
                     wf%n_J,                &
                     one,                   &
                     L_kc_J,                &
                     (wf%n_o)*(wf%n_v),     &
                     L_db_J,                &
                     (wf%n_v)*batch_length, &
                     zero,                  &
                     g_kc_db,               &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
         call deallocator(L_db_J, (wf%n_v)*batch_length, wf%n_J)
!
!        Reorder to g_cd_kb = g_kc_db = g_kcbd 
!
         call allocator(g_cd_kb, (wf%n_v)**2, (wf%n_o)*batch_length)
         g_cd_kb = zero
!
         do b = 1, batch_length
            do k = 1, wf%n_o
!
               kb = index_two(k, b, wf%n_o)
!
               do d = 1, wf%n_v
!
                  db = index_two(d, b, wf%n_v)
!
                  do c = 1, wf%n_v
!
                     kc = index_two(k, c, wf%n_o)
                     cd = index_two(c, d, wf%n_v)
!
                     g_cd_kb(cd, kb) = g_kc_db(kc, db) ! g_kcbd 
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_kc_db, (wf%n_o)*(wf%n_v), (wf%n_v)*batch_length)
!
!        Order amplitudes as t_ij_cd = t_ij^cd 
!
         call allocator(t_ij_cd, (wf%n_o)**2, (wf%n_v)**2)
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
         call allocator(X_ij_kb, (wf%n_o)**2, (wf%n_o)*batch_length)
!
         call dgemm('N','N',                &
                     (wf%n_o)**2,           &
                     (wf%n_o)*batch_length, &
                     (wf%n_v)**2,           &
                     one,                   &
                     t_ij_cd,               &
                     (wf%n_o)**2,           &
                     g_cd_kb,               &
                     (wf%n_v)**2,           &
                     zero,                  &
                     X_ij_kb,               &
                     (wf%n_o)**2)
!
         call deallocator(t_ij_cd, (wf%n_o)**2, (wf%n_v)**2)
!
!        sum_kcd g_kcbd t_ij^cd c_ak = sum_k X_ij_kb c_ak 
!        Reorder to X_k_ijb = X_ij_kb 
!
         call allocator(X_k_ijb, (wf%n_o), batch_length*(wf%n_o)**2)
         X_k_ijb = zero
!
         do b = 1, batch_length
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
         call deallocator(X_ij_kb, (wf%n_o)**2, (wf%n_o)*batch_length)
!
!        Form rho_a_ijb = sum_k c_ak X_k_ijb = sum_k c_a_i(a,k) X_k_ijb(k, ijb)
!
         call allocator(rho_a_ijb, wf%n_v, batch_length*(wf%n_o)**2)
!
         call dgemm('N','N',                   &
                     wf%n_v,                   &
                     batch_length*(wf%n_o)**2, &
                     wf%n_o,                   &
                     -one,                     & 
                     c_a_i,                    &
                     wf%n_v,                   &
                     X_k_ijb,                  &
                     wf%n_o,                   &
                     zero,                     &
                     rho_a_ijb,                &
                     wf%n_v)
!
         call deallocator(X_k_ijb, wf%n_o, batch_length*(wf%n_o)**2)
!
!        Add rho_a_ijb (batch over b) to rho_ai_bj (full space)
!
         do b = 1, batch_length ! Loop over restricted space 
            do j = 1, wf%n_o
!
               Bj = index_two(b + b_begin - 1, j, wf%n_v) ! b in full space 
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
         call deallocator(rho_a_ijb, wf%n_v, batch_length*(wf%n_o)**2)
!
!
!        :: Term 2. sum_kcd g_kcbd t_kj^ad c_ci ::
!
!        Form the intermediate X_i_dkb = sum_c g_kcbd c_ci 
!                                      = sum_c c_ci g_cd_kb 
!                                      = sum_c c_a_i^T(i,c) g_cd_kb(c, dkb)
!
!        Note: g_cd_kb is interpreted as g_c_dkb in the matrix multiplication.
!
         call allocator(X_i_dkb, wf%n_o, (wf%n_v)*(wf%n_o)*batch_length)
!
         call dgemm('T','N',                         &
                     wf%n_o,                         &
                     (wf%n_v)*(wf%n_o)*batch_length, &
                     wf%n_v,                         &
                     one,                            &
                     c_a_i,                          &
                     wf%n_v,                         &
                     g_cd_kb,                        & ! "g_c_dkb" 
                     wf%n_v,                         &
                     zero,                           &
                     X_i_dkb,                        &
                     wf%n_o)
!
!        sum_kcd g_kcbd t_kj^ad c_ci = sum_kd (sum_c c_ci g_kcbd) t_kj^ad
!                                    = sum_kd X_i_dkb t_kj^ad 
!                                    = sum_kd X_ib_dk t_dk_aj
!
!        Reorder to X_ib_dk = X_i_dkb 
!
         call allocator(X_ib_dk, (wf%n_o)*batch_length, (wf%n_v)*(wf%n_o))
!
         do k = 1, wf%n_o
            do d = 1, wf%n_v
!
               dk = index_two(d, k, wf%n_v)
!
               do b = 1, batch_length
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
         call deallocator(X_i_dkb, wf%n_o, (wf%n_o)*(wf%n_v)*batch_length)
!
!        Order the amplitudes as t_dk_aj = t_kj^ad 
!
         call allocator(t_dk_aj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
!        Calculate rho_ib_aj = sum_kcd g_kcbd t_kj^ad c_ci
!                            = sum_dk X_ib_dk t_dk_aj
!
         call allocator(rho_ib_aj, (wf%n_o)*batch_length, (wf%n_o)*(wf%n_v))
!
         call dgemm('N','N',                &
                     (wf%n_o)*batch_length, &
                     (wf%n_o)*(wf%n_v),     &
                     (wf%n_o)*(wf%n_v),     &
                     -one,                  & 
                     X_ib_dk,               &
                     (wf%n_o)*batch_length, &
                     t_dk_aj,               &
                     (wf%n_o)*(wf%n_v),     &
                     zero,                  &
                     rho_ib_aj,             &
                     (wf%n_o)*batch_length)
!
         call deallocator(X_ib_dk, (wf%n_o)*batch_length, (wf%n_o)*(wf%n_v))
         call deallocator(t_dk_aj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Add rho_ib_aj (batch over b) ro rho_ai_bj (full space)
!
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               aj = index_two(a, j, wf%n_v)
!
               do b = 1, batch_length
!
                  Bj = index_two(b + b_begin - 1, j, wf%n_v) ! b is full space index 
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
         call deallocator(rho_ib_aj, (wf%n_o)*batch_length, (wf%n_v)*(wf%n_o))
!
!
!        :: Term 3. sum_kcd g_kcbd t_ik^ca c_dj ::
!
!        sum_d g_kcbd c_dj = sum_d g_cd_kb c_dj
!
!        Reorder integrals to g_cd_kb to g_ckb_d
!
         call allocator(g_ckb_d, (wf%n_o)*(wf%n_v)*batch_length, wf%n_v)
         g_ckb_d = zero
!
         do d = 1, wf%n_v
            do b = 1, batch_length
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
         call deallocator(g_cd_kb, (wf%n_v)**2, (wf%n_o)*batch_length)
!
!        Form the intermediate X_ckb_j = sum_d g_kcbd c_dj = sum_d g_ckb_d c_d_j 
!
         call allocator(X_ckb_j, (wf%n_v)*(wf%n_o)*batch_length, wf%n_o)
!
         call dgemm('N','N',                         &
                     (wf%n_v)*(wf%n_o)*batch_length, &
                     wf%n_o,                         &
                     wf%n_v,                         &
                     one,                            &
                     g_ckb_d,                        &
                     (wf%n_v)*(wf%n_o)*batch_length, &
                     c_a_i,                          &
                     wf%n_v,                         &
                     zero,                           &
                     X_ckb_j,                        &
                     (wf%n_v)*(wf%n_o)*batch_length)
!
!        Order amplitudes as t_ai_ck = t_ik^ca
!
         call allocator(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
!        Form rho_aib_j = sum_kcd g_kcbd t_ik^ca c_dj = sum_ck t_ai_ck X_ckb_j 
!
!        Note: X_ckb_j is interpreted as X_ck_bj in the matrix multiplication.
!        Note: rho_aib_j is interpreted as rho_ai_bj in the matrix multiplication.
!
         call allocator(rho_aib_j, (wf%n_v)*(wf%n_o)*batch_length, wf%n_o)
!
         call dgemm('N','N',                &
                     (wf%n_v)*(wf%n_o),     &
                     (wf%n_o)*batch_length, &
                     (wf%n_v)*(wf%n_o),     &
                     -one,                  & 
                     t_ai_ck,               &
                     (wf%n_v)*(wf%n_o),     &
                     X_ckb_j,               & ! "X_ck_bj"
                     (wf%n_v)*(wf%n_o),     &
                     zero,                  &
                     rho_aib_j,             & ! "rho_ai_bj"
                     (wf%n_v)*(wf%n_o))
!
         call deallocator(X_ckb_j, (wf%n_v)*(wf%n_o)*batch_length, wf%n_o)
         call deallocator(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!        Add rho_aib_j to rho_ai_bj 
!
         do j = 1, wf%n_o
            do b = 1, batch_length
!
               Bj = index_two(b + b_begin - 1, j, wf%n_v)
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
         call deallocator(rho_aib_j, (wf%n_o)*(wf%n_v)*batch_length, wf%n_o)
!
!
!        :: Term 4. - sum_kcd L_kcbd t_ik^ac c_dj :: 
!
!        sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj 
!
!        Form L_ckb_d = L_kcbd = 2 * g_kcbd - g_kdbc = 2 * g_ckb_d(ckb, d) - g_ckb_d(dkb, c)
!
         call allocator(L_ckb_d, (wf%n_o)*(wf%n_v)*batch_length, wf%n_v)
         L_ckb_d = zero
!
         do d = 1, wf%n_v
            do b = 1, batch_length
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
         call deallocator(g_ckb_d, (wf%n_o)*(wf%n_v)*batch_length, wf%n_v)
!
!        Form the intermediate Y_ckb_j = sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj 
!
         call allocator(Y_ckb_j, (wf%n_v)*(wf%n_o)*batch_length, wf%n_o)
!
         call dgemm('N','N',                         &
                     (wf%n_v)*(wf%n_o)*batch_length, & 
                     wf%n_o,                         &
                     wf%n_v,                         &
                     one,                            &
                     L_ckb_d,                        &
                     (wf%n_v)*(wf%n_o)*batch_length, &
                     c_a_i,                          &
                     wf%n_v,                         &
                     zero,                           &
                     Y_ckb_j,                        &
                     (wf%n_v)*(wf%n_o)*batch_length)
!
!        Order amplitudes as t_ai_ck = t_ik^ac = t2am(aick, 1)
!
         call allocator(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_ai_ck = zero 
!
         call squareup(wf%t2am, t_ai_ck, (wf%n_o)*(wf%n_v))
!
!        Form rho_aib_j = - sum_ck t_ai_ck Y_ckb_j 
!
!        Note: we interpret Y_ckb_j as Y_ck_bj in the matrix multiplication
!        Note: we interpret rho_aib_j as rho_ai_bj in the matrix multiplication
!
         call allocator(rho_aib_j, (wf%n_o)*(wf%n_v)*batch_length, wf%n_o)
!
         call dgemm('N','N',                &
                     (wf%n_o)*(wf%n_v),     &
                     (wf%n_o)*batch_length, &
                     (wf%n_o)*(wf%n_v),     &
                     one,                   & 
                     t_ai_ck,               &
                     (wf%n_o)*(wf%n_v),     &
                     Y_ckb_j,               & ! "Y_ck_bj"
                     (wf%n_o)*(wf%n_v),     &
                     zero,                  &
                     rho_aib_j,             & ! "rho_ai_bj"
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(Y_ckb_j, (wf%n_v)*(wf%n_o)*batch_length, wf%n_o)
         call deallocator(t_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!        Add rho_aib_j to rho_ai_bj 
!
         do j = 1, wf%n_o
            do b = 1, batch_length
!
               Bj = index_two(b + b_begin - 1, j, wf%n_v)
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
         call deallocator(rho_aib_j, (wf%n_v)*(wf%n_o)*batch_length, wf%n_o)
!
!
!        :: Term 5. - sum_kcd L_kcbd t_ij^ad c_ck ::
!
!        Form the intermediate X_bd = sum_ck c_ck L_kcbd = sum_ck c_ck L_ckb_d
!
!        Note: L_ckb_d is interpreted as L_ck_bd in the matrix multiplication 
!        Note: c_a_i is interpreted as c_ai in the matrix multiplication 
!
         call allocator(X_bd, 1, batch_length*(wf%n_v))
!  
         call dgemm('N','N',                &
                     1,                     &
                     batch_length*(wf%n_v), &
                     (wf%n_o)*(wf%n_v),     &
                     one,                   &
                     c_a_i,                 & ! "c_ai"
                     1,                     &
                     L_ckb_d,               & ! "L_ck_bd"
                     (wf%n_o)*(wf%n_v),     &
                     zero,                  &
                     X_bd,                  &
                     1)
!
!        Order amplitudes as t_d_aij = t_ij^ad 
!
         call allocator(t_d_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
!        Form rho_b_aij = - sum_kcd L_kcbd t_ij^ad c_ck
!                       = - sum_d X_bd t_d_aij
!
!        Note: X_bd is interpreted as X_b_d in the matrix multiplication
!
         call allocator(rho_b_aij, batch_length, (wf%n_v)*(wf%n_o)**2)
!
         call dgemm('N','N',               &
                     batch_length,         &
                     (wf%n_v)*(wf%n_o)**2, &
                     wf%n_v,               &
                     one,                  & 
                     X_bd,                 & ! "X_b_d"
                     batch_length,         &
                     t_d_aij,              &
                     wf%n_v,               &
                     zero,                 &
                     rho_b_aij,            &
                     batch_length)
!
!
         call deallocator(X_bd, 1, batch_length*(wf%n_v))
         call deallocator(t_d_aij, wf%n_v, (wf%n_v)*(wf%n_o)**2)
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
                  do b = 1, batch_length
!
                     Bj = index_two(b + b_begin - 1, j, wf%n_v)
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
         call deallocator(rho_b_aij, batch_length, (wf%n_v)*(wf%n_o)**2)
         call deallocator(L_ckb_d, (wf%n_v)*(wf%n_o)*batch_length, wf%n_v)
!
      enddo ! End of batches over b 
!
!     Destroy amplitudes from memory
!
      call wf%destruct_amplitudes
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
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))  :: c_ai_ck
!
      real(dp), dimension(:,:), allocatable :: t_dl_bj
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_kc_ld
      real(dp), dimension(:,:), allocatable :: L_ck_dl
      real(dp), dimension(:,:), allocatable :: X_ck_bj
!
      integer(i15) :: c = 0, d = 0, k = 0, l = 0
      integer(i15) :: ck = 0, dl = 0, kc = 0, kd = 0, lc = 0, ld = 0
!
!     Read T2 amplitudes from disk
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      t_dl_bj = zero
!
      call squareup(wf%t2am, t_dl_bj, (wf%n_o)*(wf%n_v))
!
      call wf%destruct_amplitudes
!
!     Construct g_kcld = sum_J L_kc_J * L_ld_J
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ia_J = zero
      call wf%get_cholesky_ia(L_ia_J)
!
      call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'T',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Construct L_kc,ld ordered as L_ck_dl
!
      call allocator(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Intermediate X_ck_bj = sum_dl L_ck_dl * t_dl_bj
!
      call allocator(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call deallocator(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call deallocator(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
      module subroutine jacobian_ccsd_f2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD F2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^F2 =  P_(ij)^(ab) (- sum_ckld t_ai,ck * L_kc,ld * c_bl,dj 
!!                                    - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                                    - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj ) 
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
         real(dp), dimension(:,:), allocatable :: L_d_clk
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
         integer(i15) :: aij = 0, aib = 0, clk = 0, ckd = 0
!
         integer(i15) :: bldj = 0, aidj = 0, bkcl = 0, aibl = 0
!
!        :: Construct L_kc_ld ::
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!   
         call allocator(L_ck_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         L_ck_dl = zero
!
!        Construct L_kc_dl ordered as L_ck_dl
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
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        :: Term 1: - sum_ckld t_ai,ck * L_kc,ld * c_bl,dj ::
!
!
!        Reorder c_bl_dj as c_dl_bj
!
         call allocator(c_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call allocator(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call dgemm('N', 'N', &
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
         call deallocator(c_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_ai_ck = zero
!
         call squareup(wf%t2am, t_ai_ck, (wf%n_o)*(wf%n_v))
         call wf%destruct_amplitudes
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
         call deallocator(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        :: Term 2: - sum_ckdl t_ai,dj * L_kc,ld * c_bk,cl
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Reorder L_ck_dl to L_d_clk
!
         call allocator(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2))
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
!        Reorder c_bl,ck as c_clk_b
!        
         call allocator(c_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         c_clk_b = zero
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
                     c_clk_b(clk,b) = c_ai_bj(bl,ck)
!
                 enddo
               enddo
            enddo
         enddo
!
!        Y_d_b = sum_clk L_d_clk * c_clk_b 
!
         call allocator(Y_d_b, wf%n_v, wf%n_v)
!
         call dgemm('N', 'N',              &
                     wf%n_v,               &
                     wf%n_v,               &
                     ((wf%n_o)**2)*wf%n_v, &
                     one,                  &
                     L_d_clk,              &
                     wf%n_v,               &
                     c_clk_b,              &
                     ((wf%n_o)**2)*wf%n_v, &
                     zero,                 &
                     Y_d_b,                &
                     wf%n_v)
!
         call deallocator(c_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         call deallocator(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2)) 
!
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         t_aij_d = zero
!
!        Reorder T2 amplitudes
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
                     aij = index_three(a, i, j, wf%n_v, wf%n_o)
                     aidj = index_packed(ai,dj)
!
                     t_aij_d(aij, d) = wf%t2am(aidj,1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_amplitudes
!
         call allocator(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
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
         call deallocator(t_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         call deallocator(Y_d_b, wf%n_v, wf%n_v)
!
!        Adding term 2 to rho_ai_bj
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
         call deallocator(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!        :: Term 3: ::
!
!        :: Construct L_kc_ld ::
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!   
         call allocator(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2))
         L_l_ckd = zero
!
!        Construct L_kc_dl ordered as L_ck_dl (L_l_ckd ?)
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
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder c_ck,dj to c_ckd_j 
!
         call allocator(c_ckd_j, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
         c_ckd_j = zero
!
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do c = 1, wf%n_v
!  
                  ck = index_two(c, k, wf%n_v)
!
                  do d = 1, wf%n_v
!
                     dj = index_two(d, j, wf%n_v)
!
                     ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
                     c_ckd_j(ckd,j) = c_ai_bj(ck,dj)
!
                 enddo
               enddo
            enddo
         enddo
!
         call allocator(Z_l_j, wf%n_o, wf%n_o)
!
         call dgemm('N', 'N',              &
                     wf%n_o,               &
                     wf%n_o,               &
                     ((wf%n_v)**2)*wf%n_o, &
                     one,                  &
                     L_l_ckd,              &
                     wf%n_o,               &
                     c_ckd_j,              &
                     ((wf%n_v)**2)*wf%n_o, &
                     zero,                 &
                     Z_l_j,                &
                     wf%n_o)
!
         call deallocator(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2)) 
         call deallocator(c_ckd_j, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
!
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_aib_l, (wf%n_o)*((wf%n_v)**2), wf%n_o)
         t_aib_l = zero
!
!        Reorder T2 amplitudes
!
         do l = 1, wf%n_o
            do i = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bl = index_two(b, l, wf%n_v)
!
                  do a = 1, wf%n_v
!
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
                     ai = index_two(a, i, wf%n_v)
!
                     aibl = index_packed(ai, bl)
!
                     t_aib_l(aib, l) = wf%t2am(aibl,1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_amplitudes
         call allocator(rho_aib_j, (wf%n_o)*((wf%n_v)**2), wf%n_o)
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
                     zero,                   &
                     rho_aib_j,              &
                     ((wf%n_v)**2)*(wf%n_o))
!
         call deallocator(t_aib_l, (wf%n_o)*((wf%n_v)**2), wf%n_o)
         call deallocator(Z_l_j, wf%n_o, wf%n_o)
!
!        Adding term 3 to rho_ai_bj
!
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_aib_j(aib, j) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(rho_aib_j, (wf%n_o)*((wf%n_v)**2), wf%n_o)
! 
      end subroutine jacobian_ccsd_f2_ccsd
!
      module subroutine jacobian_ccsd_g2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD G2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^G2 =  P_(ij)^(ab) (- sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck 
!!                                    - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj 
!!                                    - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl )
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
!        :: Term 1: - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck  ::
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!   
         call allocator(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         L_ck_dl = zero
!
!        Construct L_kc_dl ordered as L_ck_dl
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
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder t_bl_dj as t_dl_bj
!
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
                     t_dl_bj(dl,bj) = wf%t2am(bldj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_amplitudes
!
         call allocator(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call dgemm('N', 'N', &
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
         call deallocator(t_dl_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(L_ck_dl,(wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call deallocator(X_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        :: Term 2: - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Reorder L_ck_dl to L_d_clk
!
         call allocator(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2))
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
!        Reorder t_ck,bl as t_clk_b
!        
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
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
                     ckbl = index_packed(ck, bl)
!
                     t_clk_b(clk,b) = wf%t2am(ckbl, 1)
!
                 enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_amplitudes
!
!        Y_d_b = sum_clk L_d_clk * c_clk_b 
!
         call allocator(Y_d_b, wf%n_v, wf%n_v)
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
         call deallocator(t_clk_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         call deallocator(L_d_clk, wf%n_v, (wf%n_v)*((wf%n_o)**2)) 
!
         call allocator(c_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         c_aij_d = zero
!
!        Reorder c_ai_dl to c_aij_d
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
                     aij = index_three(a, i, j, wf%n_v, wf%n_o)
!
                     c_aij_d(aij, d) = c_ai_bj(ai, dj)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_amplitudes
!
         call allocator(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
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
         call deallocator(c_aij_d, (wf%n_v)*((wf%n_o)**2), wf%n_v)
         call deallocator(Y_d_b, wf%n_v, wf%n_v)
!
!        Adding term 2 to rho_ai_bj
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
         call deallocator(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!        :: Term 3: - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl ::
!
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!   
         call allocator(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2))
         L_l_ckd = zero
!
!        Construct L_kc_dl ordered as L_ck_dl (L_l_ckd?)
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
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder t_ck,dj to t_ckd_j 
!
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_ckd_j, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
         t_ckd_j = zero
!
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do c = 1, wf%n_v
!  
                  ck = index_two(c, k, wf%n_v)
!
                  do d = 1, wf%n_v
!
                     dj = index_two(d, j, wf%n_v)
                     ckd = index_three(c, k, d, wf%n_v, wf%n_o)
!
                     ckdj = index_packed(ck, dj)
!
                     t_ckd_j(ckd,j) = wf%t2am(ckdj,1)
!
                 enddo
               enddo
            enddo
         enddo
!
         call wf%destruct_amplitudes
!
         call allocator(Z_l_j, wf%n_o, wf%n_o)
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
         call deallocator(L_l_ckd,(wf%n_o), (wf%n_o)*((wf%n_v)**2)) 
         call deallocator(t_ckd_j, ((wf%n_v)**2)*(wf%n_o), wf%n_o)
!
         call allocator(c_aib_l, (wf%n_o)*((wf%n_v)**2), wf%n_o)
         c_aib_l = zero
!
!        Reorder c_ai_bl to c_aib_l
!
         do l = 1, wf%n_o
            do i = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bl = index_two(b, l, wf%n_v)
!
                  do a = 1, wf%n_v
!
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     ai = index_two(a, i, wf%n_v)
!
                     c_aib_l(aib, l) = c_ai_bj(ai, bl)
!
                  enddo
               enddo
            enddo
         enddo
!
         call allocator(rho_aib_j, (wf%n_o)*((wf%n_v)**2), wf%n_o)
!
         call dgemm('N','N',                 &
                     ((wf%n_v)**2)*(wf%n_o), &
                     wf%n_o,                 &
                     wf%n_o,                 &
                     -one,                   &
                     c_aib_l,                &
                     ((wf%n_v)**2)*(wf%n_o), &
                     Z_l_j,                  &
                     wf%n_o,                 &
                     zero,                   &
                     rho_aib_j,              &
                     ((wf%n_v)**2)*(wf%n_o))
!
         call deallocator(c_aib_l, (wf%n_o)*((wf%n_v)**2), wf%n_o)
         call deallocator(Z_l_j, wf%n_o, wf%n_o)
!
!        Adding term 3 to rho_ai_bj
!
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_aib_j(aib, j) 
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(rho_aib_j, (wf%n_o)*((wf%n_v)**2), wf%n_o)
! 
      end subroutine jacobian_ccsd_g2_ccsd
!
      module subroutine jacobian_ccsd_h2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD H2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^H2 =  P_(ij)^(ab) ( sum_ckld t_ci,ak * g_kc,ld * c_bl,dj 
!!                                   + sum_ckdl t_cj,al * g_kc,ld * c_bk,di)
!!                
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
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        t_ak,ci ordered as t_ai_kc
!  
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_ai_kc, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
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
         call wf%destruct_amplitudes
!  
         call allocator(X_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call deallocator(t_ai_kc, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call allocator(c_ld_bj, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
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
         call dgemm('N','N',            &
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
         call deallocator(c_ld_bj, (wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
         call deallocator(X_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        :: Term 2: sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!
!        Construct g_kc_ld
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Reorder g_kc_ld to g_lc_kd 
!
         call allocator(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        t_al,cj ordered as t_aj_lc
!  
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_aj_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         t_aj_lc = zero
!
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               aj = index_two(a, j, wf%n_v)
!
               do c = 1, wf%n_v
                  do l = 1, wf%n_o
!
                     cj = index_two(c, j, wf%n_v)
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
         call wf%destruct_amplitudes
!
         call allocator(Y_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call deallocator(g_lc_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(t_aj_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!  
!        Reorder c_bk,di as c_kd_bi
!
         call allocator(c_kd_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call allocator(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call deallocator(c_kd_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(Y_aj_kd, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call deallocator(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      end subroutine jacobian_ccsd_h2_ccsd
!
!
      module subroutine jacobian_ccsd_i2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD I2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^I2 =  P_ij^ab ( sum_c F_bc * c_ai,cj - sum_k F_jk * c_ai,bk
!!                               + sum_ck L_bj,kc * c_ai,ck
!!                               - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj ) )
!!                
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
         real(dp), dimension(:,:), allocatable :: L_bj_J
         real(dp), dimension(:,:), allocatable :: L_kc_J
         real(dp), dimension(:,:), allocatable :: L_kj_J
         real(dp), dimension(:,:), allocatable :: L_bc_J
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
         integer(i15) :: offset = 0
!
         integer(i15) :: required = 0, available = 0, max_batch_length = 0, n_batch = 0, batch_dimension = 0
         integer(i15) :: c_batch = 0, c_first = 0, c_last = 0, c_length = 0
!
!        :: sum_c F_bc * c_ai,cj ::
!
!
!        Reorder c_ai,cj to c_aij_c
!
         call allocator(c_aij_c, (wf%n_v)*((wf%n_o)**2), wf%n_v)
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
         call allocator(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!        sum_c F_bc * c_ai,cj = sum_c c_aij_c(aij,c) F_ab(b,c) = sum_c c_aij_c(aij,c) F_ab^T(c,b)
!
         call dgemm('N','T',                 & ! Second factor is transposed here (E, 28 May)
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
         call deallocator(c_aij_c, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!        Reorder rho_aij_b into rho_ai_bj
!
         do i = 1, wf%n_o
            do j = 1, wf%n_o
               do a = 1, wf%n_v
!  
                  ai  = index_two(a, i, wf%n_v)
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
         call deallocator(rho_aij_b, (wf%n_v)*((wf%n_o)**2), wf%n_v)
!
!       ::  - sum_k F_jk * c_ai,bk  ::
!
!        Reorder c_ai,bk to c_aib_k
!
         call allocator(c_aib_k, (wf%n_o)*((wf%n_v)**2), wf%n_o)
         c_aib_k = zero
!        
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               do k = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bk  = index_two(b, k, wf%n_v)
                     aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                     c_aib_k(aib, k) = c_ai_bj(ai, bk)
!
                  enddo
               enddo
            enddo
         enddo
!
         call allocator(rho_aib_j, (wf%n_o)*((wf%n_v)**2), wf%n_o)
!
         call dgemm('N','T',                &  
                     (wf%n_o)*((wf%n_v)**2),&
                     wf%n_o,                & 
                     wf%n_o,                & 
                     -one,                  & 
                     c_aib_k,               & 
                     (wf%n_o)*((wf%n_v)**2),&
                     wf%fock_ij,            & 
                     wf%n_o,                & 
                     zero,                  & 
                     rho_aib_j,             & 
                     (wf%n_o)*((wf%n_v)**2))
!
         call deallocator(c_aib_k, (wf%n_o)*((wf%n_v)**2), wf%n_o)
!
!        Reorder rho_aib_j into rho_ai_bj
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!  
               ai  = index_two(a, i, wf%n_v)
!                 
               do b = 1, wf%n_v
!
                  aib = index_three(a, i, b, wf%n_v, wf%n_o)
!
                  do j = 1, wf%n_o
!
                     bj = index_two(b, j, wf%n_v)
!
                     rho_ai_bj(ai, bj) = rho_ai_bj(ai, bj) + rho_aib_j(aib, j)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(rho_aib_j, (wf%n_o)*((wf%n_v)**2), wf%n_o)
!
!        ::   sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj ) ::
!            
!         sum_ck ( g_bj,kc*(2*c_ai,ck - c_ak,ci) - g_bc,kj*c_ai,ck - g_ki,bc*c_ak,cj ) 
!
!        Construct g_bj,kc 
!
         call allocator(L_bj_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_bj_J = zero
         call wf%get_cholesky_ai(L_bj_J)
!
         call allocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_kc_J = zero
         call wf%get_cholesky_ia(L_kc_J)
!
         call allocator(g_bj_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call dgemm('N', 'T',          & 
                     (wf%n_o)*(wf%n_v),&
                     (wf%n_o)*(wf%n_v),&
                     wf%n_J,           &
                     one,              &
                     L_bj_J,           &
                     (wf%n_o)*(wf%n_v),&
                     L_kc_J,           &
                     (wf%n_o)*(wf%n_v),&
                     zero,             &
                     g_bj_kc,          &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(L_bj_J, (wf%n_o)*(wf%n_v), wf%n_J)
         call deallocator(L_kc_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Reordering g_bj_kc to g_ck_bj
!
         call allocator(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         g_ck_bj = zero
!
         do j = 1, wf%n_o
            do b = 1, wf%n_v
!
               bj = index_two(b, j, wf%n_v) ! E: bug fix 28 may, was n_o
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
         call deallocator(g_bj_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        sum_ck 2*c_ai_ck * g_ck_bj
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
!        Reorder c_ak,ci to c_ai_ck
!
         call allocator(c_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
!        - sum_ck g_ck_bj*c_ai_ck
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
         call deallocator(c_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Construct g_bc_kj       
!
         call allocator(g_bc_kj, (wf%n_v)**2, (wf%n_o)**2)
         g_bc_kj = zero
!
!        Start batching over c
! 
         required = 2*(wf%n_J)*((wf%n_v)**2) &
                  + 4*(wf%n_J)*(wf%n_v)*(wf%n_o) &
                  + 2*(wf%n_J)*((wf%n_o)**2)
!     
         required = 4*required         ! In words
         available = get_available()
!
         batch_dimension  = wf%n_v ! Batch over the virtual index a
         max_batch_length = 0      ! Initilization of unset variables 
         n_batch          = 0
!
         call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!        Loop over the number of a batches 
!
         do c_batch = 1, n_batch
!
!           For each batch, get the limits for the a index 
!
            call batch_limits(c_first, c_last, c_batch, max_batch_length, batch_dimension)
            c_length = c_last - c_first + 1
!
            call allocator(L_bc_J, (wf%n_v)*c_length, wf%n_J)
            L_bc_J = zero
!
            call wf%get_cholesky_ab(L_bc_J, c_first, c_last, c_length*(wf%n_v), .false.)
!
            call allocator(L_kj_J, (wf%n_o)**2, wf%n_J)
            L_kj_J = zero
!
            call wf%get_cholesky_ij(L_kj_J)
!
            offset = index_two(1, c_first, wf%n_v)
!
            call dgemm('N', 'T',            & 
                        (wf%n_v)*c_length,  &
                        (wf%n_o)**2,        &
                        wf%n_J,             &
                        one,                &
                        L_bc_J,             &               
                        (wf%n_v)*c_length,  &                        
                        L_kj_J,             &
                        (wf%n_o)**2,        &                                
                        one,                &                        
                        g_bc_kj(offset, 1), &   
                        (wf%n_v)**2)
!
            call deallocator(L_bc_J, (wf%n_v)*c_length, wf%n_J)
            call deallocator(L_kj_J, (wf%n_o)**2, wf%n_J)
!
         enddo
!
!        Reorder g_bc_kj
!
         call allocator(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         g_ck_bj = zero
!
         do c = 1, wf%n_v
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
                     ck = index_two(c, k, wf%n_v)
 !    
                     g_ck_bj(ck, bj) = g_bc_kj(bc, kj)
 !
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_bc_kj, (wf%n_v)**2, (wf%n_o)**2)
!
!         - sum_ck c_ai_ck * g_ck_bj       
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
!        Reorder  c_ak,cj to c_aj_ck
!
         call allocator(c_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
         call allocator(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call dgemm('N', 'N',           &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     (wf%n_o)*(wf%n_v), &
                     -one,              &
                     c_aj_ck,           &
                     (wf%n_o)*(wf%n_v), &
                     g_ck_bj,           &   
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     rho_aj_bi,         &
                     (wf%n_o)*(wf%n_v))
!
         call deallocator(g_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
         call deallocator(c_aj_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reorder rho_aj_bi into rho_ai_bj
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
         call deallocator(rho_aj_bi, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      end subroutine jacobian_ccsd_i2_ccsd
!

      module subroutine jacobian_ccsd_j2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD J2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl 
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!                
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
!        Constructing g_kc_ld
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
         call allocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call wf%get_cholesky_ia(L_ia_J)
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_kc_ld,           &
                  (wf%n_o)*(wf%n_v))
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
         call allocator(g_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
         g_kl_cd = zero
!
!        Reorder g_kc_ld to g_kl_cd
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
         call deallocator(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Reordered T2 amplitudes        
!
         call wf%initialize_amplitudes
         call wf%read_double_amplitudes
!
         call allocator(t_ab_ij, (wf%n_v)**2, (wf%n_o)**2 )
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

         call wf%destruct_amplitudes
!
!        X_kl_ij = g_kl_cd * t_cd_ij
!
         call allocator(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
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
!        rho_ab_ij = c_ab_kl * X_kl_ij
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
!        X_kl_ij = g_kl_cd * c_cd_ij
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
!        rho_ab_ij = t_ab_kl * X_kl_ij
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
         call deallocator(X_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
         call deallocator(g_kl_cd, (wf%n_o)**2, (wf%n_v)**2)
         call deallocator(t_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      end subroutine jacobian_ccsd_j2_ccsd
!
!
      module subroutine jacobian_ccsd_k2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD K2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^K2 =    sum_kl g_ki,lj * c_ak,bl 
!!                       + sum_cd g_ac,bd * c_ci,dj
!!                
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: c_ab_ij
!
         real(dp), dimension(:,:), allocatable :: L_ij_J
         real(dp), dimension(:,:), allocatable :: L_ca_J
         real(dp), dimension(:,:), allocatable :: L_db_J
!
         real(dp), dimension(:,:), allocatable :: g_ki_lj
         real(dp), dimension(:,:), allocatable :: g_kl_ij
         real(dp), dimension(:,:), allocatable :: g_ca_db
         real(dp), dimension(:,:), allocatable :: g_ab_cd
!
         real(dp), dimension(:,:), allocatable :: rho_batch_ab_ij
!
         integer(i15) :: i = 0, j = 0, k = 0, l = 0  
         integer(i15) :: a = 0, b = 0, c = 0, d = 0  
!
         integer(i15) :: ab = 0, db = 0, ca = 0, cd = 0, full_ab = 0
         integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0
!
!        Batching and memory handling variables
!
         integer(i15) :: a_n_batch = 0, a_first = 0, a_last = 0, a_length = 0, a_max_length = 0, a_batch = 0
         integer(i15) :: b_n_batch = 0, b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0
!
         integer(i15) :: required = 0, available = 0
!
         call allocator(L_ij_J, (wf%n_o)**2, wf%n_J)
         L_ij_J = zero
!
         call wf%get_cholesky_ij(L_ij_J)
!
         call allocator(g_ki_lj, (wf%n_o)**2, (wf%n_o)**2)
!
         call dgemm('N', 'T',     &
                     (wf%n_o)**2, &
                     (wf%n_o)**2, &
                     wf%n_J,      &
                     one,         &
                     L_ij_J,      &
                     (wf%n_o)**2, &
                     L_ij_J,      &
                     (wf%n_o)**2, &
                     zero,        &
                     g_ki_lj,     &
                     (wf%n_o)**2)
!
         call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J)
!
!        Reorder g_ki_lj to g_kl_ij
!
         call allocator(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
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
         call deallocator(g_ki_lj, (wf%n_o)**2, (wf%n_o)**2)
!
!        sum_kl g_ki,lj * c_ak,bl = sum_kl c_ab_ij(ab,kl) g_kl_ij(kl,ij)
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
         call deallocator(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!        Prepare for batching over a and b
!
!        ::  sum_cd g_ac,bd * c_ci,dj ::
!
         required = max(3*(wf%n_v)**2*(wf%n_J) + 2*(wf%n_v)*(wf%n_o)*(wf%n_J),      & ! Needed to get  L_db_J
                     (wf%n_v)**4 + 2*(wf%n_v)**2*(wf%n_J))                            ! Needed to get g_ac_bd
!
         required = required*4  ! Words

         available = get_available()
!
         a_max_length = 0
         call num_two_batch(required, available, a_max_length, a_n_batch, wf%n_v)
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
            call batch_limits(a_first, a_last, a_batch, a_max_length, wf%n_v)
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
               call batch_limits(b_first ,b_last ,b_batch, b_max_length, wf%n_v)
               b_length = b_last - b_first + 1 
!
!              Get Cholesky vectors L_ac^J ordered as L_ca_J
!
               call allocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
               L_ca_J = zero
!
               call wf%get_cholesky_ab(L_ca_J, a_first, a_last, (wf%n_v)*a_length, .true.)
!
!              Get Cholesky vectors L_bd^J ordered as L_db_J
!
               call allocator(L_db_J, (wf%n_v)*b_length, wf%n_J)
               L_db_J = zero
!  
               call wf%get_cholesky_ab(L_db_J, b_first, b_last, (wf%n_v)*b_length, .true.)
!
!              Allocate g_ca_db = g_acbd
!
               call allocator(g_ca_db, (wf%n_v)*a_length, (wf%n_v)*b_length)
               g_ca_db = zero
!
!              g_ca_db = sum_J L_ca_J*L_db_J
!     
               call dgemm('N','T',            &
                           (wf%n_v)*a_length, &
                           (wf%n_v)*b_length, &
                           wf%n_J,            &
                           one,               &
                           L_ca_J,            &
                           (wf%n_v)*a_length, &  
                           L_db_J,            &
                           (wf%n_v)*b_length, &
                           zero,              &
                           g_ca_db,           &
                           (wf%n_v)*a_length)
!
!              Deallocate L_db_J 
!
               call deallocator(L_db_J, (wf%n_v)*b_length, wf%n_J)
!
!              Deallocate L_ca_J
!
               call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J) 
!
!              sum_cd g_ac,bd * c_ci,dj = sum_cd g_ac,bd c_cd,ij = sum_cd g_ab_cd c_cd_ij 
!
!              Reorder g_ca_db into g_ab_cd 
!              (Here, g_ab_cd = g_acbd = g_ca_db.)
!
               call allocator(g_ab_cd, a_length*b_length, (wf%n_v)**2) 
               g_ab_cd = zero
!
               do b = 1, b_length
                  do a = 1, a_length
!
                     ab = index_two(a, b, a_length)
!
                     do d = 1, wf%n_v
!
                        db = index_two(d, b, wf%n_v)
!
                        do c = 1, wf%n_v
!
                           ca = index_two(c, a, wf%n_v)
                           cd = index_two(c, d, wf%n_v)
!
                           g_ab_cd(ab, cd) = g_ca_db(ca, db) ! = g_acbd 
!
                        enddo
                     enddo
                  enddo
               enddo
!
               call deallocator(g_ca_db, (wf%n_v)*a_length, (wf%n_v)*b_length) 
!
               call allocator(rho_batch_ab_ij,  a_length*b_length, (wf%n_o)**2)
!
!              rho_ab_ij = sum_cd g_ac,bd * c_ci,dj = sum_cd g_ab_cd(ab, cd) c_ab_ij(cd, ij)
!
               call dgemm('N', 'N',            &  
                            a_length*b_length, &
                            (wf%n_o)**2,       &  
                            (wf%n_v)**2,       &  
                            one,               &  
                            g_ab_cd,           &
                            a_length*b_length, &
                            c_ab_ij,           & ! E: Should take care of batching here, eventually.
                            (wf%n_v)**2,       &  
                            zero,              &
                            rho_batch_ab_ij,   &
                            a_length*b_length)
!               
               call deallocator(g_ab_cd, a_length*b_length, (wf%n_v)**2)
!
!              Reorder into rho_ab_ij
!
               do b = 1, b_length
                  do a = 1, a_length
!
                     ab = index_two(a, b, a_length)
                     full_ab = index_two(a + a_first - 1, b + b_first - 1, wf%n_v)
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
               call deallocator(rho_batch_ab_ij,  a_length*b_length, (wf%n_o)**2) 
!
            enddo ! End batches of b 
         enddo ! End batches of a 
!         
      end subroutine jacobian_ccsd_k2_ccsd
!
!
end submodule jacobian