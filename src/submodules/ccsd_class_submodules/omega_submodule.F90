submodule (ccsd_class) omega
!
!!
!!    Omega submodule (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    initialize_omega: allocates the projection vector (omega1, omega2),
!!                      if it is not allocated, and sets it to zero.
!!
!!    destruct_omega:   deallocates the projection vector (omega1, omega2), 
!!                      if it is allocated.
!!
!!    construct_omega:  constructs the projection vector (omega1, omega2) 
!!                      for the current amplitudes (t1am, t2am) for the
!!                      wavefunction object wf. The routine assumes that
!!                      the projection vector is allocated.
!!
!!    The construct omega routine adds the CCS as well as the CCSD contributions,
!!    which are defined in the following routines (defined below):
!!
!!    omega_ccsd_a1: adds A1 term to omega1
!!    omega_ccsd_b1: adds B1 term to omega1
!!    omega_ccsd_c1: adds C1 term to omega1
!!
!!    omega_ccsd_a2: adds A2 term to omega2
!!    omega_ccsd_b2: adds B2 term to omega2
!!    omega_ccsd_c2: adds C2 term to omega2
!!    omega_ccsd_d2: adds D2 term to omega2
!!    omega_ccsd_e2: adds E2 term to omega2
!!
!
   implicit none 
!
   logical :: debug = .false.
!
   real(dp) :: t0 
   real(dp) :: t1
!
   character(len=40) :: integral_type ! Here: electronic repulsion, 1/r_ij, g_pqrs 
!
contains
!
!
   module subroutine initialize_omega_ccsd(wf)
!!
!!    Initialize omega (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Allocates the projection vector (omega1, omega2) and sets it
!!    to zero.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      if (.not. allocated(wf%omega1)) call wf%mem%alloc(wf%omega1, wf%n_v, wf%n_o)
      wf%omega1 = zero
!
      if (.not. allocated(wf%omega2)) call wf%mem%alloc(wf%omega2, wf%n_t2am, 1)
      wf%omega2 = zero
!
   end subroutine initialize_omega_ccsd
!
!
   module subroutine destruct_omega_ccsd(wf)
!!
!!    Destruct omega (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the projection vector.
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%omega1)) call wf%mem%dealloc(wf%omega1, wf%n_v, wf%n_o)
      if (allocated(wf%omega2)) call wf%mem%dealloc(wf%omega2, wf%n_t2am, 1)
!
   end subroutine destruct_omega_ccsd
!
!
   module subroutine construct_omega_ccsd(wf)
!!
!!     Construct omega (CCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!     for the current amplitudes of the object wfn 
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp) :: begin_timer, end_timer
!
      real(dp) :: ccsd_a1_time
      real(dp) :: ccsd_b1_time 
      real(dp) :: ccsd_c1_time
      real(dp) :: ccs_a1_time 
!
      real(dp) :: ccsd_a2_time
      real(dp) :: ccsd_b2_time 
      real(dp) :: ccsd_c2_time 
      real(dp) :: ccsd_d2_time 
      real(dp) :: ccsd_e2_time  
!
!     Set the omega vector to zero 
!
      wf%omega1 = zero
      wf%omega2 = zero
!
!     Construct singles contributions 
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_a1
      call cpu_time(end_timer)
      ccsd_a1_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_b1
      call cpu_time(end_timer)
      ccsd_b1_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_c1
      call cpu_time(end_timer)
      ccsd_c1_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccs_a1
      call cpu_time(end_timer)
      ccs_a1_time = end_timer - begin_timer
!
!     Construct doubles contributions 
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_a2
      call cpu_time(end_timer)
      ccsd_a2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_b2
      call cpu_time(end_timer)
      ccsd_b2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_c2
      call cpu_time(end_timer)
      ccsd_c2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_d2
      call cpu_time(end_timer)
      ccsd_d2_time = end_timer - begin_timer
!
      call cpu_time(begin_timer)
      call wf%omega_ccsd_e2
      call cpu_time(end_timer)
      ccsd_e2_time = end_timer - begin_timer
!
!     Print timings
!
      if (wf%settings%print_level == 'developer') then 
!
         write(unit_output,'(/t3,a)') 'Breakdown of CCSD omega timings:'
!
         write(unit_output,'(/t6,a26,f14.8)') 'Time in CCSD A1 (seconds):', ccsd_a1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD B1 (seconds):', ccsd_b1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD C1 (seconds):', ccsd_c1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCS  A1 (seconds):', ccs_a1_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD A2 (seconds):', ccsd_a2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD B2 (seconds):', ccsd_b2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD C2 (seconds):', ccsd_c2_time
         write(unit_output,'(t6,a26,f14.8)')  'Time in CCSD D2 (seconds):', ccsd_d2_time
         write(unit_output,'(t6,a26,f14.8/)') 'Time in CCSD E2 (seconds):', ccsd_e2_time
!
         flush(unit_output)
!
      endif
!
   end subroutine construct_omega_ccsd
!
!
   module subroutine omega_ccsd_a1_ccsd(wf)
!!
!!    Omega A1 term
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!  
!!    Calculates the A1 term, 
!!  
!!       A1: sum_ckd g_adkc * u_ki^cd,
!!  
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf.
!
      implicit none
!
      class(ccsd) :: wf
!
!     Batching variables 
!
      integer(i15) :: required        = 0
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a 
!
      integer(i15) :: ad_dim = 0
!
      real(dp), dimension(:,:), allocatable :: g_ad_kc ! g_adkc; a is being batched over
      real(dp), dimension(:,:), allocatable :: u_dk_ci ! u_ki^cd
      real(dp), dimension(:,:), allocatable :: t_ck_di ! t_ki^cd 
      real(dp), dimension(:,:), allocatable :: t_dk_ci ! t_ki^cd 
!
!     Form u_ck_di = u_ki^cd 
!
!     Let t_dk_ci = t_ki^cd 
!     Then u_dk_ci = 2 * t_ki^cd - t_ik^cd = 2 * t_dk_ci - t_di_ck
!
      call wf%mem%alloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup_and_sort_1234_to_3214(wf%t2am, t_dk_ci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%mem%alloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call add_1432_to_1234(-one, t_dk_ci, zero, u_dk_ci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_dk_ci, 1, u_dk_ci, 1)
!
      call wf%mem%dealloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Prepare for batching 
!
!     Estimated memory required to construct g_adkc
!
      required = wf%get_vvov_required_mem()
!
!     Initialization of the batching variable
!
      call batch_a%init(wf%n_v)                ! Initialize batching index a 
      call wf%mem%num_batch(batch_a, required) ! Determine batching information   
!
!     Loop over the number of a batches 
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index 
!
         call batch_a%determine_limits(current_a_batch)
!
         ad_dim = (batch_a%length)*(wf%n_v)
!
!        Form g_ad_kc = g_adkc
!
         call wf%mem%alloc(g_ad_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_ov(integral_type, & 
                           g_ad_kc,       &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Calculate the A1 term, 
!
!           sum_ckd g_adkc * u_ki^cd = sum_ckd g_ad_kc u_dkc_i,
! 
!        and add it to the omega vector
! 
         call dgemm('N','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     g_ad_kc,                    & ! g_a_dkc
                     batch_a%length,             &
                     u_dk_ci,                    & ! u_dkc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     wf%omega1(batch_a%first,1), & 
                     wf%n_v)
!
         call wf%mem%dealloc(g_ad_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
      enddo ! End of batches of the index a 
!
!     Deallocate vectors 
!
      call wf%mem%dealloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_a1_ccsd
!
!
   module subroutine omega_ccsd_b1_ccsd(wf)
!!
!!    Omega B1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Calculates the B1 term, 
!!
!!        B1:  - sum_ckl u_kl^ac * g_kilc,
!! 
!!    and adds it to the singles projection vector (omeg1) of
!!    the wavefunction object wf
!!
      implicit none
!
      class(ccsd) :: wf 
!
      integer(i15) :: a = 0, c = 0, k = 0, l = 0, ckl = 0, ki = 0
      integer(i15) :: ak = 0, akcl = 0, al = 0, alck = 0, ck = 0, ai = 0
      integer(i15) :: cl = 0, lc = 0, i = 0, j = 0
!
      real(dp), dimension(:,:), allocatable :: g_lc_ki ! g_kilc 
      real(dp), dimension(:,:), allocatable :: t_al_ck ! t_kl^ac
      real(dp), dimension(:,:), allocatable :: u_al_ck ! u_kl^ac = 2 t_kl^ac - t_lk^ac
!
!     Get g_lc_ki = g_kilc = g_lcki 
!
      call wf%mem%alloc(g_lc_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_oo(integral_type, g_lc_ki)
!
!     Form u_al_ck = u_kl^ac = 2 * t_kl^ac - t_lk^ac
!
!     Squareup amplitudes and reorder: t_ak_cl to t_al_ck
!     u_al_ck = 2 * t_kl^ac - t_lk^ac = 2 * t_al_ck(al,ck) - t_al_ck(ak,cl) 
!     
      call wf%mem%alloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call squareup_and_sort_1234_to_1432(wf%t2am, t_al_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%mem%alloc(u_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call add_1432_to_1234(-one, t_al_ck, zero, u_al_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_al_ck, 1, u_al_ck, 1)
!
      call wf%mem%dealloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Calculate the B1 term, - sum_ckl u_a_ckl g_ckl_i
!
      call dgemm('N','N',                 &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  u_al_ck,                & ! u_a_lck 
                  wf%n_v,                 &
                  g_lc_ki,                & ! g_lck_i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  wf%omega1,              &
                  wf%n_v) 
!
!     Deallocate remaining vectors 
!
      call wf%mem%dealloc(u_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%mem%dealloc(g_lc_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
   end subroutine omega_ccsd_b1_ccsd
!
!
   module subroutine omega_ccsd_c1_ccsd(wf)        
!!  
!!    Omega C1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!  
!!    Calculates the C1 term of omega,
!!  
!!       C1: sum_ck F_kc*u_ai_ck,
!!  
!!    and adds it to the projection vector (omega1) of    
!!    the wavefunction object wf                           
!!  
!!    u_ai_kc = 2*t_ck_ai - t_ci_ak
!! 
      implicit none
!
      class(ccsd) :: wf 
!  
      real(dp), dimension(:,:), allocatable :: F_c_k ! F_kc 
!
      real(dp), dimension(:,:), allocatable :: u_ai_ck  
      real(dp), dimension(:,:), allocatable :: t_ai_ck 
!
      integer(i15) :: k, c 
!
!     Form u_ai_ck = u_ik^ac = 2*t_ik^ac - t_ki^ac
!
      call wf%mem%alloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2am, t_ai_ck, (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))  
!
      call add_1432_to_1234(-one, t_ai_ck, zero, u_ai_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ai_ck, 1, u_ai_ck, 1)
!
!     Reorder the Fock matrix F_ck = F_kc 
!
      call wf%mem%alloc(F_c_k, wf%n_v, wf%n_o)
!
      call sort_12_to_21(wf%fock_ia, F_c_k, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  F_c_k,             & ! Equal to F_kc 
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  wf%omega1,         &
                  (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(F_c_k, wf%n_v, wf%n_o)
!
   end subroutine omega_ccsd_c1_ccsd
!
!
   module subroutine omega_ccsd_a2_ccsd(wf)
!!
!!    Omega A2 term
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!!
!!    A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!!
!!    Structure: Batching over both a and b for A2.2.
!!                t^+_ci_dj = t_ci_dj + t_di_cj
!!                t^-_ci_dj = t_ci_dj - t_di_cj
!!                g^+_ac_bd = g_ac_bd + g_bc_ad 
!!                g^-_ac_bd = g_ac_bd - g_bc_ad 
!! 
!!                omega_A2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!!                omega_A2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj 
      real(dp), dimension(:,:), allocatable :: g_aibj 
      real(dp), dimension(:,:), allocatable :: g_ac_bd 
      real(dp), dimension(:,:), allocatable :: g_p_ab_cd
      real(dp), dimension(:,:), allocatable :: g_m_ab_cd
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_p_cd_ij
      real(dp), dimension(:,:), allocatable :: t_m_cd_ij
!
!     Reordered omega 2 
!  
      real(dp), dimension(:,:), allocatable :: omega2_p_ab_ij
      real(dp), dimension(:,:), allocatable :: omega2_m_ab_ij
! 
!     Indices
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ab = 0, ca = 0, bc = 0, cd = 0, ad = 0, db = 0, ac = 0, bd = 0 
      integer(i15) :: ai = 0, aj = 0, bj = 0, bi = 0, ci = 0, cj = 0, dj = 0, di = 0
      integer(i15) :: ij = 0
!
      integer(i15) :: aibj = 0, biaj = 0, cidj = 0, dicj = 0
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
      integer(i15) :: threads = 0
      integer(i15) :: omp_get_num_threads
!
!     Timing variables
!
      real(dp) :: time_non_integral_part
      real(dp) :: a2_begin_time, a2_end_time
      real(dp) :: a2_begin_time_1, a2_end_time_1
!
      time_non_integral_part = zero
!
      call cpu_time(a2_begin_time)
!
!     :: Calculate the A2.1 term of omega ::
!
!     Create g_ai_bj
!
      call wf%mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_vo(integral_type, g_ai_bj)
!
!     Add A2.1 to Omega 2
!
      call add_to_packed(wf%omega2, g_ai_bj, (wf%n_o)*(wf%n_v))

      call wf%mem%dealloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!    ::  Calculate the A2.2 term  of omega ::
!
      required = wf%get_vvvv_required_mem() + 4*(wf%n_o**2)*(wf%n_v**2) + 2*(wf%n_v**4)
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
!        Determine the limits for the a-batch 
!
         call batch_a%determine_limits(current_a_batch)
!
!        Start looping over b-batches
!
         do current_b_batch = 1, batch_b%num_batches 
!  
!           Determine the limits for the b-batch 
!
            call batch_b%determine_limits(current_b_batch)
! 
!           Allocate g_ca_db
!
            call wf%mem%alloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!           Get g_ac_bd
!
            call cpu_time(a2_begin_time_1)
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
            call cpu_time(a2_end_time_1)
            time_non_integral_part = time_non_integral_part & 
                                   - a2_end_time_1 + a2_begin_time_1 
!
            if (current_b_batch .eq. current_a_batch) then
!
!              Allocate for +-g, +-t
!
               call wf%mem%alloc(g_p_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
               call wf%mem%alloc(g_m_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
               call wf%mem%alloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
               call wf%mem%alloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
               g_p_ab_cd = zero
               g_m_ab_cd = zero
               t_p_cd_ij = zero
               t_m_cd_ij = zero
!
!              Reorder g_ca_db to g_ab_cd and t_ci_dj to t_cd_ij
! 
!$omp parallel do schedule(static) private(d,b,a,ac,cd,bd,bc,ab,ad,i,j,ij,ci,dj,cj,di,cidj,dicj)
               do c = 1, wf%n_v
!
                  do d = 1, c
!
                     cd = index_packed(c, d)
!
                     do a = 1, batch_a%length
!
                        ac = index_two(a, c, batch_a%length)
                        ad = index_two(a, d, batch_a%length)
!
                        do  b = 1, batch_b%length
                           if ((a+batch_a%first-1) .ge. (b+batch_b%first-1)) then
!                            
                              bd = index_two(b, d, batch_b%length)
                              bc = index_two(b, c, batch_b%length)
!
                              ab = index_packed(a, b)
! 
                              g_p_ab_cd(ab, cd) = g_ac_bd(ac, bd) + g_ac_bd(ad, bc)
                              g_m_ab_cd(ab, cd) = g_ac_bd(ac, bd) - g_ac_bd(ad, bc)
!
                             if(c .ne. d) then
                               g_p_ab_cd(ab, cd) = two*g_p_ab_cd(ab, cd)
                               g_m_ab_cd(ab, cd) = two*g_m_ab_cd(ab, cd)
                             endif
!                             
                           endif
                        enddo
                     enddo
!
                    do i = 1, wf%n_o
                       do j = 1, i
!    
                          ij = index_packed(i, j)
!    
                           ci = index_two(c, i, wf%n_v)
                           dj = index_two(d, j, wf%n_v)
                           cj = index_two(c, j, wf%n_v)
                           di = index_two(d, i, wf%n_v)
!
                           cidj = index_packed(ci, dj)
                           dicj = index_packed(cj, di)
! 
                           t_p_cd_ij(cd, ij) = wf%t2am(cidj, 1) + wf%t2am(dicj, 1)
                           t_m_cd_ij(cd, ij) = wf%t2am(cidj, 1) - wf%t2am(dicj, 1)  
!
                       enddo
                    enddo
                 enddo
              enddo
!$omp end parallel do
!
!              Dellocate g_ac_bd 
!
               call wf%mem%dealloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!              Allocate omega +-
!
              call wf%mem%alloc(omega2_p_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
              call wf%mem%alloc(omega2_m_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
! 
!              omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
! 
              call dgemm('N','N',                      & 
                          packed_size(batch_a%length), &
                          packed_size(wf%n_o),         &
                          packed_size(wf%n_v),         &
                          one/four,                    &
                          g_p_ab_cd,                   &
                          packed_size(batch_a%length), &
                          t_p_cd_ij,                   &
                          packed_size(wf%n_v),         &
                          zero,                        &
                          omega2_p_ab_ij,              &
                          packed_size(batch_a%length))
!
              call dgemm('N','N',                      & 
                          packed_size(batch_a%length), &
                          packed_size(wf%n_o),         &
                          packed_size(wf%n_v),         &
                          one/four,                    &
                          g_m_ab_cd,                   &
                          packed_size(batch_a%length), &
                          t_m_cd_ij,                   &
                          packed_size(wf%n_v),         &
                          zero,                        &
                          omega2_m_ab_ij,              &
                          packed_size(batch_a%length) )
!
!             Deallocate +-g, +-t
! 
              call wf%mem%dealloc(g_p_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
              call wf%mem%dealloc(g_m_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
              call wf%mem%dealloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
              call wf%mem%dealloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
              do i = 1, wf%n_o
                 do j = 1, i
!
                    ij = index_packed(i, j)
!
                    do a = 1, batch_a%length
!
                       ai = index_two(a + batch_a%first - 1, i, wf%n_v) ! A is full-CCSD-space a index
                       aj = index_two(a + batch_a%first - 1, j, wf%n_v) ! A is full-CCSD-space a index
!
                       do b = 1, batch_b%length
!                
                          if ((a+batch_a%first-1) .ge. (b+batch_b%first-1)) then
                             bj = index_two(b + batch_b%first - 1, j, wf%n_v) ! B is full-CCSD-space b index
                             bi = index_two(b + batch_b%first - 1, i, wf%n_v) ! B is full-CCSD-space b index
!
                             ab = index_packed(a, b)
!    
                             aibj = index_packed(ai, bj)
                             biaj = index_packed(bi, aj)
!                         
!                            Reorder into omega2_aibj
! 
                             wf%omega2(aibj,1) = wf%omega2(aibj, 1) &
                                                   + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                             if (aibj .ne. biaj) then
                                wf%omega2(biaj,1) = wf%omega2(biaj, 1) &
                                                   + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
                             endif  
                          endif 
!    
                       enddo
                    enddo
                 enddo
              enddo
!
!              Deallocate omega +-
!
               call wf%mem%dealloc(omega2_p_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
               call wf%mem%dealloc(omega2_m_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
            else
!
!              Allocate for +-g, +-t
!
               call wf%mem%alloc(g_p_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call wf%mem%alloc(g_m_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call wf%mem%alloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
               call wf%mem%alloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
               g_p_ab_cd = zero
               g_m_ab_cd = zero
               t_p_cd_ij = zero
               t_m_cd_ij = zero
! 
!$omp parallel do schedule(static) private(d,b,a,ac,cd,bd,bc,ab,ad,i,j,ij,ci,dj,cj,di,cidj,dicj)
               do c = 1, wf%n_v 
                  do d = 1, c
!
                     cd = index_packed(c, d)
!
                     do a = 1, batch_a%length
!
                        ac = index_two(a, c, batch_a%length)
                        ad = index_two(a, d, batch_a%length)
!
                        do  b = 1, batch_b%length
!                            
                             bd = index_two(b, d, batch_b%length)
                             bc = index_two(b, c, batch_b%length)
!
                             ab = index_two(a, b, batch_a%length)
! 
                             g_p_ab_cd(ab, cd) = g_ac_bd(ac, bd) + g_ac_bd(ad, bc)
                             g_m_ab_cd(ab, cd) = g_ac_bd(ac, bd) - g_ac_bd(ad, bc)
!
                            if(c .ne. d) then
                              g_p_ab_cd(ab, cd) = two*g_p_ab_cd(ab, cd)
                              g_m_ab_cd(ab, cd) = two*g_m_ab_cd(ab, cd)
                            endif
!                            
                       enddo
                    enddo
!
                    do i = 1, wf%n_o
                       do j = 1, i
!    
                           ij = index_packed(i, j)
!    
                           ci = index_two(c, i, wf%n_v)
                           dj = index_two(d, j, wf%n_v)
                           cj = index_two(c, j, wf%n_v)
                           di = index_two(d, i, wf%n_v)
!
                           cidj = index_packed(ci, dj)
                           dicj = index_packed(cj, di)
! 
                           t_p_cd_ij(cd, ij) = wf%t2am(cidj, 1) + wf%t2am(dicj, 1)
                           t_m_cd_ij(cd, ij) = wf%t2am(cidj, 1) - wf%t2am(dicj, 1)  
!
                       enddo
                    enddo
                 enddo
              enddo
!$omp end parallel do
!
!              Dellocate g_ac_bd 
!
               call wf%mem%dealloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!              Allocate omega +-
!
               call wf%mem%alloc(omega2_p_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
               call wf%mem%alloc(omega2_m_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
!  
!               omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
! 
              call dgemm('N','N',                            & 
                          (batch_a%length)*(batch_b%length), &
                          packed_size(wf%n_o),               &
                          packed_size(wf%n_v),               &
                          one/four,                          &
                          g_p_ab_cd,                         &
                          (batch_a%length)*(batch_b%length), &
                          t_p_cd_ij,                         &
                          packed_size(wf%n_v),               &
                          zero,                              &
                          omega2_p_ab_ij,                    &
                          (batch_a%length)*(batch_b%length))
!
              call dgemm('N','N',                            & 
                          (batch_a%length)*(batch_b%length), &
                          packed_size(wf%n_o),               &
                          packed_size(wf%n_v),               &
                          one/four,                          &
                          g_m_ab_cd,                         &
                          (batch_a%length)*(batch_b%length), &
                          t_m_cd_ij,                         &
                          packed_size(wf%n_v),               &
                          zero,                              &
                          omega2_m_ab_ij,                    &
                          (batch_a%length)*(batch_b%length))
!
!              Deallocate +-g, +-t
! 
               call wf%mem%dealloc(g_p_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call wf%mem%dealloc(g_m_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call wf%mem%dealloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
               call wf%mem%dealloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
               do i = 1, wf%n_o
                  do j = 1, i
!   
                     ij = index_packed(i, j)
!
                     do a = 1, batch_a%length
!
                        ai = index_two(a + batch_a%first - 1, i, wf%n_v) ! A is full-space a index
                        aj = index_two(a + batch_a%first - 1, j, wf%n_v) ! A is full-space a index
!
                        do b = 1, batch_b%length
!                 
                              bj = index_two(b + batch_b%first - 1, j, wf%n_v) ! B is full-space b index
                              bi = index_two(b + batch_b%first - 1, i, wf%n_v) ! B is full-space b index
!
                              ab = index_two(a, b, batch_a%length)
!     
                              aibj = index_packed(ai, bj)
                              biaj = index_packed(bi, aj)
!                          
!                             Reorder into omega2_aibj
!  
                              wf%omega2(aibj,1) = wf%omega2(aibj, 1) &
                                          + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                              if (aibj .ne. biaj) then
!
                                 wf%omega2(biaj,1) = wf%omega2(biaj, 1) &
                                          + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
!
                              endif   
!     
                        enddo
                     enddo
                  enddo
               enddo
!
!              Deallocate omega +-
!
               call wf%mem%dealloc(omega2_p_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
               call wf%mem%dealloc(omega2_m_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
!
            endif
!
         enddo ! End batching over b
      enddo ! End batching over a
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine omega_ccsd_b2_ccsd(wf)
!!
!!    Omega B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 11 Mar 2017
!! 
!!    Omega B2 = sum_(kl) t_ak_bl*(g_kilj + sum_(cd) t_ci_dj * g_kc_ld)
!!
!!    Structure: g_kilj is constructed first and reordered as g_kl_ij. 
!!    Then the contraction over cd is performed, and the results added to g_kl_ij.
!!    t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!!
      implicit none
!
      class(ccsd) :: wf 
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kc_ld    
      real(dp), dimension(:,:), allocatable :: g_kl_cd    
      real(dp), dimension(:,:), allocatable :: g_kl_ij    
      real(dp), dimension(:,:), allocatable :: g_ki_lj 
!
!     Reordered T2 apmlitudes
!   
      real(dp), dimension(:,:), allocatable :: t_cd_ij    
      real(dp), dimension(:,:), allocatable :: t_ab_kl 
      real(dp), dimension(:,:), allocatable :: t_ci_dj  
!
!     Intermediate for matrix multiplication
! 
      real(dp), dimension(:,:), allocatable :: X_kl_ij 
!
!     Reordered omega
!   
      real(dp), dimension(:,:), allocatable :: omega_ab_ij
      real(dp), dimension(:,:), allocatable :: omega_ai_bj
!
!     Indices
!   
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ab = 0, cd = 0
      integer(i15) :: ai = 0, ak = 0, bj = 0, bl = 0, ci = 0, dj = 0
      integer(i15) :: kc = 0, ld = 0
      integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0
!
      integer(i15) :: aibj = 0, akbl = 0, cidj = 0 
!
!     Allocate and construct g_ki_lj
!
      call wf%mem%alloc(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o)) 
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_oo(integral_type,g_ki_lj)
!
      call wf%mem%alloc(g_kl_ij, (wf%n_o)*(wf%n_o),(wf%n_o)*(wf%n_o))
!
      call sort_1234_to_1324(g_ki_lj, g_kl_ij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%mem%dealloc(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
!
!     Allocate and construct g_kc_ld
!
      call wf%mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!  
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kc_ld)
!
!     Reorder g_kc_ld as g_kl_cd
!
      call wf%mem%alloc(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
      call sort_1234_to_1324(g_kc_ld, g_kl_cd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ci_dj as t_cd_ij
!
      call wf%mem%alloc(t_cd_ij, (wf%n_v)**2, (wf%n_o)**2)
      call squareup_and_sort_1234_to_1324(wf%t2am, t_cd_ij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_cd,     &
                  (wf%n_o)**2, &
                  t_cd_ij,     &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_ij,     &
                  (wf%n_o)**2)
!
!     Deallocate t_cd_ij and g_kl_cd
!
      call wf%mem%dealloc(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     omega_ab_ij = sum_(kl) t_ab_kl*X_kl_ij
!
      call wf%mem%alloc(omega_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_cd_ij,     & ! t_ab_kl
                  (wf%n_v)**2, &
                  g_kl_ij,     &
                  (wf%n_o)**2, &
                  zero,        &
                  omega_ab_ij, &
                  (wf%n_v)**2)
!
      call wf%mem%dealloc(t_cd_ij, (wf%n_v)**2, (wf%n_o)**2)
      call wf%mem%dealloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Reorder into omega2
!
      call wf%mem%alloc(omega_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_1324(omega_ab_ij, omega_ai_bj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%mem%dealloc(omega_ab_ij, (wf%n_v)**2, (wf%n_o)**2) 
!
      call add_to_packed(wf%omega2, omega_ai_bj, (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(omega_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   module subroutine omega_ccsd_c2_ccsd(wf)
!!
!!    Omega C2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!    
!!    Omega C2 = -1/2 * sum_(ck) t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!!                    - sum_(ck) t_bk_ci*(g_kj_ac - sum_(dl)t_al_dj * g_kd_lc)
!!                    - 1/2 * sum_ck u_jk^bc g_acki
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kd_lc
      real(dp), dimension(:,:), allocatable :: g_dl_ck
      real(dp), dimension(:,:), allocatable :: g_ki_ac
      real(dp), dimension(:,:), allocatable :: g_ai_ck
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_al_di
      real(dp), dimension(:,:), allocatable :: t_ai_dl
      real(dp), dimension(:,:), allocatable :: t_ck_bj
      real(dp), dimension(:,:), allocatable :: t_bk_cj
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:), allocatable :: X_ai_ck
      real(dp), dimension(:,:), allocatable :: Y_ai_bj
!
!     Reordered U2 amplitudes 
!
      real(dp), dimension(:,:), allocatable :: u_ck_bj
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! Holds term temporarily
!  
!     Indices
!     
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, cj = 0, ck = 0, cl = 0, di = 0, dk = 0, dl = 0
      integer(i15) :: kd = 0, lc = 0, ca = 0, ac = 0
      integer(i15) :: ki = 0, ai_offset = 0
!
      integer(i15) :: aldi = 0, aibj = 0, cldk = 0, bkcj = 0, bjck = 0
!
!     Batching and memory handling
!
      integer(i15) :: required = 0
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a 
!
!     Sort t_al_di = t_li^ad as t_ai_dl (1234 to 1432)
!
      call wf%mem%alloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call squareup_and_sort_1234_to_1432(wf%t2am, t_ai_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Allocate and construct g_kd_lc 
!
      call wf%mem%alloc(g_kd_lc,(wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_kd_lc)
!
!     Sort g_kd_lc to g_dl_ck (1234 to 2341)
!
      call wf%mem%alloc(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call sort_1234_to_2341(g_kd_lc, g_dl_ck, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%mem%dealloc(g_kd_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     -1/2*sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call wf%mem%alloc(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -half,             &
                  t_ai_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  g_dl_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ia_J and g_dl_ck, 
!
      call wf%mem%dealloc(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%mem%dealloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form u_ck_bj = u_jk^bc = u_kj^cb = 2 * t_jk^bc - t_kj^bc 
!
      call wf%mem%alloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2am, t_ck_bj, (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(u_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call add_1432_to_1234(-one, t_ck_bj, zero, u_ck_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ck_bj, 1, u_ck_bj, 1)
!
      call wf%mem%dealloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate a holder for - 1/2 * sum_ck u_jk^bc g_acki,
!     constructed in batches over the a index below
!
      call wf%mem%alloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      omega2_ai_bj = zero
!
!     Constructing g_ki_ac
!
      required = wf%get_vvoo_required_mem()
!
!     Initialize batching variable 
!
      call batch_a%init(wf%n_v)
      call wf%mem%num_batch(batch_a, required)
!
!     Start looping over a-batches 
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Determine batch limits for the a-batch 
!
         call batch_a%determine_limits(current_a_batch)
!
!        Allocate and construct g_ki_ac
!
         call wf%mem%alloc(g_ki_ac, (wf%n_o)**2, (batch_a%length)*(wf%n_v))
!
         integral_type = 'electronic_repulsion'
         call wf%get_oo_vv(integral_type, & 
                           g_ki_ac,       &
                           1,             & 
                           wf%n_o,        &
                           1,             & 
                           wf%n_o,        &
                           batch_a%first, & 
                           batch_a%last,  &
                           1,             &
                           wf%n_v)
!
!        X_ai_ck = X_ai_ck + g_ki_ac
!  
         do i = 1, wf%n_o
            do k = 1, wf%n_o
!
               ki = index_two(k, i, wf%n_o)
!
               do c = 1, wf%n_v
!
                  ck = index_two(c, k, wf%n_v)
!
                  do a = 1, batch_a%length
!
                     ai = index_two(a + batch_a%first -1, i, wf%n_v)
                     ac = index_two(a, c, batch_a%length)
!
                     X_ai_ck(ai, ck) = X_ai_ck(ai, ck) + g_ki_ac(ki, ac)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Calculate the contribution to the term
!
!           omega_ai_bj = - 1/2 * sum_ck u_jk^bc g_acki
!
!        Reorder g_ki_ac to g_ai_ck (1234 to 3241)
!
         call wf%mem%alloc(g_ai_ck, (batch_a%length)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
         call sort_1234_to_3241(g_ki_ac, g_ai_ck, wf%n_o, wf%n_o, batch_a%length, wf%n_v)
!
         call wf%mem%dealloc(g_ki_ac, (wf%n_o)**2, (batch_a%length)*(wf%n_v))
!
!        - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj
!
         ai_offset = index_two(batch_a%first, 1, wf%n_v)
!
         call dgemm('N','N',                 &
                  (wf%n_o)*(batch_a%length), &
                  (wf%n_o)*(wf%n_v),         &
                  (wf%n_o)*(wf%n_v),         &
                  -one/two,                  &
                  g_ai_ck,                   &
                  (wf%n_o)*(batch_a%length), &
                  u_ck_bj,                   &
                  (wf%n_o)*(wf%n_v),         &
                  zero,                      &
                  omega2_ai_bj(ai_offset,1), &
                  (wf%n_o)*(wf%n_v))
!
         call wf%mem%dealloc(g_ai_ck, (batch_a%length)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
      enddo ! End of batching
!
!     Deallocate reordered u_ck_bj vector 
!
      call wf%mem%dealloc(u_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the - 1/2 * sum_ck u_jk^bc g_acki term to omega  
!
!     Omega2(aibj, 1) =+ omega2_ai_bj(ai,bj) + omega2_bj_ai(bj,ai)
!
      call symmetrize_and_add_to_packed(wf%omega2, omega2_ai_bj, (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
!
!     Reorder t_bk_cj = t_kj^bc as t_ck_bj
!
      call wf%mem%alloc(t_bk_cj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2am, t_bk_cj, (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_3214(t_bk_cj, t_ck_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%mem%dealloc(t_bk_cj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate intermediate Y_ai_bj
!
      call wf%mem%alloc(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Y_ai_bj = - sum_(ck) X_ai_ck*t_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  t_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ai_bj,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the X intermediate
!
      call wf%mem%dealloc(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Deallocate t_ck_bj
!
      call wf%mem%dealloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Omega_aibj,1 = P_ai_bj ( 1/2*Y_ai_bj + Y_aj_bi )
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
                  if (ai .ge. bj) then
!
                     aj = index_two(a, j, wf%n_v)
                     bi = index_two(b, i, wf%n_v)
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + half*Y_ai_bj(ai, bj) + Y_ai_bj(aj, bi) &
                                                               + half*Y_ai_bj(bj, ai) + Y_ai_bj(bi, aj)
!
                  endif
!  
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate intermediate Y_ai_bj
!
      call wf%mem%dealloc(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_c2_ccsd
!
!
   module subroutine omega_ccsd_d2_ccsd(wf)
!!
!!    Omega D2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Calculates the D2 term,
!!
!!      D2: sum_ck u_jk^bc g_aikc 
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!    where 
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!    The first and second terms are referred to as D2.1 and D2.2. 
!!    All terms are added to the omega vector of the wavefunction object wf.
!!
!!    The routine adds the terms in the following order: D2.2, D2.1
!!
      implicit none 
!
      class(ccsd) :: wf 
!
!     Indices 
!
      integer(i15) :: ai = 0, aidl = 0, al = 0, aldi = 0, a = 0, i = 0, b = 0, ca = 0, ac = 0
      integer(i15) :: j = 0, c = 0, d = 0, di = 0, dl = 0, k = 0, kc = 0, kd = 0, l = 0, ki = 0
      integer(i15) :: lc = 0, ld = 0, aibj = 0, bj = 0, bjck = 0, bk = 0, bkcj = 0, cj = 0, ck = 0
!
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing D2.3, D2.2 & D2.1
!
!     Vectors for D2.2 term 
!
      real(dp), dimension(:,:), allocatable :: g_ld_kc ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: L_ld_kc ! L_ldkc = 2 * g_ldkc - g_lckd 
      real(dp), dimension(:,:), allocatable :: t_ai_dl ! t_il^ad
      real(dp), dimension(:,:), allocatable :: u_ai_dl ! u_il^ad = 2 * t_il^ad - t_li^ad
      real(dp), dimension(:,:), allocatable :: u_ai_ld ! u_il^ad = 2 * t_il^ad - t_li^ad 
      real(dp), dimension(:,:), allocatable :: Z_ai_kc ! An intermediate, see below
!
      real(dp), dimension(:,:), allocatable :: g_ai_kc ! g_aikc 
      real(dp), dimension(:,:), allocatable :: u_kc_bj ! u_jk^bc
!
!     Vectors for D2.1 term 
!
      real(dp), dimension(:,:), allocatable :: g_ai_ck ! g_acki
      real(dp), dimension(:,:), allocatable :: g_ac_ki ! g_acki; a is batched over 
      real(dp), dimension(:,:), allocatable :: u_ck_bj ! u_jk^bc
!
!     :: Calculate the D2.2 term of omega ::
!
!     Get g_ld_kc = g_ldkc 
!
      call wf%mem%alloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_ld_kc)
!
!     Form L_ld_kc = L_ldkc = 2*g_ld_kc(ld,kc) - g_ld_kc(lc,kd)    
!
      call wf%mem%alloc(L_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call add_1432_to_1234(-one, g_ld_kc, zero, L_ld_kc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, g_ld_kc, 1, L_ld_kc, 1)
!
!     Deallocate g_ld_kc and L_kc_J
!
      call wf%mem%dealloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
! 
!     Form u_ai_ld = u_il^ad = 2 * t_il^ad - t_li^ad 
!
      call wf%mem%alloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2am, t_ai_dl, (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(u_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call add_1432_to_1234(-one, t_ai_dl, zero, u_ai_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ai_dl, 1, u_ai_dl, 1)
!
      call wf%mem%dealloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%mem%alloc(u_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_1243(u_ai_dl, u_ai_ld, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%mem%dealloc(u_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc and set it to zero
!
      call wf%mem%alloc(Z_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  L_ld_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Z_ai_kc,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ld_kc
!
      call wf%mem%dealloc(L_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the D2.2 term, 1/4 sum_kc Z_ai_kc u_kc_bj = 1/4 sum_kc Z_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call wf%mem%alloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one/four,          &
                  Z_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), & 
                  zero,              &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
!     Some mathematical justification for the above matrix multiplication. We have 
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_ai_ld(ai,ld) = u_il^ad, 
!     which means that u_ai_ld(bj,kc)^T = u_ai_ld(kc,bj) = u_kj^cb = u_jk^bc.
!
!     Deallocate the Z_ai_kc intermediate 
!
      call wf%mem%dealloc(Z_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Calculate the D2.1 term of omega :: 
!
!     Form g_ai_kc = g_aikc
!
      call wf%mem%alloc(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_ov(integral_type,g_ai_kc)
!
!     Calculate the D2.1 term, sum_ck u_jk^bc g_aikc = sum_ck g_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(u_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
      call wf%mem%dealloc(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the D2.1 term to the omega vector 
!
      call symmetrize_and_add_to_packed(wf%omega2, omega2_ai_bj, (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_d2_ccsd
!
! 
   module subroutine omega_ccsd_e2_ccsd(wf)
!
!     Omega E2
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Calculates the E2 term,
!
!      E2: sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) 
!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!
!     where
!
!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!
!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!
!     Both are permuted added to the projection vector element omega2(ai,bj) of
!     the wavefunction object wf.
!
      implicit none 
!
      class(ccsd) :: wf 
!
!     Indices 
!
      integer(i15) :: aib = 0, aibk = 0, bk = 0, bja = 0, ibj = 0, aibj = 0, dlck = 0
      integer(i15) :: b = 0, c = 0, k = 0, d = 0, ck = 0, ckdl = 0, cl = 0, cldk = 0
      integer(i15) :: dk = 0, dl = 0, kc = 0, ldk = 0, l = 0, ld = 0, a = 0, ai = 0
      integer(i15) :: bj = 0, aicj = 0, cj = 0, i = 0, j = 0, jai = 0, cld = 0, dkcl = 0
!
!     Vectors for E2.1 term 
!
      real(dp), dimension(:,:), allocatable :: omega2_bj_ai ! For storing the E2.1 term temporarily
      real(dp), dimension(:,:), allocatable :: g_ld_kc      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: u_bl_dk      ! u_kl^bd 
      real(dp), dimension(:,:), allocatable :: t_bk_dl      ! t_kl^bd 
      real(dp), dimension(:,:), allocatable :: u_bk_dl      ! u_kl^bd 
      real(dp), dimension(:,:), allocatable :: X_b_c        ! An intermediate, see below for definition
      real(dp), dimension(:,:), allocatable :: t_cj_ai      ! t_ij^ac
!
!     Vectors for E2.2 term 
! 
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing the E2.2 term temporarily
      real(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition 
!
!     :: Calculate the E2.1 term of omega ::
!
!     Form u_bl_dk = u_kl^bd = u_bk_dl 
!
!     = 2 * t_kl^bd - t_lk^bd = 2 * t_bk_dl(bk,dl) - t_bk_dl(bl, dk)
! 
      call wf%mem%alloc(t_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call squareup(wf%t2am, t_bk_dl, (wf%n_v)*(wf%n_o))
!
      call wf%mem%alloc(u_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call add_1432_to_1234(-one, t_bk_dl, zero, u_bk_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_bk_dl, 1, u_bk_dl, 1)
!
      call wf%mem%alloc(u_bl_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_1432(u_bk_dl, u_bl_dk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%mem%dealloc(u_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form g_ld_kc = g_ldkc 
!
      call wf%mem%alloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_ld_kc)
!
!     Allocate the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call wf%mem%alloc(X_b_c, wf%n_v, wf%n_v)
!
!     Copy the virtual-virtual Fock matrix into the intermediate 
!
      call dcopy((wf%n_v)**2, wf%fock_ab, 1, X_b_c, 1) ! X_b_c = F_bc 
!
!     Add the second contribution, 
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  u_bl_dk,              & ! u_b_ldk 
                  wf%n_v,               &
                  g_ld_kc,              & ! g_ldk_c
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_b_c,                &
                  wf%n_v)
!
!     Form t_cj_ai = t_ij^ac = t_ji^ca 
!
      call wf%mem%alloc(t_cj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call squareup(wf%t2am, t_cj_ai, (wf%n_v)*(wf%n_o))
!
!     Form the E2.1 term 
!
      call wf%mem%alloc(omega2_bj_ai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  X_b_c,                &
                  wf%n_v,               &
                  t_cj_ai,              & ! t_c_jai
                  wf%n_v,               &
                  zero,                 &
                  omega2_bj_ai,         & ! omega2_b_jai 
                  wf%n_v)
!
!     Add the E2.1 term to the omega vector 
!
      call symmetrize_and_add_to_packed(wf%omega2, omega2_bj_ai, (wf%n_o)*(wf%n_v))
!
!     Deallocate the E2.1 term, the X intermediate, keep the amplitudes 
!
      call wf%mem%dealloc(omega2_bj_ai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
      call wf%mem%dealloc(X_b_c, wf%n_v, wf%n_v)
!
!     :: Calculate E2.2 term of omega ::
!
!     Allocate the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc 
!                                     = F_kj  + sum_cdl u_lj^dc g_kcld 
!                                     = F_k_j + sum_cdl g_k_cld * u_cld_j
!
      call wf%mem%alloc(Y_k_j, wf%n_o, wf%n_o)
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj 
!
      call dcopy((wf%n_o)**2, wf%fock_ij, 1, Y_k_j, 1)
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that 
!     Y_k_j = F_k_j + sum_cdl g_k_cld u_cld_j
!
!     Note:
!
!     g_ld_kc(kc,ld) = g_kcld              -> pretend that this is g_k_cld 
!     u_b_ldk(c,ldj) = u_jl^cd (= u_lj^dc) -> pretend that this is u_cld_j
!
      call dgemm('N','N',              &
                 wf%n_o,               &
                 wf%n_o,               &
                 (wf%n_o)*(wf%n_v)**2, &
                 one,                  &
                 g_ld_kc,              & ! g_k_cld
                 wf%n_o,               &
                 u_bl_dk,              & ! u_cld_j
                 (wf%n_o)*(wf%n_v)**2, &
                 one,                  &
                 Y_k_j,                &
                 wf%n_o)
!
!     Deallocate u_b_ldk and g_ld_kc
!
      call wf%mem%dealloc(u_bl_dk, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
      call wf%mem%dealloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!    Allocate the E2.2 term and set to zero 
!
     call wf%mem%alloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!    Calculate the E2.2 term, 
!    - sum_k t_aib_k Y_k_j = - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
! 
!    Note: t_cj_ai = t_ji^ca => t_cj_ai(ai,bk) = t_ik^ab;
!    thus, we can treat t_cj_ai as t_aib_k = t_ik^ab.
!
      call dgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one,                 &
                 t_cj_ai,              & ! t_aib_k 
                 (wf%n_o)*(wf%n_v)**2, &
                 Y_k_j,                &
                 wf%n_o,               &
                 zero,                 &
                 omega2_ai_bj,         & ! omega2_aib_j 
                 (wf%n_o)*(wf%n_v)**2)
!
!     Deallocate Y_k_j and the amplitudes 
!
      call wf%mem%dealloc(Y_k_j, wf%n_o, wf%n_o)
      call wf%mem%dealloc(t_cj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Add the E2.2 term to the omega vector 
!
      call symmetrize_and_add_to_packed(wf%omega2, omega2_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Deallocate the E2.2 term 
!
      call wf%mem%dealloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_e2_ccsd
!
! 
end submodule omega
