submodule (cc2_class) omega
!
!!
!!    Omega submodule (CC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the CC2 class:
!!
!!    initialize_omega: allocates the projection vector omega1
!!                      and sets it to zero.
!!
!!    construct_omega: constructs the projection vector omega1
!!                     for the current amplitudes t1am for the
!!                     wavefunction object wf.
!!                     The routine assumes that the projection
!!                     vector is allocated.
!!
!!    omega_a1:  adds A1 term to omega1
!!    omega_b1:  adds B1 term to omega1
!!    omega_c1:  adds C1 term to omega1
!!
!
   implicit none
!
   logical :: debug = .false.
   logical :: timings = .false.
!
   character(len=40) :: integral_type
!
!
contains
!
    module subroutine construct_omega_cc2(wf)
!
!     Construct Omega (CC2)
!     Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!
!     The routine also sets up timing variables.
!
      implicit none
!
      class(cc2) :: wf
!
!     Timing variables
!
      real(dp) :: omega_start = zero
      real(dp) :: omega_end   = zero
!
!     Start timing of omega
!
      call cpu_time(omega_start)
!
!     Set the omega vector to zero
!
      wf%omega1 = zero
!
!     :: Calculate CCS omega contributions ::
!
      call wf%omega_ccs_a1
!
!     :: Calculate CC2 omega contributions ::
!
      call wf%omega_cc2_a1
!
      call wf%omega_cc2_b1
!
!     Timings
!
      call cpu_time(omega_end)
      if (timings) write(unit_output,*)'Time in omega:', omega_end-omega_start
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf)
!!
!!     Omega A1
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!     Calculates the A1 term of omega for the active space,
!!
!!     A1: sum_bcj g_ab,jc * u_ij^bc,
!!
!!     and adds it to the projection vector (omega1) of
!!     the wavefunction object wf
!!
!!     u_ij^bc = 2*s_ij^bc - s_ij^cb
!!
!!    Batching over A and c
!!
!!
      implicit none
!
      class(cc2)   :: wf
!
!     Batching variables
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: A_n_batch = 0, A_first = 0, A_last = 0, A_batch = 0, A_length = 0
      integer(i15) :: c_n_batch = 0, c_first = 0, c_last = 0, c_batch = 0, c_length = 0
!
!     Indices
!
      integer(i15) :: i = 0, j = 0, a = 0, c = 0, b = 0
!
      integer(i15) :: ba = 0, ab = 0, bi = 0
      integer(i15) :: ic = 0, ib = 0, jc = 0, jb = 0
      integer(i15) :: bjc = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: g_ib_jc
      real(dp), dimension(:,:), allocatable :: s_ib_jc
      real(dp), dimension(:,:), allocatable :: u_bjc_i
      real(dp), dimension(:,:), allocatable :: g_Ab_jc
      real(dp), dimension(:,:), allocatable :: L_Ab_J  ! L_Ab^J; A is batched over
!
      logical :: reorder  ! To get L_ab_J reordered, for batching over a
!
!     Prepare for batching ocer c and A
!
      required = wf%n_v*(wf%n_v)*(wf%n_J) &
               + wf%n_v*(wf%n_o)*(wf%n_J) &
               + (wf%n_v**2)*(wf%n_o)*wf%n_o &
               + ((wf%n_o)**2)*((wf%n_v)**2)
!
!
      required = required*4  ! Words

      max_length = 0
!
      call num_two_batch(required, wf%mem%available, max_length, c_n_batch, wf%n_v)
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
         call batch_limits(c_first ,c_last ,c_batch, max_length, wf%n_v)
!
!        Length of batch
!
         c_length = c_last - c_first + 1
!
!        :: Construct u_ib_jc ::
!
!        u_ij^bc = 2*s_ij^bc - s_ij^cb =  (2*g_ij^bc - g_ij^cb)/ε_ij^cb
!
         call wf%mem%alloc(s_ib_jc, (wf%n_o)*(wf%n_v), (wf%n_o)*c_length )
         call wf%get_s2am(s_ib_jc, c_first, c_length)
!
         call wf%mem%alloc(u_bjc_i, wf%n_v*wf%n_o*c_length, wf%n_o)
!
         do b = 1, wf%n_v
            do i = 1, wf%n_o
!
               ib = index_two(i, b, wf%n_o)
!
               do c = 1, c_length
                  do j = 1, wf%n_o
!
                     jc = index_two(j, c, wf%n_o)
                     jb = index_two(j, b, wf%n_o)
                     ic = index_two(i, c, wf%n_o)
!
                     bjc = index_three(b, j, c, wf%n_v, wf%n_o)
!
                     u_bjc_i(bjc,i) = (two*s_ib_jc(ib,jc)-s_ib_jc(jb, ic))
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(s_ib_jc, (wf%n_o)*(wf%n_v), (wf%n_o)*c_length)
!
!        Prepare for batching over A
!
         A_first  = 0
         A_last   = 0
         A_length = 0
         A_n_batch = c_n_batch
!
!        Start looping over a-batches
!
         do A_batch = 1, A_n_batch
!
            call batch_limits(A_first ,A_last ,A_batch, max_length, wf%n_v)
            A_length = A_last - A_first + 1
!
!           :: Construct integral g_ab,jc ::
!
            call wf%mem%alloc(g_ab_jc, wf%n_v*a_length, wf%n_o*c_length)
!
            integral_type = 'electronic_repulsion'
            call wf%get_vv_ov(integral_type, g_ab_jc, a_first, a_last, 1, wf%n_v, 1, wf%n_o, c_first, c_last)
!
!
!           :: Add contributions to omega ::
!
            call dgemm('N', 'N',                      &
                        A_length,                     &
                        wf%n_o,                       &
                        (wf%n_v)*(wf%n_o)*c_length,   &
                        one,                          &
                        g_ab_jc,                      &
                        A_length,                     &
                        u_bjc_i,                      &
                        (wf%n_v)*(wf%n_o)*c_length,   &
                        one,                          &
                        wf%omega1(A_first,1),         &
                        wf%n_v)
!
            call wf%mem%dealloc(g_ab_jc, (wf%n_v)*a_length, wf%n_o*c_length)
!
         enddo ! Batching over a
!
         call wf%mem%dealloc(u_bjc_i, (wf%n_v)*(wf%n_o)*c_length, wf%n_o)
!
      enddo ! Batching over c

!
   end subroutine omega_cc2_a1_cc2
!
!
    module subroutine omega_cc2_b1_cc2(wf)
!!
!!    Omega B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the B1 term of omega,
!!
!!    B1: - sum_bjk u_jk^ab*g_kbjI + sum_bj u_ij^ab F_jb,
!!
!!    with u_ij^ab = 2*s_ij^ab - s_ij^ba.
!!
!!    Batching over b.
!!
      implicit none
!
      class(cc2)   :: wf
!
!     Batching
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: b_n_batch = 0, b_first = 0, b_last = 0, b_batch = 0, b_length = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: s_ja_kb
      real(dp), dimension(:,:), allocatable :: u_a_kbj
      real(dp), dimension(:,:), allocatable :: g_kb_ji
!
!     looping indices
!
      integer(i15) :: i = 0, j = 0, k = 0, a  = 0, b = 0
      integer(i15) :: I_full = 0, J_full = 0, A_full  = 0, B_full = 0
      integer(i15) :: ja = 0, kb = 0, ka  = 0, jb = 0, aj = 0
      integer(i15) :: kbj = 0, jbi = 0
!
!     Prepare for batching
!
      required = ((wf%n_o)**3)*(wf%n_v) + ((wf%n_o)**2)*(wf%n_J) &
               + (wf%n_o)*(wf%n_v)*(wf%n_J) + ((wf%n_o)**2)*((wf%n_v)**2)
!
!
      required = required*4  ! Words

      max_length = 0
      call num_batch(required, wf%mem%available, max_length, b_n_batch, wf%n_v)
!
!     Initialize some variables for batching
!
      b_first  = 0
      b_last   = 0
      b_length = 0
!
!     Start looping over a-batches
!
      do b_batch = 1, b_n_batch
!
         call batch_limits(b_first ,b_last ,b_batch, max_length, wf%n_v)
!
!        b is active index, and thus b_first and b_last must be displaced
!
         b_length = b_last - b_first + 1
!
!        :: Construct u_jk^ab ::
!
!        u_jk^ab = 2*s_jk^ab - s_jk^ba  (place in u_a_jkb)
!
         call wf%mem%alloc(s_ja_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
         call wf%get_s2am(s_ja_kb, b_first, b_length)
!
         call wf%mem%alloc(u_a_kbj, wf%n_v, (wf%n_o**2)*b_length)

         do k = 1, wf%n_o
            do b = 1, b_length
!
                kb = index_two(k, b, wf%n_o)
!
               do j = 1, wf%n_o
!
                  jb = index_two(j, b, wf%n_o)
                  kbj = index_three(k, b, j, wf%n_o, b_length)
!
                  do a = 1, wf%n_v
!
                     ja = index_two(j, a, wf%n_o)
                     ka = index_two(k, a, wf%n_o)

!
                     u_a_kbj(a,kbj) = (two*s_ja_kb(ja,kb)-s_ja_kb(ka, jb))
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(s_ja_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(b_length))
!
!        :: - sum_bjk u_ja_kb * g_kb_jI ::
!
         call wf%mem%alloc(g_kb_ji, (wf%n_o)*b_length, (wf%n_o)**2 )
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_oo(integral_type, g_kb_ji, 1, wf%n_o, b_first, b_last, 1, wf%n_o, 1, wf%n_o)
!
!        Add contributions to omega
!
         call dgemm('N', 'N',                &
                     wf%n_v,                 &
                     (wf%n_o),               &
                     b_length*((wf%n_o)**2), &
                     -one,                   &
                     u_a_kbj,                &
                     wf%n_v,                 &
                     g_kb_ji,                &
                     b_length*((wf%n_o)**2), &
                     one,                    &
                     wf%omega1,              &
                     (wf%n_v))
!
         call wf%mem%dealloc(g_kb_ji, wf%n_o*b_length, wf%n_o*(wf%n_o))
!
!        :: sum_jb F_jb u_ij^ab ::
!
         do i = 1, wf%n_o
!
            do a = 1, wf%n_v
!
               do j = 1, wf%n_o
!
                  do b = 1, b_length
!
                     B_full = b + b_first - 1
!
                     jbi = index_three(j, b, i, wf%n_o, b_length)

                     wf%omega1(a, i) = wf%omega1(a, i) + u_a_kbj(a, jbi)*wf%fock_ia(j, B_full)
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(u_a_kbj, wf%n_v, (wf%n_o**2)*b_length)
!
      enddo
!
   end subroutine omega_cc2_b1_cc2
!
   module subroutine get_s2am_cc2(wf, s_ia_jb, b_first, b_length)
!!
!!    Get S_2 amplitudes,
!!    Written by Sarai D. Folkestad, July 2017
!!
!!    Construct
!!
!!       s_ai_bj = - 1/ε_ij^ab * g_aibj,
!!
!!    while batching over b.
!!
      implicit none
!
      class(cc2) :: wf
!
      integer(i15) :: b_first, b_length
      real(dp), dimension((wf%n_v)*(wf%n_o), b_length*(wf%n_o)) :: s_ia_jb
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, ia = 0, jb = 0, ai = 0, bj = 0
!
!
      call wf%mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_vo(integral_type, g_ai_bj, 1, wf%n_v, 1, wf%n_o, b_first, b_first + b_length - 1, 1, wf%n_o)
!
         do a = 1, wf%n_v
            do i = 1, wf%n_o
!
               ia = index_two(i, a, wf%n_o)
               ai = index_two(a, i, wf%n_v)
!
               do b = 1, b_length
                  do j = 1, wf%n_o
!
                     jb = index_two(j, b, wf%n_o)
                     bj = index_two(b, j, wf%n_v)
!
                     s_ia_jb(ia, jb) = g_ai_bj(ai, bj)/(wf%fock_diagonal(i,1)&
                                           + wf%fock_diagonal(j,1) &
                                           - wf%fock_diagonal(wf%n_o + b + b_first - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + a ,1))
!
                  enddo
               enddo
            enddo
         enddo
!
      call wf%mem%dealloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
!
   end subroutine get_s2am_cc2
!
!
end submodule
