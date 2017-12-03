submodule (mlcc2_class) omega
!
!!
!!    Omega submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!
!!    construct_omega: constructs the projection vector omega1
!!                     for the current amplitudes t1am for the
!!                     wavefunction object wf. 
!!                     The routine assumes that the projection
!!                     vector is allocated.
!!
!!    omega_a1:  adds A1 term to omega1
!!    omega_b1:  adds B1 term to omega1
!!
!!
!!    Upper case indices are general indices, lower case indices are restricted
!!    to the CC2 orbital space.
!! 
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
  module subroutine construct_omega_mlcc2(wf)
!! 
!!    Construct Omega (MLCC2)
!!    Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!!
!!    Constructs the MlCC2 omega.
!! 
!!    s2-amplitudes are constructed on the fly, according to the CC2
!!    expression for the doubles amplitudes. 
!!
!!    Calculated by looping over active spaces, 
!!    Adding the omega contribution from each active space in turn.
!! 
      implicit none 
!
      class(mlcc2) :: wf
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
      call wf%omega_mlcc2_a1
!
      call wf%omega_mlcc2_b1
!
!     Timings
!
      call cpu_time(omega_end)
      if (timings) write(unit_output,*)'Time in omega:', omega_end-omega_start    
!
   end subroutine construct_omega_mlcc2
!
    module subroutine omega_mlcc2_a1_mlcc2(wf)
!! 
!!     Omega A1
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!     Calculates the A1 term of omega for the active space, 
!!   
!!     A1: sum_bcj g_Abjc * u_ij^bc,
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
      class(mlcc2)   :: wf
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
      integer(i15) :: ci = 0, bi = 0, cj = 0, bj = 0 
      integer(i15) :: bjc = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: g_ib_jc 
      real(dp), dimension(:,:), allocatable :: s_bi_cj 
      real(dp), dimension(:,:), allocatable :: u_bjc_i 
      real(dp), dimension(:,:), allocatable :: g_Ab_jc 
      real(dp), dimension(:,:), allocatable :: L_Ab_J  ! L_Ab^J; A is batched over
!
      logical :: reorder  ! To get L_ab_J reordered, for batching over a
!
!     Active space variables
!
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
      integer(i15) :: n_active_o    ! number of active occupied
      integer(i15) :: n_active_v    ! number of active virual
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_v, first_active_o)
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     Prepare for batching ocer c and A     
!
      required = n_active_v*(wf%n_v)*(wf%n_J) &
               + n_active_v*(n_active_o)*(wf%n_J) &
               + (n_active_v**2)*(wf%n_o)*n_active_o &
               + ((n_active_o)**2)*((n_active_v)**2)
!
!
      required = required*4  ! Words

      max_length = 0
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
         call batch_limits(c_first ,c_last ,c_batch, max_length, n_active_v)
!
!        c is active index, and thus c_first and c_last must be displaced
!
         c_first  = c_first + (first_active_v - 1)
         c_last   = c_last  + (first_active_v - 1)
!
         if (c_last .gt. last_active_v) c_last = last_active_v
!
!        Length of batch
!
         c_length = c_last - c_first + 1 
!
!      :: Construct u_ib_jc ::
!
!        u_ij^bc = 2*s_ij^bc - s_ij^cb =  (2*g_ij^bc - g_ij^cb)/ε_ij^cb
!
         call wf%mem%alloc(s_bi_cj, (n_active_o)*(n_active_v), (n_active_o)*c_length )
         call wf%get_s2am(s_bi_cj, c_first, c_length)
!
         call wf%mem%alloc(u_bjc_i, n_active_v*n_active_o*c_length, n_active_o)
!
         do b = 1, n_active_v
            do i = 1, n_active_o
!
               bi = index_two(b, i, n_active_v)
!
               do c = 1, c_length
                  do j = 1, n_active_o
!
                     cj = index_two(c, j, c_length)
                     bj = index_two(b, j, n_active_v)
                     ci = index_two(c, i, c_length)
!
                     bjc = index_three(b, j, c, n_active_v, n_active_o)
!
                     u_bjc_i(bjc,i) = (two*s_bi_cj(bi,cj)-s_bi_cj(bj, ci))
!
                  enddo
               enddo
            enddo
         enddo
!         
         call wf%mem%dealloc(s_bi_cj, (n_active_o)*(n_active_v), (n_active_o)*c_length)
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
!           :: Construct integral g_Ab,jc ::
!
            call wf%mem%alloc(g_Ab_jc, n_active_v*a_length, n_active_o*c_length)
! 
            integral_type = 'electronic_repulsion'
            call wf%get_vv_ov(integral_type, g_Ab_jc, &
                              A_first, A_last, &
                              first_active_v, last_active_v, &     
                              first_active_o, last_active_o, &     
                              first_active_v, last_active_v)     
!
!           :: Add contributions to omega ::
!
            call dgemm('N', 'N',                            &
                        A_length,                           &
                        n_active_o,                         &
                        (n_active_v)*(n_active_o)*c_length, &
                        one,                                &
                        g_Ab_jc,                            &
                        A_length,                           &
                        u_bjc_i,                            &
                        (n_active_v)*(n_active_o)*c_length, &
                        one,                                &
                        wf%omega1(A_first,1),               &
                        wf%n_v)
! 
            call wf%mem%dealloc(g_Ab_jc, n_active_v*a_length, n_active_o*c_length)
!
         enddo ! Batching over a
!
         call wf%mem%dealloc(u_bjc_i, n_active_v*n_active_o*c_length, n_active_o)
!
         if (c_last .eq. last_active_v) exit ! exit loop over c; This is necessary because n_active_v may be less than n_v
!
      enddo ! Batching over c

!
   end subroutine omega_mlcc2_a1_mlcc2
!
!
    module subroutine omega_mlcc2_b1_mlcc2(wf)
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
      class(mlcc2)   :: wf
!
!     Batching 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: b_n_batch = 0, b_first = 0, b_last = 0, b_batch = 0, b_length = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: s_aj_bk 
      real(dp), dimension(:,:), allocatable :: u_a_kbj 
      real(dp), dimension(:,:), allocatable :: g_kb_ji 
!
!     looping indices
!
      integer(i15) :: i = 0, j = 0, k = 0, a  = 0, b = 0
      integer(i15) :: I_full = 0, J_full = 0, A_full  = 0, B_full = 0
      integer(i15) :: aj = 0, bk = 0, ak  = 0, bj = 0
      integer(i15) :: kbj = 0, jbi = 0
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
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
!     Prepare for batching
!
      required = ((n_active_o)**3)*(n_active_v) + ((n_active_o)**2)*(wf%n_J) &
               + (n_active_o)*(n_active_v)*(wf%n_J) + ((n_active_o)**2)*((n_active_v)**2)  
!
!
      required = required*4  ! Words

      max_length = 0
      call num_batch(required, wf%mem%available, max_length, b_n_batch, n_active_v)
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
         call batch_limits(b_first ,b_last ,b_batch, max_length, n_active_v)
!
!        b is active index, and thus b_first and b_last must be displaced
!
         b_first  = b_first + (first_active_v - 1)
         b_last   = b_last  + (first_active_v - 1)
!
         b_length = b_last - b_first + 1 
!
!        :: Construct u_jk^ab ::
!  
!        u_jk^ab = 2*s_jk^ab - s_jk^ba  (place in u_a_jkb)        
!  
         call wf%mem%alloc(s_aj_bk, (n_active_o)*n_active_v, (n_active_o)*b_length)
         call wf%get_s2am(s_aj_bk, b_first, b_length)
!
         call wf%mem%alloc(u_a_kbj, n_active_v, (n_active_o**2)*b_length)

         do k = 1, n_active_o
            do b = 1, b_length         
!  
                bk = index_two(b, k, b_length)
!  
               do j = 1, n_active_o
!
                  bj = index_two(b, j, b_length)
                  kbj = index_three(k, b, j, n_active_o, b_length)
!
                  do a = 1, n_active_v             
!  
                     aj = index_two(a, j, n_active_v)
                     ak = index_two(a, k, n_active_v)
                     
!  
                     u_a_kbj(a,kbj) = (two*s_aj_bk(aj,bk)-s_aj_bk(ak, bj))
!
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(s_aj_bk, (n_active_o)*(n_active_v), (n_active_o)*(b_length))
!
!        :: - sum_bjk u_ja_kb * g_kb_jI ::
!
         call wf%mem%alloc(g_kb_jI, n_active_o*b_length, n_active_o*(wf%n_o) )
!
         integral_type = 'electronic_repulsion'
         call wf%get_ov_oo(integral_type, g_kb_jI,        &
                           first_active_o, last_active_o, &
                           b_first, b_last,               &
                           first_active_o, last_active_o, &
                           1, wf%n_o)
!
!        Add contributions to omega
!
         call dgemm('N', 'N',                        &
                     n_active_v,                     &
                     (wf%n_o),                       &
                     b_length*((n_active_o)**2),     &
                     -one,                           &
                     u_a_kbj,                        &
                     n_active_v,                     &
                     g_kb_jI,                        &
                     b_length*((n_active_o)**2),     &
                     one,                            &
                     wf%omega1(first_active_v, 1),   &
                     (wf%n_v))
!
         call wf%mem%dealloc(g_kb_jI, n_active_o*b_length, n_active_o*(wf%n_o))
!
!        :: sum_jb F_jb u_ij^ab ::
!
         do i = 1, n_active_o
!
            I_full = i + first_active_o - 1
!
            do a = 1, n_active_v
!
               A_full = a + first_active_v - 1 
!          
               do j = 1, n_active_o
!
                  J_full = j + first_active_o - 1
!
                  do b = 1, b_length
! 
                     B_full = b + b_first - 1
! 
                     jbi = index_three(j, b, i, n_active_o, b_length)
                    
                     wf%omega1(A_full, I_full) = wf%omega1(A_full, I_full) + u_a_kbj(a, jbi)*wf%fock_ia(J_full, B_full)
!
                  enddo
               enddo
            enddo
         enddo
! 
         call wf%mem%dealloc(u_a_kbj, n_active_v, (n_active_o**2)*b_length)
!
      enddo 
!      
   end subroutine omega_mlcc2_b1_mlcc2
!
   module subroutine get_s2am_mlcc2(wf, s_ai_bj, b_first, b_length)
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
      class(mlcc2) :: wf
! 
      integer(i15) :: b_first, b_length
      real(dp), dimension((wf%n_CC2_v)*(wf%n_CC2_o), b_length*(wf%n_CC2_o)) :: s_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
!
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
      integer(i15) :: offset
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, ia = 0, jb = 0, ai = 0, bj = 0
!
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
      call wf%mem%alloc(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*b_length)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_vo(integral_type, g_ai_bj,          &
                        first_active_v, last_active_v,   &
                        first_active_o, last_active_o,   &
                        b_first, b_first + b_length - 1, &
                        first_active_o, last_active_o)
!
!
         do a = 1, n_active_v
            do i = 1, n_active_o
!
               ai = index_two(a, i, n_active_v)
!
               do b = 1, b_length
                  do j = 1, n_active_o
!
                     bj = index_two(b, j, b_length)
!
                     s_ai_bj(ai, bj) = g_ai_bj(ai, bj)/(wf%fock_diagonal(i + first_active_o - 1,1)&
                                           + wf%fock_diagonal(j+ first_active_o - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + b + b_first - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + a + first_active_v - 1,1))
!
                  enddo
               enddo
            enddo
         enddo
!
      call wf%mem%dealloc(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*b_length)
!
   end subroutine get_s2am_mlcc2
!
end submodule