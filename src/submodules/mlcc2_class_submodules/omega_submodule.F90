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
      integer(i15) :: ic = 0, ib = 0, jc = 0, jb = 0 
      integer(i15) :: bjc = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: L_ib_J 
      real(dp), dimension(:,:), allocatable :: L_bi_J 
      real(dp), dimension(:,:), allocatable :: L_jc_J
      real(dp), dimension(:,:), allocatable :: g_ib_jc 
      real(dp), dimension(:,:), allocatable :: s_ib_jc 
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

      available = get_available()

      max_length = 0
      call num_two_batch(required, available, max_length, c_n_batch, wf%n_v)
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
         call allocator(s_ib_jc, (n_active_o)*(n_active_v), (n_active_o)*c_length )
         call wf%get_s2am(s_ib_jc, c_first, c_length)
!
         call allocator(u_bjc_i, n_active_v*n_active_o*c_length, n_active_o)
!
         do b = 1, n_active_v
            do i = 1, n_active_o
!
               ib = index_two(i, b, n_active_o)
!
               do c = 1, c_length
                  do j = 1, n_active_o
!
                     jc = index_two(j, c, n_active_o)
                     jb = index_two(j, b, n_active_o)
                     ic = index_two(i, c, n_active_o)
!
                     bjc = index_three(b, j, c, n_active_v, n_active_o)
!
                     u_bjc_i(bjc,i) = (two*s_ib_jc(ib,jc)-s_ib_jc(ic, jb))
!
                  enddo
               enddo
            enddo
         enddo
!         
         call deallocator(s_ib_jc, (n_active_o)*(n_active_v), (n_active_o)*c_length)
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
!           Get cholesky vectors
!
            call allocator(L_jc_J, n_active_o*c_length, wf%n_J)
            L_jc_J = zero
            call wf%get_cholesky_ia(L_jc_J, first_active_o, last_active_o, c_first, c_last)
!
            call allocator(L_Ab_J, (n_active_v)*a_length, wf%n_J) 
            L_Ab_J = zero
!
            call wf%get_cholesky_ab(L_Ab_J, a_first, a_last, first_active_v, last_active_v)
!
            call allocator(g_Ab_jc, n_active_v*a_length, n_active_o*c_length)      
!
!           g_Ab,jc = sum_J L_Ab^J*L_jc^J
!
            call dgemm('N', 'T',                  &
                        (n_active_v)*a_length,    &
                        (n_active_o)*c_length,    &
                        (wf%n_J),                 &
                        one,                      &
                        L_Ab_J,                   &
                        (n_active_v)*a_length,    &
                        L_jc_J,                   &
                        (n_active_o)*c_length,    &
                        zero,                     &
                        g_Ab_jc,                  &
                        (n_active_v)*a_length) 
!
            call deallocator(L_jc_J, (n_active_o)*(c_length), wf%n_J)
            call deallocator(L_Ab_J, (n_active_v)*(a_length), wf%n_J)
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
            call deallocator(g_Ab_jc, n_active_v*a_length, n_active_o*c_length)
!
         enddo ! Batching over a
!
         call deallocator(u_bjc_i, n_active_v*n_active_o*c_length, n_active_o)
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
      real(dp), dimension(:,:), allocatable :: L_aj_J 
      real(dp), dimension(:,:), allocatable :: L_ja_J 
      real(dp), dimension(:,:), allocatable :: L_kb_J 
      real(dp), dimension(:,:), allocatable :: L_ji_J 
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

      available = get_available()

      max_length = 0
      call num_batch(required, available, max_length, b_n_batch, n_active_v)
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
         call allocator(s_ja_kb, (n_active_o)*n_active_v, (n_active_o)*b_length)
         call wf%get_s2am(s_ja_kb, b_first, b_length)
!
         call allocator(u_a_kbj, n_active_v, (n_active_o**2)*b_length)

         do k = 1, n_active_o
            do b = 1, b_length         
!  
                kb = index_two(k, b, n_active_o)
!  
               do j = 1, n_active_o
!
                  jb = index_two(j, b, n_active_o)
                  kbj = index_three(k, b, j, n_active_o, b_length)
!
                  do a = 1, n_active_v             
!  
                     ja = index_two(j, a, n_active_o)
                     ka = index_two(k, a, n_active_o)
                     
!  
                     u_a_kbj(a,kbj) = (two*s_ja_kb(ja,kb)-s_ja_kb(jb, ka))
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(s_ja_kb, (n_active_o)*(n_active_v), (n_active_o)*(b_length))
!
!        :: - sum_bjk u_ja_kb * g_kb_jI ::
!
!        Get cholesky vectors
! 
         call allocator(L_jI_J, n_active_o*(wf%n_o), wf%n_J)
!
         call wf%get_cholesky_ij(L_jI_J, first_active_o, last_active_o, 1, wf%n_o)

!
         call allocator(L_kb_J, n_active_o*b_length, wf%n_J)
!
         call wf%get_cholesky_ia(L_kb_J, first_active_o, last_active_o, b_first, b_last)
!
         call allocator(g_kb_jI, n_active_o*b_length, n_active_o*(wf%n_o) )
!
!        g_kb_jI = sum_J L_kb^J*L_jI_J
!
         call dgemm('N', 'T',                    &
                     (n_active_o)*b_length,      &
                     (n_active_o)*(wf%n_o),      &
                     (wf%n_J),                   &
                     one,                        &
                     L_kb_J,                     &
                     (n_active_o)*b_length,      &
                     L_ji_J,                     &
                     (n_active_o)*(wf%n_o),      &
                     zero,                       &
                     g_kb_jI,                    &
                     (n_active_o)*b_length)
!
         call deallocator(L_jI_J, n_active_o*(wf%n_o), wf%n_J)
         call deallocator(L_kb_J, n_active_o*b_length, wf%n_J)
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
         call deallocator(g_kb_jI, n_active_o*b_length, n_active_o*(wf%n_o))
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
         call deallocator(u_a_kbj, n_active_v, (n_active_o**2)*b_length)
!
      enddo 
!      
   end subroutine omega_mlcc2_b1_mlcc2
!
   module subroutine get_s2am_mlcc2(wf, s_ia_jb, b_first, b_length)
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
      real(dp), dimension((wf%n_CC2_v)*(wf%n_CC2_o), b_length*(wf%n_CC2_o)) :: s_ia_jb
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: g_ia_jb
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
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, ia = 0, jb = 0, ai = 0
!
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
      call allocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
      L_ai_J = zero
!
      call wf%get_cholesky_ai(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
      call allocator(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!     Reorder L_bi_J to L_ib_J
!
      do a = 1, n_active_v
         do i = 1, n_active_o
!
            ia = index_two(i, a, n_active_o)
            ai = index_two(a, i, n_active_v)
!
            do J = 1, wf%n_J
!
               L_ia_J(ia, J) = L_ai_J(ai, J) 
!
            enddo
         enddo
      enddo
!
      call deallocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
!
      call allocator(g_ia_jb, (n_active_o)*(n_active_v), (n_active_o)*b_length)
!
      offset = index_two(1, b_first, n_active_o)
!
!     g_ib_jc = g_bi,cj = sum_J L_bj^J*L_ci^J
!
      call dgemm('N', 'T',                    &
                  (n_active_o)*(n_active_v),  &
                  (n_active_o)*b_length,      &
                  (wf%n_J),                   &
                  one,                        &
                  L_ia_J,                     &
                  (n_active_o)*(n_active_v),  &
                  L_ia_J(offset,1),           &
                  (n_active_o)*(n_active_v),  &
                  zero,                       &
                  g_ia_jb,                    &
                  (n_active_o)*(n_active_v))
!
      call deallocator(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!
         do a = 1, n_active_v
            do i = 1, n_active_o
!
               ia = index_two(i, a, n_active_o)
!
               do b = 1, b_length
                  do j = 1, n_active_o
!
                     jb = index_two(j, b, n_active_o)
!
                     s_ia_jb(ia, jb) = g_ia_jb(ia, jb)/(wf%fock_diagonal(i + first_active_o - 1,1)&
                                           + wf%fock_diagonal(j+ first_active_o - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + b + b_first - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + a + first_active_v - 1,1))
!
                  enddo
               enddo
            enddo
         enddo
!
      call deallocator(g_ia_jb, (n_active_o)*(n_active_v), (n_active_o)*b_length)
!
   end subroutine get_s2am_mlcc2
!
end submodule