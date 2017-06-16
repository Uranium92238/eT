submodule (mlcc2_class) omega
!
!!
!!    Omega submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
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
!
!
contains
  module subroutine construct_omega_mlcc2(wf)
!  
!     Construct Omega (CC2)
!     Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!  
!     Constructs t2-amplitudes on the fly, according to the CC2
!     expression for the doubles amplitudes 
!  
      implicit none 
!
      class(mlcc2) :: wf
!     
!     Looping variables
!
      integer(i15) :: active_space
!
!     Timing variables
!
      real(dp) :: omega_start = zero
      real(dp) :: omega_end   = zero
      if (debug) write(unit_output,*)'In omega_constructor'
      flush(unit_output)
!
!     Start timing of omega
!
      call cpu_time(omega_start)
!
!     Set the omega vector to zero 
!
      wf%omega1 = zero
      if (debug) write(unit_output,*)'omega ccs a1'
      flush(unit_output)
      call wf%omega_ccs_a1
!
!     Loop over active spaces
!
      do active_space = 1, wf%n_active_spaces
!
!        :: Calculate omega contributions ::
!
         call wf%omega_mlcc2_a1(active_space)
         call wf%omega_mlcc2_b1(active_space)
!
      enddo
!
!       Timings
!
      call cpu_time(omega_end)
      if (debug) write(unit_output,*)'Time in omega:', omega_end-omega_start    
!
   end subroutine construct_omega_mlcc2
!
    module subroutine omega_mlcc2_a1_mlcc2(wf, active_space)
! 
!     Omega A1
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!     Calculates the A1 term of omega, 
!   
!     A1: sum_bcj g_Abjc * u_ij^bc,
!  
!     and adds it to the projection vector (omega1) of
!     the wavefunction object wf
! 
!     u_ij^bc = 2*s_ij^bc - s_ij^cb 
! 
      implicit none
!
      class(mlcc2)   :: wf
      integer(i15)   :: active_space
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: A_n_batch = 0, A_first = 0, A_last = 0, A_batch = 0, A_length = 0
      integer(i15) :: c_n_batch = 0, c_first = 0, c_last = 0, c_batch = 0, c_length = 0
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, a = 0, c = 0, b = 0
!
      integer(i15) :: ba = 0, ab = 0, bi = 0
      integer(i15) :: ic = 0, ib = 0, jc = 0, jb = 0 
      integer(i15) :: bjc = 0
!
!
      real(dp), dimension(:,:), allocatable :: L_ib_J 
      real(dp), dimension(:,:), allocatable :: L_bi_J 
      real(dp), dimension(:,:), allocatable :: L_jc_J
      real(dp), dimension(:,:), allocatable :: s_ib_jc 
      real(dp), dimension(:,:), allocatable :: u_bjc_i 
      real(dp), dimension(:,:), allocatable :: g_Ab_jc 
      real(dp), dimension(:,:), allocatable :: L_Ab_J  ! L_Ab^J; A is batched over
!
      logical :: reorder  ! To get L_ab_J reordered, for batching over a
!
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_v, first_active_o, active_space)
!
      n_active_o = wf%n_CC2_o(active_space,1) 
      n_active_v = wf%n_CC2_v(active_space,1)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     ::  Calculate the A1 term  of omega ::
!
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
      call num_two_batch(required, available, max_length, A_n_batch, wf%n_v)
!
!     Initialize some variables for batching
!
      A_first  = 0
      A_last   = 0
      A_length = 0
!
!     Start looping over a-batches
!
      do A_batch = 1, A_n_batch
!   
         call batch_limits(A_first ,A_last ,A_batch, max_length, wf%n_v)
         A_length = A_last - A_first + 1   
!
!        Start looping over batches of c
!
         c_first  = 0
         c_last   = 0
         c_length = 0
         c_n_batch = a_n_batch
!
         do c_batch = 1, c_n_batch
!
            call batch_limits(c_first ,c_last ,c_batch, max_length, n_active_v)
!
!           c is active index, and thus c_first and c_last must be displaced
!
            c_first  = c_first + (first_active_v - 1)
            c_last   = c_last  + (first_active_v - 1)
!
            if (c_last .gt. last_active_v) c_last = last_active_v
!
!           Length of batch
!
            c_length = c_last - c_first + 1 
!
            call allocator(L_bi_J, n_active_o*n_active_v, wf%n_J)
            L_bi_J = zero
!
            call wf%get_cholesky_ai(L_bi_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
            call allocator(L_ib_J, (n_active_o)*(n_active_v), wf%n_J)
!
!           reorder and constrain L_bi_J
!
            do b = 1, n_active_v
               do i = 1, n_active_o
!
                  ib = index_two(i, b, n_active_o)
                  bi = index_two(b, i, n_active_v)
!
                  do J = 1, wf%n_J
!
                     L_ib_J(ib, J) = L_bi_J(bi, J) 
!
                  enddo
               enddo
            enddo
!
            call deallocator(L_bi_J, n_active_o*n_active_v, wf%n_J)
!
            call allocator(s_ib_jc, (n_active_o)*(n_active_v), (n_active_o)*c_length)
!
            offset = index_two(1, c_first, n_active_o)
!
            call dgemm('N', 'T',                    &
                        (n_active_o)*(n_active_v),  &
                        (n_active_o)*c_length,      &
                        (wf%n_J),                   &
                        one,                        &
                        L_ib_J,                     &
                        (n_active_o)*(n_active_v),  &
                        L_ib_J(offset,1),           &
                        (n_active_o)*(n_active_v),  &
                        zero,                       &
                        s_ib_jc,                    &
                        (n_active_o)*(n_active_v))        
!
            call deallocator(L_ib_J, (n_active_o)*(n_active_v), wf%n_J)
!
!           u_ij^bc = 2*s_ij^bc - s_ij^cb  (place in s_bjc_i)        
!
            call allocator(u_bjc_i, n_active_v*n_active_o*c_length, n_active_o)
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
                        bjc = index_three(b, j, c, n_active_v, n_active_o)
!
                        u_bjc_i(bjc,i) = (two*s_ib_jc(ib,jc)-s_ib_jc(ic, jb))/(wf%fock_diagonal(i + first_active_o - 1,1)&
                                              + wf%fock_diagonal(j+ first_active_o - 1,1) &
                                              - wf%fock_diagonal(wf%n_o + c + first_active_v - 1,1) &
                                              - wf%fock_diagonal(wf%n_o + b + first_active_v - 1,1))
!
                     enddo
                  enddo
               enddo
            enddo
!            
            call deallocator(s_ib_jc, (n_active_o)*(n_active_v), (n_active_o)*c_length)
!
!           Construct integral g_Ab,jc batching over A
!
            call allocator(L_Ab_J, (n_active_v)*a_length, wf%n_J) 
            L_Ab_J = zero
!
            reorder = .false.
            call wf%get_cholesky_ab(L_Ab_J, a_first, a_last, reorder, first_active_v, last_active_v)
!
!
            call allocator(L_jc_J, n_active_o*c_length, wf%n_J)
            call wf%get_cholesky_ia(L_jc_J, first_active_o, last_active_o, c_first, c_last)

!
            call allocator(g_Ab_jc, n_active_v*a_length, n_active_o*c_length)      
!
!        
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
            call deallocator(u_bjc_i, n_active_v*n_active_o*c_length, n_active_o)
!
            if (c_batch*max_length .ge. n_active_v) exit ! exit loop over c; This is necessary because n_active_v may be less than n_v
!
         enddo ! Batching over c
      enddo ! Batching over a

!
   end subroutine omega_mlcc2_a1_mlcc2
    module subroutine omega_mlcc2_b1_mlcc2(wf, active_space)
! 
!     Omega B1
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!     Calculates the B1 term of omega, 
!   
!     B1: - sum_bjk u_jk^ab*g_kbji + sum_bj u_ij^ab F_jb
!
!     Batching over a
!
      implicit none
!
      class(mlcc2)   :: wf
      integer(i15)   :: active_space 
!
!     Batching 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: a_n_batch = 0, a_first = 0, a_last = 0, a_batch = 0, a_length = 0
!
!     Arrays
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
!     ML variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_active_v, first_active_o, active_space)
!
      n_active_o = wf%n_CC2_o(active_space,1) 
      n_active_v = wf%n_CC2_v(active_space,1)
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
      call num_batch(required, available, max_length, a_n_batch, n_active_v)
!
!     Initialize some variables for batching
!
      a_first  = 0
      a_last   = 0
      a_length = 0
!
!     Start looping over a-batches
!
      do a_batch = 1, a_n_batch
!   
         call batch_limits(a_first ,a_last ,a_batch, max_length, wf%n_v)
!
!        a is active index, and thus a_first and a_last must be displaced
!
         a_first  = a_first + (first_active_v - 1)
         a_last   = a_last  + (first_active_v - 1)
!
         if (a_last .gt. last_active_v) a_last = last_active_v
!
         a_length = a_last - a_first + 1 
!
         call allocator(L_aj_J, n_active_o*n_active_v, wf%n_J)
         L_aj_J = zero
!
         call wf%get_cholesky_ai(L_aj_J, first_active_v, last_active_v, first_active_o, last_active_o)
!

         call allocator(s_ja_kb, (n_active_o)*(a_length), (n_active_o)*(n_active_v))
!
         call allocator(L_ja_J, (n_active_o)*(n_active_v), wf%n_J)
!
!        Reorder L_aj_J to L_ja_J
!
            do a = 1, n_active_v
               do k = 1, n_active_o
!
                  ja = index_two(k, a, n_active_o)
                  aj = index_two(a, k, n_active_v)
!
                  do J = 1, wf%n_J
!
                     L_ja_J(ja, J) = L_aj_J(aj, J) 
!
                  enddo
               enddo
            enddo
!
         call deallocator(L_aj_J, n_active_o*n_active_v, wf%n_J)
!  
         offset = index_two(1, a_first, n_active_o)
         call dgemm('N', 'T',                    &
                     (n_active_o)*(a_length),    &
                     (n_active_o)*(n_active_v),  &
                     (wf%n_J),                   &
                     one,                        &
                     L_ja_J(offset,1),           &
                     (n_active_o)*(n_active_v),  &
                     L_ja_J,                     &
                     (n_active_o)*(n_active_v),  &
                     zero,                       &
                     s_ja_kb,                    &
                     (n_active_o)*(a_length))

!  
         call deallocator(L_ja_J, (n_active_o)*(n_active_v), wf%n_J)
!  
!        u_jk^ab = 2*s_jk^ab - s_jk^ba  (place in s_a_jkb)        
!  
         call allocator(u_a_kbj, a_length, (n_active_o**2)*(n_active_v))

         do k = 1, n_active_o
            do b = 1, n_active_v         
!  
                kb = index_two(k, b, n_active_o)
!  
               do j = 1, n_active_o
!
                  jb = index_two(j, b, n_active_o)
                  kbj = index_three(k, b, j, n_active_o, n_active_v)
!
                  do a = 1, a_length             
!  
                     ja = index_two(j, a, n_active_o)
                     ka = index_two(k, a, n_active_o)
                     
!  
                     u_a_kbj(a,kbj) = (two*s_ja_kb(ja,kb)-s_ja_kb(jb, ka))/(wf%fock_diagonal(j + first_active_o - 1,1)&
                                              + wf%fock_diagonal(k+ first_active_o - 1,1) &
                                              - wf%fock_diagonal(wf%n_o + a + first_active_v - 1,1) &
                                              - wf%fock_diagonal(wf%n_o + b + first_active_v - 1,1))
!
!  
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(s_ja_kb, (n_active_o)*(a_length), (n_active_o)*(n_active_v))
!
!        Construct g_kb_ji
!
         call allocator(L_jI_J, n_active_o*(wf%n_o), wf%n_J)
!
         call wf%get_cholesky_ij(L_jI_J, first_active_o, last_active_o, 1, wf%n_o)
!
         call allocator(L_kb_J, n_active_o*n_active_v, wf%n_J)
!
         call wf%get_cholesky_ia(L_kb_J, first_active_o, last_active_o, first_active_v, last_active_v)
!
         call allocator(g_kb_jI, n_active_o*n_active_v, n_active_o*(wf%n_o) )
!
         call dgemm('N', 'T',                    &
                     (n_active_o)*(n_active_v),  &
                     (n_active_o)*(wf%n_o),      &
                     (wf%n_J),                   &
                     one,                        &
                     L_kb_J,                     &
                     (n_active_o)*(n_active_v),  &
                     L_ji_J,                     &
                     (n_active_o)*(wf%n_o),      &
                     zero,                       &
                     g_kb_ji,                    &
                     (n_active_o)*(n_active_v)) 
!
         call deallocator(L_jI_J, n_active_o*(wf%n_o), wf%n_J)
         call deallocator(L_kb_J, n_active_o*n_active_v, wf%n_J)

!
         call dgemm('N', 'N',                        &
                     a_length,                       &
                     (wf%n_o),                       &
                     (n_active_v)*((n_active_o)**2), &
                     -one,                           &
                     u_a_kbj,                        &
                     a_length,                       &
                     g_kb_ji,                        &
                     (n_active_v)*((n_active_o)**2), &
                     one,                            &
                     wf%omega1(a_first, 1),          &
                     (wf%n_v))
!
         call deallocator(g_kb_jI, n_active_o*n_active_v, n_active_o*(wf%n_o))
!
!        u_a_kbj = u_a_jbi (in s_a_kbj)
!
        do i = 1, n_active_o
!
           I_full = i + first_active_o - 1
!
           do a = 1, a_length
!
              A_full = a + a_first - 1 
!          
              do j = 1, n_active_o
!
                 J_full = j + first_active_o - 1
!
                 do b = 1, n_active_v
! 
                    B_full = b + first_active_v - 1
! 
                    jbi = index_three(j, b, i, n_active_o, n_active_v)
                    
                    wf%omega1(A_full, I_full) = wf%omega1(A_full, I_full) + u_a_kbj(a, jbi)*wf%fock_ia(J_full, B_full)
!
                 enddo
              enddo
           enddo
        enddo
! 
         call deallocator(u_a_kbj, a_length, (n_active_o**2)*(n_active_v))
!
      enddo 
!      
   end subroutine omega_mlcc2_b1_mlcc2
!
end submodule