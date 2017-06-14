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
!
    module subroutine omega_mlcc2_a1_mlcc2(wf, active_space)
! 
!     Omega A1
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!     Calculates the A1 term of omega, 
!   
!     A1: sum_ckd g_adkc * u_ki^cd,
!  
!     and adds it to the projection vector (omega1) of
!     the wavefunction object wf
! 
!     u_ki^cd = 2*s_ki^cd - s_ik^cd 
! 
      implicit none
!
      class(mlcc2)   :: wf
      integer(i15)   :: active_space
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: n_batch = 0, a_first = 0, a_last = 0, a_batch = 0, a_length = 0
      integer(i15) :: c_first = 0, c_last = 0, c_batch = 0, c_length = 0
!
      integer(i15) :: n_a_o = 0, n_a_v = 0
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, k = 0, a = 0, c = 0, d = 0, b = 0
!
      integer(i15) :: ad = 0, da = 0, ba = 0, ab = 0
      integer(i15) :: ci = 0, ck = 0, di = 0, dk = 0, bi = 0
      integer(i15) :: kc = 0, ic = 0, ib = 0, jc = 0, jc_a = 0 
!
      integer(i15) :: dkc = 0 
      integer(i15) :: cidk = 0, ckdi = 0
!
      integer(i15) :: ad_dim = 0   
!
      real(dp), dimension(:,:), allocatable :: L_kc_J 
      real(dp), dimension(:,:), allocatable :: L_ib_J 
      real(dp), dimension(:,:), allocatable :: L_BI_J 
      real(dp), dimension(:,:), allocatable :: L_JC_J 
      real(dp), dimension(:,:), allocatable :: L_jc_J_a 
      real(dp), dimension(:,:), allocatable :: s_ib_jc 
      real(dp), dimension(:,:), allocatable :: g_Ab_jc 
      real(dp), dimension(:,:), allocatable :: L_BA_J  ! L_AB^J; A is batched over
      real(dp), dimension(:,:), allocatable :: L_Ab_J  ! L_Ab^J; A is batched over
!
      logical :: reorder = .true. ! To get L_ab_J reordered, for batching over a
!
 !   n_a_o = wf%n_CC2_o(active_space,1) 
 !   n_a_v = wf%n_CC2_v(active_space,1) 
!
!!   ::  Calculate the A1 term  of omega ::
!
!
 !   required = n_a_v*(wf%n_v)*(wf%n_J) &
 !            + n_a_v*(n_a_o)*(wf%n_J) &
 !            + (n_a_v**2)*(wf%n_o)*n_a_o &
 !            + ((n_a_o)**2)*((n_a_v)**2)
!
!
 !   required = required*4  ! Words

 !   available = get_available()

 !   max_length = 0
 !   call num_two_batch(required, available, max_length, n_batch, wf%n_v)
!
!!   Initialize some variables for batching
!
 !   a_first  = 0
 !   a_last   = 0
 !   a_length = 0
!
!!   Start looping over a-batches
!
 !   do a_batch = 1, n_batch
!! 
 !      call batch_limits(a_first ,a_last ,a_batch, max_length, wf%n_v)
 !      a_length = a_last - a_first + 1     
!
!!      Start looping over batches of c
!
 !      c_first  = 0
 !      c_last   = 0
 !      c_length = 0
!
 !      do c_batch = 1, n_batch
!
 !         call batch_limits(c_first ,c_last ,c_batch, max_length, n_a_v)
 !         c_length = c_last - c_first + 1 
!
 !         call allocator(L_BI_J, (wf%n_o)*(wf%n_v), wf%n_J)
 !         L_BI_J = zero
!
 !         call wf%get_cholesky_ai(L_BI_J)
!
 !         call allocator(L_ib_J, (n_a_o)*(n_a_v), wf%n_J)
!
!!         reorder and constrain L_bi_J
!
 !         do b = 1, n_a_v
 !            do i = 1, n_a_o
!
 !               ib = index_two(i, b, n_a_o)
 !               BI = index_two(b, i, wf%n_v)
!
 !               do J = 1, wf%n_J
!
 !                  L_ib_J(ib, J) = L_BI_J(BI, J) 
!
 !               enddo
 !            enddo
 !         enddo
!
 !         call deallocator(L_BI_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
 !         call allocator(s_ib_jc, (n_a_o)*(n_a_v), (n_a_o)*c_length)
!
 !         offset = index_two(1, c_first, n_a_o)
 !         call dgemm('N', 'T', &
 !                     (n_a_o)*(n_a_v),  &
 !                     (n_a_o)*c_length, &
 !                     (wf%n_J),         &
 !                     one,              &
 !                     L_ib_J,           &
 !                     (n_a_o)*(n_a_v),  &
 !                     L_ib_J(offset,1), &
 !                     (n_a_o)*(n_a_v),  &
 !                     zero,             &
 !                     s_ib_jc,          &
 !                     (n_a_o)*(n_a_v))        
!
 !         call deallocator(L_ib_J, (n_a_o)*(n_a_v), wf%n_J)
!
 !         do b = 1, n_a_v
 !            do i = 1, n_a_o
 !               ib = index_two(i, b, n_a_o)
 !               do c = 1, c_length
 !                  do j = 1, n_a_o
 !                     jc = index_two(j, c, n_a_o)
!
 !                     s_ib_jc(ib,jc) = s_ib_jc(ib,jc)/(wf%fock_diagonal(i,1)+wf%fock_diagonal(j,1) &
 !                                           - wf%fock_diagonal(wf%n_o + a,1) - wf%fock_diagonal(wf%n_o + b,1))
!
 !                  enddo
 !               enddo
 !            enddo
 !         enddo
!
!!         Construct integral g_Ab,jc batching over A
!
 !         call allocator(L_BA_J, (wf%n_v)*a_length, wf%n_J) 
 !         L_BA_J = zero
!
 !         call wf%get_cholesky_ab(L_BA_J, a_first, a_last, wf%n_v*a_length, reorder)
!
!!         restrict b index
!!  
 !         call allocator(L_Ab_J, (n_a_v)*a_length, wf%n_J) 
 !         L_Ab_J = zero
!
 !         do a = 1, a_length
 !            do b = 1, n_a_v
 !               ab   = index_two(a, b, a_length)
 !               ba   = index_two(b, a, wf%n_v)
 !               do J = 1, wf%n_J
 !                  L_Ab_J(ab, J) = L_BA_J(ba, J)
 !               enddo
 !            enddo
 !         enddo
!
 !         call deallocator(L_BA_J, (wf%n_v)*a_length, wf%n_J)  
!
 !         call allocator(L_JC_J, (wf%n_v)*(wf%n_o), wf%n_J)
 !         call wf%get_cholesky_ia(L_JC_J)
!
 !         call allocator(L_jc_J_a, c_length*(n_a_o), wf%n_J)
!
 !         do c = 1, c_length
 !            do k = 1, n_a_o
!
 !               jc_a = index_two(k, c, n_a_o)
 !               jc   = index_two(k, c + c_first - 1, wf%n_o)
!
 !               do J = 1, wf%n_J
 !                  L_jc_J_a(jc_a, J) = L_JC_J(jc, J)
 !               enddo
 !            enddo
 !         enddo
 !         call deallocator(L_JC_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
 !         call allocator(g_Ab_jc, n_a_v*a_length, n_a_o*c_length)          
!
 !         offset = index_two(1, c_first, n_a_o)
 !         call dgemm('N', 'T', &
 !                     (n_a_v)*a_length,    &
 !                     (n_a_o)*c_length,    &
 !                     (wf%n_J),            &
 !                     one,                 &
 !                     L_Ab_J,              &
 !                     (n_a_v)*a_length,    &
 !                     L_jc_J_a(offset,1),  &
 !                     (n_a_o)*c_length,    &
 !                     zero,                &
 !                     g_Ab_jc,             &
 !                     (n_a_v)*a_length) 
!
 !         call deallocator(L_Ab_J, (n_a_v)*a_length, wf%n_J)
!
 !         call dgemm('N', 'T',                  &
 !                     a_length,                 &
 !                     n_a_o,                    &
 !                     (n_a_v)*(n_a_o)*c_length, &
 !                     one,                      &
 !                     g_Ab_jc,                  &
 !                     a_length,                 &
 !                     s_ib_jc,                  &
 !                     n_a_o,                    &
 !                     one,                      &
 !                     wf%omega1(a_first,1),     &
 !                     wf%n_v)
!!
 !         call deallocator(g_Ab_jc, n_a_v*a_length, n_a_o*c_length)

 !         call deallocator(s_ib_jc, (n_a_o)*(n_a_v), (n_a_o)*c_length)
!
 !      enddo ! Batching over c
 !   enddo ! Batching over a

!
   end subroutine omega_mlcc2_a1_mlcc2

!
end submodule