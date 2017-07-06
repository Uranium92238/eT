submodule (mlccsd_class) omega
!
   implicit none 
!
   logical :: debug = .false.
   logical :: timings = .false. 
!
!
contains
!
   module subroutine initialize_omega_mlccsd(wf)
!
!      Initialize Omega (MLCCSD)
!      Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!      Allocates the projection vector (omega1, omega2) and sets it
!      to zero.
!
      implicit none 
!
      class(mlccsd) :: wf
!
      call allocator(wf%omega1, wf%n_v, wf%n_o)
      wf%omega1 = zero
!
      call allocator(wf%omega2, wf%n_t2am, 1)
      wf%omega2 = zero
!
   end subroutine initialize_omega_mlccsd
!
    module subroutine omega_mlcc2_a1_mlccsd(wf, active_space)
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
!!    Batching over A
!!
!! 
      implicit none
!
      class(mlccsd)   :: wf
      integer(i15)   :: active_space ! Current active space IS ALWAYS ONE FOR MLCCSD !
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: A_n_batch = 0, A_first = 0, A_last = 0, A_batch = 0, A_length = 0
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, a = 0, c = 0, b = 0
      integer(i15) :: ib = 0, ic = 0, jb = 0, jc = 0
      integer(i15) :: bjc = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: s_ib_jc
      real(dp), dimension(:,:), allocatable :: u_bjc_i
      real(dp), dimension(:,:), allocatable :: g_Ab_jc
      real(dp), dimension(:,:), allocatable :: L_Ab_J
      real(dp), dimension(:,:), allocatable :: L_jc_J
!
      logical :: reorder  ! To get L_ab_J reordered, for batching over a
!
!     Active space variables
!
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      first_active_v = 1
      first_active_o = 1
!
      last_active_o = first_active_o + wf%n_total_active_o - 1
      last_active_v = first_active_v + wf%n_total_active_v - 1 
!
!     :: Construct u_ib_jc ::
!
!     u_ij^bc = 2*s_ij^bc - s_ij^cb =  (2*g_ij^bc - g_ij^cb)/ε_ij^cb
!
      call allocator(s_ib_jc, (wf%n_total_active_o)*(wf%n_total_active_v), (wf%n_total_active_o)*(wf%n_total_active_v))
      s_ib_jc = zero
      call wf%get_mlccsd_s2am(s_ib_jc)
!
      call allocator(u_bjc_i, (wf%n_total_active_v**2)*(wf%n_total_active_o), (wf%n_total_active_o))
!
      do b = 1, wf%n_total_active_v
         do i = 1, wf%n_total_active_o
!
            ib = index_two(i, b, wf%n_total_active_o)
!
            do c = 1, wf%n_total_active_v
               do j = 1, wf%n_total_active_o
!
                  jc = index_two(j, c, wf%n_total_active_o)
                  jb = index_two(j, b, wf%n_total_active_o)
                  ic = index_two(i, c, wf%n_total_active_o)
!
                  bjc = index_three(b, j, c, wf%n_total_active_v, wf%n_total_active_o)
!
                  u_bjc_i(bjc,i) = (two*s_ib_jc(ib,jc)-s_ib_jc(ic, jb))
!
               enddo
            enddo
         enddo
      enddo
!      
      call deallocator(s_ib_jc,(wf%n_total_active_o)*(wf%n_total_active_v), (wf%n_total_active_o)*(wf%n_total_active_v))
!
!     Prepare for batching over A
!
      required = (wf%n_total_active_v)*(wf%n_v)*(wf%n_J) &
               + (wf%n_total_active_v)*(wf%n_total_active_o)*(wf%n_J) &
               + ((wf%n_total_active_v)**2)*(wf%n_o)*(wf%n_total_active_o) &
               + ((wf%n_total_active_o)**2)*((wf%n_total_active_v)**2)
!
!
      required = required*4  ! Words

      available = get_available()

      max_length = 0
      call num_batch(required, available, max_length, A_n_batch, wf%n_v)
!
      A_first  = 0
      A_last   = 0
      A_length = 0
!
!     Start looping over a-batches
!
      do A_batch = 1, A_n_batch
!   
         call batch_limits(A_first, A_last, A_batch, max_length, wf%n_v)
         A_length = A_last - A_first + 1   
!
!        :: Construct integral g_Ab,jc ::
!
!        Get cholesky vectors
!
         call allocator(L_jc_J, (wf%n_total_active_o)*(wf%n_total_active_v), wf%n_J)
         L_jc_J = zero
         call wf%get_cholesky_ia(L_jc_J, first_active_o, last_active_o, first_active_v, last_active_v)
!
         call allocator(L_Ab_J, (wf%n_total_active_v)*a_length, wf%n_J) 
         L_Ab_J = zero
!
         call wf%get_cholesky_ab(L_Ab_J, a_first, a_last, first_active_v, last_active_v)
!
         call allocator(g_Ab_jc, (wf%n_total_active_v)*a_length, (wf%n_total_active_o)*(wf%n_total_active_v))      
!
!        g_Ab,jc = sum_J L_Ab^J*L_jc^J
!
         call dgemm('N', 'T',                                        &
                     (wf%n_total_active_v)*a_length,                 &
                     (wf%n_total_active_o)*(wf%n_total_active_v),    &
                     (wf%n_J),                 &
                     one,                      &
                     L_Ab_J,                   &
                     (wf%n_total_active_v)*a_length,    &
                     L_jc_J,                   &
                     (wf%n_total_active_o)*(wf%n_total_active_v),    &
                     zero,                     &
                     g_Ab_jc,                  &
                     (wf%n_total_active_v)*a_length) 
!
         call deallocator(L_jc_J, (wf%n_total_active_o)*(wf%n_total_active_v), wf%n_J)
         call deallocator(L_Ab_J, (wf%n_total_active_v)*a_length, wf%n_J) 
!
!        :: Add contributions to omega ::
!
         call dgemm('N', 'N',                            &
                     A_length,                           &
                     wf%n_total_active_o,                &
                    (wf%n_total_active_o)*(wf%n_total_active_v**2), &
                     one,                                &
                     g_Ab_jc,                            &
                     A_length,                           &
                     u_bjc_i,                            &
                     (wf%n_total_active_o)*(wf%n_total_active_v**2), &
                     one,                                &
                     wf%omega1(A_first,1),               &
                     wf%n_v)
! 
         call deallocator(g_Ab_jc, (wf%n_total_active_v)*a_length, (wf%n_total_active_o)*(wf%n_total_active_v))      
!
      enddo ! Batching over a
!
      call deallocator(u_bjc_i, (wf%n_total_active_v**2)*(wf%n_total_active_o), (wf%n_total_active_o))
!
   end subroutine omega_mlcc2_a1_mlccsd
!
!
    module subroutine omega_mlcc2_b1_mlccsd(wf, active_space)
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
      implicit none
!
      class(mlccsd)   :: wf
      integer(i15)   :: active_space 
!
!     Batching 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: s_ja_kb 
      real(dp), dimension(:,:), allocatable :: u_a_kbj
      real(dp), dimension(:,:), allocatable :: g_kb_jI
      real(dp), dimension(:,:), allocatable :: L_kb_J
      real(dp), dimension(:,:), allocatable :: L_jI_J
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
      first_active_v = 1
      first_active_o = 1
!
      n_active_o = wf%n_total_active_o 
      n_active_v = wf%n_total_active_v
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
!     :: Construct u_jk^ab ::
!  
!     u_jk^ab = 2*s_jk^ab - s_jk^ba  (place in u_a_jkb)        
!  
      call allocator(s_ja_kb, (n_active_o)*n_active_v, (n_active_o)*(n_active_v))
      call wf%get_mlccsd_s2am(s_ja_kb)
!
      call allocator(u_a_kbj, n_active_v, (n_active_o**2)*n_active_v)

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
         call deallocator(s_ja_kb, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
!        :: - sum_bjk u_ja_kb * g_kb_jI ::
!
!        Get cholesky vectors
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
!        g_kb_jI = sum_J L_kb^J*L_jI_J
!
         call dgemm('N', 'T',                    &
                     (n_active_o)*n_active_v,    &
                     (n_active_o)*(wf%n_o),      &
                     (wf%n_J),                   &
                     one,                        &
                     L_kb_J,                     &
                     (n_active_o)*n_active_v,    &
                     L_ji_J,                     &
                     (n_active_o)*(wf%n_o),      &
                     zero,                       &
                     g_kb_jI,                    &
                     (n_active_o)*n_active_v) 
!
         call deallocator(L_jI_J, n_active_o*(wf%n_o), wf%n_J)
         call deallocator(L_kb_J, n_active_o*n_active_v, wf%n_J)
!
!        Add contributions to omega
!
         call dgemm('N', 'N',                        &
                     n_active_v,                     &
                     (wf%n_o),                       &
                     n_active_v*((n_active_o)**2),   &
                     -one,                           &
                     u_a_kbj,                        &
                     n_active_v,                     &
                     g_kb_jI,                        &
                     n_active_v*((n_active_o)**2),   &
                     one,                            &
                     wf%omega1(first_active_v, 1),   &
                     (wf%n_v))
!
         call deallocator(g_kb_jI, n_active_o*n_active_v, n_active_o*(wf%n_o))
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
         call deallocator(u_a_kbj, n_active_v, (n_active_o**2)*n_active_v)
!      
   end subroutine omega_mlcc2_b1_mlccsd
!
!
       module subroutine omega_mlccsd_a1_mlccsd(wf, active_space)
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
!!     u_ij^bc = 2*t_ij^bc - t_ij^cb 
!!
!!    Batching over A
!!
!! 
      implicit none
!
      class(mlccsd)   :: wf
      integer(i15)   :: active_space
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
      integer(i15) :: A_n_batch = 0, A_first = 0, A_last = 0, A_batch = 0, A_length = 0
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, a = 0, c = 0, b = 0
      integer(i15) :: bi = 0, ci = 0, bj = 0, cj = 0
      integer(i15) :: bjc = 0
      integer(i15) :: bicj = 0, cibj = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: u_bjc_i
      real(dp), dimension(:,:), allocatable :: g_Ab_jc
      real(dp), dimension(:,:), allocatable :: L_Ab_J
      real(dp), dimension(:,:), allocatable :: L_jc_J
!
      logical :: reorder  ! To get L_ab_J reordered, for batching over a
!
!     Active space variables
!
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
      integer(i15) :: n_active_o    ! number of active occupied index 
      integer(i15) :: n_active_v    ! number of active virtual index
!
      integer(i15) :: offset = 0
!
!     Calculate first/last indeces
! 
      call wf%get_CCSD_active_indices(first_active_o, first_active_v, active_space)
!
      n_active_o = wf%n_CCSD_o(active_space, 1)
      n_active_v = wf%n_CCSD_v(active_space, 1)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     :: Construct u_ib_jc ::
!
!     u_ij^bc = 2*t_ij^bc - t_ij^cb
!
      call allocator(u_bjc_i, (n_active_v**2)*(n_active_o), (n_active_o))
!
      offset = wf%active_space_t2am_offset(active_space)
!
      do b = 1, n_active_v
         do i = 1, n_active_o
!
            bi = index_two(b, i, n_active_v)
!
            do c = 1, n_active_v
               do j = 1, n_active_o
!
                  cj = index_two(c, j, n_active_v)
                  bj = index_two(b, j, n_active_v)
                  ci = index_two(c, i, n_active_v)
!
                  cibj = index_packed(ci, bj)
                  bicj = index_packed(bi, cj)
!
                  bjc = index_three(b, j, c, n_active_v, n_active_o)
!
                  u_bjc_i(bjc,i) = (two*wf%t2am(bicj + offset,1)-wf%t2am(cibj + offset,1))
!
               enddo
            enddo
         enddo
      enddo
!
!     Prepare for batching over A
!
      required = (n_active_v)*(wf%n_v)*(wf%n_J) &
               + (n_active_v)*(n_active_o)*(wf%n_J) &
               + ((n_active_v)**2)*(wf%n_o)*(n_active_o) &
               + ((n_active_o)**2)*((n_active_v)**2)
!
!
      required = required*4  ! Words

      available = get_available()

      max_length = 0
!
      call num_batch(required, available, max_length, A_n_batch, wf%n_v)
!
      A_first  = 0
      A_last   = 0
      A_length = 0
!
!     Start looping over a-batches
!
      do A_batch = 1, A_n_batch
!   
         call batch_limits(A_first, A_last, A_batch, max_length, wf%n_v)
         A_length = A_last - A_first + 1   
!
!        :: Construct integral g_Ab,jc ::
!
!        Get cholesky vectors
!
         call allocator(L_jc_J, (n_active_o)*(n_active_v), wf%n_J)
         L_jc_J = zero
         call wf%get_cholesky_ia(L_jc_J, first_active_o, last_active_o, first_active_v, last_active_v)
!
         call allocator(L_Ab_J, (n_active_v)*a_length, wf%n_J) 
         L_Ab_J = zero
!
         call wf%get_cholesky_ab(L_Ab_J, a_first, a_last, first_active_v, last_active_v)
!
         call allocator(g_Ab_jc, (n_active_v)*a_length, (n_active_o)*(n_active_v))      
!
!        g_Ab,jc = sum_J L_Ab^J*L_jc^J
!
         call dgemm('N', 'T',                  &
                     (n_active_v)*a_length,    &
                     (n_active_o)*(n_active_v),&
                     (wf%n_J),                 &
                     one,                      &
                     L_Ab_J,                   &
                     (n_active_v)*A_length,    &
                     L_jc_J,                   &
                     (n_active_o)*(n_active_v),&
                     zero,                     &
                     g_Ab_jc,                  &
                     (n_active_v)*A_length) 
!
         call deallocator(L_jc_J, (wf%n_total_active_o)*(wf%n_total_active_v), wf%n_J)
         call deallocator(L_Ab_J, (wf%n_total_active_v)*a_length, wf%n_J) 
!
!        :: Add contributions to omega ::
!
         call dgemm('N', 'N',                            &
                     A_length,                           &
                     n_active_o,                &
                    (n_active_o)*(n_active_v**2), &
                     one,                                &
                     g_Ab_jc,                            &
                     A_length,                           &
                     u_bjc_i,                            &
                     (n_active_o)*(n_active_v**2), &
                     one,                                &
                     wf%omega1(A_first,1),               &
                     wf%n_v)
! 
         call deallocator(g_Ab_jc, (n_active_v)*A_length, (n_active_o)*(n_active_v))      
!
      enddo ! Batching over a
!
      call deallocator(u_bjc_i, (n_active_v**2)*(n_active_o), (n_active_o))
!
   end subroutine omega_mlccsd_a1_mlccsd
!
!
    module subroutine omega_mlccsd_b1_mlccsd(wf, active_space)
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
      implicit none
!
      class(mlccsd)   :: wf
      integer(i15)   :: active_space 
!
!     Batching 
!
      integer(i15) :: required = 0, available = 0, max_length = 0, batch_dimension = 0, offset = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: u_a_kbj
      real(dp), dimension(:,:), allocatable :: g_kb_jI
      real(dp), dimension(:,:), allocatable :: L_kb_J
      real(dp), dimension(:,:), allocatable :: L_jI_J
!
!     looping indices
!
      integer(i15) :: i = 0, j = 0, k = 0, a  = 0, b = 0
      integer(i15) :: I_full = 0, J_full = 0, A_full  = 0, B_full = 0
      integer(i15) :: ja = 0, bk = 0, ak  = 0, bj = 0, aj = 0
      integer(i15) :: kbj = 0, jbi = 0
      integer(i15) :: ajbk = 0, bjak = 0
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
      call wf%get_CCSD_active_indices(first_active_o, first_active_v, active_space)
!
      n_active_o = wf%n_CCSD_o(active_space, 1)
      n_active_v = wf%n_CCSD_v(active_space, 1)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     :: Construct u_jk^ab ::
!  
!     u_jk^ab = 2*s_jk^ab - s_jk^ba  (place in u_a_jkb)        
!  
!
      call allocator(u_a_kbj, n_active_v, (n_active_o**2)*n_active_v)

      do k = 1, n_active_o
         do b = 1, n_active_v         
!  
             bk = index_two(b, k, n_active_v)
!  
            do j = 1, n_active_o
!
               bj = index_two(b, j, n_active_v)
               kbj = index_three(k, b, j, n_active_o, n_active_v)
!
               do a = 1, n_active_v             
!  
                  aj = index_two(a, j, n_active_v)
                  ak = index_two(a, k, n_active_v)
!
                  ajbk = index_packed(aj, bk) 
                  bjak = index_packed(bj, ak)                 
!  
                  u_a_kbj(a,kbj) = (two*wf%t2am(ajbk + offset, 1)-wf%t2am(bjak + offset, 1))
!
               enddo
            enddo
         enddo
      enddo
!
!     :: - sum_bjk u_ja_kb * g_kb_jI ::
!
!     Get cholesky vectors
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
!     g_kb_jI = sum_J L_kb^J*L_jI_J
!
      call dgemm('N', 'T',                    &
                  (n_active_o)*n_active_v,    &
                  (n_active_o)*(wf%n_o),      &
                  (wf%n_J),                   &
                  one,                        &
                  L_kb_J,                     &
                  (n_active_o)*n_active_v,    &
                  L_ji_J,                     &
                  (n_active_o)*(wf%n_o),      &
                  zero,                       &
                  g_kb_jI,                    &
                  (n_active_o)*n_active_v) 
!
      call deallocator(L_jI_J, n_active_o*(wf%n_o), wf%n_J)
      call deallocator(L_kb_J, n_active_o*n_active_v, wf%n_J)
!
!     Add contributions to omega
!
      call dgemm('N', 'N',                        &
                  n_active_v,                     &
                  (wf%n_o),                       &
                  n_active_v*((n_active_o)**2),   &
                  -one,                           &
                  u_a_kbj,                        &
                  n_active_v,                     &
                  g_kb_jI,                        &
                  n_active_v*((n_active_o)**2),   &
                  one,                            &
                  wf%omega1(first_active_v, 1),   &
                  (wf%n_v))
!
      call deallocator(g_kb_jI, n_active_o*n_active_v, n_active_o*(wf%n_o))
!
!     :: sum_jb F_jb u_ij^ab ::
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
      call deallocator(u_a_kbj, n_active_v, (n_active_o**2)*n_active_v)
!      
   end subroutine omega_mlccsd_b1_mlccsd
!
!
  module subroutine omega_ccsd_a2_mlccsd(wf, active_space)
!
!     Omega A2 term: Omega A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!     Structure: Batching over both a and b for A2.2.
!                t^+_ci_dj = t_ci_dj + t_di_cj
!                t^-_ci_dj = t_ci_dj - t_di_cj
!                g^+_ac_bd = g_ac_bd + g_bc_ad 
!                g^-_ac_bd = g_ac_bd - g_bc_ad 
! 
!                omega_A2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!                omega_A2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!
      implicit none
!
      class(ccsd)  :: wf
      integer(i15) :: active_space
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj 
      real(dp), dimension(:,:), allocatable :: g_ca_db 
      real(dp), dimension(:,:), allocatable :: g_p_ab_cd
      real(dp), dimension(:,:), allocatable :: g_m_ab_cd
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_ca_J 
      real(dp), dimension(:,:), allocatable :: L_ac_J 
      real(dp), dimension(:,:), allocatable :: L_db_J 
      real(dp), dimension(:,:), allocatable :: L_bd_J 
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
      integer(i15) :: ab = 0, ca = 0, cb = 0, cd = 0, da = 0, db = 0, ac = 0, bd = 0 
      integer(i15) :: ai = 0, aj = 0, bj = 0, bi = 0, ci = 0, cj = 0, dj = 0, di = 0
      integer(i15) :: ij = 0
!
      integer(i15) :: aibj = 0, biaj = 0, cidj = 0, cjdi = 0
!
!     Batching and memory handling variables
!
      integer(i15) :: a_n_batch = 0, a_first = 0, a_last = 0, a_length = 0, a_max_length = 0, a_batch = 0
      integer(i15) :: b_n_batch = 0, b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0

      integer(i15) :: required = 0, available = 0
!
!     Logical for reordering in L_ab_J when batching over the last index
!
      logical :: reorder
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
      call wf%get_CCSD_active_indices(first_active_o, first_active_v, active_space)
!
      n_active_o = wf%n_CCSD_o(active_space, 1)
      n_active_v = wf%n_CCSD_v(active_space, 1)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
      offset = wf%active_space_t2am_offset(active_space)
!
!     ::  Calculate the A2.1 term of omega ::
!
!     Create g_ai_bj
!  
      call allocator(g_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
      call allocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
!
      call wf%get_cholesky_ai(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
!     g_ai_bj = sum_J L_ai_J*L_bj_J
!     
      call dgemm('N','T',            &
                  n_active_o*n_active_v, &
                  n_active_o*n_active_v, &
                  wf%n_J,            &
                  one,               & 
                  L_ai_J,            &
                  n_active_o*n_active_v, &
                  L_ai_J,            &
                  n_active_o*n_active_v, &
                  zero,              &
                  g_ai_bj,           &
                  n_active_o*n_active_v)
!
!
      call deallocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
!
!     Add A2.1 to Omega 2
!     
      do i = 1, n_active_o
         do a = 1, n_active_v
!
            ai = index_two(a, i, n_active_v)
!
            do j = 1, n_active_o
               do b = 1, n_active_v
!
                  bj = index_two(b, j, n_active_v)
!
                  if(ai .ge. bj) then
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj + offset, 1) = wf%omega2(aibj + offset, 1) + g_ai_bj(ai, bj)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ai_bj, n_active_o*n_active_v, n_active_o*n_active_v)
!
!     ::  Calculate the A2.2 term  of omega ::
!
!
      required = max(3*(n_active_v)**2*(wf%n_J) + 2*(n_active_v)*(n_active_o)*(wf%n_J),      & ! Needed to get  L_db_J
                     (n_active_v)**4 + 2*(n_active_v)**2*(wf%n_J), &                       ! Needed to get g_ac_bd
                     (n_active_v)**4 + 2*(packed_size(n_active_v))*(packed_size(n_active_v)) & ! Needed to get g+- and t+-
                     + 2*(packed_size(n_active_v))*(packed_size(n_active_o)), &            !
                       2*(packed_size(n_active_v))*(packed_size(n_active_v)) &             ! Needed for g+- and t+- and Omega+-
                     + 2*(packed_size(n_active_v))*(packed_size(n_active_o)) &             !
                     + 2*(n_active_v)**2*(packed_size(n_active_v)))                        !
!
      required = required*4  ! Words

      available=get_available()
!
      a_max_length = 0
      call num_two_batch(required, available, a_max_length, a_n_batch, n_active_v)
!
!     Initialize some variables for batching
!
      a_first  = 0
      a_last   = 0
      a_length = 0
!
!     Start looping over a-batches
!
      do a_batch = 1,a_n_batch
!   
         call batch_limits(a_first ,a_last ,a_batch, a_max_length, n_active_v)
         a_length = a_last - a_first + 1     
!
!        Start looping over batches of b
!
         b_first  = 0
         b_last   = 0
         b_length = 0
!
         b_max_length = a_max_length
!
         do b_batch = 1, a_batch
!
            call batch_limits(b_first ,b_last ,b_batch, b_max_length, n_active_v)
            b_length = b_last - b_first + 1 
!
!           Get ab-cholesky vectors for the batch, L_ac^J, then reorder from L_ac_J to L_ca_J
!
            call allocator(L_ac_J, (n_active_v)*a_length, wf%n_J)
            L_ac_J = zero
!
            call wf%get_cholesky_ab(L_ac_J, a_first, a_last, first_active_v,last_active_v)
!
!
!
!          Get ab-cholesky vectors for the batch, L_bd^J, then reorder from L_bd_J to L_db_J
!
           call allocator(L_bd_J, (wf%n_v)*b_length, wf%n_J)
           L_bd_J = zero
!
           call wf%get_cholesky_ab(L_bd_J, b_first, b_last, first_active_v, last_active_v)         
! 
!          Allocate g_ca_db
!
           call allocator(g_ac_bd, (n_active_v)*a_length, (n_active_v)*b_length)
           g_ac_bd = zero
!
!          g_ca_db = sum_J L_ca_J*L_db_J
!    
           call dgemm('N','T',            &
                       (wf%n_v)*a_length, &
                       (wf%n_v)*b_length, &
                       wf%n_J,            &
                       one,               &
                       L_ac_J,            &
                       (wf%n_v)*a_length, &
                       L_bd_J,            &
                       (wf%n_v)*b_length, &
                       zero,              &
                       g_ac_bd,           &
                       (wf%n_v)*a_length)
!
            call deallocator(L_ac_J, (wf%n_v)*a_length, wf%n_J)
            call deallocator(L_bd_J, (wf%n_v)*b_length, wf%n_J) 
!
            if (b_batch .eq. a_batch) then
!
!
!           Allocate for +-g, +-t
!
               call allocator(g_p_ab_cd, packed_size(a_length), packed_size(n_active_v))
               call allocator(g_m_ab_cd, packed_size(a_length), packed_size(n_active_v))
               call allocator(t_p_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
               call allocator(t_m_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
!
               g_p_ab_cd = zero
               g_m_ab_cd = zero
               t_p_cd_ij = zero
               t_m_cd_ij = zero
!
!              Reorder g_ca_db to g_ab_cd and t_ci_dj to t_cd_ij
! 
               do c = 1, n_active_v 
                  do d = 1, c
!
                     cd = index_packed(c, d)
!
                     do a = 1, a_length
!
                        ac = index_two(a, c, a_length)
                        ad = index_two(a, d, a_length)
!
                        do  b = 1, b_length
                           if ((a+a_first-1) .ge. (b+b_first-1)) then
!                            
                              bd = index_two(b, d, b_length)
                              bc = index_two(b, c, b_length)
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
                    do i = 1, n_active_o
                       do j = 1, i
!    
                          ij = index_packed(i, j)
!    
                          ci = index_two(c, i, n_active_v)
                          dj = index_two(d, j, n_active_v)
                          cj = index_two(c, j, n_active_v)
                          di = index_two(d, i, n_active_v)
! 
                          cidj = index_packed(ci, dj)
                          cjdi = index_packed(cj, di)
! 
                          t_p_cd_ij(cd, ij) = wf%t2am(cidj + offset, 1) + wf%t2am(cjdi + offset, 1)
                          t_m_cd_ij(cd, ij) = wf%t2am(cidj + offset, 1) - wf%t2am(cjdi + offset, 1)  
!
                       enddo
                    enddo
                 enddo
              enddo
!
!              Dellocate g_ac_bd 
!
               call deallocator(g_ac_bd, (n_active_v)*a_length, (n_active_v)*b_length)
!
!              Allocate omega +-
!
              call allocator(omega2_p_ab_ij, packed_size(a_length), packed_size(n_active_o))
              call allocator(omega2_m_ab_ij, packed_size(a_length), packed_size(n_active_o))
! 
!              omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
! 
              call dgemm('N','N',                & 
                          packed_size(a_length), &
                          packed_size(n_active_o),   &
                          packed_size(n_active_v),   &
                          one/four,              &
                          g_p_ab_cd,             &
                          packed_size(a_length), &
                          t_p_cd_ij,             &
                          packed_size(n_active_v),   &
                          zero,                  &
                          omega2_p_ab_ij,        &
                          packed_size(a_length))
!
              call dgemm('N','N',                & 
                          packed_size(a_length), &
                          packed_size(n_active_o),   &
                          packed_size(n_active_v),   &
                          one/four,              &
                          g_m_ab_cd,             &
                          packed_size(a_length), &
                          t_m_cd_ij,             &
                          packed_size(n_active_v),   &
                          zero,                  &
                          omega2_m_ab_ij,        &
                          packed_size(a_length) )
!
!             Deallocate +-g, +-t
! 
              call deallocator(g_p_ab_cd, packed_size(a_length), packed_size(n_active_v))
              call deallocator(g_m_ab_cd, packed_size(a_length), packed_size(n_active_v))
              call deallocator(t_p_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
              call deallocator(t_m_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
!
              do i = 1, n_active_o
                 do j = 1, i
!
                    
                    ij = index_packed(i, j)
!
                    do a = 1, a_length
!
                       Ai = index_two(a + a_first - 1, i, n_active_v) ! A is full-space a index
                       Aj = index_two(a + a_first - 1, j, n_active_v) ! A is full-space a index
!
                       do b = 1, b_length
!                
                          if ((a+a_first-1) .ge. (b+b_first-1)) then
                             Bj = index_two(b + b_first - 1, j, n_active_v) ! B is full-space b index
                             Bi = index_two(b + b_first - 1, i, n_active_v) ! B is full-space b index
!
!
                             ab = index_packed(a, b)
!    
                             AiBj = index_packed(Ai, Bj)
                             BiAj = index_packed(Bi, Aj)
!                         
!                            Reorder into omega2_aibj
! 
                             wf%omega2(AiBj + offset,1) = wf%omega2(AiBjoffset + , 1) + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                             if (AiBj .ne. BiAj) then
                                wf%omega2(BiAj + offset,1) = wf%omega2(BiAj + offset, 1) + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
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
               call deallocator(omega2_p_ab_ij, packed_size(a_length), packed_size(n_active_o))
               call deallocator(omega2_m_ab_ij, packed_size(a_length), packed_size(n_active_o))
            else
!
!              Allocate for +-g, +-t
!
               call allocator(g_p_ab_cd, a_length*b_length, packed_size(n_active_v))
               call allocator(g_m_ab_cd, a_length*b_length, packed_size(n_active_v))
               call allocator(t_p_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
               call allocator(t_m_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
!
               g_p_ab_cd = zero
               g_m_ab_cd = zero
               t_p_cd_ij = zero
               t_m_cd_ij = zero
! 
               do c = 1, n_active_v 
                  do d = 1, c
!
                     cd = index_packed(c, d)
!
                     do a = 1, a_length
!
                        ac = index_two(a, c, a_length)
                        ad = index_two(a, d, a_length)
!
                        do  b = 1, b_length
!                            
                             bd = index_two(b, d, b_length)
                             bc = index_two(b, c, b_length)
!
                             ab = index_two(a, b, a_length)
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
                    do i = 1, n_active_o
                       do j = 1, i
!    
                          ij = index_packed(i, j)
!    
                          ci = index_two(c, i, n_active_v)
                          dj = index_two(d, j, n_active_v)
                          cj = index_two(c, j, n_active_v)
                          di = index_two(d, i, n_active_v)
! 
                          cidj = index_packed(ci, dj)
                          cjdi = index_packed(cj, di)
! 
                          t_p_cd_ij(cd, ij) = wf%t2am(cidj + offset, 1) + wf%t2am(cjdi + offset, 1)
                          t_m_cd_ij(cd, ij) = wf%t2am(cidj + offset, 1) - wf%t2am(cjdi + offset, 1)  
!
                       enddo
                    enddo
                 enddo
              enddo
!
!              Dellocate g_ac_bd 
!
               call deallocator(g_ac_bd, (n_active_v)*a_length, (n_active_v)*b_length)
!
!              Allocate omega +-
!

               call allocator(omega2_p_ab_ij, b_length*a_length, packed_size(n_active_o))
               call allocator(omega2_m_ab_ij, b_length*a_length, packed_size(n_active_o))
!  
!               omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
! 
              call dgemm('N','N',                & 
                          b_length*a_length,     &
                          packed_size(n_active_o),   &
                          packed_size(n_active_v),   &
                          one/four,              &
                          g_p_ab_cd,             &
                          b_length*a_length,     &
                          t_p_cd_ij,             &
                          packed_size(n_active_v),   &
                          zero,                  &
                          omega2_p_ab_ij,        &
                          b_length*a_length)
!
              call dgemm('N','N',                & 
                          b_length*a_length,     &
                          packed_size(n_active_o),   &
                          packed_size(n_active_v),   &
                          one/four,              &
                          g_m_ab_cd,             &
                          b_length*a_length,     &
                          t_m_cd_ij,             &
                          packed_size(n_active_v),   &
                          zero,                  &
                          omega2_m_ab_ij,        &
                          b_length*a_length)
!
!          Deallocate +-g, +-t
! 
              call deallocator(g_p_ab_cd, b_length*a_length, packed_size(n_active_v))
              call deallocator(g_m_ab_cd, b_length*a_length, packed_size(n_active_v))
              call deallocator(t_p_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
              call deallocator(t_m_cd_ij, packed_size(n_active_v), packed_size(n_active_o))
!
               do i = 1, n_active_o
                  do j = 1, i
!
                     
                     ij = index_packed(i, j)
!
                     do a = 1, a_length
!
                        Ai = index_two(a + a_first - 1, i, n_active_v) ! A is full-space a index
                        Aj = index_two(a + a_first - 1, j, n_active_v) ! A is full-space a index
!
                        do b = 1, b_length
!                 
                              Bj = index_two(b + b_first - 1, j, n_active_v) ! B is full-space b index
                              Bi = index_two(b + b_first - 1, i, n_active_v) ! B is full-space b index
!
!
                              ab = index_two(a, b, a_length)

!     
                              AiBj = index_packed(Ai, Bj)
                              BiAj = index_packed(Bi, Aj)
!                          
!                             Reorder into omega2_aibj
!  
                              wf%omega2(AiBj + offset,1) = wf%omega2(AiBj + offset, 1) + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                              if (AiBj .ne. BiAj) then
                                 wf%omega2(BiAj + offset,1) = wf%omega2(BiAj + offset, 1) + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
                              endif   
!     
                        enddo
                     enddo
                  enddo
               enddo
!
!              Deallocate omega +-
!
               call deallocator(omega2_p_ab_ij, b_length*a_length, packed_size(n_active_o))
               call deallocator(omega2_m_ab_ij, b_length*a_length, packed_size(n_active_o))
            endif
!
         enddo
!
      enddo
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine get_mlccsd_s2am_mlccsd(wf, s_ia_jb)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:) :: s_ia_jb
!
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: s_ai_bj_CC2
      real(dp), dimension(:,:), allocatable :: X1 ! Intermediato of transformation
      real(dp), dimension(:,:), allocatable :: X2 ! Intermediato of transformation
      real(dp), dimension(:,:), allocatable :: X3 ! Intermediato of transformation
      real(dp), dimension(:,:), allocatable :: X4 ! Intermediato of transformation
!
!     Active space variables
!
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
      integer(i15) :: offset
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, ia = 0, bj = 0, ai = 0, aIB, jB
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: CCSD_active_space, CC2_active_space
!
      do CC2_active_space = 1, wf%n_active_spaces
         if (wf%n_CC2_o(CC2_active_space, 1) .gt. 0) then
!
            n_active_o = wf%n_CC2_o(CC2_active_space, 1) + wf%n_CCSD_o(CC2_active_space, 1)
            n_active_v = wf%n_CC2_v(CC2_active_space, 1) + wf%n_CCSD_v(CC2_active_space, 1)
!  
            call wf%get_CCSD_active_indices(first_active_v, first_active_o, CC2_active_space)
!
            last_active_o = first_active_o + n_active_o - 1
            last_active_v = first_active_v + n_active_v - 1
!
            call allocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
            L_ai_J = zero
!
            call wf%get_cholesky_ai_for_cc2_amplitudes(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
!
            call allocator(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
!           g_ib_jc = g_bi,cj = sum_J L_bj^J*L_ci^J
!
            call dgemm('N', 'T',                                &
                        (n_active_o)*(n_active_v),  &
                        (n_active_o)*(n_active_v),  &
                        (wf%n_J),                               &
                        one,                                    &
                        L_ai_J,                                 &
                        (n_active_o)*(n_active_v),  &
                        L_ai_J,                                 &
                        (n_active_o)*(n_active_v),  &
                        zero,                                   &
                        g_ai_bj,                                &
                        (n_active_o)*(n_active_v))        
!
            call deallocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
!
            call allocator(s_AI_BJ_CC2, n_active_o*n_active_v, n_active_o*n_active_v)
!
            do A = 1, n_active_v
               do I = 1, n_active_o
!
                  AI = index_two(A, I, n_active_v)
!
                  do B = 1, n_active_v
                     do J = 1, n_active_o
!
                        BJ = index_two(B, J, n_active_v)
!
                        s_AI_BJ_CC2(AI, BJ) = g_ai_bj(AI, BJ)/(wf%fock_diagonal(i + first_active_o - 1,1)&
                                           + wf%fock_diagonal(j + first_active_o - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + b + first_active_v - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + a + first_active_v - 1,1))
!
                     enddo
                  enddo
               enddo
            enddo
!
!           Change basis
!
!           Transform index A : X1_a_IBJ = sum_A T_aA*s_A_IBJ_CC2
!
            call allocator(X1, wf%n_total_active_v, (n_active_o**2)*n_active_v)
!
            call dgemm('N','N', &
                     wf%n_total_active_v,    &
                     (n_active_o**2)*n_active_v, &
                     n_active_v, &
                     one, &
                     wf%T(wf%n_total_active_o + 1, wf%n_total_active_o + first_active_v), & !!! FEIL
                     wf%n_total_active_o + wf%n_total_active_v, &
                     s_AI_BJ_CC2, &
                     n_active_v, &
                     zero, &
                     X1, &
                     wf%n_total_active_v)
!
            call deallocator(s_AI_BJ_CC2, n_active_o*n_active_v, n_active_o*n_active_v)
!
!           Transform index J : X2_aIB_j = sum_J X1_aIB_J*T_jJ 
!
            call allocator(X2, wf%n_total_active_o*n_active_v*n_active_o, wf%n_total_active_v)
!
            call dgemm('N', 'T', &
                     wf%n_total_active_v*n_active_v*n_active_o, &
                     wf%n_total_active_o, &
                     n_active_o, &
                     one, &
                     X1, &
                     wf%n_total_active_v*n_active_v*n_active_o, &
                     wf%T(1, first_active_o), &
                     wf%n_total_active_o + wf%n_total_active_v, &
                     zero, &
                     X2, &
                     wf%n_total_active_v*n_active_v*n_active_o)
!
            call deallocator(X1, wf%n_total_active_v, (n_active_o**2)*n_active_v)
!
!           Reorder X2_aIB_j to X3_Ia_jB
!
            call allocator(X3, n_active_o*(wf%n_total_active_v), n_active_v*(wf%n_total_active_o))
            X3 = zero
!
            do a = 1, wf%n_total_active_v
               do j = 1, wf%n_total_active_o
                  do B = 1, n_active_v
!
                     jB = index_two(j, B, wf%n_total_active_o)
!
                     do I = 1, n_active_o
!
                        aIB = index_three(a, I, B, wf%n_total_active_v, n_active_o)
                        Ia = index_two(I, a, n_active_o)
!
                        X3(Ia, jB) = X2(aIB, j)
!
                     enddo
                  enddo
               enddo
            enddo
            call deallocator(X2, wf%n_total_active_o*n_active_v*n_active_o, wf%n_total_active_v)
!
!           Transform index I : X4_i_ajB = sum_I T_iI*X_3_Ia_jB 
!
            call allocator(X4, wf%n_total_active_o, (wf%n_total_active_o)*(wf%n_total_active_v)*n_active_v)
!
            call dgemm('N', 'N', &
                     wf%n_total_active_o, &
                     (wf%n_total_active_o)*(wf%n_total_active_v)*n_active_v, &
                     n_active_o, &
                     one, &
                     wf%T(1, first_active_o), &
                     wf%n_total_active_o + wf%n_total_active_v, &
                     X3, &
                     n_active_o, &
                     zero, &
                     X4, &
                     wf%n_total_active_o)
!
            call deallocator(X3, n_active_o*(wf%n_total_active_v), n_active_v*(wf%n_total_active_o))
!
!           Transform index B : s_ia_jb = sum_B X4_i_ajB*T_bB
!
            call dgemm('N', 'T', &
                     (wf%n_total_active_o**2)*(wf%n_total_active_v), &
                     wf%n_total_active_v, &
                     n_active_v, &
                     one, &
                     X4, &
                     (wf%n_total_active_o**2)*(wf%n_total_active_v), &
                     wf%T(wf%n_total_active_o + 1, wf%n_total_active_o + first_active_v), &
                     wf%n_total_active_o + wf%n_total_active_v, &
                     zero, &
                     s_ia_jb, &
                     (wf%n_total_active_o**2)*(wf%n_total_active_v))
!
            call deallocator(X4, wf%n_total_active_o, (wf%n_total_active_o)*(wf%n_total_active_v)*n_active_v)
!
         endif
      enddo
!
!     Zero out terms in CCSD area
!
      do CCSD_active_space = 1, wf%n_active_spaces
!
         n_active_o =  wf%n_CCSD_o(CCSD_active_space, 1)
         n_active_v =  wf%n_CCSD_v(CCSD_active_space, 1)
!
         call wf%get_CCSD_active_indices(first_active_v, first_active_o, CCSD_active_space)
!
         last_active_o = first_active_o + n_active_o - 1
         last_active_v = first_active_v + n_active_v - 1
!
         do i = first_active_o, last_active_o
            do j = first_active_o, last_active_o
               do a = first_active_v, last_active_v
                  do b = first_active_v, last_active_v
!
                     ia = index_two(i, a, wf%n_total_active_o)
                     jb = index_two(j, b, wf%n_total_active_o)
!
                     s_ia_jb(ia, jb) = zero
!
                  enddo
               enddo
            enddo
         enddo
!
      enddo        
!
   end subroutine get_mlccsd_s2am_mlccsd
!
end submodule omega