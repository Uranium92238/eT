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
   subroutine initialize_omega_mlccsd(wf)
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
    subroutine omega_mlcc2_a1_mlccsd(wf, active_space)
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
    subroutine omega_mlcc2_b1_mlccsd(wf, active_space)
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
       subroutine omega_mlccsd_a1_mlccsd(wf, active_space)
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
    subroutine omega_mlccsd_b1_mlccsd(wf, active_space)
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

  subroutine omega_mlccsd_a2_mlccsd(wf, active_space)
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
      integer(i15) :: offset
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
   end subroutine omega_mlccsd_a2_mlccsd
!
!
   subroutine omega_mlccsd_b2_mlccsd(wf, active_space)
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
      real(dp), dimension(:,:), allocatable :: L_kc_J     
      real(dp), dimension(:,:), allocatable :: L_ij_J  
      real(dp), dimension(:,:), allocatable :: g_kc_ld    
      real(dp), dimension(:,:), allocatable :: g_kl_cd    
      real(dp), dimension(:,:), allocatable :: g_kl_ij    
      real(dp), dimension(:,:), allocatable :: g_ki_lj 
!
!     Reordered T2 apmlitudes
!   
      real(dp), dimension(:,:), allocatable :: t_cd_ij    
      real(dp), dimension(:,:), allocatable :: t_ab_kl   
!
!     Intermediate for matrix multiplication
! 
      real(dp), dimension(:,:), allocatable :: X_kl_ij 
!
!     Reordered omega
!   
      real(dp), dimension(:,:), allocatable :: omega_ab_ij
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
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
      integer(i15) :: offset ! first active virtual index
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
!     Read Cholesky vector of type L_ij_J
!
      call allocator(L_ij_J, (n_active_o)*(n_active_o), wf%n_J)
!
      call wf%get_cholesky_ij(L_ij_J, first_active_o, last_active_o, first_active_o, last_active_o)
!
!     Create g_ki_lj = sum_J L_li_J*L_lj_J
!
      call allocator(g_ki_lj, (n_active_o)*(n_active_o), (n_active_o)*(n_active_o)) 
!
      call dgemm('N','T',            &
                  (n_active_o)*(n_active_o), &
                  (n_active_o)*(n_active_o), &
                  wf%n_J,            &
                  one,               &
                  L_ij_J,            &
                  (n_active_o)*(n_active_o), &
                  L_ij_J,            &
                  (n_active_o)*(n_active_o), &
                  zero,              &
                  g_ki_lj,           &
                  (n_active_o)*(n_active_o))
!
!
      call deallocator(L_ij_J, (n_active_o)*(n_active_o), wf%n_J)
!
      call allocator(g_kl_ij, (n_active_o)*(n_active_o),(n_active_o)*(n_active_o))
!
      do k = 1, n_active_o
         do l = 1, n_active_o
            do i = 1, n_active_o
               do j=1, n_active_o
!
!                 Calculate compound indices
!
                  ki = index_two(k, i, n_active_o)
                  lj = index_two(l, j, n_active_o)
                  kl = index_two(k, l, n_active_o)
                  ij = index_two(i, j, n_active_o)
!
!                 Reordering g_ki_lj to g_kl_ij
!
                  g_kl_ij(kl, ij) = g_ki_lj(ki, lj)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ki_lj, (n_active_o)*(n_active_o), (n_active_o)*(n_active_o))
!
!     Read Cholesky vectors of ia-type into L_kc_J
!
      call allocator(L_kc_J, (n_active_o)*(n_active_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_kc_J, first_active_o, last_active_o, first_active_v, last_active_v)
!
!     Create g_ck_ld = sum_(J) L_kc_J*L_ld_J
!
      call allocator(g_kc_ld, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!  
      call dgemm('N','T',            &
                  (n_active_o)*(n_active_v), &
                  (n_active_o)*(n_active_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (n_active_o)*(n_active_v), &
                  L_kc_J,            &
                  (n_active_o)*(n_active_v), &
                  zero,              &
                  g_kc_ld,           &
                  (n_active_o)*(n_active_v))
!
!     Deallocate cholesky vectors L_ck_J
!
      call deallocator(L_kc_J, (n_active_o)*(n_active_v), wf%n_J)
!
!     Reorder g_kc_ld as g_kl_cd, also reordering t_ci_dj as t_cd_ij
!
      call allocator(t_cd_ij, (n_active_v)**2, (n_active_o)**2)
      call allocator(g_kl_cd, (n_active_o)**2, (n_active_v)**2)
!
      do d = 1, n_active_v
         do c = 1, n_active_v
!
            cd = index_two(c, d, n_active_v)
!
            do l = 1, n_active_o
!  
               ld = index_two(l, d, n_active_o)
               dj = index_two(d, l, n_active_v)
!
               do k = 1, wf%n_o              
!
                  kl = index_two(k, l, n_active_o)
                  kc = index_two(k, c, n_active_o)
                  ci = index_two(c, k, n_active_v)
                  ij = kl
!
                  cidj = index_packed(ci, dj)
!
                  g_kl_cd(kl, cd) = g_kc_ld(kc, ld)
                  t_cd_ij(cd, ij) = wf%t2am(cidj + offset, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_kc_ld
!
      call deallocator(g_kc_ld, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      call dgemm('N','N',      &
                  (n_active_o)**2, &
                  (n_active_o)**2, &
                  (n_active_v)**2, &
                  one,         &
                  g_kl_cd,     &
                  (n_active_o)**2, &
                  t_cd_ij,     &
                  (n_active_v)**2, &
                  one,         &
                  g_kl_ij,     &
                  (n_active_o)**2)
!
!     Deallocate t_cd_ij and g_kl_cd
!
      call deallocator(t_cd_ij, (n_active_v)**2, (n_active_o)**2)
      call deallocator(g_kl_cd, (n_active_o)**2, (n_active_v)**2)
!
!     Reorder t_ak_bl to t_ab_kl
!
      call allocator(t_ab_kl, (n_active_v)**2, (n_active_o)**2)
!
      do l=1, n_active_o
         do k = 1, n_active_o
!
            kl = index_two(k, l, n_active_o)
!
            do b = 1, n_active_v
!
               bl = index_two(b, l, n_active_v)
!
               do a = 1, n_active_v
!
                  ak = index_two(a, k, n_active_v)
                  ab = index_two(a, b, n_active_v)
                  
!
                  akbl = index_packed(ak, bl)
!
                  t_ab_kl(ab, kl) = wf%t2am(akbl + offset, 1)
!
               enddo
            enddo
         enddo
      enddo
!
!     omega_ab_ij = sum_(kl) t_ab_kl*X_kl_ij
!
      call allocator(omega_ab_ij, (n_active_v)**2, (n_active_o)**2)
!
      call dgemm('N','N',      &
                  (n_active_v)**2, &
                  (n_active_o)**2, &
                  (n_active_o)**2, &
                  one,         &
                  t_ab_kl,     &
                  (n_active_v)**2, &
                  g_kl_ij,     &
                  (n_active_o)**2, &
                  zero,        &
                  omega_ab_ij, &
                  (n_active_v)**2)
!
      call deallocator(t_ab_kl, (n_active_v)**2, (n_active_o)**2)
      call deallocator(g_kl_ij, (n_active_o)**2, (n_active_o)**2)
!
!     Reorder into omega2
!
      do i = 1, n_active_o
         do j = 1, n_active_o
!
            ij = index_two(i, j, n_active_o)
!
            do b = 1, n_active_v
!
               bj = index_two(b, j, n_active_v)
!
               do a = 1, n_active_v
!
                  ai = index_two(a, i, n_active_v)                 
!
                  if (ai .ge. bj) then
!
                     ab = index_two(a, b, n_active_v)
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj + offset, 1) = wf%omega2(aibj + offset, 1) + omega_ab_ij(ab, ij)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ab_ij,n_active_v**2,n_active_o**2)  
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   subroutine omega_ccsd_c2_ccsd(wf)
!!
!!    Omega C2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!    
!!    Omega C2 = -1/2* sum_(ck)t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!!                                  - sum_(ck) t_bk_ci (g_kj_ac-sum_(dl)t_al_dj*g_kd_lc)
!!    
      implicit none
!
      class(ccsd) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_ki_J 
      real(dp), dimension(:,:), allocatable :: L_ca_J 
      real(dp), dimension(:,:), allocatable :: L_ac_J 
      real(dp), dimension(:,:), allocatable :: g_kd_lc
      real(dp), dimension(:,:), allocatable :: g_dl_ck
      real(dp), dimension(:,:), allocatable :: g_ki_ca
      real(dp), dimension(:,:), allocatable :: g_ai_ck
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_ai_dl
      real(dp), dimension(:,:), allocatable :: t_ck_bj
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:), allocatable :: X_ai_ck
      real(dp), dimension(:,:), allocatable :: Y_ai_bj
!  
!     Indices
!     
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ca = 0
      integer(i15) :: ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, cj = 0, ck = 0, cl = 0, di = 0, dk = 0, dl = 0
      integer(i15) :: kd = 0, lc = 0, ca = 0, ac = 0
      integer(i15) :: ki = 0
!
      integer(i15) :: aldi = 0, aibj = 0, cldk = 0, bkcj = 0
!
!     Batching and memory handling
!
      integer(i15) :: required = 0, available = 0
!
      integer(i15) :: n_batch = 0, max_batch_length = 0
      integer(i15) :: a_batch = 0, a_start = 0, a_end = 0, a_length = 0 
!
!     Logical for reordering L_ab_J when batching over last index 
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
      integer(i15) :: offset ! first active virtual index
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
!     Allocate L_ia_J
!
      call allocator(L_ia_J,(n_active_o)*(n_active_v),(wf%n_J))
!
!     Get L_ia_J
!
      call wf%get_cholesky_ia(L_ia_J, first_active_o, last_active_o, first_active_v, last_active_v)
!
!     Create g_kd_lc = sum_J L_kd_J * L_lc_J
!
      call allocator(g_kd_lc, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      call dgemm('N','T',            &
                  (n_active_o)*(n_active_v), &
                  (n_active_o)*(n_active_v), &
                  wf%n_J,            &   
                  one,               &
                  L_ia_J,            &
                  (n_active_o)*(n_active_v), &
                  L_ia_J,            &
                  (n_active_o)*(n_active_v), &
                  zero,              &
                  g_kd_lc,           &
                  (n_active_o)*(n_active_v))
!
!     Deallocate L_ia_J
!
      call deallocator(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!     Reorder g_kd_lc as g_dl_ck and t_al_di as t_ai_dl
!
      call allocator(g_dl_ck, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
      call allocator(t_ai_dl, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      do k = 1, n_active_o        
         do d = 1, n_active_v
!
            kd = index_two(k, d, n_active_o)
            dk = index_two(d, k, n_active_v)
!
            do l = 1, n_active_o
!
               dl = index_two(d, l, n_active_v)
!
               do c = 1, n_active_v               
!  
                  lc = index_two(l, c, n_active_o)
                  ck = index_two(c, k, n_active_v)
                  cl = index_two(c, l, n_active_v)
!
                  cldk = index_packed(cl, dk)
!
                  g_dl_ck(dl, ck) = g_kd_lc(kd, lc)
                  t_ai_dl(ck, dl) = wf%t2am(cldk + offset, 1)
!
               enddo
            enddo
         enddo
      enddo

      call deallocator(g_kd_lc, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
!     -1/2*sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call allocator(X_ai_ck, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      call dgemm('N','N',            &
                  (n_active_o)*(n_active_v), &
                  (n_active_o)*(n_active_v), &
                  (n_active_o)*(n_active_v), &
                  -half,             &
                  t_ai_dl,           &
                  (n_active_o)*(n_active_v), &
                  g_dl_ck,           &
                  (n_active_o)*(n_active_v), &
                  zero,              &
                  X_ai_ck,           &
                  (n_active_o)*(n_active_v))
!
!     Deallocate L_ia_J and g_dl_ck, 
!
      call deallocator(g_dl_ck, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
      call deallocator(t_ai_dl, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
!     Constructing g_ki_ac ordered as g_ki_ca
!
!     Allocate g_ki_ca
!
      call allocator(g_ki_ca, (n_active_o)**2, (n_active_v)**2)
      g_ki_ca = zero
!
!     Allocate L_ki_J
!
      call allocator(L_ki_J, (n_active_o)**2, wf%n_J) 
!
!     Get cholesky vectors of ij-type
!
      call wf%get_cholesky_ij(L_ki_J, first_active_o, last_active_o, first_active_o, last_active_o)
!
!     Prepare batching over a 
!
!     Setup of variables needed for batching
!
      available = get_available()
   !  required = 2*(n_active_v**2)*(wf%n_J) + 2*(n_active_v)*(wf%n_o)*(wf%n_J)
   !  required = 4*required
   !  call num_batch(required, available, max_batch_length, n_batch, wf%n_v)
!
   !  a_start  = 1
   !  a_end    = 0
   !  a_length = 0
!
!  !  Start looping over batches
!
   !  do a_batch = 1,n_batch
!
!  !     Get batch limits  and  length of batch
!
   !     call batch_limits(a_start, a_end, a_batch, max_batch_length, wf%n_v)
   !     a_length = a_end - a_start + 1
!
!  !     Get ab-cholesky vectors for the batch, L_ac^J, then reorder from L_ac_J to L_ca_J
!
   !     call allocator(L_ac_J,(wf%n_v)*a_length, wf%n_J)
   !     L_ac_J = zero
   !     call wf%get_cholesky_ab(L_ac_J, a_start, a_end, 1, wf%n_v)
!
   !     call allocator(L_ca_J,(wf%n_v)*a_length, wf%n_J)
   !     L_ca_J = zero
   !     do a = 1, a_length
   !        do c = 1, wf%n_v
   !           ac = index_two(a, c, a_length) 
   !           ca = index_two(c, a, wf%n_v) 
   !           do J = 1, wf%n_J
   !             L_ca_J(ca, J) = L_ac_J(ac, J)
   !           enddo
   !        enddo
   !     enddo

   !     call deallocator(L_ac_J,(wf%n_v)*a_length, wf%n_J)
!
!  !     g_ki_ca = sum_J L_ki_J*L_ca_J
!
   !     call dgemm('N','T',                                   &
   !                 (wf%n_o)*(wf%n_o),                        &
   !                 (wf%n_v)*a_length,                        &
   !                 wf%n_J,                                   &   
   !                 one,                                      &   
   !                 L_ki_J,                                   &
   !                 (wf%n_o)*(wf%n_o),                        &
   !                 L_ca_J,                                   &
   !                 (wf%n_v)*a_length,                        &
   !                 one,                                      &
   !                 g_ki_ca(1,index_two(1, a_start, wf%n_v)), &
   !                 (wf%n_o)*(wf%n_o))
!
!  !     Deallocate L_ca_J
!
   !     call deallocator(L_ca_J, (wf%n_v)*a_length, wf%n_J)
!
   !  enddo ! End of batching
!
!  !  Deallocate L_ki_J
!
   !  call deallocator(L_ki_J, (wf%n_o)*(wf%n_o), wf%n_J)
!
!  !  Reorder g_ki_ca to g_ai_ck
!
   !  call allocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   !  do i = 1, wf%n_o
   !     do k = 1, wf%n_o
!
   !        ki = index_two(k, i, wf%n_o)
!
   !        do a = 1, wf%n_v
!
   !           ai = index_two(a, i, wf%n_v)
!
   !           do c = 1, wf%n_v
!
   !              ca = index_two(c, a, wf%n_v)
   !              ck = index_two(c, k, wf%n_v)
!
   !              g_ai_ck(ai, ck) = g_ki_ca(ki, ca)
!
   !           enddo
   !        enddo
   !     enddo
   !  enddo
!
   !  call deallocator(g_ki_ca, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!  !  X_ai_ck = X_ai_ck + g_ai_ck
!
   !  call daxpy((wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v), one, g_ai_ck, 1, X_ai_ck, 1)
!
!  !  Deallocate g_ai_kc
!
   !  call deallocator(g_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!  !  Reorder t_bkcj_1 as t_ck_bj
!
   !  call allocator(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   !  do j = 1, wf%n_o
   !     do k = 1, wf%n_o
   !        do b = 1, wf%n_v
!
   !           bk = index_two(b, k, wf%n_v)
   !           bj = index_two(b, j, wf%n_v)
!
   !           do c = 1, wf%n_v
!
   !              cj = index_two(c, j, wf%n_v)
   !              ck = index_two(c, k, wf%n_v)
!
   !              bkcj = index_packed(bk, cj)
!
   !              t_ck_bj(ck, bj) = wf%t2am(bkcj, 1)
!
   !           enddo
   !        enddo
   !     enddo
   !  enddo
!
!  !  Allocate intermediate Y_ai_bj
!
   !  call allocator(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!  !  Y_ai_bj = - sum_(ck) X_ai_ck*t_ck_bj
!
   !  call dgemm('N','N',            &
   !              (wf%n_o)*(wf%n_v), &
   !              (wf%n_o)*(wf%n_v), &
   !              (wf%n_o)*(wf%n_v), &
   !              -one,              &
   !              X_ai_ck,           &
   !              (wf%n_o)*(wf%n_v), &
   !              t_ck_bj,           &
   !              (wf%n_o)*(wf%n_v), &
   !              zero,              &
   !              Y_ai_bj,           &
   !              (wf%n_o)*(wf%n_v))
!
!  !  Deallocate the X intermediate
!
   !  call deallocator(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!  !  Deallocate t_ck_bj
!
   !  call deallocator(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!  !  Omega_aibj,1 = P_ai_bj ( 1/2*Y_ai_bj + Y_aj_bi )
!
   !     do i = 1, wf%n_o
   !        do a = 1, wf%n_v
!
   !           ai = index_two(a, i, wf%n_v)
!
   !           do j = 1, wf%n_o    
   !              do b = 1, wf%n_v
!
   !              bj = index_two(b, j, wf%n_v)
!
   !              if (ai .ge. bj) then
   !                 aj = index_two(a, j, wf%n_v)
   !                 bi = index_two(b, i, wf%n_v)
!
   !                 aibj=index_packed(ai, bj)
!
   !                 wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + half*Y_ai_bj(ai, bj) + Y_ai_bj(aj, bi) &
   !                                                           + half*Y_ai_bj(bj, ai) + Y_ai_bj(bi, aj)
!
   !              endif
!  
   !           enddo
   !        enddo
   !     enddo
   !  enddo
!
!  !  Deallocate intermediate Y_ai_bj
!
   !  call deallocator(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_c2_ccsd
   subroutine get_mlccsd_s2am_mlccsd(wf, s_ia_jb)

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
!           :: Change to local basis ::
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