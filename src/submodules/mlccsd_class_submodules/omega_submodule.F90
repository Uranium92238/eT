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
  module subroutine construct_omega_mlccsd(wf)
!! 
!!    Construct Omega (MLCCSD)
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
      class(mlccsd) :: wf
!     
!     Looping variables
!
      integer(i15) :: active_space
!
!     Timing variables
!
      real(dp) :: omega_start = zero
      real(dp) :: omega_end   = zero
!
      real(dp), dimension(:,:), allocatable :: x_IA_JB
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
!
!     Start timing of omega
!
      call cpu_time(omega_start)
!
!     Set the omega vector to zero 
!
      wf%omega1 = zero
      wf%omega2 = zero
!
!     :: Calculate CCS omega contributions ::
!
      call wf%omega_ccs_a1
!
!     OBS! T2 could be deleted here to be allocated and read from file
!           at the end of this routine for memory savings.
!
!     :: Calculate CCSD omega contributions ::
!
!     Construct x2 amplitude: x^AB_IJ = s^AB_IJ + t^AB_IJ  
!     s amplitudes are CC2 amplitudes, t amplitudes are CCSD corrections.      
!     x is given in entire CC2 space (upper case letters), however
!     t is zero for external/semi-external and s is zero for internal.
!
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
      call allocator(x_IA_JB, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
      call wf%get_mlccsd_x2am(x_IA_JB)
!
!     Omega 1 contributions
!
      call wf%omega_mlccsd_a1(x_IA_JB)
      call wf%omega_mlccsd_b1(x_IA_JB)
!
!     Omega 2 contributions
!
      call wf%omega_mlccsd_a2(x_IA_JB)
      call wf%omega_mlccsd_b2(x_IA_JB)
      call wf%omega_mlccsd_c2(x_IA_JB)
      call wf%omega_mlccsd_d2(x_IA_JB)
      call wf%omega_mlccsd_e2(x_IA_JB)
!
      call deallocator(x_IA_JB, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      call cpu_time(omega_end)
      if (timings) write(unit_output,*)'Time in omega:', omega_end-omega_start    
!
   end subroutine construct_omega_mlccsd
!
!
   module subroutine get_mlccsd_x2am_mlccsd(wf, x_ia_jb)
!!
!!    Constructs x_ia_jb amplitudes for current active space 
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:) :: x_ia_jb
!
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: s_ai_bj_CC2
      real(dp), dimension(:,:), allocatable :: I1_a_IBJ ! Intermediato of transformation
      real(dp), dimension(:,:), allocatable :: I2_aIB_j ! Intermediato of transformation
      real(dp), dimension(:,:), allocatable :: I3_Ia_jB ! Intermediato of transformation
      real(dp), dimension(:,:), allocatable :: I4_i_ajB ! Intermediato of transformation
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, ia = 0, bj = 0, ai = 0, aIB = 0, jB = 0, aibj = 0
!
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: first_CCSD_o ! first active occupied index 
      integer(i15) :: first_CCSD_v ! first active virtual index
!
      integer(i15) :: last_CC2_o ! first active occupied index 
      integer(i15) :: last_CC2_v ! first active virtual index
!
!     Calculate first/last indeces
!
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)

      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      if(wf%mlcc_settings%CC2) then
!
!        Construct s2-amplitudes in CC2/CCS block diagonal basis
!
         call allocator(L_ai_J, n_CC2_o*n_CC2_v, wf%n_J)
         L_ai_J = zero
!
         call wf%get_cholesky_ai_for_cc2_amplitudes(L_ai_J, first_CC2_v, last_CC2_v, first_CC2_o, last_CC2_o)
!
         call allocator(g_ai_bj, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        g_ib_jc = g_bi,cj = sum_J L_bj^J*L_ci^J
!
         call dgemm('N', 'T',       &
              (n_CC2_o)*(n_CC2_v),  &
              (n_CC2_o)*(n_CC2_v),  &
              (wf%n_J),             &
              one,                  &
              L_ai_J,               &
              (n_CC2_o)*(n_CC2_v),  &
              L_ai_J,               &
              (n_CC2_o)*(n_CC2_v),  &
              zero,                 &
              g_ai_bj,              &
              (n_CC2_o)*(n_CC2_v)) 
!
         call deallocator(L_ai_J, n_CC2_o*n_CC2_v, wf%n_J)
!
         call allocator(s_AI_BJ_CC2, n_CC2_o*n_CC2_v, n_CC2_o*n_CC2_v)
!
         do A = 1, n_CC2_v
            do I = 1, n_CC2_o
!
               AI = index_two(A, I, n_CC2_v)
!
               do B = 1, n_CC2_v
                  do J = 1, n_CC2_o
!
                     BJ = index_two(B, J, n_CC2_v)
!
                     s_AI_BJ_CC2(AI, BJ) = g_ai_bj(AI, BJ)/(wf%fock_diagonal_cc2_ccs(I + first_CC2_o - 1,1)&
                                        + wf%fock_diagonal_cc2_ccs(J + first_CC2_o - 1,1) &
                                        - wf%fock_diagonal_cc2_ccs(wf%n_o + B + first_CC2_v - 1,1) &
                                        - wf%fock_diagonal_cc2_ccs(wf%n_o + A + first_CC2_v - 1,1))
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ai_bj, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!        Transform to CCSD/CC2/CCS block diagonal basis (Cholesky or cnto)
!
!        Transform index A : X1_a_IBJ = sum_A T_aA*s_A_IBJ_CC2
!
         call allocator(I1_a_IBJ, n_CC2_v, (n_CC2_o**2)*n_CC2_v)
!
         call dgemm('N','N', &
                  n_CC2_v,&
                  (n_CC2_o**2)*n_CC2_v, &
                  n_CC2_v, &
                  one, &
                  wf%T_v(first_CC2_v, first_CC2_v), &
                  wf%n_v, &
                  s_AI_BJ_CC2, &
                  n_CC2_v, &
                  zero, &
                  I1_a_IBJ, &
                  n_CC2_v)
!
         call deallocator(s_AI_BJ_CC2, n_CC2_o*n_CC2_v, n_CC2_o*n_CC2_v)
!
!        Transform index J : X2_aIB_j = sum_J X1_aIB_J*T_jJ 
!
         call allocator(I2_aIB_j, (n_CC2_v**2)*n_CC2_o, n_CC2_o)
!
         call dgemm('N', 'T', &
                  (n_CC2_v**2)*n_CC2_o, &
                  n_CC2_o, &
                  n_CC2_o, &
                  one, &
                  I1_a_IBJ, &
                  (n_CC2_v**2)*n_CC2_o, &
                  wf%T_o(first_CC2_o, first_CC2_o), &
                  wf%n_o, &
                  zero, &
                  I2_aIB_j, &
                  (n_CC2_v**2)*n_CC2_o)
!
         call deallocator(I1_A_IBJ, n_CC2_v, (n_CC2_o**2)*n_CC2_v)
!
!        Reorder X2_aIB_j to X3_Ia_jB
!
         call allocator(I3_Ia_jB, (n_CC2_o)*(n_CC2_v), (n_CC2_v)*(n_CC2_o))
         I3_Ia_jB = zero
!
         do a = 1, n_CC2_v
           do j = 1, n_CC2_o
              do B = 1, n_CC2_v
!
                 jB = index_two(j, B, n_CC2_o)
!
                 do I = 1, n_CC2_o
!
                    aIB = index_three(a, I, B, n_CC2_v, n_CC2_o)
                    Ia = index_two(I, a, n_CC2_o)
!
                    I3_Ia_jB(Ia, jB) = I2_aIB_j(aIB, j)
!
                 enddo
              enddo
           enddo
         enddo
!
         call deallocator(I2_aIB_j, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
!
!        Transform index I : I4_i_ajB = sum_I T_iI*I3_Ia_jB 
!
         call allocator(I4_i_ajB, (n_CC2_o), (n_CC2_o)*(n_CC2_v)**2)
!
         call dgemm('N', 'N', &
                  n_CC2_o, &
                  (n_CC2_o)*(n_CC2_v**2), &
                  n_CC2_o, &
                  one, &
                  wf%T_o(first_CC2_o, first_CC2_o), &
                  wf%n_o, &
                  I3_Ia_jB, &
                  n_CC2_o, &
                  zero, &
                  I4_i_ajB, &
                  n_CC2_o)
!
         call deallocator(I3_Ia_jB, (n_CC2_o)*(n_CC2_v), (n_CC2_v)*(n_CC2_o))
!
!        Transform index B : s_ia_jb = sum_B I4_i_ajB*T_bB
!
         call dgemm('N', 'T', &
                  (n_CC2_o**2)*(n_CC2_v), &
                  n_CC2_v, &
                  n_CC2_v, &
                  one, &
                  I4_i_ajB, &
                  (n_CC2_o**2)*(n_CC2_v), &
                  wf%T_v(first_CC2_v, first_CC2_v), &
                  wf%n_v, &
                  zero, &
                  x_ia_jb, &
                  (n_CC2_o**2)*(n_CC2_v))
!
         call deallocator(I4_i_ajB, n_CC2_o, (n_CC2_o)*(n_CC2_v**2))
!
!        Replace internals with CCSD amplitudes
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
                  do b = 1, n_CCSD_v
!
                     ia = index_two(i, a, n_CC2_o)
                     jb = index_two(j, b, n_CC2_o)
                     ai = index_two(a, i, n_CCSD_v)
                     bj = index_two(b, j, n_CCSD_v)
!
                     aibj = index_packed(ai, bj)
!
                     x_ia_jb(ia, jb) =  wf%t2am(aibj, 1)
!
                  enddo
               enddo
            enddo
         enddo
!
      else ! No CC2 amplitudes
!
         do i = 1, n_CCSD_o
            do j = 1, n_CCSD_o
               do a = 1, n_CCSD_v
                  do b = 1, n_CCSD_v
!
                     ia = index_two(i, a, n_CC2_o)
                     jb = index_two(j, b, n_CC2_o)
                     ai = index_two(a, i, n_CCSD_v)
                     bj = index_two(b, j, n_CCSD_v)
!
                     aibj = index_packed(ai, bj)
!
                     x_ia_jb(ia, jb) = wf%t2am(aibj, 1)
!
                  enddo
               enddo
            enddo
         enddo
      endif
!
   end subroutine get_mlccsd_x2am_mlccsd
!
!
    module subroutine omega_mlccsd_a1_mlccsd(wf, x_ib_jc)
!! 
!!    Omega A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!    Calculates the A1 term of omega for the active space, 
!!   
!!    A1: sum_bcj g_Abjc * u_ij^bc,
!!
!!    where upper case letters indicate CCS space, and 
!!    lower case letters are the combined CC2/CCSD spaces.
!!  
!!    A1 is added to the projection vector (omega1) of
!!    the wavefunction object wf.
!! 
!!    u_ij^bc = 2*x_ij^bc - x_ij^cb 
!!
!!    Batching over A.
!!
!! 
      implicit none
!
      class(mlccsd)            :: wf
      real(dp), dimension(:,:) :: x_ib_jc
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
      real(dp), dimension(:,:), allocatable :: u_bjc_i
      real(dp), dimension(:,:), allocatable :: g_Ab_jc
      real(dp), dimension(:,:), allocatable :: L_Ab_J
      real(dp), dimension(:,:), allocatable :: L_jc_J
!
      logical :: reorder  ! To get L_ab_J reordered, for batching over a
!
!     Active space variables
!
      integer(i15) :: first_CC2_o   ! first active occupied index 
      integer(i15) :: first_CC2_v   ! first active virtual index
      integer(i15) :: last_CC2_o    ! last active occupied index 
      integer(i15) :: last_CC2_v    ! last active virtual index
      integer(i15) :: n_CC2_o       ! number of active occupied index 
      integer(i15) :: n_CC2_v       ! number of active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
!     :: Construct u_ib_jc ::
!
!     u_ij^bc = 2*s_ij^bc - s_ij^cb =  (2*g_ij^bc - g_ij^cb)/ε_ij^cb
!
!
      call allocator(u_bjc_i, (n_CC2_v**2)*(n_CC2_o), (n_CC2_o))
!
      do b = 1, n_CC2_v
         do i = 1, n_CC2_o
!
            ib = index_two(i, b, n_CC2_o)
!
            do c = 1, n_CC2_v
               do j = 1, n_CC2_o
!
                  jc = index_two(j, c, n_CC2_o)
                  jb = index_two(j, b, n_CC2_o)
                  ic = index_two(i, c, n_CC2_o)
!
                  bjc = index_three(b, j, c, n_CC2_v, n_CC2_o)
!
                  u_bjc_i(bjc,i) = (two*x_ib_jc(ib,jc)-x_ib_jc(ic, jb))
!
               enddo
            enddo
         enddo
      enddo
!
!     Prepare for batching over A
!
      required = (n_CC2_v)*(wf%n_v)*(wf%n_J) &
               + (n_CC2_v)*(n_CC2_o)*(wf%n_J) &
               + ((n_CC2_v)**2)*(wf%n_o)*(n_CC2_o) &
               + ((n_CC2_o)**2)*((n_CC2_v)**2)
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
         call allocator(L_jc_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
         L_jc_J = zero
         call wf%get_cholesky_ia(L_jc_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
         call allocator(L_Ab_J, (n_CC2_v)*a_length, wf%n_J) 
         L_Ab_J = zero
!
         call wf%get_cholesky_ab(L_Ab_J, a_first, a_last, first_CC2_v, last_CC2_v)
!
         call allocator(g_Ab_jc, (n_CC2_v)*a_length, (n_CC2_o)*(n_CC2_v))      
!
!        g_Ab,jc = sum_J L_Ab^J*L_jc^J
!
         call dgemm('N', 'T',             &
                     (n_CC2_v)*a_length,  &
                     (n_CC2_o)*(n_CC2_v), &
                     (wf%n_J),            &
                     one,                 &
                     L_Ab_J,              &
                     (n_CC2_v)*a_length,  &
                     L_jc_J,              &
                     (n_CC2_o)*(n_CC2_v), &
                     zero,                &
                     g_Ab_jc,             &
                     (n_CC2_v)*a_length) 
!
         call deallocator(L_jc_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
         call deallocator(L_Ab_J, (n_CC2_v)*a_length, wf%n_J) 
!
!        :: Add contributions to omega ::
!
         call dgemm('N', 'N',                &
                     A_length,               &
                     n_CC2_o,                &
                     (n_CC2_o)*(n_CC2_v**2), &
                     one,                    &
                     g_Ab_jc,                &
                     A_length,               &
                     u_bjc_i,                &
                     (n_CC2_o)*(n_CC2_v**2), &
                     one,                    &
                     wf%omega1(A_first,1),   &
                     wf%n_v)
! 
         call deallocator(g_Ab_jc, (n_CC2_v)*a_length, (n_CC2_o)*(n_CC2_v))      
!
      enddo ! Batching over a
!
      call deallocator(u_bjc_i, (n_CC2_v**2)*(n_CC2_o), (n_CC2_o))
!
   end subroutine omega_mlccsd_a1_mlccsd
!
!
    module subroutine omega_mlccsd_b1_mlccsd(wf, x_ja_kb)
!! 
!!    Omega B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!    Calculates the B1 term of omega, 
!!   
!!    B1: - sum_bjk u_jk^ab*g_kbjI + sum_bj u_ij^ab F_jb,
!!
!!    with u_ij^ab = 2*x_ij^ab - x_ij^ba. 
!!
      implicit none
!
      class(mlccsd)            :: wf
      real(dp), dimension(:,:) :: x_ja_kb  
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
      integer(i15) :: ja = 0, kb = 0, ka  = 0, jb = 0, aj = 0
      integer(i15) :: kbj = 0, jbi = 0
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: last_CC2_o ! first active occupied index 
      integer(i15) :: last_CC2_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1
!
!     :: Construct u_jk^ab ::
!  
!     u_jk^ab = 2*s_jk^ab - s_jk^ba  (place in u_a_jkb)        
!  
!
      call allocator(u_a_kbj, n_CC2_v, (n_CC2_o**2)*n_CC2_v)

      do k = 1, n_CC2_o
         do b = 1, n_CC2_v         
!  
             kb = index_two(k, b, n_CC2_o)
!  
            do j = 1, n_CC2_o
!
               jb = index_two(j, b, n_CC2_o)
               kbj = index_three(k, b, j, n_CC2_o, n_CC2_v)
!
               do a = 1, n_CC2_v             
!  
                  ja = index_two(j, a, n_CC2_o)
                  ka = index_two(k, a, n_CC2_o)
                  
!  
                  u_a_kbj(a,kbj) = (two*x_ja_kb(ja,kb)-x_ja_kb(jb, ka))
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
      call allocator(L_jI_J, n_CC2_o*(wf%n_o), wf%n_J)
!
      call wf%get_cholesky_ij(L_jI_J, first_CC2_o, last_CC2_o, 1, wf%n_o)
!
      call allocator(L_kb_J, n_CC2_o*n_CC2_v, wf%n_J)
!
      call wf%get_cholesky_ia(L_kb_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
      call allocator(g_kb_jI, n_CC2_o*n_CC2_v, n_CC2_o*(wf%n_o) )
!
!     g_kb_jI = sum_J L_kb^J*L_jI_J
!
      call dgemm('N', 'T',             &
                  (n_CC2_o)*(n_CC2_v), &
                  (n_CC2_o)*(wf%n_o),  &
                  (wf%n_J),            &
                  one,                 &
                  L_kb_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  L_ji_J,              &
                  (n_CC2_o)*(wf%n_o),  &
                  zero,                &
                  g_kb_jI,             &
                  (n_CC2_o)*(n_CC2_v)) 
!
      call deallocator(L_jI_J, n_CC2_o*(wf%n_o), wf%n_J)
      call deallocator(L_kb_J, n_CC2_o*n_CC2_v, wf%n_J)
!
!     Add contributions to omega
!
      call dgemm('N', 'N',                   &
                  n_CC2_v,                   &
                  (wf%n_o),                  &
                  n_CC2_v*((n_CC2_o)**2),    &
                  -one,                      &
                  u_a_kbj,                   &
                  n_CC2_v,                   &
                  g_kb_jI,                   &
                  n_CC2_v*((n_CC2_o)**2),    &
                  one,                       &
                  wf%omega1(first_CC2_v, 1),  &
                  (wf%n_v))
!
      call deallocator(g_kb_jI, n_CC2_o*n_CC2_v, n_CC2_o*(wf%n_o))
!
!     :: sum_jb F_jb u_ij^ab ::
!
      do i = 1, n_CC2_o
!
         I_full = i + first_CC2_o - 1
!
         do a = 1, n_CC2_v
!
            A_full = a + first_CC2_v - 1 
!       
            do j = 1, n_CC2_o
!
               J_full = j + first_CC2_o - 1
!
               do b = 1, n_CC2_v
! 
                  B_full = b + first_CC2_v - 1
! 
                  jbi = index_three(j, b, i, n_CC2_o, n_CC2_v)
                 
                  wf%omega1(A_full, I_full) = wf%omega1(A_full, I_full) + u_a_kbj(a, jbi)*wf%fock_ia(J_full, B_full)
!
               enddo
            enddo
         enddo
      enddo
! 
      call deallocator(u_a_kbj, n_CC2_v, (n_CC2_o**2)*n_CC2_v)
!      
   end subroutine omega_mlccsd_b1_mlccsd
!
!
  module subroutine omega_mlccsd_a2_mlccsd(wf, x_IC_JD)
!
!     Omega A2 term: Omega A2 = sum_(cd)g_aC_bD * x_Ci_Dj
!
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!     Structure: Batching over both a and b for A2.2.
!                x^+_Ci_Dj = x_Ci_Dj + x_Di_Cj
!                x^-_Ci_Dj = x_Ci_Dj - x_Di_Cj
!                g^+_aC_bD = g_aC_bD + g_bC_aD 
!                g^-_aC_bD = g_aC_bD - g_bC_aD 
! 
!                omega_A2_ai_bj = 1/4*(g^+_aC_bD*x^+_Ci_Dj + g^-_aC_bD*x^-_Ci_Dj)
!                omega_A2_aj_bi = 1/4*(g^+_aC_bD*x^+_Ci_Dj - g^-_aC_bD*x^-_Ci_Dj)
!
      implicit none
!
      class(mlccsd)  :: wf
      real(dp), dimension(:,:) :: x_IC_JD
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_aC_bD 
      real(dp), dimension(:,:), allocatable :: g_ai_bj 
      real(dp), dimension(:,:), allocatable :: g_p_ab_CD
      real(dp), dimension(:,:), allocatable :: g_m_ab_CD
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_Ca_J 
      real(dp), dimension(:,:), allocatable :: L_aC_J 
      real(dp), dimension(:,:), allocatable :: L_Db_J 
      real(dp), dimension(:,:), allocatable :: L_bD_J 
!
!     Reordered T2 amplitudes
!
      
      real(dp), dimension(:,:), allocatable :: x_p_CD_ij
      real(dp), dimension(:,:), allocatable :: x_m_CD_ij
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
      integer(i15) :: ab = 0, ca = 0, cb = 0, cd = 0, da = 0, db = 0, ac = 0, bd = 0, ad = 0, bc = 0
      integer(i15) :: ai = 0, aj = 0, bj = 0, bi = 0, ic = 0, jc = 0, jd = 0, id = 0
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
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o ! first active occupied index 
      integer(i15) :: first_CC2_v ! first active virtual index
      integer(i15) :: first_CCSD_o ! first active occupied index 
      integer(i15) :: first_CCSD_v ! first active virtual index
!
      integer(i15) :: last_CC2_o ! first active occupied index 
      integer(i15) :: last_CC2_v ! first active virtual index
!
      integer(i15) :: last_CCSD_o ! first active occupied index 
      integer(i15) :: last_CCSD_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
!
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!
!     ::  Calculate the A2.1 term of omega ::
!
!     Create g_ai_bj
!  
      call allocator(g_ai_bj, n_CCSD_o*n_CCSD_v, n_CCSD_o*n_CCSD_v)
      call allocator(L_ai_J, n_CCSD_o*n_CCSD_v, wf%n_J)
!
      call wf%get_cholesky_ai(L_ai_J, first_CCSD_v, last_CCSD_v, first_CCSD_o, last_CCSD_o)
!
!     g_ai_bj = sum_J L_ai_J*L_bj_J
!     
      call dgemm('N','T',            &
                  n_CCSD_o*n_CCSD_v, &
                  n_CCSD_o*n_CCSD_v, &
                  wf%n_J,            &
                  one,               & 
                  L_ai_J,            &
                  n_CCSD_o*n_CCSD_v, &
                  L_ai_J,            &
                  n_CCSD_o*n_CCSD_v, &
                  zero,              &
                  g_ai_bj,           &
                  n_CCSD_o*n_CCSD_v)
!
!
      call deallocator(L_ai_J, n_CCSD_o*n_CCSD_v, wf%n_J)
!
!     Add A2.1 to Omega 2
!     
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CCSD_v)
!
            do j = 1, n_CCSD_o
               do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CCSD_v)
!
                  if(ai .ge. bj) then
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + g_ai_bj(ai, bj)
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ai_bj, n_CCSD_o*n_CCSD_v, n_CCSD_o*n_CCSD_v)
!
!     ::  Calculate the A2.2 term  of omega ::
!
!     Calc of required should be fixed, is now conservative.
!
      required = max(3*((n_CC2_v)**2)*(wf%n_J) + 2*(n_CC2_v)*(n_CC2_o)*(wf%n_J),  &      ! Needed to get  L_db_J
                     (n_CC2_v)**4 + 2*(n_CC2_v)**2*(wf%n_J), &                         ! Needed to get g_ac_bd
                     (n_CC2_v)**4 + 2*((n_CC2_v)**2)*(packed_size(n_CC2_v))      &  ! Needed to get g+- and t+-
                     + 2*(packed_size(n_CC2_v))*(packed_size(n_CC2_o)), &              !
                       2*(packed_size(n_CC2_o))*(packed_size(n_CC2_v)) &               ! Needed for g+- and t+- and Omega+-
                     + 2*(packed_size(n_CC2_v))*(packed_size(n_CC2_o)) &               !
                     + 2*(n_CC2_v)**2*(packed_size(n_CC2_v)))                          !
!
      required = required*4  ! Words
      write(unit_output,*)'Required in A2', required

      available=get_available()
!
      a_max_length = 0
      call num_two_batch(required, available, a_max_length, a_n_batch, n_CCSD_v)
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
         call batch_limits(a_first, a_last ,a_batch, a_max_length, n_CCSD_v)
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
            call batch_limits(b_first ,b_last ,b_batch, b_max_length, n_CCSD_v)
            b_length = b_last - b_first + 1 
!
!           Get ab-cholesky vectors for the batch, L_ac^J, then reorder from L_ac_J to L_ca_J
!
            call allocator(L_aC_J, (n_CC2_v)*a_length, wf%n_J)
            L_aC_J = zero
!
            call wf%get_cholesky_ab(L_aC_J, a_first, a_last, first_CC2_v, last_CC2_v)
!
!          Get ab-cholesky vectors for the batch, L_bd^J, then reorder from L_bd_J to L_db_J
!
           call allocator(L_bD_J, (n_CC2_v)*b_length, wf%n_J)
           L_bD_J = zero
!
           call wf%get_cholesky_ab(L_bD_J, b_first, b_last, first_CC2_v, last_CC2_v)         
! 
!          Allocate g_ca_db
!
           call allocator(g_aC_bD, (n_CC2_v)*a_length, (n_CC2_v)*b_length)
           g_aC_bD = zero
!
!          g_ca_db = sum_J L_ca_J*L_db_J
!     
           call dgemm('N','T',                  &
                       (n_CC2_v)*a_length,      &
                       (n_CC2_v)*b_length,      &
                       wf%n_J,                  &
                       one,                     &
                       L_aC_J,                  &
                       (n_CC2_v)*a_length,      &
                       L_bD_J,                  &
                       (n_CC2_v)*b_length,      &
                       zero,                    &
                       g_aC_bD,                 &
                       (n_CC2_v)*a_length)
!
            call deallocator(L_aC_J, (n_CC2_v)*a_length, wf%n_J)
            call deallocator(L_bD_J, (n_CC2_v)*b_length, wf%n_J) 
!
            if (b_batch .eq. a_batch) then
!
!
!           Allocate for +-g, +-t
!
               call allocator(g_p_ab_CD, packed_size(a_length), packed_size(n_CC2_v))
               call allocator(g_m_ab_CD, packed_size(a_length), packed_size(n_CC2_v))
               call allocator(x_p_CD_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
               call allocator(x_m_CD_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
!
               g_p_ab_CD = zero
               g_m_ab_CD = zero
               x_p_CD_ij = zero
               x_m_CD_ij = zero
!
!              Reorder g_ca_db to g_ab_cd and t_ci_dj to t_cd_ij
! 
               do C = 1, n_CC2_v 
                  do D = 1, C
!
                     CD = index_packed(c, d)
!
                     do a = 1, a_length
!
                        ac = index_two(a, C, a_length)
                        ad = index_two(a, D, a_length)
!
                        do  b = 1, b_length
                           if ((a+a_first-1) .ge. (b+b_first-1)) then
!                            
                              bD = index_two(b, D, b_length)
                              bC = index_two(b, C, b_length)
!
                              ab = index_packed(a, b)
! 
                              g_p_ab_cd(ab, cd) = g_aC_bD(aC, bD) + g_aC_bD(aD, bC)
                              g_m_ab_cd(ab, cd) = g_aC_bD(aC, bD) - g_aC_bD(aD, bC)
!
                             if(C .ne. D) then
                               g_p_ab_CD(ab, CD) = two*g_p_ab_CD(ab, CD)
                               g_m_ab_CD(ab, CD) = two*g_m_ab_CD(ab, CD)
                             endif
!                             
                           endif
                        enddo
                     enddo
!
                    do i = 1, n_CCSD_o
                       do j = 1, i
!    
                          ij = index_packed(i, j)
!    
                          iC = index_two(i, C, n_CC2_o)
                          jD = index_two(j, D, n_CC2_o)
                          jC = index_two(j, C, n_CC2_o)
                          iD = index_two(i, D, n_CC2_o)
! 
                          x_p_CD_ij(CD, ij) = x_IC_JD(iC, jD) + x_IC_JD(iD, jC)
                          x_m_CD_ij(CD, ij) = x_IC_JD(iC, jD) - x_IC_JD(iD, jC)  
!
                       enddo
                    enddo
                 enddo
              enddo
!
!              Dellocate g_ac_bd 
!
               call deallocator(g_aC_bD, (n_CC2_v)*a_length, (n_CC2_v)*b_length)
!
!              Allocate omega +-
!
              call allocator(omega2_p_ab_ij, packed_size(a_length), packed_size(n_CCSD_o))
              call allocator(omega2_m_ab_ij, packed_size(a_length), packed_size(n_CCSD_o))
! 
!              omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
! 
              call dgemm('N','N',                & 
                          packed_size(a_length), &
                          packed_size(n_CCSD_o), &
                          packed_size(n_CC2_v),  &
                          one/four,              &
                          g_p_ab_CD,             &
                          packed_size(a_length), &
                          x_p_CD_ij,             &
                          packed_size(n_CC2_v),  &
                          zero,                  &
                          omega2_p_ab_ij,        &
                          packed_size(a_length))
!
              call dgemm('N','N',                & 
                          packed_size(a_length), &
                          packed_size(n_CCSD_o), &
                          packed_size(n_CC2_v),  &
                          one/four,              &
                          g_m_ab_CD,             &
                          packed_size(a_length), &
                          x_m_CD_ij,             &
                          packed_size(n_CC2_v),  &
                          zero,                  &
                          omega2_m_ab_ij,        &
                          packed_size(a_length) )
!
!             Deallocate +-g, +-t
! 
              call deallocator(g_p_ab_CD, packed_size(a_length), packed_size(n_CC2_v))
              call deallocator(g_m_ab_CD, packed_size(a_length), packed_size(n_CC2_v))
              call deallocator(x_p_CD_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
              call deallocator(x_m_CD_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
!
              do i = 1, n_CCSD_o
                 do j = 1, i
!
                    
                    ij = index_packed(i, j)
!
                    do a = 1, a_length
!
                       Ai = index_two(a + a_first - 1, i, n_CCSD_v) ! A is full-CCSD-space a index
                       Aj = index_two(a + a_first - 1, j, n_CCSD_v) ! A is full-CCSD-space a index
!
                       do b = 1, b_length
!                
                          if ((a+a_first-1) .ge. (b+b_first-1)) then
                             Bj = index_two(b + b_first - 1, j, n_CCSD_v) ! B is full-CCSD-space b index
                             Bi = index_two(b + b_first - 1, i, n_CCSD_v) ! B is full-CCSD-space b index
!
!
                             ab = index_packed(a, b)
!    
                             AiBj = index_packed(Ai, Bj)
                             BiAj = index_packed(Bi, Aj)
!                         
!                            Reorder into omega2_aibj
! 
                             wf%omega2(AiBj,1) = wf%omega2(AiBj, 1) &
                                                   + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                             if (AiBj .ne. BiAj) then
                                wf%omega2(BiAj,1) = wf%omega2(BiAj, 1) &
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
               call deallocator(omega2_p_ab_ij, packed_size(a_length), packed_size(n_CCSD_o))
               call deallocator(omega2_m_ab_ij, packed_size(a_length), packed_size(n_CCSD_o))
            else
!
!              Allocate for +-g, +-t
!
               call allocator(g_p_ab_CD, a_length*b_length, packed_size(n_CC2_v))
               call allocator(g_m_ab_CD, a_length*b_length, packed_size(n_CC2_v))
               call allocator(x_p_CD_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
               call allocator(x_m_CD_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
!
               g_p_ab_CD = zero
               g_m_ab_CD = zero
               x_p_CD_ij = zero
               x_m_CD_ij = zero
! 
               do C = 1, n_CC2_v 
                  do D = 1, C
!
                     CD = index_packed(C, D)
!
                     do a = 1, a_length
!
                        ac = index_two(a, C, a_length)
                        ad = index_two(a, D, a_length)
!
                        do  b = 1, b_length
!                            
                             bD = index_two(b, D, b_length)
                             bC = index_two(b, C, b_length)
!
                             ab = index_two(a, b, a_length)
! 
                             g_p_ab_CD(ab, CD) = g_aC_bD(aC, bD) + g_aC_bD(aD, bC)
                             g_m_ab_CD(ab, CD) = g_aC_bD(aC, bD) - g_aC_bD(aD, bC)
!
                            if(c .ne. d) then
                              g_p_ab_CD(ab, CD) = two*g_p_ab_CD(ab, CD)
                              g_m_ab_CD(ab, CD) = two*g_m_ab_CD(ab, CD)
                            endif
!                            
                       enddo
                    enddo
!
                    do i = 1, n_CCSD_o
                       do j = 1, i
!    
                          ij = index_packed(i, j)
!    
                          iC = index_two(i, C, n_CC2_o)
                          jD = index_two(j, D, n_CC2_o)
                          jC = index_two(j, C, n_CC2_o)
                          iD = index_two(i, D, n_CC2_o)
! 
                          x_p_CD_ij(CD, ij) = x_IC_JD(iC, jD) + x_IC_JD(iD, jC)
                          x_m_CD_ij(CD, ij) = x_IC_JD(iC, jD) - x_IC_JD(iD, jC)  
!
                       enddo
                    enddo
                 enddo
              enddo
!
!              Dellocate g_ac_bd 
!
               call deallocator(g_aC_bD, (n_CC2_v)*a_length, (n_CC2_v)*b_length)
!
!              Allocate omega +-
!

               call allocator(omega2_p_ab_ij, b_length*a_length, packed_size(n_CCSD_o))
               call allocator(omega2_m_ab_ij, b_length*a_length, packed_size(n_CCSD_o))
!  
!               omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
! 
              call dgemm('N','N',                  & 
                          b_length*a_length,       &
                          packed_size(n_CCSD_o),   &
                          packed_size(n_CC2_v),    &
                          one/four,                &
                          g_p_ab_CD,               &
                          b_length*a_length,       &
                          x_p_CD_ij,               &
                          packed_size(n_CC2_v),    &
                          zero,                    &
                          omega2_p_ab_ij,          &
                          b_length*a_length)
!
              call dgemm('N','N',                  & 
                          b_length*a_length,       &
                          packed_size(n_CCSD_o),   &
                          packed_size(n_CC2_v),    &
                          one/four,                &
                          g_m_ab_CD,               &
                          b_length*a_length,       &
                          x_m_CD_ij,               &
                          packed_size(n_CC2_v),    &
                          zero,                    &
                          omega2_m_ab_ij,          &
                          b_length*a_length)
!
!          Deallocate +-g, +-t
! 
              call deallocator(g_p_ab_cd, b_length*a_length, packed_size(n_CC2_v))
              call deallocator(g_m_ab_cd, b_length*a_length, packed_size(n_CC2_v))
              call deallocator(x_p_cd_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
              call deallocator(x_m_cd_ij, packed_size(n_CC2_v), packed_size(n_CCSD_o))
!
               do i = 1, n_CCSD_o
                  do j = 1, i
!
                     
                     ij = index_packed(i, j)
!
                     do a = 1, a_length
!
                        Ai = index_two(a + a_first - 1, i, n_CCSD_v) ! A is full-space a index
                        Aj = index_two(a + a_first - 1, j, n_CCSD_v) ! A is full-space a index
!
                        do b = 1, b_length
!                 
                              Bj = index_two(b + b_first - 1, j, n_CCSD_v) ! B is full-space b index
                              Bi = index_two(b + b_first - 1, i, n_CCSD_v) ! B is full-space b index
!
!
                              ab = index_two(a, b, a_length)

!     
                              AiBj = index_packed(Ai, Bj)
                              BiAj = index_packed(Bi, Aj)
!                          
!                             Reorder into omega2_aibj
!  
                              wf%omega2(AiBj,1) = wf%omega2(AiBj, 1) &
                                          + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                              if (AiBj .ne. BiAj) then
                                 wf%omega2(BiAj,1) = wf%omega2(BiAj, 1) &
                                          + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
                              endif   
!     
                        enddo
                     enddo
                  enddo
               enddo
!
!              Deallocate omega +-
!
               call deallocator(omega2_p_ab_ij, b_length*a_length, packed_size(n_CCSD_o))
               call deallocator(omega2_m_ab_ij, b_length*a_length, packed_size(n_CCSD_o))
            endif
!
         enddo ! End batching over b
      enddo ! End batching over a
!
   end subroutine omega_mlccsd_a2_mlccsd
!
!
   module subroutine omega_mlccsd_b2_mlccsd(wf, x_kc_ld)
!!
!!    Omega B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 11 Mar 2017
!! 
!!    Omega B2 = sum_(kl) x_ak_bl*(g_kilj + sum_(cd) x_ci_dj * g_kc_ld)
!!
!!    Structure: g_kilj is constructed first and reordered as g_kl_ij. 
!!    Then the contraction over cd is performed, and the results added to g_kl_ij.
!!    x_ka_lb is then reordered as x_ab_kl and the contraction over kl is performed.
!!
      implicit none
!
      class(mlccsd) :: wf 
      real(dp), dimension(:,:) :: x_kc_ld 
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_kc_J     
      real(dp), dimension(:,:), allocatable :: L_Ki_J  
      real(dp), dimension(:,:), allocatable :: g_kc_ld    
      real(dp), dimension(:,:), allocatable :: g_kl_cd     
      real(dp), dimension(:,:), allocatable :: g_Ki_Lj 
      real(dp), dimension(:,:), allocatable :: I_KL_ij 
!
!     Reordered T2 apmlitudes
!   
      real(dp), dimension(:,:), allocatable :: x_cd_ij     
      real(dp), dimension(:,:), allocatable :: x_ab_kl    
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
      integer(i15) :: kc = 0, ld = 0, ka = 0, lb = 0
      integer(i15) :: ij = 0, ki = 0, kl = 0, lj = 0
!
      integer(i15) :: aibj = 0, akbl = 0, cidj = 0 
!
!     Active space variables
!
      integer(i15) :: n_CC2_o  = 0, n_CC2_v  = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o 
      integer(i15) :: first_CC2_v 
      integer(i15) :: first_CCSD_o 
      integer(i15) :: first_CCSD_v 
!
      integer(i15) :: last_CC2_o 
      integer(i15) :: last_CC2_v 
      integer(i15) :: last_CCSD_o  
      integer(i15) :: last_CCSD_v 
!
!     Calculate first/last indeces
! 
!     CC2
!  
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
!     CCSD
!
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1  
!
!     Read Cholesky vector L_Ki_J
!
      call allocator(L_Ki_J, (n_CC2_o)*(n_CCSD_o), wf%n_J)
!
      call wf%get_cholesky_ij(L_Ki_J, first_CC2_o, last_CC2_o, first_CCSD_o, last_CCSD_o)
!
!     Create g_Ki_Lj = sum_J L_Ki_J*L_Lj_J
!
      call allocator(g_Ki_Lj, (n_CC2_o)*(n_CCSD_o), (n_CC2_o)*(n_CCSD_o)) 
!
      call dgemm('N','T',                 &
                  (n_CC2_o)*(n_CCSD_o),   &
                  (n_CC2_o)*(n_CCSD_o),   &
                  wf%n_J,                 &
                  one,                    &
                  L_Ki_J,                 &
                  (n_CC2_o)*(n_CCSD_o),   &
                  L_Ki_J,                 &
                  (n_CC2_o)*(n_CCSD_o),   &
                  zero,                   &
                  g_Ki_Lj,                &
                  (n_CC2_o)*(n_CCSD_o))
!
!
      call deallocator(L_Ki_J, (n_CC2_o)*(n_CCSD_o), wf%n_J)
!
!     Reorder g_Ki_Lj to I_KL_ij
!
      call allocator(I_KL_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
!
      do K = 1, n_CC2_o
         do L = 1, n_CC2_o
            do i = 1, n_CCSD_o
               do j=1, n_CCSD_o
!
!                 Calculate compound indices
!
                  Ki = index_two(K, i, n_CC2_o)
                  Lj = index_two(L, j, n_CC2_o)
                  KL = index_two(K, L, n_CC2_o)
                  ij = index_two(i, j, n_CCSD_o)
!
!                 Reordering g_ki_lj to g_kl_ij
!
                  I_KL_ij(KL, ij) = g_Ki_Lj(Ki, Lj)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_Ki_Lj, (n_CC2_o)*(n_CCSD_o), (n_CC2_o)*(n_CCSD_o))
!
!     Create g_ck_ld = sum_J L_kc_J*L_ld_J
!
!     Read Cholesky vectors L_kc^J
!
      call allocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_KC_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
      call allocator(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!  
      call dgemm('N','T',              &
                  (n_CC2_o)*(n_CC2_v), &
                  (n_CC2_o)*(n_CC2_v), &
                  wf%n_J,              &
                  one,                 &
                  L_kc_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  L_kc_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  zero,                &
                  g_kc_ld,             &
                  (n_CC2_o)*(n_CC2_v))
!
      call deallocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     :: s2 contribution to I_kl_ij intermediate ::
!
      call allocator(x_CD_ij, (n_CC2_v)**2, (n_CCSD_o)**2)
      call allocator(g_KL_CD, (n_CC2_o)**2, (n_CC2_v)**2)
!
      do D = 1, n_CC2_v
         do C= 1, n_CC2_v
!
            CD = index_two(C, D, n_CC2_v)
!
            do L = 1, n_CC2_o
!  
               LD = index_two(L, D, n_CC2_o)
!
               do K = 1, n_CC2_o              
!
                  KL = index_two(K, L, n_CC2_o)
!
                  KC = index_two(K, C, n_CC2_o)
!
                  g_KL_CD(KL, CD) = g_KC_LD(KC, LD)
!
                  if ((K .le. n_CCSD_o) .and. (L .le. n_CCSD_o)) then
!
                     ij = index_two(k, l, n_CCSD_o)
                     x_cd_ij(cd, ij) = x_KC_LD(KC, LD)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_KC_LD, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      call dgemm('N','N',        &
                  (n_CC2_o)**2,  &
                  (n_CCSD_o)**2, &
                  (n_CC2_v)**2,  &
                  one,           &
                  g_kl_cd,       &
                  (n_CC2_o)**2,  &
                  x_CD_ij,       &
                  (n_CC2_v)**2,  &
                  one,           &
                  I_kl_ij,       &
                  (n_CC2_o)**2)
!
      call deallocator(x_CD_ij, (n_CC2_v)**2, (n_CCSD_o)**2)
      call deallocator(g_KL_CD, (n_CC2_o)**2, (n_CC2_v)**2)
!
!     Reorder s_KC_LD into s_ab_KL
!
      call allocator(x_ab_KL, (n_CCSD_v)**2, (n_CC2_o)**2)
!
      do L = 1, n_CC2_o
         do K = 1, n_CC2_o
!
            KL = index_two(K, L, n_CC2_o)
!
            do b = 1, n_CCSD_v
!
               Lb = index_two( L, b, n_CC2_o)
!
               do a = 1, n_CCSD_v
!
                  Ka = index_two(K, a, n_CC2_o)
                  ab = index_two(a, b, n_CCSD_v)
!
                  x_ab_KL(ab, KL) = x_KC_LD(Ka, Lb)
!
               enddo
            enddo
         enddo
      enddo
!
!     omega_ab_ij = sum_(kl) s_ab_kl*I_kl_ij 
!                 = sum_(kl) s_ab_kl*(g_kilj + sum_(cd)(s_ci_dj + t_ci_dj)*g_kc_ld)
!
      call allocator(omega_ab_ij, (n_CCSD_v)**2, (n_CCSD_o)**2)
!
      call dgemm('N','N',        &
                  (n_CCSD_v)**2, &
                  (n_CCSD_o)**2, &
                  (n_CC2_o)**2,  &
                  one,           &
                  x_ab_KL,       &
                  (n_CCSD_v)**2, &
                  I_KL_ij,       &
                  (n_CC2_o)**2,  &
                  zero,          &
                  omega_ab_ij,   &
                  (n_CCSD_v)**2)
!
      call deallocator(x_ab_KL, (n_CCSD_v)**2, (n_CC2_o)**2)
      call deallocator(I_KL_ij, (n_CC2_o)**2, (n_CCSD_o)**2)
!
!     Reorder into omega_aibj
!
      do i = 1, n_CCSD_o
         do j = 1, n_CCSD_o
!
            ij = index_two(i, j, n_CCSD_o)
!
            do b = 1, n_CCSD_v
!
               bj = index_two(b, j, n_CCSD_v)
!
               do a = 1, n_CCSD_v
!
                  ai = index_two(a, i, n_CCSD_v)                 
!
                  if (ai .ge. bj) then
!
                     ab = index_two(a, b, n_CCSD_v)
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega_ab_ij(ab, ij)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ab_ij, (n_CCSD_v)**2, (n_CCSD_o)**2)
!
   end subroutine omega_mlccsd_b2_mlccsd
!
!
  module subroutine omega_mlccsd_c2_mlccsd(wf, x_lc_kd)
!!
!!    Omega C2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!    
!!    Z_ai_bj = - sum_(ck) x^bc_kj*( g_ki,ac - 1/2 * sum_(dl) x^ad_li*g_kd,cl )
!!
!!    Omega_ai_bj = P_ij^ab (1/2 Z_ai_bj + Z_aj_bi) = 1/2 Z_ai_bj + 1/2 Z_bj_ai + Z_aj_bi+ Z_bi_aj
!!
!!    
      implicit none
!
      class(mlccsd)            :: wf
      real(dp), dimension(:,:) :: x_lc_kd
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_Ki_J 
      real(dp), dimension(:,:), allocatable :: L_Ca_J 
      real(dp), dimension(:,:), allocatable :: L_aC_J 
      real(dp), dimension(:,:), allocatable :: L_KD_J 
      real(dp), dimension(:,:), allocatable :: g_kd_lc
      real(dp), dimension(:,:), allocatable :: g_dl_ck
      real(dp), dimension(:,:), allocatable :: g_ki_ca
      real(dp), dimension(:,:), allocatable :: g_ai_ck
!
!     Reordered X2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: x_ai_dl
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:), allocatable :: I_ai_ck
      real(dp), dimension(:,:), allocatable :: Z_ai_bj
!  
!     Indices
!     
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0, k = 0, l = 0
!
      integer(i15) :: ai = 0, aj = 0, al = 0, bi = 0, bj = 0, bk = 0, cj = 0
      integer(i15) :: ck = 0, cl = 0, di = 0, dk = 0, dl = 0, ck_restricted = 0
      integer(i15) :: kd = 0, lc = 0, ca = 0, ac = 0, id = 0, la = 0
      integer(i15) :: ki = 0
!
      integer(i15) :: aldi = 0, aibj = 0, cldk = 0, bkcj = 0, cjbk = 0
!
!     Batching and memory handling
!
      integer(i15) :: required = 0, available = 0
!
      integer(i15) :: n_batch = 0, max_batch_length = 0
      integer(i15) :: a_batch = 0, a_start = 0, a_end = 0, a_length = 0 
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o
      integer(i15) :: first_CC2_v
      integer(i15) :: first_CCSD_o
      integer(i15) :: first_CCSD_v
!
      integer(i15) :: last_CC2_o 
      integer(i15) :: last_CC2_v 
      integer(i15) :: last_CCSD_o 
      integer(i15) :: last_CCSD_v 
!
!     Calculate first/last indeces
!
!     CC2
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
!     CCSD
!
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1
!
!     Construct g_Ki,aC ordered as g_Ki_Ca
!
!     Allocate g_ki_ca
!
      call allocator(g_Ki_Ca, (n_CC2_o)*(n_CCSD_o), (n_CC2_v)*(n_CCSD_v))
      g_Ki_Ca = zero
!
!     Allocate L_ki_J
!
      call allocator(L_Ki_J, (n_CC2_o)*(n_CCSD_o), wf%n_J) 
!
!     Get cholesky vectors of ij-type
!
      call wf%get_cholesky_ij(L_Ki_J, first_CC2_o, last_CC2_o, first_CCSD_o, last_CCSD_o)
!
!     Prepare batching over a 
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = 2*(n_CC2_v**2)*(wf%n_J) + 2*(n_CC2_v)*(n_CC2_o)*(wf%n_J)
      required = 4*required
      call num_batch(required, available, max_batch_length, n_batch, n_CCSD_v)
!
      a_start  = 1
      a_end    = 0
      a_length = 0
!
!     Start looping over batches
!
      do a_batch = 1,n_batch
!
!        Get batch limits  and  length of batch
!
         call batch_limits(a_start, a_end, a_batch, max_batch_length, n_CCSD_v)
         a_length = a_end - a_start + 1
!
!        Get ab-cholesky vectors for the batch, L_ac^J, then reorder from L_ac_J to L_ca_J
!
         call allocator(L_aC_J, (n_CC2_v)*a_length, wf%n_J)
         L_ac_J = zero
         call wf%get_cholesky_ab(L_aC_J, a_start, a_end, first_CC2_v, last_CC2_v)
!
         call allocator(L_ca_J, (n_CC2_v)*a_length, wf%n_J)
         L_ca_J = zero
!
         do a = 1, a_length
            do c = 1, n_CC2_v
!
               ac = index_two(a, c, a_length) 
               ca = index_two(c, a, n_CC2_v) 
!
               do J = 1, wf%n_J
!
                 L_ca_J(ca, J) = L_ac_J(ac, J)
!
               enddo
            enddo
         enddo

        call deallocator(L_ac_J,(n_CC2_v)*a_length, wf%n_J)
!
!       g_Ki_Ca = sum_J L_Ki_J*L_Ca_J
!
        call dgemm('N','T',                                    &
                    n_CC2_o*n_CCSD_o,                          &
                    (n_CC2_v)*a_length,                        &
                    wf%n_J,                                    &   
                    one,                                       &   
                    L_ki_J,                                    &
                    n_CC2_o*n_CCSD_o,                          &
                    L_ca_J,                                    &
                    (n_CC2_v)*a_length,                        &
                    one,                                       &
                    g_Ki_Ca(1,index_two(1, a_start, n_CC2_v)), &
                    n_CC2_o*n_CCSD_o)
!
!        Deallocate L_ca_J
!
         call deallocator(L_Ca_J, (n_CC2_v)*a_length, wf%n_J)
!
      enddo ! End of batching
!
!     Deallocate L_Ki_J
!
      call deallocator(L_Ki_J, n_CC2_o*n_CCSD_o, wf%n_J)
!
!
!     Reorder g_ki_ca to g_ai_ck
!
      call allocator(I_ai_CK, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
      do i = 1, n_CCSD_o
         do K = 1, n_CC2_o
!
            Ki = index_two(k, i, n_CC2_o)
!
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CCSD_v)
!
               do C = 1, n_CC2_v
!
                  Ca = index_two(C, a, n_CC2_v)
                  CK = index_two(C, K, n_CC2_v)
!
                  I_ai_CK(ai, CK) = g_ki_ca(Ki, Ca)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_Ki_Ca, n_CC2_o*n_CCSD_o, n_CC2_v*n_CCSD_v)
!
!
!     Reorder s_la_id to s_ai_dl
!
      call allocator(x_ai_DL, n_CCSD_o*n_CCSD_v, n_CC2_o*n_CC2_v)
!
      do L = 1, n_CC2_o        
         do D = 1, n_CC2_v
!
            DL = index_two(D, L, n_CC2_v)
!
            do i = 1, n_CCSD_o
!
               do a = 1, n_CCSD_v               
!
                  La = index_two(L, a, n_CC2_o)
                  ai = index_two(a, i, n_CCSD_v)
                  iD = index_two(i, d, n_CC2_o)
!
                  x_ai_DL(ai, DL) = x_lc_kd(La, iD)
!             
               enddo
            enddo
         enddo
      enddo
!
!     Construct and reorder g_kd,lc
!
      call allocator(L_KD_J, (n_CC2_o)*(n_CC2_v),(wf%n_J))
!
!     Get L_ia_J
!
      call wf%get_cholesky_ia(L_KD_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
!     Create g_kd_lc = sum_J L_kd_J * L_lc_J
!
      call allocator(g_KD_LC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      call dgemm('N','T',              &
                  (n_CC2_o)*(n_CC2_v), &
                  (n_CC2_o)*(n_CC2_v), &
                  wf%n_J,              &   
                  one,                 &
                  L_KD_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  L_KD_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  zero,                &
                  g_KD_LC,             &
                  (n_CC2_o)*(n_CC2_v))
!
!     Deallocate L_ia_J
!
      call deallocator(L_KD_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
      call allocator(g_DL_CK, n_CC2_v*n_CC2_o, n_CC2_v*n_CC2_o)
!
      do L = 1, n_CC2_o        
         do D = 1, n_CC2_v
!
            DL = index_two(D, L, n_CC2_v)
!
            do K = 1, n_CC2_o
!
               KD = index_two(K, D, n_CC2_o)
!
               do C = 1, n_CC2_v               
!
                  LC = index_two(L, C, n_CC2_o)
                  CK = index_two(C, K, n_CC2_v)
!
                  g_DL_CK(DL, CK) = g_KD_LC(KD,LC)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_KD_LC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      call dgemm('N','N',                 &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  (n_CC2_o)*(n_CC2_v),    &
                  -half,                  &
                  x_ai_DL,                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  g_DL_CK,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  one,                   &
                  I_ai_CK,                &
                  (n_CCSD_o)*(n_CCSD_v))
!
      call deallocator(g_DL_CK, n_CC2_v*n_CC2_o, n_CC2_v*n_CC2_o)
!
!     Create Z_ai_bj = -sum(ck) I_ai_ck*s_ck_bj (s^bc_kj)
!
      call allocator(Z_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
      call dgemm('N','T',                 &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  -one,                   &
                  I_ai_ck,                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  x_ai_dl,                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  zero,                   &
                  Z_ai_bj,                &
                  (n_CCSD_o)*(n_CCSD_v))
!
      call deallocator(I_ai_CK, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
      call deallocator(x_ai_DL, n_CCSD_o*n_CCSD_v, n_CC2_o*n_CC2_v)   
!
!     Omega_aibj,1 = P_ai_bj ( 1/2*Z_ai_bj + Z_aj_bi )
!
         do i = 1, n_CCSD_o
            do a = 1, n_CCSD_v
!
               ai = index_two(a, i, n_CCSD_v)
!
               do j = 1, n_CCSD_o    
                  do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CCSD_v)
!
                  if (ai .ge. bj) then
                     aj = index_two(a, j, n_CCSD_v)
                     bi = index_two(b, i, n_CCSD_v)
!
                     aibj=index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + half*Z_ai_bj(ai, bj) + Z_ai_bj(aj, bi) &
                                                               + half*Z_ai_bj(bj, ai) + Z_ai_bj(bi, aj)
!
                  endif
!  
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(Z_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
   end subroutine omega_mlccsd_c2_mlccsd
!
!
  module subroutine omega_mlccsd_d2_mlccsd(wf, x_KC_LD)
!!
!!     Omega D2 
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Calculates the D2 term,
!!
!!      D2: sum_ck u_jk^bc g_aikc 
!!        - 1/2 * sum_ck u_jk^bc g_acki 
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!     where 
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!     The first, second, and third terms are referred to as D2.1, D2.2, and D2.3, 
!!     and comes out ordered as (ai,bj). All terms are added to the omega vector of the 
!!     wavefunction object wf.
!
!     The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
      implicit none 
!
      class(mlccsd) :: wf 
      real(dp), dimension(:,:) :: x_KC_LD
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, batch_dimension = 0, n_batch = 0
      integer(i15) :: a_begin = 0, a_end = 0, a_batch = 0, batch_length = 0, a_full = 0, ac_dim = 0 
!
!     Indices 
!
      integer(i15) :: ai = 0, aidl = 0, al = 0, aldi = 0, a = 0, i = 0, b = 0, ca = 0, ac = 0
      integer(i15) :: j = 0, c = 0, d = 0, di = 0, dl = 0, k = 0, kc = 0, kd = 0, l = 0, ki = 0
      integer(i15) :: lc = 0, ld = 0, aibj = 0, bj = 0, bjck = 0, bk = 0, bkcj = 0, cj = 0, ck = 0
      integer(i15) :: ia = 0, kb = 0, jb = 0, jc = 0, la = 0, id = 0
!
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing D2.3, D2.2 & D2.1
!
!     Vectors for D2.3 term 
!
      real(dp), dimension(:,:), allocatable :: L_kc_J  ! L_kc^J 
      real(dp), dimension(:,:), allocatable :: g_ld_kc ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: L_ld_kc ! L_ldkc = 2 * g_ldkc - g_lckd 
      real(dp), dimension(:,:), allocatable :: u_ai_ld ! u_il^ad = 2 * t_il^ad - t_li^ad 
      
      real(dp), dimension(:,:), allocatable :: Z_ai_kc ! An intermediate, see below
!
!     Vectors for D2.2 term 
!
      real(dp), dimension(:,:), allocatable :: g_ai_kc ! g_aikc 
      real(dp), dimension(:,:), allocatable :: u_kc_bj ! u_jk^bc
      real(dp), dimension(:,:), allocatable :: L_ai_J  ! L_ai^J 
!
!     Vectors for D2.1 term 
!
      real(dp), dimension(:,:), allocatable :: g_ai_ck ! g_acki
      real(dp), dimension(:,:), allocatable :: g_ac_ki ! g_acki; a is batched over 
      real(dp), dimension(:,:), allocatable :: L_ac_J  ! L_ac^J; a is batched over 
      real(dp), dimension(:,:), allocatable :: L_ki_J  ! L_ki^J 
      real(dp), dimension(:,:), allocatable :: u_ck_bj ! u_jk^bc
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o  
      integer(i15) :: first_CC2_v 
!
      integer(i15) :: first_CCSD_o  
      integer(i15) :: first_CCSD_v
! 
      integer(i15) :: last_CC2_o
      integer(i15) :: last_CC2_v
! 
      integer(i15) :: last_CCSD_o
      integer(i15) :: last_CCSD_v
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1 
!
!     :: Calculate the D2.3 term of omega ::
!
!     Allocate the Cholesky vector L_kc_J = L_kc^J and set to zero 
!
      call allocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     Read the Cholesky vector L_kc_J from file
!
      call wf%get_cholesky_ia(L_KC_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
!     Allocate g_ld_kc = g_ldkc and set to zero 
!
      call allocator(g_LD_KC,(n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     Calculate g_LD_KC = g_LDKC = sum_J L_LD^J L__KC^J
!
      call dgemm('N','T',              &
                  (n_CC2_o)*(n_CC2_v), &
                  (n_CC2_o)*(n_CC2_v), &
                  wf%n_J,              &
                  one,                 &
                  L_kc_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  L_kc_J,              &
                  (n_CC2_o)*(n_CC2_v), &
                  zero,                &
                  g_ld_kc,             &
                  (n_CC2_o)*(n_CC2_v))
!
      call deallocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     Allocate L_ld_kc = L_ldkc and set to zero    
!
      call allocator(L_LD_KC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     Determine L_ld_kc = L_ldkc from g_ld_kc = g_ldkc 
!
      do L = 1, n_CC2_o
         do D = 1, n_CC2_v
            do K = 1, n_CC2_o
               do C = 1, n_CC2_v
!
!                 Calculate the necessary indices 
!
                  LD = index_two(L, D, n_CC2_o)
                  KC = index_two(K, C, n_CC2_o)
!
                  LC = index_two(L, C, n_CC2_o)
                  KD = index_two(K, D, n_CC2_o)
!
!                 Set the value of L_ld_kc = 2 * g_ldkc - g_lckd 
!
                  L_LD_KC(LD, KC) = two*g_LD_KC(LD, KC) - g_LD_KC(LC, KD)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_LD_KC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
! 
!     Determine u_ai_ld = u_il^ad = 2 * x_il^ad - x_li^ad 
!
!
      call allocator(u_ai_LD, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
      do i = 1, n_CCSD_o
         do l = 1, n_CC2_o
            do a = 1, n_CCSD_v
!
               ai   = index_two(a, i, n_CCSD_v)
               ia   = index_two(i, a, n_CC2_o)
               la   = index_two(l, a, n_CC2_o)
!
               do d = 1, n_CC2_v
!
                  ld   = index_two(l, d, n_CC2_o)
                  id   = index_two(i, d, n_CC2_o)
!
                  u_ai_LD(ai, LD) = two*(x_KC_LD(IA, LD)) - x_KC_LD(LA, ID)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc and set it to zero
!
      call allocator(Z_ai_KC, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!     Form the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc
!
      call dgemm('N','N',                 &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  (n_CC2_o)*(n_CC2_v),    &
                  one,                    &
                  u_ai_LD,                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  L_LD_KC,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  zero,                   &
                  Z_ai_KC,                &
                  (n_CCSD_o)*(n_CCSD_v))
!
!     Deallocate L_ld_kc
!
      call deallocator(L_ld_kc, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     Allocate the D2.3 term omega2_ai_bj and set it to zero
!
      call allocator(omega2_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
!     Form the D2.3 term, 1/4 sum_kc Z_ai_kc u_kc_bj = 1/4 sum_kc Z_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call dgemm('N','T',                 &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  one/four,               &
                  Z_ai_kc,                &
                  (n_CCSD_o)*(n_CCSD_v),  & 
                  u_ai_ld,                &
                  (n_CCSD_o)*(n_CCSD_v),  & 
                  zero,                   &
                  omega2_ai_bj,           &
                  (n_CCSD_o)*(n_CCSD_v))
!
!     Some mathematical justification for the above matrix multiplication. We have 
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_ai_ld(ai,ld) = u_il^ad, 
!     which means that u_ai_ld(bj,kc)^T = u_ai_ld(kc,bj) = u_kj^cb = u_jk^bc.
!
!
!     Deallocate the Z_ai_kc intermediate 
!
      call deallocator(Z_ai_kc, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!     Add the D2.3 term to the omega vector 
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai   = index_two(a, i, n_CCSD_v)
!
            do b = 1, n_CCSD_v
               do j = 1, n_CCSD_o
!
                  bj   = index_two(b, j, n_CCSD_v)
!
                  aibj = index_packed(ai, bj)
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega2_ai_bj(ai, bj) & 
                                                               + omega2_ai_bj(bj, ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the omega2_ai_bj and u_ai_ld(ai,ld) = u_il^ad vector
!
      call deallocator(omega2_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
      call deallocator(u_ai_LD, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v)) 
!
!     :: Calculate the D2.1 term of omega :: 
!           
!
!     Allocate the L_ai_J and L_KC_J terms and set them to zero 
!
      call allocator(L_ai_J, (n_CCSD_o)*(n_CCSD_v), wf%n_J)
      call allocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     Read the Cholesky vectors from file 
!
      call wf%get_cholesky_ai(L_ai_J, first_CCSD_v, last_CCSD_v, first_CCSD_o, last_CCSD_o)
      call wf%get_cholesky_ia(L_KC_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
!     Allocate g_ai_kc = g_aikc and set it zero
!
      call allocator(g_ai_KC, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
!
!     Form the g_ai_kc integrals from L_ai_J and L_kc_J
!  
     call dgemm('N','T',                  &
                 (n_CCSD_o)*(n_CCSD_v),   &
                 (n_CC2_o)*(n_CC2_v),     &   
                 wf%n_J,                  &
                 one,                     &
                 L_ai_J,                  &
                 (n_CCSD_o)*(n_CCSD_v),   & 
                 L_KC_J,                  &
                 (n_CC2_o)*(n_CC2_v),     &
                 zero,                    &
                 g_ai_KC,                 &
                 (n_CCSD_o)*(n_CCSD_v))
!
!     Deallocate the Cholesky vectors L_ai_J and L_kc_J
!
      call deallocator(L_ai_J, (n_CCSD_o)*(n_CCSD_v), wf%n_J)
      call deallocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     Allocate u_kc_bj and set it to zero 
!
      call allocator(u_KC_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
!     Determine u_kc_bj = u_jk^bc = 2 * t_jk^bc - t_kj^bc 
!
      do K = 1, n_CC2_o
         do C = 1, n_CC2_v
!
            KC = index_two(K, C, n_CC2_o)
!
            do j = 1, n_CCSD_o
!
               JC   = index_two(J, C, n_CC2_o)
!
               do b = 1, n_CCSD_v
!                 
                  bj   = index_two(b, j, n_CCSD_v)
                  KB   = index_two(k, b, n_CC2_o)            
                  JB   = index_two(j, b, n_CC2_o)            
!    
                  u_KC_bj(kc, bj) = two*(x_KC_LD(KC, JB)) - x_KC_LD(KB, JC)
!
               enddo
            enddo
         enddo
      enddo 
!
!     Allocate omega2_ai_bj and set it to zero 
!
      call allocator(omega2_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v)) 
!
!     Calculate the D2.1 term sum_ck u_jk^bc g_aikc = sum_ck g_ai_kc u_kc_bj
!
      call dgemm('N','N',                 &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  one,                    &
                  g_ai_kc,                &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  u_kc_bj,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  zero,                   &
                  omega2_ai_bj,           &
                  (n_CCSD_o)*(n_CCSD_v))
!
!     Add the D2.1 term to the omega vector 
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai   = index_two(a, i, n_CCSD_v)
!
            do j = 1, n_CCSD_o
               do b = 1, n_CCSD_v
!                 
                  bj   = index_two(b, j, n_CCSD_v)
!
                  aibj = index_packed(ai, bj)
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj,1) = wf%omega2(aibj,1) + omega2_ai_bj(ai,bj) & 
                                                            + omega2_ai_bj(bj,ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!    Deallocate g_ai_kc, u_kc_bj, and the omega2_ai_bj vectors 
!
     call deallocator(g_ai_KC, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
     call deallocator(u_KC_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
!
     call deallocator(omega2_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSd_o)*(n_CCSD_v))
!
!    :: Calculate D2.2 term of Omega ::
!
!    - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj 
!
!    Allocate L_ki_J and set it to zero
!
     call allocator(L_Ki_J, (n_CCSD_o)*(n_CC2_o), wf%n_J)
!
!    Read the Cholesky vector L_ki_J from file 
!
     call wf%get_cholesky_ij(L_Ki_J, first_CC2_o, last_CC2_o, first_CCSD_o, last_CCSD_o)
!
!    Allocate the full g_ai_ck = g_acki and set it to zero 
!
     call allocator(g_ai_CK, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
     g_ai_CK = zero
!
!    Prepare for batching over the index a to calculate g_ai_ck = g_acki
!
!    To calculate this term, we need to first create L_ac^J, then hold L_ac^J and g_acki
!    in memory simultaneously 
!
     required = (wf%n_J)*(n_CC2_v)*(n_CCSD_v)  ! Holding L_ac^J
!
     required = required &                                             
                 + max( ((n_CC2_v)**2)*((n_CC2_o)**2), &                     ! Testing if it is more demanding 
                  (wf%n_J)*(n_CC2_v)**2 + 2*(wf%n_J)*(n_CC2_o)*(n_CC2_v)) ! to hold g_acki or to create L_ac^J
!
     required = 4*required ! In words
!
     available = get_available()
     batch_dimension = n_CCSD_v
!
!     Determine the batching variables 
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension) 
!
!     Determine g_ai_ck = g_acki successively in batches over a 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index        
!
         call batch_limits(a_begin, a_end, a_batch, max_batch_length, batch_dimension)
         batch_length = a_end - a_begin + 1 
!
!        Get ab-cholesky vectors for the batch, L_ac^J, then reorder from L_ac_J to L_ca_J
!
         aC_dim = batch_length*(n_CC2_v) ! Dimension of ac for the batch over index a 
!
         call allocator(L_aC_J, ac_dim, wf%n_J)
         L_ac_J = zero
         call wf%get_cholesky_ab(L_ac_J, a_begin, a_end, first_CC2_v, last_CC2_v)
!
!       Allocate the integral g_ca_ki = g_acki and set to zero 
!
        call allocator(g_aC_Ki, aC_dim, (n_CC2_o)*(n_CCSD_o))
!
!       Calculate g_ac_ki = g_acki from L_ca_J = L_ac^J and L_ki_J = L_ki^J
!
       call dgemm('N','T',                &
                   aC_dim,                &
                   (n_CC2_o)*(n_CCSD_o),  &
                   wf%n_J,                &
                   one,                   &
                   L_ac_J,                &
                   ac_dim,                &
                   L_ki_J,                &
                   (n_CC2_o)*(n_CCSD_o),  &
                   zero,                  & 
                   g_ac_ki,               &
                   ac_dim)
!
       call deallocator(L_ac_J, aC_dim, wf%n_J)
!
!        Reorder the integrals g_ca_ki (reduced a) = g_acki = g_ai_ck (full a)
!
        do a = 1, batch_length
!
           a_full = a - 1 + a_begin ! The full matrix index a
!
           do i = 1, n_CCSD_o
!
              ai = index_two(a_full, i, n_CCSD_v)
!
              do k = 1, n_CC2_o
!
                 ki = index_two(k, i, n_CC2_o)
!
                 do c = 1, n_CC2_v
!
                    ac = index_two(a, c, batch_length)
                    ck = index_two(c, k, n_CC2_v)
!
                    g_ai_ck(ai, ck) = g_ac_ki(ac, ki)
!
                 enddo
              enddo
           enddo
        enddo
!
!       Deallocate the g_ca_ki and L_ca_J vectors
!
        call deallocator(g_ac_ki, ac_dim, (n_CC2_o)*(n_CCSD_o))
!
      enddo ! End of loop over batches of a
!     Allocate the u_ck_bj = u_jk^bc vector and set it to zero 
!
      call allocator(u_CK_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
      u_CK_bj = zero
!
!     Determine u_ck_bj = u_jk^bc = 2 * x_jk^bc - x_kj^bc 
!
      do k = 1, n_CC2_o   
         do j = 1, n_CCSD_o
            do c = 1, n_CC2_v
!
               CK = index_two(C, K, n_CC2_v)
               KC = index_two(K, C, n_CC2_o)
               JC = index_two(J, C, n_CC2_o)
!
               do b = 1, n_CCSD_v
!
                  bj = index_two(b, j, n_CCSD_v)
                  JB = index_two(J, B, n_CC2_o)
                  KB = index_two(K, B, n_CC2_o)
!
                  u_CK_bj(CK, bj) = two*(x_KC_LD(JB, KC)) - x_KC_LD(KB, JC)
!
               enddo
            enddo
         enddo
      enddo
!
!    Allocate the D2.2 term and set it to zero 
!
     call allocator(omega2_ai_bj, (n_CCSD_o)*(n_CCSD_v), (n_CCSD_o)*(n_CCSD_v))
!
!    Calculate the D2.2 term, - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj
!
      call dgemm('N','N',                 &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  (n_CC2_o)*(n_CC2_v),    &
                  -one/two,               &
                   g_ai_ck,               &
                  (n_CCSD_o)*(n_CCSD_v),  &
                  u_ck_bj,                &
                  (n_CC2_o)*(n_CC2_v),    &
                  zero,                   &
                  omega2_ai_bj,           &
                  (n_CCSD_o)*(n_CCSD_v))
!
!    Add the D2.2 term to the omega vector 
!
     do i = 1, n_CCSD_o
        do a = 1, n_CCSD_v
!
           ai = index_two(a, i, n_CCSD_v)
!
           do j = 1, n_CCSD_o
              do b = 1, n_CCSD_v              
!
                 bj = index_two(b, j, n_CCSD_v)
!
                 aibj = index_packed(ai, bj)
!
                 if (ai .ge. bj) then
!
                    wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega2_ai_bj(ai,bj) &
                                                              + omega2_ai_bj(bj,ai)
!
                 endif
!
              enddo
           enddo
        enddo
     enddo  
!
!    Deallocations
!
     call deallocator(g_ai_ck, (n_CCSD_o)*(n_CCSD_v), (n_CC2_o)*(n_CC2_v))
     call deallocator(u_ck_bj, (n_CC2_o)*(n_CC2_v), (n_CCSD_o)*(n_CCSD_v))
     call deallocator(L_ki_J, (n_CC2_o)*(n_CCSD_o), wf%n_J)
!
     call deallocator(omega2_ai_bj,  (n_CCSD_o)*(n_CCSD_v),  (n_CCSD_o)*(n_CCSD_v))
!
   end subroutine omega_mlccsd_d2_mlccsd
!
!
   module subroutine omega_mlccsd_e2_mlccsd(wf, x_kc_ld)
!!
!!     Omega E2
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Calculates the E2 term,
!!
!!      E2: sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) 
!!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!     where
!!
!!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!!
!!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!!
!!     Both are permuted added to the projection vector element omega2(ai,bj) of
!!     the wavefunction object wf.
!!
       implicit none 
!
      class(mlccsd) :: wf 
      real(dp), dimension(:,:) :: x_kc_ld
!
!     Indices 
!
      integer(i15) :: aib = 0, aibk = 0, bk = 0, bja = 0, ibj = 0, aibj = 0, dlck = 0
      integer(i15) :: b = 0, c = 0, k = 0, d = 0, ck = 0, ckdl = 0, cl = 0, cldk = 0
      integer(i15) :: dk = 0, dl = 0, kc = 0, kdl = 0, l = 0, ld = 0, a = 0, ai = 0
      integer(i15) :: ia = 0, kd = 0, lc = 0, jc = 0, kb = 0 
      integer(i15) :: bj = 0, aicj = 0, cj = 0, i = 0, j = 0, jai = 0, dlc = 0, dkcl = 0
!
!     Vectors for E2.1 term 
!
      real(dp), dimension(:,:), allocatable :: omega2_b_jai ! For storing the E2.1 term temporarily
      real(dp), dimension(:,:), allocatable :: L_kc_J       ! L_kc^J
      real(dp), dimension(:,:), allocatable :: g_ld_kc      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: g_kdl_c      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: u_b_kdl      ! u_kl^bd 
      real(dp), dimension(:,:), allocatable :: X_b_c        ! An intermediate, see below for definition
      real(dp), dimension(:,:), allocatable :: x_c_jai      ! t_ij^ac 
      
!
!     Vectors for E2.2 term 
!
      real(dp), dimension(:,:), allocatable :: g_k_dlc      ! g_ldkc 
      real(dp), dimension(:,:), allocatable :: u_dlc_j      ! u_lj^dc 
      real(dp), dimension(:,:), allocatable :: omega2_aib_j ! For storing the E2.2 term temporarily
      real(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition 
      real(dp), dimension(:,:), allocatable :: x_aib_k      ! t_ik^ab 
!
!     Active space variables
!
      integer(i15) :: n_CC2_o = 0, n_CC2_v = 0
      integer(i15) :: n_CCSD_o = 0, n_CCSD_v = 0
!
      integer(i15) :: first_CC2_o  
      integer(i15) :: first_CC2_v 
!
      integer(i15) :: first_CCSD_o  
      integer(i15) :: first_CCSD_v
! 
      integer(i15) :: last_CC2_o
      integer(i15) :: last_CC2_v
! 
      integer(i15) :: last_CCSD_o
      integer(i15) :: last_CCSD_v
!
      integer(i15) :: offset 
!
!     Calculate first/last indeces
! 
      call wf%get_CC2_active_indices(first_CC2_o, first_CC2_v)
      call wf%get_CC2_n_active(n_CC2_o, n_CC2_v)
!
      call wf%get_CCSD_active_indices(first_CCSD_o, first_CCSD_v)
      call wf%get_CCSD_n_active(n_CCSD_o, n_CCSD_v)
!
      last_CC2_o = first_CC2_o + n_CC2_o - 1
      last_CC2_v = first_CC2_v + n_CC2_v - 1 
!
      last_CCSD_o = first_CCSD_o + n_CCSD_o - 1
      last_CCSD_v = first_CCSD_v + n_CCSD_v - 1
!
!     :: Calculate the E2.1 term of omega ::
!
!     Form the Cholesky vector L_kc_J = L_kc^J 
!
      call allocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_KC_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
!     Form g_ld_kc = g_ldkc 
!
      call allocator(g_LD_KC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      call dgemm('N','T',            & 
                  (n_CC2_o)*(n_CC2_v), &
                  (n_CC2_o)*(n_CC2_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (n_CC2_o)*(n_CC2_v), &
                  L_kc_J,            &
                  (n_CC2_o)*(n_CC2_v), &
                  zero,              &
                  g_ld_kc,           &
                  (n_CC2_o)*(n_CC2_v))
!
!     Deallocate the Cholesky vector, L_kc_J
!
      call deallocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     Allocate u_b_kdl = u_kl^bd 
!
      call allocator(u_b_kdl, n_CCSD_v, (n_CC2_v)*((n_CC2_o)**2))
!
!     Allocate g_kdl_c = g_ldkc 
!
      call allocator(g_kdl_c, (n_CC2_v)*((n_CC2_o)**2), n_CC2_v)
!
!     Determine u_b_kdl = u_kl^bd and g_kdl_c = g_ldkc
! 
!
      do l = 1, n_CC2_o
         do d = 1, n_CC2_v
!
            LD = index_two(L, D, n_CC2_o)
!
            do k = 1, n_CC2_o
!
               KD  = index_two(K, D, n_CC2_o)
               KDL = index_three(K, D, L, n_CC2_o, n_CC2_v)
!
               do c = 1, n_CC2_v ! Use as though "b" for u_b_kdl term
!                  
                  KC   = index_two(K, C, n_CC2_o)
                  LC   = index_two(L, C, n_CC2_o)
                  
!
!                 Set the values of u_b_kdl and g_kdl_c
!
                  if (C .le. n_CCSD_v) then ! C as b
                     u_b_KDL(c, KDL) = two*(x_KC_LD(KC, LD)) - x_KC_LD(LC, KD)
                  endif
!
                  g_KDL_C(KDL, C) = g_LD_KC(LD, KC) ! g_ldkc 
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the unordered integrals g_ld_kc = g_ldkc
!
!        Note: It might be better to reorder g_kdl_c to g_k_dlc in the 
!        calculation of the E2.2 term. For now (8 Mar 2017), it remains 
!        simple & stupid.
!
      call deallocator(g_LD_KC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     Allocate the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call allocator(X_b_C, n_CCSD_v, n_CC2_v)
!
!     Copy the virtual-virtual Fock matrix into the intermediate 
!
      do b = 1, n_CCSD_v
         do c = 1, n_CC2_v
             X_b_c(b, C) = wf%fock_ab(b, C) 
         enddo
      enddo
!
!     Add the second contribution, 
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',                 &
                  n_CCSD_v,               &
                  n_CC2_v,                &
                  (n_CC2_v)*(n_CC2_o)**2, &
                  -one,                   &
                  u_b_KDL,                &
                  n_CCSD_v,               &
                  g_KDL_C,                &
                  (n_CC2_v)*(n_cc2_o)**2, &
                  one,                    &
                  X_b_C,                  &
                  n_CCSD_v)
!
!     Deallocate u_b_kdl and g_kdl_c
!
      call deallocator(u_b_KDL, n_CCSD_v, (n_CC2_v)*(n_CC2_o)**2)
      call deallocator(g_KDL_C, (n_CC2_v)*(n_CC2_o)**2, n_CC2_v)
!
!     Form the reordered t_c_jai = t_ij^ac
!
      call allocator(x_c_jai, n_CC2_v, (n_CCSD_v)*(n_CCSD_o)**2)
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ia = index_two(i, a, n_CC2_o)
!
            do j = 1, n_CCSD_o
!
               jai = index_three(j, a, i, n_CCSD_o, n_CCSD_v)
!
               do c = 1, n_CC2_v
!               
                  jc   = index_two(j, c, n_CC2_o)
!
                  x_c_jai(c, jai) = x_KC_LD(IA, JC)
!
               enddo
            enddo
         enddo
      enddo
!
!     Form the E2.1 term 
!
      call allocator(omega2_b_jai, n_CCSD_v, (n_CCSD_v)*(n_CCSD_o)**2)
!
      call dgemm('N','N',                    &
                  n_CCSD_v,                  &
                  (n_CCSD_v)*(n_CCSD_o)**2,  &
                  n_CC2_v,                   &
                  one,                       &
                  X_b_c,                     &
                  n_CCSD_v,                  &
                  x_c_jai,                   &
                  n_CC2_v,                   &
                  zero,                      &
                  omega2_b_jai,              &
                  n_CCSD_v)
!
!     Add the E2.1 term to the omega vector 
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CCSD_v)
!
            do j = 1, n_CCSD_o
!
               jai = index_three(j, a, i, n_CCSD_o, n_CCSD_v)
!
               do b = 1, n_CCSD_v
!
                  ibj  = index_three(i, b, j, n_CCSD_o, n_CCSD_v)
!             
                  bj   = index_two(b, j, n_CCSD_v)
!
                  aibj = index_packed(ai,bj)
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj,1) = wf%omega2(aibj,1) + omega2_b_jai(b,jai) &
                                                            + omega2_b_jai(a,ibj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the E2.1 term, the X intermediate, and the reordered amplitudes 
!
      call deallocator(omega2_b_jai, n_CCSD_v, (n_CCSD_v)*(n_CCSD_o)**2)
      call deallocator(X_b_c, n_CCSD_v, n_CC2_v)
      call deallocator(x_c_jai, n_CC2_v, (n_CCSD_v)*(n_CCSD_o)**2)
!
!     :: Calculate E.2.2 term of omega ::
!
!
!     Form the Cholesky vector L_kc_J = L_kc^J 
!
      call allocator(L_KC_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_KC_J, first_CC2_o, last_CC2_o, first_CC2_v, last_CC2_v)
!
!     Form g_ld_kc = g_ldkc = sum_J L_ld^J L_kc^J 
!
      call allocator(g_LD_KC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
      call dgemm('N','T',            &
                  (n_CC2_o)*(n_CC2_v), &
                  (n_CC2_o)*(n_CC2_v), &
                  wf%n_J,            &
                  one,               &
                  L_kc_J,            &
                  (n_CC2_o)*(n_CC2_v), &
                  L_kc_J,            &
                  (n_CC2_o)*(n_CC2_v), &
                  zero,              &
                  g_ld_kc,           &
                  (n_CC2_o)*(n_CC2_v))
!
!     Deallocate the Cholesky vector, L_kc_J
!
      call deallocator(L_kc_J, (n_CC2_o)*(n_CC2_v), wf%n_J)
!
!     Allocate g_k_dlc = g_ldkc
!
      call allocator(g_k_dlc, n_CC2_o, (n_CC2_o)*((n_CC2_v)**2))
!
!     Allocate u_dlc_j = u_lj^dc
!
      call allocator(u_dlc_j, (n_CC2_o)*((n_CC2_v)**2), n_CCSD_o)
!
!     Determine g_k_dlc = g_ldkc and u_dlc_j = u_lj^dc 
!
      do k = 1, n_CC2_o ! Use as though "j" for u_dlc_j term 
         do c = 1, n_CC2_v
!
            KC = index_two(K, C, n_CC2_o)
!
            do L = 1, n_CC2_o
!
               LC = index_two(L, C, n_CC2_o)
!
               do D = 1, n_CC2_v
!
                  DLC  = index_three(D, L, C, n_CC2_v, n_CC2_o)
!
                  LD = index_two(L, D, n_CC2_o)
                  KD = index_two(K, D, n_CC2_o)
!
!                 Set the value of g_k_dlc and u_dlc_j 
!
                  g_k_dlc(K, DLC) = g_ld_kc(ld, kc)   ! g_ldkc
                  if (k .le. n_CCSD_o) then !k = j
!
                     u_dlc_j(dlc, k) = two*(x_KC_LD(LD, KC)) - x_KC_LD(KD, LC) ! u_lk^dc = 2 * t_lk^dc - t_kl^dc 
!                 
                  endif
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the integrals g_ld_kc = g_ldkc 
!
      call deallocator(g_LD_KC, (n_CC2_o)*(n_CC2_v), (n_CC2_o)*(n_CC2_v))
!
!     Allocate the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc 
!                                     = F_k_j + sum_cdl g_k_dlc * u_dlc_j
!
      call allocator(Y_k_j, n_CC2_o, n_CCSD_o)
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj 
!
      do j = 1, n_CCSD_o
         do k = 1, n_CC2_o
            Y_k_j(k, j) = wf%fock_ij(k, j) 
         enddo
      enddo
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that 
!     Y_k_j = F_k_j + sum_cdl g_k_dlc u_dlc_j
!
      call dgemm('N','N',                    &
                  n_CC2_o,                   &
                  n_CCSD_o,                  &
                  (n_CC2_o)*((n_CC2_v)**2),  &
                  one,                       &
                  g_K_DLC,                   &
                  n_CC2_o,                   &
                  u_dlc_j,                   &
                  (n_CC2_o)*((n_CC2_v)**2),  &
                  one,                       &
                  Y_k_j,                     &
                  n_CC2_o)
!
!     Deallocate u_dlc_j and g_k_dlc 
!
      call deallocator(u_dlc_j, (n_CC2_o)*((n_CC2_v)**2), n_CCSD_o)
      call deallocator(g_k_dlc, n_CC2_o, (n_CC2_o)*((n_CC2_v)**2))
!
!     Allocate t_aib_k = t_ik^ab and set it to zero 
!
      call allocator(x_aib_k, (n_CCSD_o)*((n_CCSD_v)**2), n_CC2_o)
!
!     Determine t_aib_k = t_ik^ab 
!
      do k = 1, n_CC2_o
         do b = 1, n_CCSD_v
!
            KB = index_two(K, B, n_CC2_o)
!
            do i = 1, n_CCSD_o
               do a = 1, n_CCSD_v
!
                  IA = index_two(I, A, n_CC2_o)
                  aib = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
!
                  x_aib_k(aib, k) = x_KC_LD(ia, kb)
!
               enddo
            enddo
         enddo
      enddo
!
!     Allocate the E2.2 term and set to zero 
!
      call allocator(omega2_aib_j, (n_CCSD_o)*((n_CCSD_v)**2), n_CCSD_o)
!
!     Calculate the E2.2 term, 
!     - sum_k t_aib_k Y_k_j = - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
!
      call dgemm('N','N',                       &
                  (n_CCSD_o)*((n_CCSD_v)**2),   &
                  n_CCSD_o,                     &
                  n_CC2_o,                      &
                  -one,                         &
                  x_aib_k,                      &
                  (n_CCSD_o)*((n_CCSD_v)**2),   &
                  Y_k_j,                        &
                  n_CC2_o,                      &
                  zero,                         &
                  omega2_aib_j,                 &
                  (n_CCSD_o)*((n_CCSD_v)**2))
!
!     Deallocate t_aib_k and Y_k_j 
!
      call deallocator(x_aib_k, (n_CCSD_o)*((n_CCSD_v)**2), n_CC2_o)
      call deallocator(Y_k_j, n_CC2_o, n_CCSD_o)
!
!     Add the E2.2 term to the omega vector 
!
      do i = 1, n_CCSD_o
         do a = 1, n_CCSD_v
!
            ai = index_two(a, i, n_CCSD_v)
!
            do j = 1, n_CCSD_o
               do b = 1, n_CCSD_v
!
                  bj   = index_two(b, j, n_CCSD_v)
!
                  aibj = index_packed(ai, bj)
!
                  aib  = index_three(a, i, b, n_CCSD_v, n_CCSD_o)
                  bja  = index_three(b, j, a, n_CCSD_v, n_CCSD_o) 
!
!                 Restrict the indices to avoid adding both (ai,bj) and (bj,ai), 
!                 as they are identical in packed indices
!
                  if (ai .ge. bj) then
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) + omega2_aib_j(aib, j) & 
                                                               + omega2_aib_j(bja, i)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate the E2.2 term 
!
      call deallocator(omega2_aib_j, (n_CCSD_o)*((n_CCSD_v)**2), n_CCSD_o)
!
   end subroutine omega_mlccsd_e2_mlccsd
!
!
end submodule omega