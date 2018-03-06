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
!!     A1: sum_bcj L_Abjc * s_ij^bc,
!!  
!!     and adds it to the projection vector (omega1) of
!!     the wavefunction object wf
!! 
!! 
      implicit none
!
      class(mlcc2)   :: wf
!
!     Indices 
!
      integer(i15) :: i = 0, j = 0, a = 0, c = 0, b = 0
!
      integer(i15) :: Ac = 0, Ab = 0
      integer(i15) :: jc = 0, jb = 0 
      integer(i15) :: cjb = 0
!
!     Allocatables
! 
      real(dp), dimension(:,:), allocatable :: s_cj_bi 
      real(dp), dimension(:,:), allocatable :: L_A_cjb 
      real(dp), dimension(:,:), allocatable :: g_Ab_jc 
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
!     :: Construct L_Abjc
!
      call wf%mem%alloc(g_Ab_jc, n_active_v*wf%n_v, n_active_o*n_active_v)
! 
      integral_type = 'electronic_repulsion'
      call wf%get_vv_ov(integral_type, g_Ab_jc,          &
                        1, wf%n_v,                       &
                        first_active_v, last_active_v,   &     
                        first_active_o, last_active_o,   &     
                        first_active_v, last_active_v) 
!
!     L_Ab_jc = 2*g_Ab_jc - g_Ac_jb (ordered as L_A_cjb)
!
      call wf%mem%alloc(L_A_cjb, wf%n_v, n_active_o*(n_active_v**2))
!
      do A = 1, wf%n_v
         do b = 1, n_active_v
!
            Ab = index_two(A, b, wf%n_v)
!
            do c = 1, n_active_v
!
               Ac = index_two(A, c, wf%n_v)
!
               do j = 1, n_active_o
!
                  jb = index_two(j, b, n_active_o)
                  jc = index_two(j, c, n_active_o)
                  cjb = index_three(c, j, b, n_active_v, n_active_o)
!
                  L_A_cjb(A, cjb) = two*g_Ab_jc(Ab, jc) - g_Ab_jc(Ac, jb)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_Ab_jc, n_active_v*wf%n_v, n_active_o*n_active_v)
!
      call wf%mem%alloc(s_cj_bi, (n_active_o)*(n_active_v), (n_active_o)*n_active_v )
      call wf%get_s2am(s_cj_bi)
!
!     :: Add contributions to omega ::
!
      call dgemm('N', 'N',                                  &
                  wf%n_v,                                   &
                  n_active_o,                               &
                  (n_active_v**2)*(n_active_o),             &
                  one,                                      &
                  L_A_cjb,                                  &
                  wf%n_v,                                   &
                  s_cj_bi,                                  &
                  (n_active_v**2)*(n_active_o),             &
                  one,                                      &
                  wf%omega1(1,first_active_o),              &
                  wf%n_v)
! 
      call wf%mem%dealloc(L_A_cjb, n_active_v*wf%n_v, n_active_o*n_active_v)
!
      call wf%mem%dealloc(s_cj_bi, (n_active_o)*(n_active_v), (n_active_o)*n_active_v )
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
!!    B1: - sum_bjk s_jk^ab*L_kbjI + sum_bj u_ij^ab F_jb,
!!
!!    with u_ij^ab = 2*s_ij^ab - s_ij^ba. 
!!
!!    Batching over b.
!!
      implicit none
!
      class(mlcc2)   :: wf
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: s_aj_bk 
      real(dp), dimension(:,:), allocatable :: g_kb_jI 
      real(dp), dimension(:,:), allocatable :: L_jbk_I 
      real(dp), dimension(:,:), allocatable :: g_kI_jb 
!
!     looping indices
!
      integer(i15) :: i = 0, j = 0, k = 0, a  = 0, b = 0
      integer(i15) :: I_full = 0, J_full = 0, A_full  = 0, B_full = 0
      integer(i15) :: aj = 0, ai = 0, bj = 0, bi = 0, jb = 0, jI = 0, kb = 0, kI = 0
      integer(i15) :: jbk = 0
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
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1
!
!     :: Construct L_kb_jI
!
      call wf%mem%alloc(g_kb_jI, n_active_o*n_active_v, n_active_o*(wf%n_o) )
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_oo(integral_type, g_kb_jI,        &
                        first_active_o, last_active_o, &
                        first_active_v, last_active_v, &
                        first_active_o, last_active_o, &
                        1, wf%n_o)
!
      call wf%mem%alloc(g_kI_jb, n_active_o*(wf%n_o), n_active_o*n_active_v )
!
      integral_type = 'electronic_repulsion'             ! Do not think we need this can use g_kb_jI(jb, kI)
      call wf%get_oo_ov(integral_type, g_kI_jb,        &
                        first_active_o, last_active_o, &
                        1, wf%n_o,                     &
                        first_active_o, last_active_o, &
                        first_active_v, last_active_v)
!
!     L_kb_jI = 2*g_kb_jI - g_kI_jb (ordered as L_jbk_I)
!
      call wf%mem%alloc(L_jbk_I, (n_active_o**2)*n_active_v, (wf%n_o) )
!
      do I = 1, wf%n_o
!
         do k = 1, n_active_o
!
            kI = index_two(k, I, n_active_o)
!
            do j = 1, n_active_o
!
               jI = index_two(j, I, n_active_o)
!
               do b = 1, n_active_v
!
                  kb = index_two(k, b, n_active_o)
                  jb = index_two(j, b, n_active_o)
                  jbk = index_three(j, b, k, n_active_o, n_active_v)
!
                  L_jbk_I(jbk, I) = two*g_kb_jI(kb, jI) - g_kI_jb(kI, jb)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_kI_jb, n_active_o*(wf%n_o), n_active_o*n_active_v )
      call wf%mem%dealloc(g_kb_jI, n_active_o*n_active_v, n_active_o*(wf%n_o) )     
!  
      call wf%mem%alloc(s_aj_bk, (n_active_o)*n_active_v, (n_active_o)*n_active_v)
      call wf%get_s2am(s_aj_bk)
!
!     - sum_bjk s_jk^ab*L_kbjI 
!
      call dgemm('N', 'N',                        &
                  n_active_v,                     &
                  (wf%n_o),                       &
                  n_active_v*((n_active_o)**2),   &
                  -one,                           &
                  s_aj_bk,                        &
                  n_active_v,                     &
                  L_jbk_I,                        &
                  n_active_v*((n_active_o)**2),   &
                  one,                            &
                  wf%omega1(first_active_v, 1),   &
                  (wf%n_v))
!
      
      call wf%mem%dealloc(L_jbk_I, (n_active_o**2)*n_active_v, (wf%n_o) )
!
!     :: sum_jb F_jb (2*s_ij^ab-s_ji^ab) ::
!
      do i = 1, n_active_o
!
         I_full = i + first_active_o - 1
!
         do a = 1, n_active_v
!
            A_full = a + first_active_v - 1 
            ai = index_two(a, i, n_active_v)
!       
            do j = 1, n_active_o
!
               J_full = j + first_active_o - 1
               aj = index_two(a, j, n_active_v)

!
               do b = 1, n_active_v
! 
                  B_full = b + first_active_v - 1
                  bi = index_two(b, i, n_active_v)
                  bj = index_two(b, j, n_active_v)
!  
                  wf%omega1(A_full, I_full) = wf%omega1(A_full, I_full)&
                            + (2*s_aj_bk(ai, bj)-s_aj_bk(aj, bi))*wf%fock_ia(J_full, B_full)
!
               enddo
            enddo
         enddo
      enddo
!
   call wf%mem%dealloc(s_aj_bk, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
! 
!      
   end subroutine omega_mlcc2_b1_mlcc2
!
   module subroutine get_s2am_mlcc2(wf, s_ai_bj)
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
      real(dp), dimension((wf%n_CC2_v)*(wf%n_CC2_o), (wf%n_CC2_v)*(wf%n_CC2_o)) :: s_ai_bj
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
      call wf%mem%alloc(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_vo(integral_type, g_ai_bj,          &
                        first_active_v, last_active_v,   &
                        first_active_o, last_active_o,   &
                        first_active_v, last_active_v,   &
                        first_active_o, last_active_o)
!
!
         do a = 1, n_active_v
            do i = 1, n_active_o
!
               ai = index_two(a, i, n_active_v)
!
               do b = 1, n_active_v
                  do j = 1, n_active_o
!
                     bj = index_two(b, j, n_active_v)
!
                     s_ai_bj(ai, bj) = g_ai_bj(ai, bj)/(wf%fock_diagonal(i + first_active_o - 1,1)&
                                           + wf%fock_diagonal(j+ first_active_o - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + b + first_active_v - 1,1) &
                                           - wf%fock_diagonal(wf%n_o + a + first_active_v - 1,1))
!
                  enddo
               enddo
            enddo
         enddo
!
      call wf%mem%dealloc(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
   end subroutine get_s2am_mlcc2
!
end submodule