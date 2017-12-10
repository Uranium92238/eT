submodule (sccsd_class) omega
!
!!
!!    Omega submodule (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Contains the following family of procedures of the SCCSD class:
!!
!!    construct_omega: directs the calculation of the SCCSD omega vector.
!!
!
   implicit none 
!
   character(len=40) :: integral_type
!
!
contains 
!
!
   module subroutine construct_omega_sccsd(wf)
!!
!!    Construct Omega (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn 
!!
      implicit none 
!
      class(sccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_i_a ! F_ia
!
      real(dp), dimension(:,:), allocatable :: g_ac_kd ! g_ackd 
      real(dp), dimension(:,:), allocatable :: g_lc_ki ! g_lcki
!
      real(dp), dimension(:,:), allocatable :: omega_ai_bj_corr
!
      real(dp) :: norm_correction, ddot
!
      integer(i15) :: j = 0, b = 0, i = 0, a = 0, ai = 0, bj = 0, aibj = 0
!
      real(dp) :: begin_timer
      real(dp) :: end_timer
!
!     Calculate the CCSD omega vector 
!
      call cpu_time(begin_timer)
      call construct_omega_ccsd(wf)
      call cpu_time(end_timer)
!
      if (wf%settings%print_level == 'developer') then 
!
         write(unit_output,'(t6,a27,f15.8)') 'Time in CC  part (seconds):', end_timer-begin_timer
!
      endif
!
      call cpu_time(begin_timer)
!
!     Add the SCCSD singles contribution
!
      call wf%omega_sccsd_a1
!
!     The SCCSD doubles contribution to the projection vector
!     is identical to the doubles contriution to the Jacobian
!     transformation, except that the X, Y, and Z intermediates
!     are defined differently:
!
!        X_kc   = F_kc
!        Y_lcki = g_lcki
!        Z_ackd = g_ackd 
!
!     The X intermediate is already stored in memory (the wavefunction's
!     fock_ia array), and so we need only calcuate the Y and Z inter-
!     mediates.
!
!     Form Z_ackd = g_ac_kd 
!
      call allocator(g_ac_kd, (wf%n_v)**2, (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vv_ov(integral_type, g_ac_kd)
!
!     Form Y_lcki = g_lc_ki    
!
      call allocator(g_lc_ki, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_oo(integral_type, g_lc_ki)
!
!     Make X_ia = F_i_a
!
      call allocator(F_i_a, wf%n_o, wf%n_v)
      F_i_a = wf%fock_ia
!
!     Allocate correction & add the doubles terms 
!
      call allocator(omega_ai_bj_corr, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      omega_ai_bj_corr = zero
!
      call wf%jacobian_sccsd_a2(omega_ai_bj_corr, F_i_a)   ! X
      call wf%jacobian_sccsd_b2(omega_ai_bj_corr, g_lc_ki) ! Y 
      call wf%jacobian_sccsd_c2(omega_ai_bj_corr, g_ac_kd) ! Z 
!
!     Permute ai <-> bj and add to omega 
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  if (ai .ge. bj) then
!
                     aibj = index_packed(ai, bj)
!
                     wf%omega2(aibj, 1) = wf%omega2(aibj, 1) &
                                          + omega_ai_bj_corr(ai, bj) &
                                          + omega_ai_bj_corr(bj, ai)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(omega_ai_bj_corr, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call cpu_time(end_timer)
!
      if (wf%settings%print_level == 'developer') then 
!
         write(unit_output,'(t6,a27,f15.8/)') 'Time in SCC part (seconds):', end_timer-begin_timer
!
      endif
!
   end subroutine construct_omega_sccsd
!
!
   module subroutine omega_sccsd_a1_sccsd(wf)
!!
!!    Omega SCCSD A1 
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       t * P_IJK^ABC (delta_ai,AI L_JBKC - delta_ai,CI L_JBKA),
!!
!!    and adds it to omega(a,i) = omega_ai.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(:,:), allocatable :: L_jb_J
      real(dp), dimension(:,:), allocatable :: g_jb_kc
!
      real(dp), dimension(:,:), allocatable :: L ! L(JB,KC) = L_JBKC 
!
      integer(i15) :: c = 0, k = 0, j = 0, b = 0
!
      integer(i15) :: JB = 0, KC = 0, IA = 0, KC = 0, KA = 0
      integer(i15) :: JA = 0, KB = 0, IB = 0, JC = 0, IC = 0
!
!     Form L(jb, kc) = 2 * g_jb_kc(jb, kc) - g_jb_kc(jc, kb)
!
      call allocator(g_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_jb_kc)
!
      call allocator(L, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L = zero 
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            kc = index_two(k, c, wf%n_o)
!
            do b = 1, wf%n_v
!
               kb = index_two(k, b, wf%n_o)
!
               do j = 1, wf%n_o
!
                  jc = index_two(j, c, wf%n_o)
                  jb = index_two(j, b, wf%n_o)
!
                  L(jb, kc) = two*g_jb_kc(jb, kc) - g_jb_kc(jc, kb)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Scale L by the triples amplitude 
!
      call dscal(((wf%n_o)**2)*((wf%n_v)**2), wf%triples, L, 1) 
!
!     Add the contributions to the singles omega vector
!
      JB = index_two(wf%J, wf%B, wf%n_o)
      KC = index_two(wf%K, wf%C, wf%n_o)
      IA = index_two(wf%I, wf%A, wf%n_o)
      KC = index_two(wf%K, wf%C, wf%n_o)
      KA = index_two(wf%K, wf%A, wf%n_o)
      JA = index_two(wf%J, wf%A, wf%n_o)
      KB = index_two(wf%K, wf%B, wf%n_o)
      IB = index_two(wf%I, wf%B, wf%n_o)
      JC = index_two(wf%J, wf%C, wf%n_o)
      IC = index_two(wf%I, wf%C, wf%n_o)
!
      wf%omega1(wf%A, wf%I) = wf%omega1(wf%A, wf%I) + two*L(JB, KC) ! 1
      wf%omega1(wf%B, wf%J) = wf%omega1(wf%B, wf%J) + two*L(IA, KC) ! 2
      wf%omega1(wf%C, wf%K) = wf%omega1(wf%C, wf%K) + two*L(IA, JB) ! 3
      wf%omega1(wf%C, wf%I) = wf%omega1(wf%C, wf%I) - L(JB, KA)     ! 4
      wf%omega1(wf%B, wf%I) = wf%omega1(wf%B, wf%I) - L(KC, JA)     ! 5
      wf%omega1(wf%C, wf%J) = wf%omega1(wf%C, wf%J) - L(IA, KB)     ! 6
      wf%omega1(wf%A, wf%J) = wf%omega1(wf%A, wf%J) - L(KC, IB)     ! 7
      wf%omega1(wf%B, wf%K) = wf%omega1(wf%B, wf%K) - L(IA, JC)     ! 8
      wf%omega1(wf%A, wf%K) = wf%omega1(wf%A, wf%K) - L(JB, IC)     ! 9
!
   end subroutine omega_sccsd_a1_sccsd
!
!
end submodule omega