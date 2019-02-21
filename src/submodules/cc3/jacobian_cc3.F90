submodule (cc3_class) jacobian_cc3
!
!!
!!    Jacobian submodule (cc3)
!!    Written by Rolf H. Myhre and Alexander Paul
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | ν >.
!!
!
   implicit none
!
!
contains
!
!
   !module subroutine jacobian_transform_trial_vector_cc3(wf, omega, c_i)
!!!
!! !   Jacobian transform trial vector
!! !   Written by Rolf H. Myhre and Alexander Paul
!!!
   !   class(cc3), intent(in) :: wf
!!
   !   real(dp), intent(in) :: omega
!!
   !   real(dp), dimension(wf%n_amplitudes, 1) :: c_i
!!
   !   call wf%jacobian_cc3_transformation(omega, c_i)
!!
   !end subroutine jacobian_transform_trial_vector_cc3
!
!
   module subroutine jacobian_cc3_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transformation (CC3)
!!    Written by Rolf H. Myhre and Alexander Paul
!!
!!    Directs the transformation by the CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl)
!!
!!    On exit, c is overwritten by rho. That is, c_a_i = rho_a_i,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_abij
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
      type(timings) :: cc3_timer
      type(timings) :: ccsd_timer
!
      call cc3_timer%init('CC3 contribution)')
      call ccsd_timer%init('CCSD contribution)')
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      rho_ai = zero
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai, 1)
!
         enddo
      enddo
!$omp end parallel do
!
!     :: CCS contributions to the singles c vector ::
!
      call ccsd_timer%start()
!
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_ccsd_a1(rho_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c_aibj(a,i,b,j) = c(wf%n_o*wf%n_v + aibj, 1)
                     c_aibj(b,j,a,i) = c(wf%n_o*wf%n_v + aibj, 1)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!$omp parallel do schedule(static) private(a,i)
      do a = 1, wf%n_v
         do i = 1, wf%n_v
!
         c_aibj(a,i,a,i) = two*c_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
!
      call wf%jacobian_ccsd_b1(rho_ai, c_aibj)
      call wf%jacobian_ccsd_c1(rho_ai, c_aibj)
      call wf%jacobian_ccsd_d1(rho_ai, c_aibj)
!
!     :: CCSD contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      rho_aibj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_ccsd_a2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_b2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_c2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_d2(rho_aibj, c_ai)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_f2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_g2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_h2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_i2(rho_aibj, c_aibj)
!
!     Last two terms are already symmetric (J2 and K2). Perform the symmetrization
!     rho_ai_bj = P_ij^ab rho_ai_bj now, for convenience
!
      call symmetric_sum(rho_aibj, (wf%n_v)*(wf%n_o))
!
!     In preparation for last two terms, reorder
!     rho_aibj to rho_abij, and c_aibj to c_abij
!
      call mem%alloc(rho_abij,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
      call mem%alloc(c_abij,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(rho_aibj, rho_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_ccsd_j2(rho_abij, c_abij)
      call wf%jacobian_ccsd_k2(rho_abij, c_abij)
!
!     Done with reordered doubles c; deallocate
!
      call mem%dealloc(c_abij,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
!     Order rho_ab_ij back into rho_ai_bj & divide by
!     the biorthonormal factor 1 + delta_ai,bj
!
      call mem%alloc(rho_aibj,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
      call sort_1234_to_1324(rho_abij, rho_aibj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!$omp parallel do schedule(static) private(ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            rho_aibj(a,i,a,i) = half*rho_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
   !  dont want to delete rho_aibj
   !   call mem%dealloc(rho_abij,wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call ccsd_timer%freeze()
      call ccsd_timer%switch_off()
!
!     :: CC3 contributions to the transformed singles and doubles vector ::
!
      call cc3_timer%start()
      call wf%jacobian_cc3_A(rho_ai, rho_aibj)
      call cc3_timer%freeze()
      call cc3_timer%switch_off()
!
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Overwrite the incoming singles c vector for exit
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai, 1) = rho_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do

      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     Overwrite the incoming doubles c vector & pack in
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj, 1) = rho_aibj(a,i,b,j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_aibj,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
   end subroutine jacobian_cc3_transformation_cc3
!
!
   module subroutine jacobian_cc3_A_cc3(wf, omega, c_abji, rho_ai, rho_aibj)
!!
!!    CC3 jacobian terms
!!    Alex C. Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abji
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_aibj
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
      real(dp), dimension(:,:,:), allocatable :: v_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abji
!
      real(dp), dimension(:,:), allocatable :: F_kc ! Transpose the fock matrix sub block
      real(dp), dimension(:,:), allocatable :: F_kc_c1
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_klic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_iljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ilkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlkc_p => null()
!
!     L_kbjc and L_jbkc each other's transpose,
!     but utilising this makes the code more complicated and
!     error prone without any huge advantages
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_kbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_kbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_kbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_kbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
!     C1 transformed integrals
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlic_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_klic_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kljc_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_iljc_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ilkc_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlkc_c1_p => null()
!
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: current_i_batch, current_j_batch, current_k_batch
      integer :: req_0, req_1, req_2, req_3
      real(dp)     :: batch_buff = 0.0
!
!     Set up required integrals
      call wf%jacobian_cc3_integrals()
      call wf%jacobian_cc3_c1_integrals(c_ai)
!
      call mem%alloc(t_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%alloc(u_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%alloc(v_abc,wf%n_v,wf%n_v,wf%n_v)
!
      call mem%alloc(F_kc,wf%n_v,wf%n_o)
      call sort_12_to_21(wf%fock_ia,F_kc,wf%n_o,wf%n_v)
!
!     C1 transformed Fock matrix
!
      call mem%alloc(F_kc_c1, wf%n_v, wf%n_o)
      call wf%jacobian_cc3_construct_fock_ia_c1_cc3(c_ai, F_kc_c1)
!
      call mem%alloc(t_abji,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
      call squareup_and_sort_1234_to_1342(wf%t2,t_abji,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
!     Arrays for the amplitudes and intermediates
!
      t_abc = zero
      u_abc = zero
      v_abc = zero
!
      req_0 = 0
      req_1 = 4*(wf%n_v)**3
      req_2 = 4*(wf%n_o)*(wf%n_v) + 2*(wf%n_v)**2
      req_3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, batch_buff)
!
!     Allocate integral arrays and assign pointers.
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%alloc(L_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_dbic_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_jlic_c1, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
!        Ordered such that batching indices at the end
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%alloc(L_jbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_kbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_kbjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_dbic_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbjc_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbkc_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%alloc(g_jlic_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_klic_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_kljc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_iljc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_ilkc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_jlkc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
      endif
!
!
      call disk%open_file(wf%g_bdck_t,'read')
      call disk%open_file(wf%g_ljck_t,'read')
      call disk%open_file(wf%g_dbkc_t,'read')
      call disk%open_file(wf%g_jlkc_t,'read')
      call disk%open_file(wf%L_jbkc_t,'read')
!
      call disk%open_file(wf%g_bdck_c1,'read')
      call disk%open_file(wf%g_ljck_c1,'read')
      call disk%open_file(wf%g_dbkc_c1,'read')
      call disk%open_file(wf%g_jlkc_c1,'read')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call wf%jacobian_cc3_vvv_reader(batch_i, g_bdci, g_dbic, g_bdci_c1, g_dbic_c1)
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
         g_bdci_c1_p => g_bdci_c1
         g_dbic_c1_p => g_dbic_c1
!
         do current_j_batch = 1, current_i_batch
!
            call batch_j%determine_limits(current_j_batch)
!
            call wf%jacobian_cc3_ov_vv_reader(batch_j, batch_i, g_ljci, g_jlic, L_jbic, g_ljci_c1, g_jlic_c1)
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_jbic_p => L_jbic
!
            g_ljci_c1_p => g_ljci_c1
            g_jlic_c1_p => g_jlic_c1
!
            if (current_j_batch .ne. current_i_batch) then
!
               call wf%jacobian_cc3_vvv_reader(batch_j, g_bdcj, g_dbjc, g_bdcj_c1, g_dbjc_c1)
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
               g_bdcj_c1_p => g_bdcj_c1
               g_dbjc_c1_p => g_dbjc_c1
!
               call wf%jacobian_cc3_ov_vv_reader(batch_i, batch_j, g_licj, g_iljc, L_ibjc, g_licj_c1, g_iljc_c1)
               g_licj_p => g_licj
               g_iljc_p => g_iljc
               L_ibjc_p => L_ibjc
!
               g_licj_c1_p => g_licj_c1
               g_iljc_c1_p => g_iljc_c1
!
            else
!
               g_bdcj_p => g_bdci
               g_dbjc_p => g_dbic
!
               g_licj_p => g_ljci
               g_iljc_p => g_jlic
               L_ibjc_p => L_jbic
!
               g_bdcj_c1_p => g_bdci_c1
               g_dbjc_c1_p => g_dbic_c1
!
               g_licj_c1_p => g_ljci_c1
               g_iljc_c1_p => g_jlic_c1
!
            endif
!
            do current_k_batch = 1, current_j_batch
!
               call batch_k%determine_limits(current_k_batch)
!
               if (current_k_batch .ne. current_i_batch .and. current_k_batch .ne. current_j_batch) then
!
                  call wf%jacobian_cc3_vvv_reader(batch_k, g_bdck, g_dbkc, g_bdck_c1, g_dbkc_c1)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
!
                  g_bdck_c1_p => g_bdck_c1
                  g_dbkc_c1_p => g_dbkc_c1
!
                  call wf%jacobian_cc3_ov_vv_reader(batch_k, batch_i, g_lkci, g_klic, L_kbic, g_lkci_c1, g_klic_c1)
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
                  L_kbic_p => L_kbic
!
                  g_lkci_c1_p => g_lkci_c1
                  g_klic_c1_p => g_klic_c1
!
                  call wf%jacobian_cc3_ov_vv_reader(batch_i, batch_k, g_lick, g_ilkc, L_ibkc, g_lick_c1, g_ilkc_c1)
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  g_lick_c1_p => g_lick_c1
                  g_ilkc_c1_p => g_ilkc_c1
!
                  call wf%jacobian_cc3_ov_vv_reader(batch_k, batch_j, g_lkcj, g_kljc, L_kbjc, g_lkcj_c1, g_kljc_c1)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  g_lkcj_c1_p => g_lkcj_c1
                  g_kljc_c1_p => g_kljc_c1
!
                  call wf%jacobian_cc3_ov_vv_reader(batch_j, batch_k, g_ljck, g_jlkc, L_jbkc, g_ljck_c1, g_jlkc_c1)
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
!
                  g_ljck_c1_p => g_ljck_c1
                  g_jlkc_c1_p => g_jlkc_c1
!
               else if (current_k_batch .eq. current_i_batch) then
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
!
                  g_bdck_c1_p => g_bdci_c1
                  g_dbkc_c1_p => g_dbic_c1
!
                  if (current_j_batch .eq. current_i_batch) then
!
                     g_lkci_p => g_ljci
                     g_klic_p => g_jlic
                     L_kbic_p => L_jbic
!
                     g_lkci_c1_p => g_ljci_c1
                     g_klic_c1_p => g_jlic_c1
!
                     g_lick_p => g_ljci
                     g_ilkc_p => g_jlic
                     L_ibkc_p => L_jbic
!
                     g_lick_c1_p => g_ljci_c1
                     g_ilkc_c1_p => g_jlic_c1
!
                     g_lkcj_p => g_ljci
                     g_kljc_p => g_jlic
                     L_kbjc_p => L_jbic
!
                     g_lkcj_c1_p => g_ljci_c1
                     g_kljc_c1_p => g_jlic_c1
!
                     g_ljck_p => g_ljci
                     g_jlkc_p => g_jlic
                     L_jbkc_p => L_jbic
!
                     g_ljck_c1_p => g_ljci_c1
                     g_jlkc_c1_p => g_jlic_c1
!
                  else
!
                     call wf%jacobian_cc3_ov_vv_reader(batch_k, batch_i, g_lkci, g_klic, L_kbic, g_lkci_c1, g_klic_c1)
                     g_lkci_p => g_lkci
                     g_klic_p => g_klic
                     L_kbic_p => L_kbic
!
                     g_lkci_c1_p => g_lkci_c1
                     g_klic_c1_p => g_klic_c1
!
                     g_lick_p => g_lkci
                     g_ilkc_p => g_klic
                     L_ibkc_p => L_kbic
!
                     g_lick_c1_p => g_lkci_c1
                     g_ilkc_c1_p => g_klic_c1
!
                     g_lkcj_p => g_licj
                     g_kljc_p => g_iljc
                     L_kbjc_p => L_ibjc
!
                     g_lkcj_c1_p => g_licj_c1
                     g_kljc_c1_p => g_iljc_c1
!
                     g_ljck_p => g_ljci
                     g_jlkc_p => g_jlic
                     L_jbkc_p => L_jbic
!
                     g_ljck_c1_p => g_ljci_c1
                     g_jlkc_c1_p => g_jlic_c1
!
                  endif
!
               else if (current_k_batch .eq. current_j_batch) then
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
!
                  g_bdck_c1_p => g_bdcj_c1
                  g_dbkc_c1_p => g_dbjc_c1
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
                  L_kbic_p => L_jbic
!
                  g_lkci_c1_p => g_ljci_c1
                  g_klic_c1_p => g_jlic_c1
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
!
                  g_lick_c1_p => g_licj_c1
                  g_ilkc_c1_p => g_iljc_c1
!
                  call wf%jacobian_cc3_ov_vv_reader(batch_k, batch_j, g_lkcj, g_kljc, L_kbjc, g_lkcj_c1, g_kljc_c1)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  g_lkcj_c1_p => g_lkcj_c1
                  g_kljc_c1_p => g_kljc_c1
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_kbjc
!
                  g_ljck_c1_p => g_lkcj_c1
                  g_jlkc_c1_p => g_kljc_c1
!
               endif
!
               do i = batch_i%first,batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first,min(batch_j%last,i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first,min(batch_k%last,j)
!
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct C^{abc}_{ijk} for given i, j, k (t_abc)
!                       and calculate contributions to rho1 and rho2
!
                        call wf%jacobian_cc3_c3_calc(omega, i, j, k, t_abc, u_abc, t_abji, c_abji  &
                                                      g_bdci_p(:,:,:,i_rel),                       &
                                                      g_bdcj_p(:,:,:,j_rel),                       &
                                                      g_bdck_p(:,:,:,k_rel),                       &
                                                      g_ljci_p(:,:,j_rel,i_rel),                   &
                                                      g_lkci_p(:,:,k_rel,i_rel),                   &
                                                      g_lkcj_p(:,:,k_rel,j_rel),                   &
                                                      g_licj_p(:,:,i_rel,j_rel),                   &
                                                      g_lick_p(:,:,i_rel,k_rel),                   &
                                                      g_ljck_p(:,:,j_rel,k_rel),                   &
                                                      g_bdci_c1_p(:,:,:,i_rel),                    &
                                                      g_bdcj_c1_p(:,:,:,j_rel),                    &
                                                      g_bdck_c1_p(:,:,:,k_rel),                    &
                                                      g_ljci_c1_p(:,:,j_rel,i_rel),                &
                                                      g_lkci_c1_p(:,:,k_rel,i_rel),                &
                                                      g_lkcj_c1_p(:,:,k_rel,j_rel),                &
                                                      g_licj_c1_p(:,:,i_rel,j_rel),                &
                                                      g_lick_c1_p(:,:,i_rel,k_rel),                &
                                                      g_ljck_c1_p(:,:,j_rel,k_rel))
!
      !                  call wf%jacobian_cc3_eps(i, j, k, t_abc)
!
                        call wf%jacobian_cc3_rho1(i, j, k, t_abc, u_abc, rho1,   &
                                                   L_jbic_p(:,:,j_rel,i_rel),    &
                                                   L_kbic_p(:,:,k_rel,i_rel),    &
                                                   L_kbjc_p(:,:,k_rel,j_rel),    &
                                                   L_ibjc_p(:,:,i_rel,j_rel),    &
                                                   L_ibkc_p(:,:,i_rel,k_rel),    &
                                                   L_jbkc_p(:,:,j_rel,k_rel))
!
                        call wf%jacobian_cc3_rho2(i, j, k, t_abc, u_abc, v_abc, rho2, F_kc,  &
                                                      g_dbic_p(:,:,:,i_rel),                 &
                                                      g_dbjc_p(:,:,:,j_rel),                 &
                                                      g_dbkc_p(:,:,:,k_rel),                 &
                                                      g_jlic_p(:,:,j_rel,i_rel),             &
                                                      g_klic_p(:,:,k_rel,i_rel),             &
                                                      g_kljc_p(:,:,k_rel,j_rel),             &
                                                      g_iljc_p(:,:,i_rel,j_rel),             &
                                                      g_ilkc_p(:,:,i_rel,k_rel),             &
                                                      g_jlkc_p(:,:,j_rel,k_rel))
!
!                       Construct t^{abc}_{ijk} for given i, j, k (t_abc)
!                       and calculate contributions to rho1 and rho2
!
                        call wf%jacobian_cc3_t3_calc(omega, i, j, k, t_abc, u_abc, t_abji,   &
                                                      g_bdci_p(:,:,:,i_rel),                 &
                                                      g_bdcj_p(:,:,:,j_rel),                 &
                                                      g_bdck_p(:,:,:,k_rel),                 &
                                                      g_ljci_p(:,:,j_rel,i_rel),             &
                                                      g_lkci_p(:,:,k_rel,i_rel),             &
                                                      g_lkcj_p(:,:,k_rel,j_rel),             &
                                                      g_licj_p(:,:,i_rel,j_rel),             &
                                                      g_lick_p(:,:,i_rel,k_rel),             &
                                                      g_ljck_p(:,:,j_rel,k_rel))
!
      !                  call wf%jacobian_cc3_eps(i, j, k, t_abc)
!
                        call wf%jacobian_cc3_rho2(i, j, k, t_abc, u_abc, v_abc, rho2, F_kc_c1,  &
                                                      g_dbic_c1_p(:,:,:,i_rel),                 &
                                                      g_dbjc_c1_p(:,:,:,j_rel),                 &
                                                      g_dbkc_c1_p(:,:,:,k_rel),                 &
                                                      g_jlic_c1_p(:,:,j_rel,i_rel),             &
                                                      g_klic_c1_p(:,:,k_rel,i_rel),             &
                                                      g_kljc_c1_p(:,:,k_rel,j_rel),             &
                                                      g_iljc_c1_p(:,:,i_rel,j_rel),             &
                                                      g_ilkc_c1_p(:,:,i_rel,k_rel),             &
                                                      g_jlkc_c1_p(:,:,j_rel,k_rel))
!
                     enddo ! loopg over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
!     Close files
!
      call disk%close_file(wf%g_bdck_t)
      call disk%close_file(wf%g_ljck_t)
      call disk%close_file(wf%g_dbkc_t)
      call disk%close_file(wf%g_jlkc_t)
      call disk%close_file(wf%L_jbkc_t)
!
      call disk%close_file(wf%g_bdck_c1)
      call disk%close_file(wf%g_ljck_c1)
      call disk%close_file(wf%g_dbkc_c1)
      call disk%close_file(wf%g_jlkc_c1)
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(L_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%dealloc(L_jbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_kbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_kbjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_dbic_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbjc_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbkc_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_jlic_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_klic_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_kljc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_iljc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_ilkc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_jlkc_c1, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
      endif
!
      call mem%dealloc(F_kc, wf%n_v, wf%n_o)
      call mem%dealloc(F_kc_c1, wf%n_v, wf%n_o)
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine jacobian_cc3_A_cc3
!
!
   module subroutine jacobian_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Jacobian and store on disk
!!    (bd|ck) ordered as dbc,k
!!    (bd|kc) ordered as dcb,k
!!    (lj|ck) ordered as lc,jk
!!    (lj|kc) ordered as cj,lk
!!    (jb|kc) stored as L_jbkc = 2g_jbkc - g_jckb ordered as bc,jk
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs   ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs   ! Array for sorted integrals
      real(dp), dimension(:,:), allocatable     :: v2_help  ! Help array for constructing L_jbkc
!
      integer :: k, j, record
      type(batching_index) :: batch_k
!
      integer :: req_0, req_k
      integer :: current_k_batch
!
      integer :: ioerror=-1
!
      call mem%alloc(v2_help,wf%n_v,wf%n_v)
!
      call batch_k%init(wf%n_o)
!
!     (bd|ck)
!
      req_0 = wf%integrals%n_J*wf%n_v**2
      req_k = 2*wf%n_v**3 + wf%integrals%n_J*wf%n_v
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%g_bdck_t%init('g_bdck_t','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_bdck_t,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call wf%get_vvvo(g_pqrs, &
                           1,wf%n_v, &
                           1,wf%n_v, &
                           1,wf%n_v, &
                           batch_k%first,batch_k%last)
!
         call sort_1234_to_2134(g_pqrs,h_pqrs,wf%n_v,wf%n_v,wf%n_v,batch_k%length)
!
         do k = 1,batch_k%length
!
            record = batch_k%first + k -1
            write(wf%g_bdck_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write bdck_t file')
         endif
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_bdck_t,'keep')
!
!
!     (db|kc)
!     Same batching
!
      call wf%g_dbkc_t%init('g_dbkc_t','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_dbkc_t,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, batch_k%length, wf%n_v)
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call wf%get_vvov(g_pqrs, &
                           1,wf%n_v, &
                           1,wf%n_v, &
                           batch_k%first,batch_k%last, &
                           1,wf%n_v)
!
         call sort_1234_to_2413(g_pqrs,h_pqrs,wf%n_v,wf%n_v,batch_k%length,wf%n_v)
!
         do k = 1,batch_k%length
!
            record = batch_k%first + k -1
            write(wf%g_dbkc_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write dbkc_t file')
         endif
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, batch_k%length, wf%n_v)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_dbkc_t,'keep')
!
!
!     (lj|ck)
!
      req_0 = wf%integrals%n_J*wf%n_o**2
      req_k = 2*wf%n_o**2*wf%n_v + wf%integrals%n_J*wf%n_v
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%g_ljck_t%init('g_ljck_t','direct','unformatted',dp*wf%n_v*wf%n_o)
      call disk%open_file(wf%g_ljck_t,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
         call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o ,batch_k%length)
!
         call wf%get_oovo(g_pqrs, &
                           1,wf%n_o, &
                           1,wf%n_o, &
                           1,wf%n_v, &
                           batch_k%first,batch_k%last)
!
         call sort_1234_to_1324(g_pqrs,h_pqrs,wf%n_o,wf%n_o,wf%n_v,batch_k%length)
!
         do k = 1,batch_k%length
            do j = 1,wf%n_o
!
               record  = (batch_k%first + k - 2)*wf%n_o + j
               write(wf%g_ljck_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,j,k)
!
            enddo
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write ljck_t file')
         endif

         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_ljck_t,'keep')
!
!
!     (jl|kc)
!     Same batching
!
      call wf%g_jlkc_t%init('g_jlkc_t','direct','unformatted',dp*wf%n_v*wf%n_o)
      call disk%open_file(wf%g_jlkc_t,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
         call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_k%length)
!
         call wf%get_ooov(g_pqrs, &
                           1,wf%n_o, &
                           1,wf%n_o, &
                           batch_k%first,batch_k%last, &
                           1,wf%n_v)
!
         call sort_1234_to_4213(g_pqrs,h_pqrs,wf%n_o,wf%n_o,batch_k%length,wf%n_v)
!
         do k = 1,batch_k%length
            do j = 1,wf%n_o
!
               record  = (batch_k%first + k - 2)*wf%n_o + j
               write(wf%g_jlkc_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,j,k)
!
            enddo
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write jlkc_t file')
         endif

         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, batch_k%length, wf%n_v)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_jlkc_t,'keep')
!
!
!     (jb|kc)
!
      req_0 = wf%integrals%n_J*wf%n_o*wf%n_v
      req_k = 2*wf%n_v**2*wf%n_o + wf%integrals%n_J*wf%n_v
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%L_jbkc_t%init('L_jbkc_t','direct','unformatted',dp*wf%n_v**2)
      call disk%open_file(wf%L_jbkc_t,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_v, batch_k%length, wf%n_v)
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
!
         call wf%get_ovov(g_pqrs, &
                           1,wf%n_o, &
                           1,wf%n_v, &
                           batch_k%first,batch_k%last, &
                           1,wf%n_v)
!
         call sort_1234_to_2413(g_pqrs,h_pqrs,wf%n_o,wf%n_v,batch_k%length,wf%n_v)
!
         do k = 1,batch_k%length
            do j = 1,wf%n_o
!
               call sort_12_to_21(h_pqrs(:,:,j,k), v2_help, wf%n_v, wf%n_v)
!
               call dscal(wf%n_v**2, two, h_pqrs(:,:,j,k),1)
!
               call daxpy(wf%n_v**2, -one, v2_help, 1, h_pqrs(:,:,j,k), 1)
!
               record  = (batch_k%first + k - 2)*wf%n_o + j
               write(wf%L_jbkc_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,j,k)
!
            enddo
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write jbkc_t file')
         endif

         call mem%dealloc(g_pqrs, wf%n_o, wf%n_v, batch_k%length, wf%n_v)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%L_jbkc_t,'keep')
!
!
   end subroutine jacobian_cc3_integrals_cc3
!
!
   module subroutine jacobian_cc3_c1_integrals_cc3(wf, c_ai)
!!
!!    Construct c1-transformed integrals needed in CC3 jacobian
!!    from the c1-transformed Cholesky Vectors
!!
!!    Should include check if they can be safely stored in memory
!!
!!    g'_bdck = (b'd|ck) + (bd|c'k) - (bd|ck')  ordered as dbc,k
!!    g'_ljck = (lj'|ck) + (lj|ck') - (lj|c'k)  ordered as lc,jk

!!    NB: the indices d and l are part of rho_2 while j and l are summation indices in the following terms
!!    (d'b|kc) ordered as bcd,k
!!    (jl'|kc) orderd as cljk
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      real(dp), dimension(:,:), allocatable :: L_bd_J_c1, L_ck_J_c1, L_lj_J_c1 ! c1 transformed Cholesky vectors
      real(dp), dimension(:,:), allocatable :: L_db_J_c1
      real(dp), dimension(:,:), allocatable :: L_ck_J, L_bd_J, L_kc_J, L_lj_J ! Cholesky vectors
!
      integer :: k, j, record
      type(batching_index) :: batch_k
!
      integer :: req_0, req_k
      integer :: current_k_batch
!
      integer :: ioerror=-1
!
      call batch_k%init(wf%n_o)
!
!     g'_bdck = (b'd|ck) + (bd|c'k) - (bd|ck')
!
      req_0 = (wf%integrals%n_J)*(wf%n_v)**2 + (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req_k = max((wf%integrals%n_J)*(wf%n_v) + (wf%n_v)**3, 2*(wf%n_v)**3)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%g_bdck_c1%init('g_bdck_c1','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_bdck_c1,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
!        :: Term 1: g_b'dck = sum_J L_bd_J_c1 L_ck_J ::
!
         call mem%alloc(L_bd_J_c1, (wf%n_v)*(wf%n_v), wf%integrals%n_J)
!
         call wf%integrals%construct_cholesky_ab_c1(L_bd_J_c1, c_ai, 1, wf%n_v, 1, wf%n_v)
!
         call mem%alloc(L_ck_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call wf%integrals%read_cholesky_ai_t1(L_ck_J, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call dgemm('N', 'T',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_bd_J_c1,                 & ! L_bd_J  b is c1-transformed
                     (wf%n_v)**2,               &
                     L_ck_J,                    & ! L_ck_J
                     (wf%n_v)*(batch_k%length), &
                     zero,                      &
                     g_pqrs,                    & ! (b'd|ck)
                     (wf%n_v)**2)
!
         call mem%dealloc(L_bd_J_c1, (wf%n_v)*(wf%n_v), wf%integrals%n_J)
         call mem%dealloc(L_ck_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
!        :: Term 2: g_bdc'k = sum_J L_bd_J L_ck_J_c1 ::
!
         call mem%alloc(L_ck_J_c1, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call wf%integrals%construct_cholesky_ai_a_c1(L_ck_J_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call mem%alloc(L_bd_J, (wf%n_v)*(wf%n_v), wf%integrals%n_J)
!
         call wf%integrals%read_cholesky_ab_t1(L_bd_J, 1, wf%n_v, 1, wf%n_v)
!
         call dgemm('N', 'T',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_bd_J,                    & ! L_bd_J
                     (wf%n_v)**2,               &
                     L_ck_J_c1,                 & ! L_ck_J  c is c1_transformed
                     (wf%n_v)*(batch_k%length), &
                     one,                       &
                     g_pqrs,                    & ! (bd|c'k)
                     (wf%n_v)**2)
!
!        :: Term 3: g_bdck' = sum_J L_bd_J L_ck_J_c1 ::
!
!
         call wf%integrals%construct_cholesky_ai_i_c1(L_ck_J_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call dgemm('N', 'T',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     -one,                      &
                     L_bd_J,                    & ! L_bd_J
                     (wf%n_v)**2,               &
                     L_ck_J_c1,                 & ! L_ck_J  k is c1_transformed
                     (wf%n_v)*(batch_k%length), &
                     one,                       &
                     g_pqrs,                    & ! (bd|ck')
                     (wf%n_v)**2)
!
         call mem%dealloc(L_ck_J_c1, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
         call mem%dealloc(L_db_J, (wf%n_v)*(wf%n_v), wf%integrals%n_J)
!
!        Sort from g_pqrs = (b'd|ck) + (bd|c'k) - (bd|ck') to h_pqrs sorted as dbck
!
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length) ! order dbck
!
         call sort_1234_to_2134(g_pqrs, h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
!        Write to file
!        Should implement possibility to have them in mem if possible
!
         do k = 1, batch_k%length
!
            record = batch_k%first + k -1
            write(wf%g_bdck_c1%unit, rec=record, iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write bdck_c1 file')
         endif
!
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_bdck_c1,'keep')
!
!
!     (d'b|kc) same batching (req0 and req_k)
!
!
      call wf%g_dbkc_c1%init('g_dbkc_c1','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_dbkc_c1,'write')
!
      do current_k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(L_db_J_c1, (wf%n_v)*(wf%n_v), wf%integrals%n_J)
!
         call wf%integrals%construct_cholesky_ab_c1(L_db_J_c1, c_ai, 1, wf%n_v, 1, wf%n_v)
!
         call mem%alloc(L_kc_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call wf%integrals%read_cholesky_ia_t1(L_kc_J, batch_k%first, batch_k%last, 1, wf%n_v)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, batch_k%length, wf%n_v)
!
         call dgemm('N', 'T',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_db_J_c1,                 & ! L_db_J  d is c1-transformed
                     (wf%n_v)**2,               &
                     L_kc_J,                    & ! L_kc_J
                     (wf%n_v)*(batch_k%length), &
                     zero,                      &
                     g_pqrs,                    & ! (d'b|kc)
                     (wf%n_v)**2)
!
         call mem%dealloc(L_db_J_c1, (wf%n_v)*(wf%n_v), wf%integrals%n_J)
         call mem%dealloc(L_kc_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length) ! order bcd,k
!
         call sort_1234_to_2413(g_pqrs,h_pqrs,wf%n_v,wf%n_v,batch_k%length,wf%n_v)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, batch_k%length, wf%n_v)
!
!        Write to file
!        Should implement possibility to have them in mem if possible
!
         do k = 1,batch_k%length
!
            record = batch_k%first + k -1
            write(wf%g_dbkc_c1%unit, rec=record, iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write dbkc_c1 file')
         endif
!
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_dbkc_c1,'keep')
!
!
!     g'_ljck = (lj|'ck) + (lj|ck') - (lj|c'k) ordered as lc,jk
!
!
      req_0 = wf%integrals%n_J*wf%n_o**2 + 2*(wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req_k = max(2*(wf%n_v)*(wf%n_o)**2, (wf%n_v)*(wf%n_o)**2 + wf%integrals%n_J*wf%n_v)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%g_ljck_c1%init('g_ljck_c1','direct','unformatted',dp*wf%n_v*wf%n_o)
      call disk%open_file(wf%g_ljck_c1,'write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
!        :: Term 1: g_lj'ck = sum_J L_jl_J_c1 L_ck_J ::
!
         call mem%alloc(L_lj_J_c1, (wf%n_o)*(wf%n_o), wf%integrals%n_J)
!
         call wf%integrals%construct_cholesky_ij_c1(L_lj_J_c1, c_ai, 1, wf%n_o, 1, wf%n_o)
!
         call mem%alloc(L_ck_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call wf%integrals%read_cholesky_ai_t1(L_ck_J, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         call dgemm('N', 'T',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_lj_J_c1,                 & ! L_lj_J  j is c1-transformed
                     (wf%n_o)**2,               &
                     L_ck_J,                    & ! L_ck_J
                     (wf%n_v)*(batch_k%length), &
                     zero,                      &
                     g_pqrs,                    & ! (lj'|ck)
                     (wf%n_o)**2)
!
         call mem%dealloc(L_lj_J_c1, (wf%n_o)*(wf%n_o), wf%integrals%n_J)
         call mem%dealloc(L_ck_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
!        :: Term 2: g_ljck' = sum_J L_lj_J L_ck_J_c1 ::
!
         call mem%alloc(L_ck_J_c1, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call wf%integrals%construct_cholesky_ai_i_c1(L_ck_J_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call mem%alloc(L_lj_J, (wf%n_o)*(wf%n_o), wf%integrals%n_J)
!
         call wf%integrals%read_cholesky_ij_t1(L_lj_J, 1, wf%n_o, 1, wf%n_o)
!
         call dgemm('N', 'T',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_lj_J,                    & ! L_lj_J
                     (wf%n_o)**2,               &
                     L_ck_J_c1,                 & ! L_ck_J  k is c1_transformed
                     (wf%n_v)*(batch_k%length), &
                     one,                       &
                     g_pqrs,                    & ! (lj|ck')
                     (wf%n_o)**2)
!
!        :: Term 3: g_bdck' = sum_J L_bd_J L_ck_J_c1 ::
!
!
         call wf%integrals%construct_cholesky_ai_a_c1(L_ck_J_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call dgemm('N', 'T',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     -one,                      &
                     L_lj_J,                    & ! L_lj_J
                     (wf%n_o)**2,               &
                     L_ck_J_c1,                 & ! L_ck_J  c is c1_transformed
                     (wf%n_v)*(batch_k%length), &
                     one,                       &
                     g_pqrs,                    & ! (lj|c'k)
                     (wf%n_o)**2)
!
         call mem%dealloc(L_ck_J_c1, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
         call mem%dealloc(L_lj_J, (wf%n_o)*(wf%n_o), wf%integrals%n_J)
!
!        Sort from g_pqrs = (lj'|ck) + (lj|ck') - (lj|c'k) to h_pqrs
!
         call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o ,batch_k%length) ! lcjk
!
         call sort_1234_to_1324(g_pqrs,h_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         do k = 1,batch_k%length
            do j = 1,wf%n_o
!
               record  = (batch_k%first + k - 2)*wf%n_o + j
               write(wf%g_ljck_c1%unit,rec=record,iostat=ioerror) h_pqrs(:,:,j,k)
!
            enddo
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write ljck_c file')
         endif
!
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_ljck_c1,'keep')
!
!
!     (jl'|kc) same batching (req0 and req_k)
!
!
      call wf%g_jlkc_c1%init('g_jlkc_c1','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_jlkc_c1,'write')
!
      do current_k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(L_jl_J_c1, (wf%n_o)*(wf%n_o), wf%integrals%n_J)
!
         call wf%integrals%construct_cholesky_ij_c1(L_jl_J_c1, c_ai, 1, wf%n_o, 1, wf%n_o)
!
         call mem%alloc(L_kc_J, (wf%n_v)*(batch_k%length), wf%integrals%n_J)
!
         call wf%integrals%read_cholesky_ia_t1(L_kc_J, batch_k%first, batch_k%last, 1, wf%n_v)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, batch_k%length, wf%n_v)
!
         call dgemm('N', 'T',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_jl_J_c1,                 & ! L_jl_J  l is c1-transformed
                     (wf%n_o)**2,               &
                     L_kc_J,                    & ! L_kc_J
                     (wf%n_v)*(batch_k%length), &
                     zero,                      &
                     g_pqrs,                    & ! (jl'|kc)
                     (wf%n_o)**2)
!
         call mem%dealloc(L_jl_J_c1, (wf%n_o)*(wf%n_o), wf%integrals%n_J)
         call mem%dealloc(L_kc_J, (wf%n_o)*(batch_k%length), wf%integrals%n_J)
!
         call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_k%length) ! order cl,jk
!
         call sort_1234_to_2413(g_pqrs, h_pqrs, wf%n_o, wf%n_o, batch_k%length, wf%n_v)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, batch_k%length, wf%n_v)
!
!        Write to file
!        Should implement possibility to have them in mem if possible
!
         do k = 1,batch_k%length
!
            record = batch_k%first + k -1
            write(wf%g_jlkc_c1%unit, rec=record, iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write jlkc_c file')
         endif
!
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_jlkc_c1, 'keep')
!
!
   end subroutine jacobian_cc3_c1_integrals_cc3
!
!
   subroutine jacobian_cc3_construct_fock_ia_c1_cc3(wf, c_ai, F_ia_c1)
!!
!!    Calculates C1 transformed elements of the Fock matrix required for the CC3 jacobian
!!    and returns as n_v, n_o
!!
!!    F_ia_c1 = sum_j L_iajj' = sum_j 2 g_iajj' - g_ij'ja
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: F_ia_c1
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: L_ij_J_c1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajj, g_ijja
!
      integer :: i, a, j
!
!     Construct the integrals from the Cholesky Vectors
!
      call mem%alloc(L_ia_J, (wf%n_o)*(wf%n_v), wf%integrals%n_J)
!
      call wf%integrals%read_cholesky_ia_t1(L_ia_J, 1, wf%n_o, 1, wf%n_v)
!
      call mem%alloc(L_ij_J_c1, (wf%n_o)**2, wf%integrals%n_J)
!
      call wf%integrals%construct_cholesky_ij_c1(L_ij_J_c1, c_ai, 1, wf%n_o, 1, wf%n_o)
!
      call mem%alloc(g_iajj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_o)**2,         &
                  wf%integrals%n_J,    &
                  one,                 &
                  L_ia_J,              & ! L_ia_J
                  (wf%n_v)*(wf%n_o),   &
                  L_ij_J_c1,           & ! L_jj'_J
                  (wf%n_o)**2,         &
                  zero,                &
                  g_iajj,              & ! (ia|jj')
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(g_ijja, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T',             &
                  (wf%n_o)**2,         &
                  (wf%n_v)*(wf%n_o),   &
                  wf%integrals%n_J,    &
                  one,                 &
                  L_ij_J_c1,           & ! L_ij'_J
                  (wf%n_o)**2,         &
                  L_ia_J,              & ! L_ja_J
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  g_ijja,              & ! (ia|jj')
                  (wf%n_o)**2)
!
      call mem%dealloc(L_ij_J_c1, (wf%n_o)**2, wf%integrals%n_J)
      call mem%dealloc(L_ia_J, (wf%n_o)*(wf%n_v), wf%integrals%n_J)
!
!     Add contributions and resort to F_ia_c1(a,i)
!
!$omp parallel do private(a,i,j)
      do a = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               F_ia_c1(a,i) = two*g_iajj(i,a,j,j) - g_ijja(i,j,j,a)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_ijja, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(g_iajj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!
   end subroutine jacobian_cc3_construct_fock_ia_c1_cc3
!
!
   module subroutine jacobian_cc3_vvv_reader_cc3(wf, batch_x, g_bdcx, g_dbxc, g_bdcx_c1, g_dbxc_c1)
!!
!!    Read the bdck and dbkc integrals needed in the current batch
!!    both c1- and t1-transformed integrals
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_bdcx
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_dbxc
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_bdcx_c1
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_dbxc_c1
!
      integer :: ioerror
      integer :: x, x_abs
!
!     t1-transformed integrals
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         read(wf%g_bdck_t%unit,rec=x_abs, iostat=ioerror) g_bdcx(:,:,:,x)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a)') 'Failed to read bdck_t file'
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         read(wf%g_dbkc_t%unit,rec=x_abs, iostat=ioerror) g_dbxc(:,:,:,x)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a)') 'Failed to read dbkc_t file'
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!     c1-transformed integrals
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         read(wf%g_bdck_c1%unit,rec=x_abs, iostat=ioerror) g_bdcx_c1(:,:,:,x)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a)') 'Failed to read bdck_c file'
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         read(wf%g_dbkc_c1%unit,rec=x_abs, iostat=ioerror) g_dbxc_c1(:,:,:,x)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a)') 'Failed to read dbkc_c file'
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!
   end subroutine jacobian_cc3_vvv_reader_cc3
!
!
   module subroutine jacobian_cc3_ov_vv_reader_cc3(wf, batch_y, batch_x, g_lycx, g_ylxc, L_ybxc, &
                                                   g_lycx_c1, g_ylxc_c1)
!!
!!    Read the ljck, jlkc and jbkc integrals needed in the current batches
!!    both c1- and t1-transformed integrals
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      real(dp), dimension(wf%n_o,wf%n_v,batch_y%length,batch_x%length), intent(out) :: g_lycx
      real(dp), dimension(wf%n_v,wf%n_o,batch_y%length,batch_x%length), intent(out) :: g_ylxc
      real(dp), dimension(wf%n_v,wf%n_v,batch_y%length,batch_x%length), intent(out) :: L_ybxc
      real(dp), dimension(wf%n_o,wf%n_v,batch_y%length,batch_x%length), intent(out) :: g_lycx_c1
      real(dp), dimension(wf%n_v,wf%n_o,batch_y%length,batch_x%length), intent(out) :: g_ylxc_c1
!
      integer :: ioerror, record
      integer :: x, y, x_abs, y_abs
!
!     t1-transformed integrals
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         do y = 1,batch_y%length
!
            y_abs = batch_y%first + y - 1
!
            record = wf%n_o*(x_abs - 1) + y_abs
!
            read(wf%g_ljck_t%unit,rec=record, iostat=ioerror) g_lycx(:,:,y,x)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read ljck_t file'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
      enddo
!
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         do y = 1,batch_y%length
!
            y_abs = batch_y%first + y - 1
!
            record = wf%n_o*(x_abs - 1) + y_abs
!
            read(wf%g_jlkc_t%unit,rec=record, iostat=ioerror) g_ylxc(:,:,y,x)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read jlkc file'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
      enddo
!
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         do y = 1,batch_y%length
!
            y_abs = batch_y%first + y - 1
!
            record = wf%n_o*(x_abs - 1) + y_abs
!
            read(wf%L_jbkc_t%unit,rec=record, iostat=ioerror) L_ybxc(:,:,y,x)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read jbkc file'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
      enddo
!
!     c1-transformed integrals
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         do y = 1,batch_y%length
!
            y_abs = batch_y%first + y - 1
!
            record = wf%n_o*(x_abs - 1) + y_abs
!
            read(wf%g_ljck_c1%unit,rec=record, iostat=ioerror) g_lycx_c1(:,:,y,x)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read ljck_c1 file'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
      enddo
!
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         do y = 1,batch_y%length
!
            y_abs = batch_y%first + y - 1
!
            record = wf%n_o*(x_abs - 1) + y_abs
!
            read(wf%g_jlkc_c1%unit,rec=record, iostat=ioerror) g_ylxc_c1(:,:,y,x)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read jlkc_c1 file'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
      enddo
!
   end subroutine jacobian_cc3_ov_vv_reader_cc3
!
!
   module subroutine jacobian_cc3_c3_calc_cc3(wf, omega, i, j, k, c_abc, u_abc, t_abji, c_abji, &
                                                g_bdci, g_bdcj, g_bdck, g_ljci, g_lkci,         &
                                                g_lkcj, g_licj, g_lick, g_ljck,                 &
                                                g_bdci_c1, g_bdcj_c1, g_bdck_c1, g_ljci_c1,     &
                                                g_lkci_c1, g_lkcj_c1, g_licj_c1, g_lick_c1, g_ljck_c1)
!!
!!    Construct c^abc_ijk amplitudes for the fixed indices i, j, k
!!
!!    c^abc = (omega - ε^abc_ijk)^-1 * P^abc_ijk (sum_d c^ad_ij g_ckbd - sum_l c^ab_il g_cklj
!!             + sum_d t^ad_ij g'_bdck - sum_l t^ab_il g'_cklj
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) intent(in) :: t_abji
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) intent(in) :: c_abji
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdck
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_ljck
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdci_c1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdcj_c1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdck_c1
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_ljci_c1
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lkci_c1
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lkcj_c1
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_licj_c1
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lick_c1
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_ljck_c1
!
      integer :: a, b, c
!
      real(dp) :: epsilon_c, epsilon_cb, epsilon_ijk
!
!     c^ad_ij*(bd|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 c_abji(:,:,j,i),   & ! c^ad_ij
                 wf%n_v,            &
                 g_bdck,            & ! g_dbck
                 wf%n_v,            &
                 zero,              &
                 c_abc,             & ! c^abc_ijk
                 wf%n_v)
!
!     t^ad_ij*g'_bdck
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,j,i),   & ! t^ad_ij
                 wf%n_v,            &
                 g_bdck_c1,         & ! g'_dbck
                 wf%n_v,            &
                 one,               &
                 c_abc,             & ! c^abc_ijk
                 wf%n_v)
!
!     - c^ab_il*(lj|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 c_abji(:,:,:,i),   & ! c^ab_il
                 wf%n_v**2,         &
                 g_ljck,            & ! g_lcjk
                 wf%n_o,            &
                 one,               & ! c^abc_ijk
                 c_abc,             &
                 wf%n_v**2)
!
!     - t^ab_il*g'_ljck
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,i),   & ! t^ab_il
                 wf%n_v**2,         &
                 g_ljck_c1,         & ! g'_lcjk
                 wf%n_o,            &
                 one,               & ! c^abc_ijk
                 c_abc,             &
                 wf%n_v**2)
!
!     c^bd_ji*(ad|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 c_abji(:,:,i,j),   & ! c^bd_ji
                 wf%n_v,            &
                 g_bdck,            & ! g_dack
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u_bac_ijk
                 wf%n_v)
!
!     t^bd_ji*g'_adck
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,i,j),   & ! t^bd_ij
                 wf%n_v,            &
                 g_bdck_c1,         & ! g_dack
                 wf%n_v,            &
                 one,               &
                 u_abc,             & ! u_bac_jik
                 wf%n_v)
!
!     - c^ba_jl*(li|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 c_abji(:,:,:,j),   & ! c^ba_jl
                 wf%n_v**2,         &
                 g_lick,            & ! g_lcik
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u_bac_jik
                 wf%n_v**2)
!
!     - t^ba_jl*g'_lick
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,j),   & ! t^ba_jl
                 wf%n_v**2,         &
                 g_lick_c1,         & ! g'_lcik
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u_bac_jik
                 wf%n_v**2)
!
      call sort_123_to_213_and_add(u_abc, c_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     c^ad_ik*(cd|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 c_abji(:,:,k,i),   & ! c^ad_ik
                 wf%n_v,            &
                 g_bdcj,            & ! g_dcbj
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u^acb_ikj
                 wf%n_v)
!
!     t^ad_ik*g'_cdbj
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,k,i),   & ! t^ad_ik
                 wf%n_v,            &
                 g_bdcj_c1,         & ! g'_dcbj
                 wf%n_v,            &
                 one,               &
                 u_abc,             & ! u^acb_ikj
                 wf%n_v)
!
!     - c^ac_il*(lk|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 c_abji(:,:,:,i),   & ! c^ac_il
                 wf%n_v**2,         &
                 g_lkcj,            & ! g_lbkj
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^acb_ikj
                 wf%n_v**2)
!
!     - t^ac_il*g'_lkbj
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,i),   & ! t^ac_il
                 wf%n_v**2,         &
                 g_lkcj_c1,         & ! g'_lbkj
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^acb_ikj
                 wf%n_v**2)
!
      call sort_123_to_132_and_add(u_abc, c_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     c^cd_ki*(ad|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 c_abji(:,:,i,k),   & ! c^cd_ki
                 wf%n_v,            &
                 g_bdcj,            & ! g_dabj
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u_cab_kij
                 wf%n_v)
!
!     t^cd_ki*g'_adbj
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,i,k),   & ! t^cd_ki
                 wf%n_v,            &
                 g_bdcj_c1,         & ! g'_dabj
                 wf%n_v,            &
                 one,               &
                 u_abc,             & ! u_cab_kij
                 wf%n_v)
!
!     - c^ca_kl*(li|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 c_abji(:,:,:,k),   & ! c^ca_kl
                 wf%n_v**2,         &
                 g_licj,            & ! g_lbij
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^cab_kij
                 wf%n_v**2)
!
!     - t^ca_kl*g'_libj
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,k),   & ! t^ca_kl
                 wf%n_v**2,         &
                 g_licj_c1,         & ! g'_lbij
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^cab_kij
                 wf%n_v**2)
!
      call sort_123_to_231_and_add(u_abc, c_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     c^bd_jk*(cd|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 c_abji(:,:,k,j),   & ! c^bd_jk
                 wf%n_v,            &
                 g_bdci,            & ! g_dcai
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u^bca_jki
                 wf%n_v)
!
!     t^bd_jk*g'_cd|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,k,j),   & ! t^bd_jk
                 wf%n_v,            &
                 g_bdci_c1,         & ! g'_dcai
                 wf%n_v,            &
                 one,               &
                 u_abc,             & ! u^bca_jki
                 wf%n_v)
!
!     - c^bc_jl*(lk|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 c_abji(:,:,:,j),   & ! c^bc_jl
                 wf%n_v**2,         &
                 g_lkci,            & ! g_laki
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^bca_jki
                 wf%n_v**2)
!
!     - t^bc_jl*g_lkai
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,j),   & ! t^bc_jl
                 wf%n_v**2,         &
                 g_lkci_c1,         & ! g'_laki
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^bca_jki
                 wf%n_v**2)
!
      call sort_123_to_312_and_add(u_abc, c_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     c^cd_kj*(bd|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 c_abji(:,:,j,k),   & ! c^cd_kj
                 wf%n_v,            &
                 g_bdci,            & ! g_dbai
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u^cba_kji
                 wf%n_v)
!
!     t^cd_kj*g'_bdai
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,j,k),   & ! t^cd_kj
                 wf%n_v,            &
                 g_bdci_c1,         & ! g'_dbai
                 wf%n_v,            &
                 one,               &
                 u_abc,             & ! u^cba_kji
                 wf%n_v)
!
!     - c^cb_kl*(lj|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 c_abji(:,:,:,k),   & ! c^cb_kl
                 wf%n_v**2,         &
                 g_ljci,            & ! g_laji
                 wf%n_o,            &
                 one,               & ! u^cba_kji
                 u_abc,             &
                 wf%n_v**2)
!
!     - t^cb_kl*g'_ljai
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,k),   & ! t^cb_kl
                 wf%n_v**2,         &
                 g_ljci_c1,         & ! g_laji
                 wf%n_o,            &
                 one,               & ! u^cba_kji
                 u_abc,             &
                 wf%n_v**2)
!
      call sort_123_to_321_and_add(u_abc, c_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Scale by (omega - ε^abc_ijk)^-1
!
      epsilon_ijk = omega + wf%fock_diagonal(i,1) + wf%fock_diagonal(j,1) + wf%fock_diagonal(k,1)
!
!$omp parallel do schedule(static) private(a)
      do a = 1,wf%n_v
!
         c_abc(a,a,a) = zero
!
      enddo
!$omp end parallel do
!
!$omp parallel do schedule(static) private(c,b,a,epsilon_c,epsilon_cb)
      do c = 1,wf%n_v
!
         epsilon_c = epsilon_ijk - wf%fock_diagonal(wf%n_o + c, 1)
!
         do b = 1,wf%n_v
!
            epsilon_cb = epsilon_c - wf%fock_diagonal(wf%n_o + b, 1)
!
            do a = 1,wf%n_v
!
               c_abc(a,b,c) = c_abc(a,b,c)*one/(epsilon_cb - wf%fock_diagonal(wf%n_o + a, 1))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine jacobian_cc3_c3_calc_cc3
!
!
   module subroutine jacobian_cc3_t3_calc_cc3(omega, i, j, k, t_abc, u_abc, t_abji              &
                                             g_bdci_p, g_bdcj_p, g_bdck_p, g_ljci_p, g_lkci_p,  &
                                             g_lkcj_p, g_licj_p, g_lick_p, g_ljck_p
!!
!!    Construct t^abc_ijk amplitudes for the fixed indices i, j, k
!!
!!    t^abc = (omega - ε^abc_ijk)^-1 * P^abc_ijk (sum_d t^ad_ij (bd|ck) - sum_l t^ab_il (ck|lj))
!!
!!    Alexander Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) intent(in) :: t_abji
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdck
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_ljck
!
      integer :: a, b, c
!
      real(dp) :: epsilon_c, epsilon_cb, epsilon_ijk
!
!     t^ad_ij*(bd|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,j,i),   & ! t^ad_ij
                 wf%n_v,            &
                 g_bdck,            & ! g_dbck
                 wf%n_v,            &
                 zero,              &
                 t_abc,             & ! t^abc_ijk
                 wf%n_v)
!
!     - t^ab_il*(lj|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,i),   & ! t^ab_il
                 wf%n_v**2,         &
                 g_ljck,            & ! g_lcjk
                 wf%n_o,            &
                 one,               & ! t^abc_ijk
                 t_abc,             &
                 wf%n_v**2)
!
!     t^bd_ji*(ad|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,i,j),   & ! t^bd_ji
                 wf%n_v,            &
                 g_bdck,            & ! g_dack
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u_bac_ijk
                 wf%n_v)
!
!     - t^ba_jl*(li|ck)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,j),   & ! t^ba_jl
                 wf%n_v**2,         &
                 g_lick,            & ! g_lcik
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u_bac_jik
                 wf%n_v**2)
!
      call sort_123_to_213_and_add(u_abc, t_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     t^ad_ik*(cd|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,k,i),   & ! t^ad_ik
                 wf%n_v,            &
                 g_bdcj,            & ! g_dcbj
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u^acb_ikj
                 wf%n_v)
!
!     - t^ac_il*(lk|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,i),   & ! t^ac_il
                 wf%n_v**2,         &
                 g_lkcj,            & ! g_lbkj
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^acb_ikj
                 wf%n_v**2)
!
      call sort_123_to_132_and_add(u_abc, t_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     t^cd_ki*(ad|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,i,k),   & ! t^cd_ki
                 wf%n_v,            &
                 g_bdcj,            & ! g_dabj
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u_cab_kij
                 wf%n_v)
!
!     - t^ca_kl*(li|bj)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,k),   & ! t^ca_kl
                 wf%n_v**2,         &
                 g_licj_c1,         & ! g_lbij
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^cab_kij
                 wf%n_v**2)
!
      call sort_123_to_231_and_add(u_abc, t_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     t^bd_jk*(cd|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,k,j),   & ! t^bd_jk
                 wf%n_v,            &
                 g_bdci,            & ! g_dcai
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u^bca_jki
                 wf%n_v)
!
!     - t^bc_jl*(lk|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,j),   & ! t^bc_jl
                 wf%n_v**2,         &
                 g_lkci,            & ! g_laki
                 wf%n_o,            &
                 one,               &
                 u_abc,             & ! u^bca_jki
                 wf%n_v**2)
!
      call sort_123_to_312_and_add(u_abc, t_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     t^cd_kj*(bd|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 t_abji(:,:,j,k),   & ! t^cd_kj
                 wf%n_v,            &
                 g_bdci,            & ! g_dbai
                 wf%n_v,            &
                 zero,              &
                 u_abc,             & ! u^cba_kji
                 wf%n_v)
!
!     - t^cb_kl*(lj|ai)
!
      call dgemm('N', 'N',          &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 wf%n_o,            &
                 -one,              &
                 t_abji(:,:,:,k),   & ! t^cb_kl
                 wf%n_v**2,         &
                 g_ljci,            & ! g_laji
                 wf%n_o,            &
                 one,               & ! u^cba_kji
                 u_abc,             &
                 wf%n_v**2)
!
      call sort_123_to_321_and_add(u_abc, t_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Scale by (omega - ε^abc_ijk)^-1
!
      epsilon_ijk = omega + wf%fock_diagonal(i,1) + wf%fock_diagonal(j,1) + wf%fock_diagonal(k,1)
!
!$omp parallel do schedule(static) private(a)
      do a = 1,wf%n_v
!
         t_abc(a,a,a) = zero
!
      enddo
!$omp end parallel do
!
!$omp parallel do schedule(static) private(c, b, a, epsilon_c, epsilon_cb)
      do c = 1,wf%n_v
!
         epsilon_c = epsilon_ijk - wf%fock_diagonal(wf%n_o + c, 1)
!
         do b = 1,wf%n_v
!
            epsilon_cb = epsilon_c - wf%fock_diagonal(wf%n_o + b, 1)
!
            do a = 1,wf%n_v
!
               t_abc(a,b,c) = t_abc(a,b,c)*one/(epsilon_cb - wf%fock_diagonal(wf%n_o + a, 1))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine jacobian_cc3_t3_calc_cc3
!
!
   module subroutine jacobian_cc3_rho1_cc3(wf, i, j, k, c_abc, u_abc, rho1,               &
                                          L_jbic, L_kbic, L_kbjc, L_ibjc, L_ibkc, L_jbkc)
!!
!!    Calculate the triples contribution to rho1 for fixed i,j and k
!!
!!    rho_1 =+ P^abc_ijk sum_bc (C^abc_ijk - C^cba_ijk) L_jbkc
!!
!!    Alexander Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: rho1
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!
!     Construct u_abc = c_abc - c_cba
!
      call construct_123_minus_321(c_abc, u_abc, wf%n_v)
!
!     rho_ai += sum_bc (c^abc - c^cba)*L_jbkc
!
      call dgemv('N',         &
                 wf%n_v,      &
                 wf%n_v**2,   &
                 one,         &
                 u_abc,       & ! u^abc
                 wf%n_v,      &
                 L_jbkc,      & ! L_bcjk
                 1,           &
                 one,         &
                 rho1(:,i),   & ! rho_ai
                 1)
!
!     rho_ak += sum_cb (t^cba - t^abc)*L_jbic
!     Same as above if i == k, but this is never true
!
      call dgemv('N',         &
                 wf%n_v,      &
                 wf%n_v**2,   &
                 -one,        &
                 u_abc,       & ! u^cba = - u^abc
                 wf%n_v,      &
                 L_jbic,      & ! L_bcji
                 1,           &
                 one,         & ! rho_ak
                 rho1(:,k),   &
                 1)
!
!
      if (j .ne. k) then
!
!        Construct u_abc = c_acb - c_cab
!
         call construct_132_minus_312(c_abc, u_abc, wf%n_v)
!
!        omega_ai += sum_cb (c^acb - c^cab)*L_kbjc
!
         call dgemv('N',         &
                    wf%n_v,      &
                    wf%n_v**2,   &
                    one,         &
                    u_abc,       & ! u^acb
                    wf%n_v,      &
                    L_kbjc,      & ! L_bckj
                    1,           &
                    one,         &
                    rho1(:,i),   & ! rho_ai
                    1)
!
!
         if (i .ne. j) then
!
!           rho_aj += sum_cb (c^cab - c^acb)*L_kbic
!
            call dgemv('N',         &
                       wf%n_v,      &
                       wf%n_v**2,   &
                       -one,        &
                       u_abc,       & ! u^cab = -u^acb
                       wf%n_v,      &
                       L_kbic,      & ! L_bcki
                       1,           &
                       one,         &
                       rho1(:,j),   & ! rho_aj
                       1)
!
         end if
!
      end if
!
      if (i .ne. j .and. j .ne. k) then
!
!        Construct u_abc = c_bac - c_bca
!        This is zero if j == k
!
         call construct_213_minus_231(c_abc, u_abc, wf%n_v)
!
!        rho_aj += sum_cb (c^bac - c^bca)*L_ibkc
!
         call dgemv('N',         &
                    wf%n_v,      &
                    wf%n_v**2,   &
                    one,         &
                    u_abc,       & ! u^bac
                    wf%n_v,      &
                    L_ibkc,      & ! L_bcik
                    1,           &
                    one,         &
                    rho1(:,j),   & ! rho_aj
                    1)
!
!        rho_ak += sum_cb (c^bca - c^bac)*L_ibjc
!
         call dgemv('N',         &
                    wf%n_v,      &
                    wf%n_v**2,   &
                    -one,        &
                    u_abc,       & ! u^bca = - u^bac
                    wf%n_v,      &
                    L_ibjc,      & ! L_bcij
                    1,           &
                    one,         &
                    rho1(:,k),   & ! rho_ak
                    1)
!
      end if
!
   end subroutine jacobian_cc3_rho1_cc3
!
!
   module subroutine jacobian_cc3_rho2_cc3(wf, i, j, k, t_abc, u_abc, v_abc, rho2,        &
                                          F_kc, g_dbic, g_dbjc, g_dbkc,                   &
                                          g_jlic, g_klic, g_kljc, g_iljc, g_ilkc, g_jlkc)
!!
!!    Calculate the triples contribution to rho2 for fixed i,j and k
!!
!!    C^abc_ijk contributions:
!!    rho2 =+ P^abc_ijk (2 C^abc_ijk - C^acb_ijk - C^cba_ijk) (F_kc + (db|kc) - (jl|kc))
!!
!!    t^abc_ijk contributions:
!!    rho2 =+ P^abc_ijk (2 t^abc_ijk - t^acb_ijk - t^cba_ijk) (F'_kc + g'_dbkc - g'_jlkc)
!!
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_kc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbkc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_klic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jlkc
!
      real(dp) :: alpha
!
      if (i .ne. j) then
         alpha = one
      else
         alpha = half
      end if
!
!     construct u_abc = 2*t_abc - t_acb - t_cba
!
      call construct_123_min_132_min_321(t_abc, u_abc, wf%n_v)
!
!     rho_adij += sum_bc (2*t_abc - t_acb - t_cba)*g_dbkc
!
      call dgemm('N','N',        &
                 wf%n_v,         &
                 wf%n_v,         &
                 wf%n_v**2,      &
                 alpha,          &
                 u_abc,          & ! u^abc
                 wf%n_v,         &
                 g_dbkc,         & ! g_bc_dk
                 wf%n_v**2,      &
                 one,            &
                 rho2(:,:,i,j),  & ! rho^ad_ij
                 wf%n_v)
!
      if (i .ne. j) then
!
!        rho_ablj += sum_c (2*t_abc - t_acb - t_cba)*g_ilkc
!
         call dgemm('N','N',        &
                    wf%n_v**2,      &
                    wf%n_o,         &
                    wf%n_v,         &
                    -one,           &
                    u_abc,          & ! u^abc
                    wf%n_v**2,      &
                    g_ilkc,         & ! g_c_ljk
                    wf%n_v,         &
                    one,            &
                    rho2(:,:,:,j),  & ! rho^ab_lj
                    wf%n_v**2)
!
!        rho_abij += sum_c (2*t_abc - t_acb - t_cba)*F_kc
!
         call dgemv('N',            &
                    wf%n_v**2,      &
                    wf%n_v,         &
                    one,            &
                    u_abc,          & ! u^abc
                    wf%n_v**2,      &
                    F_kc(:,k),      & ! F_ck
                    1,              &
                    one,            &
                    rho2(:,:,i,j),  & ! rho^ab_lj
                    wf%n_v**2)
!
      end if
!
!     resort to u_abc = 2*t_bac - t_cab - t_bca
!
      call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     rho_adji += sum_bc (2*t_bac - t_cab - t_bca)*g_dbkc
!
      call dgemm('N','N',        &
                 wf%n_v,         &
                 wf%n_v,         &
                 wf%n_v**2,      &
                 alpha,          &
                 v_abc,          & ! v^abc
                 wf%n_v,         &
                 g_dbkc,         & ! g_bc_dk
                 wf%n_v**2,      &
                 one,            &
                 rho2(:,:,j,i),  & ! rho^ad_ji
                 wf%n_v)
!
!     rho_abli += sum_c (2*t_bac - t_cab - t_bca)*g_jlkc
!
      call dgemm('N','N',        &
                 wf%n_v**2,      &
                 wf%n_o,         &
                 wf%n_v,         &
                 -one,           &
                 v_abc,          & ! v^bac
                 wf%n_v**2,      &
                 g_jlkc,         & ! g_c_ljk
                 wf%n_v,         &
                 one,            &
                 rho2(:,:,:,i),  & ! rho^ab_li
                 wf%n_v**2)
!
!        rho_abij += sum_c (2*t_bac - t_cab - t_bca)*F_kc
!
         call dgemv('N',            &
                    wf%n_v**2,      &
                    wf%n_v,         &
                    one,            &
                    u_abc,          & ! u^bac
                    wf%n_v**2,      &
                    F_kc(:,k),      & ! F_ck
                    1,              &
                    one,            &
                    rho2(:,:,j,i),  & ! rho^ba_ji
                    wf%n_v**2)
!
!
      if (j .ne. k) then
!
!        construct u_abc = 2*t_acb - t_abc - t_cab
!
         call construct_132_min_123_min_312(t_abc, u_abc, wf%n_v)
!
         if (i .ne. k) then
            alpha = one
         else
            alpha = half
         end if
!
!        rho_adik += sum_bc (2*t_acb - t_abc - t_cab)*g_dbjc
!
         call dgemm('N','N',        &
                    wf%n_v,         &
                    wf%n_v,         &
                    wf%n_v**2,      &
                    alpha,          &
                    u_abc,          & ! u^acb
                    wf%n_v,         &
                    g_dbjc,         & ! g_bc_dj
                    wf%n_v**2,      &
                    one,            &
                    rho2(:,:,i,k),  & ! rho^ad_ik
                    wf%n_v)
!
         if (i .ne. k .and. i .ne. j) then
!
!           rho_ablk += sum_c (2*t_acb - t_abc - t_cab)*g_iljc
!
            call dgemm('N','N',        &
                       wf%n_v**2,      &
                       wf%n_o,         &
                       wf%n_v,         &
                       -one,           &
                       u_abc,          & ! u^acb
                       wf%n_v**2,      &
                       g_iljc,         & ! g_c_lij
                       wf%n_v,         &
                       one,            &
                       rho2(:,:,:,k),  &
                       wf%n_v**2)
!
!        rho_abik += sum_c (2*t_acb - t_abc - t_cab)*F_jc
!
         call dgemv('N',            &
                    wf%n_v**2,      &
                    wf%n_v,         &
                    one,            &
                    u_abc,          & ! u^acb
                    wf%n_v**2,      &
                    F_kc(:,j),      & ! F_cj
                    1,              &
                    one,            &
                    rho2(:,:,i,k),  & ! rho^ab_ik
                    wf%n_v**2)
!
         end if
!
!        resort to u_abc = 2*t_cab - t_bac - t_cba
!
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        rho_adjk += sum_bc (2*t_cab - t_bac - t_cba)*g_dbic
!
         call dgemm('N','N',        &
                    wf%n_v,         &
                    wf%n_v,         &
                    wf%n_v**2,      &
                    alpha,          &
                    v_abc,          & ! v^cab
                    wf%n_v,         &
                    g_dbic,         & ! g_bc_di
                    wf%n_v**2,      &
                    one,            &
                    rho2(:,:,j,k),  & ! rho^ad_jk
                    wf%n_v)
!
!        rho_ablk += sum_c (2*t_cab - t_bac - t_cba)*g_jlic
!
         call dgemm('N','N',        &
                    wf%n_v**2,      &
                    wf%n_o,         &
                    wf%n_v,         &
                    -one,           &
                    v_abc,          & ! v^cab
                    wf%n_v**2,      &
                    g_jlic,         & ! g_c_lji
                    wf%n_v,         &
                    one,            &
                    rho2(:,:,:,i),  & ! rho^ab_lk
                    wf%n_v**2)
!
!        rho_abjk += sum_c (2*t_cab - t_bac - t_cba)*F_ic
!
         call dgemv('N',            &
                    wf%n_v**2,      &
                    wf%n_v,         &
                    one,            &
                    u_abc,          & ! u^cab
                    wf%n_v**2,      &
                    F_kc(:,i),      & ! F_ci
                    1,              &
                    one,            &
                    rho2(:,:,j,k),  & ! rho^ab_jk
                    wf%n_v**2)
!
      end if
!
!
      if(i .ne. k .and. (j .ne. k .or. i .ne. j)) then
!
!        construct u_abc = 2*t_cba - t_bca - t_abc
!
         call construct_321_min_231_min_123(t_abc, u_abc, wf%n_v)
!
         if (i .ne. j) then
!
            if (k .ne. j) then
               alpha = one
            else
               alpha = half
            end if
!
!           rho_adkj += sum_bc (2*t_cba - t_bca - t_abc)*g_dbic
!
            call dgemm('N','N',        &
                       wf%n_v,         &
                       wf%n_v,         &
                       wf%n_v**2,      &
                       alpha,          &
                       u_abc,          & ! u^cba
                       wf%n_v,         &
                       g_dbic,         & ! g_cb_di
                       wf%n_v**2,      &
                       one,            &
                       rho2(:,:,k,j),  & ! rho^ad_kj
                       wf%n_v)
!
!           omega_ablj += sum_c (2*t_cba - t_bca - t_abc)*g_klic
!
            call dgemm('N','N',        &
                       wf%n_v**2,      &
                       wf%n_o,         &
                       wf%n_v,         &
                       -one,           &
                       u_abc,          & ! u^cba
                       wf%n_v**2,      &
                       g_klic,         & ! g_c_lki
                       wf%n_v,         &
                       one,            &
                       rho2(:,:,:,j),  & ! rho^ab_lj
                       wf%n_v**2)
!
!        rho_abjk += sum_c (2*t_cab - t_bac - t_cba)*F_ic
!
         call dgemv('N',            &
                    wf%n_v**2,      &
                    wf%n_v,         &
                    one,            &
                    u_abc,          & ! u^cab
                    wf%n_v**2,      &
                    F_kc(:,i),      & ! F_ci
                    1,              &
                    one,            &
                    rho2(:,:,k,j),  & ! rho^ba_kj
                    wf%n_v**2)
!
         end if
!
!        resort to u_abc = 2*t_bca - t_cba - t_bac
!
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
         if (i .ne. j) then
!
!           rho_adki += sum_bc (2*t_bca - t_cba - t_bac)*g_dbjc
!
            call dgemm('N','N',        &
                       wf%n_v,         &
                       wf%n_v,         &
                       wf%n_v**2,      &
                       alpha,          &
                       v_abc,          & ! v^bca
                       wf%n_v,         &
                       g_dbjc,         & ! g_bc_dj
                       wf%n_v**2,      &
                       one,            &
                       rho2(:,:,k,i),  & ! rho^ad_ki
                       wf%n_v)
!
         end if
!
         if (j .ne. k) then
!
!           rho_abli += sum_c (2*t_bca - t_cba - t_bac)*g_kljc
!
            call dgemm('N','N',        &
                       wf%n_v**2,      &
                       wf%n_o,         &
                       wf%n_v,         &
                       -one,           &
                       v_abc,          & ! v^bca
                       wf%n_v**2,      &
                       g_kljc,         & ! g_c_lkj
                       wf%n_v,         &
                       one,            &
                       rho2(:,:,:,i),  & ! rho^ab_li
                       wf%n_v**2)
!
!        rho_abik += sum_c (2*t_bca - t_cba - t_bac)*F_jc
!
         call dgemv('N',            &
                    wf%n_v**2,      &
                    wf%n_v,         &
                    one,            &
                    u_abc,          & ! u^bca
                    wf%n_v**2,      &
                    F_kc(:,j),      & ! F_cj
                    1,              &
                    one,            &
                    rho2(:,:,k,i),  & ! rho^ba_ki
                    wf%n_v**2)
!
         end if
!
      end if
!
!
   end subroutine jacobian_cc3_rho2_cc3
!
!
end submodule jacobian_cc3
