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
   module subroutine jacobian_transform_trial_vector_cc3(wf, c_i)
!!
!!    Jacobian transform trial vector
!!    Written by Rolf H. Myhre and Alexander Paul
!!
      class(cc3), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c_i
!
      call wf%jacobian_cc3_transformation(c_i)
!
   end subroutine jacobian_transform_trial_vector_cc3
!
!
   module subroutine jacobian_cc3_transformation_cc3(wf, c)
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
      real(dp), dimension(wf%n_amplitudes, 1) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_abij
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
      logical, private :: eri_c1_mem = .false.
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
   module subroutine jacobian_cc3_A_cc3(wf, rho, rho_ai, rho_aibj)
!!
!!    CC3 terms depending on the C^{abc}_{ijk} amplitudes
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
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
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
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
      call mem%alloc(F_ia_c1, wf%n_o, wf%n_v)
      call wf%construct_fock_ia_c1_cc3(wf, c_ai, F_ia_c1)
!
      call mem%alloc(t_abji,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
      call squareup_and_sort_1234_to_1342(wf%t2,t_abji,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
      t_abc = zero
      u_abc = zero
      v_abc = zero
!
      req_0 = 0
      req_1 = 0
      req_2 = 0
      req_3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, batch_buff)
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

!!    NB: the indices d and l are part of sigma_2 while j and l are summation indices in the following terms
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
      real(dp), dimension(:,:), allocatable :: L_db_J_c1, L_lj_J_c1
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
      call wf%g_bdck_c%init('g_bdck_c','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_bdck_c,'write')
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
                     integrals%n_J,             &
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
                     integrals%n_J,             &
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
                     integrals%n_J,             &
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
            write(wf%g_bdck_c%unit, rec=record, iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write bdck_c file')
         endif
!
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_bdck_c,'keep')
!
!
!     (d'b|kc) same batching (req0 and req_k)
!
!
      call wf%g_dbkc_t%init('g_dbkc_c','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_dbkc_c,'write')
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
                     integrals%n_J,             &
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
            write(wf%g_dbkc_c%unit, rec=record, iostat=ioerror) h_pqrs(:,:,:,k)
!
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write dbkc_c file')
         endif
!
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_dbkc_c,'keep')
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
      call wf%g_ljck_c%init('g_ljck_c','direct','unformatted',dp*wf%n_v*wf%n_o)
      call disk%open_file(wf%g_ljck_c,'write')
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
                     integrals%n_J,             &
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
                     integrals%n_J,             &
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
                     integrals%n_J,             &
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
               write(wf%g_ljck_c%unit,rec=record,iostat=ioerror) h_pqrs(:,:,j,k)
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
      call disk%close_file(wf%g_ljck_c,'keep')
!
!
!     (jl'|kc) same batching (req0 and req_k)
!
!
      call wf%g_jlkc_c%init('g_jlkc_c','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_jlkc_c,'write')
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
                     integrals%n_J,             &
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
            write(wf%g_jlkc_c%unit, rec=record, iostat=ioerror) h_pqrs(:,:,:,k)
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
      call disk%close_file(wf%g_jlkc_c, 'keep')
!
!
   end subroutine jacobian_cc3_c1_integrals_cc3
!
!
   subroutine construct_fock_ia_c1_cc3(wf, c_ai, F_ia_c1)
!!
!!    Calculates C1 transformed elements of the Fock matrix required for the CC3 jacobian
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
      real(dp), dimension(wf%n_o, wf%n_v), intent(out) :: F_ia_c1
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
                  integrals%n_J,       &
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
                  integrals%n_J,       &
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
!$omp parallel do private(a,i,j)
      do a = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               F_ia_c1(i,a) = two*g_iajj(i,a,j,j) - g_ijja(i,j,j,a)
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
   end subroutine construct_fock_ia_c1_cc3
!
!
end submodule jacobian_cc3
