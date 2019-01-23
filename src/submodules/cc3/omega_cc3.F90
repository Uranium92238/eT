submodule (cc3_class) omega_cc3
!
!!
!!    Omega submodule (cc3)
!!    Alex C. Paul and Rolf H. Myhre January 2019
!!
!!    Routines to construct 
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_omega_cc3(wf, omega)
!!
!!    Construct omega (CC3)
!!    Written by Alex C. Paul and Rolf H. Myhre 2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable     :: omega1
      real(dp), dimension(:), allocatable       :: omega2
      real(dp), dimension(:,:,:,:), allocatable :: omega_abij
!
      integer(i15)  :: i,j,a,b,a_end,aibj
!
      type(timings) :: cc3_timer
      type(timings) :: ccsd_timer
!
      call cc3_timer%init('CC3 contribution)')
      call ccsd_timer%init('CCSD contribution)')
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call mem%alloc(omega2, wf%n_t2)
!
!     Set the omega vector to zero
!
      omega1 = zero
      omega2 = zero
!
!     Construct CCSD singles contributions
!
      call ccsd_timer%start()
!
      call wf%omega_ccsd_a1(omega1)
      call wf%omega_ccsd_b1(omega1)
      call wf%omega_ccsd_c1(omega1)
!
      call wf%omega_ccs_a1(omega1)
!
!     Construct CCSD doubles contributions
!
      call wf%omega_ccsd_a2(omega2)
      call wf%omega_ccsd_b2(omega2)
      call wf%omega_ccsd_c2(omega2)
      call wf%omega_ccsd_d2(omega2)
      call wf%omega_ccsd_e2(omega2)
!
      call dcopy(wf%n_t2, omega2, 1, omega(wf%n_t1+1, 1), 1)
!
      call ccsd_timer%freeze()
      call ccsd_timer%switch_off()
!
      call mem%dealloc(omega2, wf%n_t2)
      call mem%alloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      omega_abij = zero
!
      call cc3_timer%start()
      call wf%omega_cc3_a(omega1,omega_abij)
      call cc3_timer%freeze()
      call cc3_timer%switch_off()
!
      call dcopy(wf%n_t1, omega1, 1, omega, 1)
!
      aibj = 0
!
      do j = 1,wf%n_o
         do b = 1,wf%n_v
            do i = 1,j
!
               if(i .ne. j) then
                  a_end = wf%n_v
               else
                  a_end = b
               end if
!
               do a = 1,a_end
!
                  aibj = aibj + 1
!
                  omega(wf%n_t1+aibj,1) = omega(wf%n_t1+aibj,1) + omega_abij(a,b,i,j) + omega_abij(b,a,j,i)
!
               end do
            end do
         end do
      end do
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
      call mem%dealloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2)
!!
!!    CC3 Omega terms
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc 
      real(dp), dimension(:,:,:), allocatable :: u_abc 
      real(dp), dimension(:,:,:), allocatable :: v_abc 
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abji
!
      real(dp), dimension(:,:), allocatable :: F_kc !Transpose the fock matrix sub block
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
      integer(i15) :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer(i15) :: i_batch, j_batch, k_batch
      integer(i15) :: req_0, req_1, req_2, req_3
      real(dp)     :: batch_buff = 0.0
!
!     Set up required integrals on disk
      call wf%omega_cc3_integrals() 
!
      call mem%alloc(t_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%alloc(u_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%alloc(v_abc,wf%n_v,wf%n_v,wf%n_v)
!
      call mem%alloc(F_kc,wf%n_v,wf%n_o)
      call sort_12_to_21(wf%fock_ia,F_kc,wf%n_o,wf%n_v)
!
      call mem%alloc(t_abji,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
      call squareup_and_sort_1234_to_1342(wf%t2,t_abji,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
      t_abc = zero
      u_abc = zero
      v_abc = zero
!
      req_0 = 0
      req_1 = 2*wf%n_v**3
      req_2 = 2*wf%n_o*wf%n_v+wf%n_v**2
      req_3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, batch_buff)
!
!
!     Allocate integral arrays and assign pointers.
!     Without pointers we'll have to use three times as much 
!     memory for the non-batching case
!
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o) 
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o) 
!
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,wf%n_o,wf%n_o) 
!
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o) 
!
         call mem%alloc(L_jbic,wf%n_v,wf%n_v,wf%n_o,wf%n_o) 
!
      else !batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
!
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
!
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_klic,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_kljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_iljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_ilkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_jlkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
!
         call mem%alloc(L_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(L_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(L_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(L_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(L_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(L_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
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
      do i_batch = 1,batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%omega_cc3_vvv_reader(batch_i,g_bdci,g_dbic)
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
         do j_batch = 1,i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%omega_cc3_ov_vv_reader(batch_j,batch_i,g_ljci,g_jlic,L_jbic)
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_jbic_p => L_jbic
!
            if (j_batch .ne. i_batch) then
!
               call wf%omega_cc3_vvv_reader(batch_j,g_bdcj,g_dbjc)
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
               call wf%omega_cc3_ov_vv_reader(batch_i,batch_j,g_licj,g_iljc,L_ibjc)
               g_licj_p => g_licj
               g_iljc_p => g_iljc
               L_ibjc_p => L_ibjc
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
            endif
!
            do k_batch = 1,j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. i_batch .and. k_batch .ne. j_batch) then
!
                  call wf%omega_cc3_vvv_reader(batch_k,g_bdck,g_dbkc)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
!
                  call wf%omega_cc3_ov_vv_reader(batch_k,batch_i,g_lkci,g_klic,L_kbic)
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
                  L_kbic_p => L_kbic
!
                  call wf%omega_cc3_ov_vv_reader(batch_i,batch_k,g_lick,g_ilkc,L_ibkc)
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  call wf%omega_cc3_ov_vv_reader(batch_k,batch_j,g_lkcj,g_kljc,L_kbjc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  call wf%omega_cc3_ov_vv_reader(batch_j,batch_k,g_ljck,g_jlkc,L_jbkc)
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
!
               else if (k_batch .eq. i_batch) then
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
!
                  if (j_batch .eq. i_batch) then
!
                     g_lkci_p => g_ljci
                     g_klic_p => g_jlic
                     L_kbic_p => L_jbic
!
                     g_lick_p => g_ljci
                     g_ilkc_p => g_jlic
                     L_ibkc_p => L_jbic
!
                     g_lkcj_p => g_ljci
                     g_kljc_p => g_jlic
                     L_kbjc_p => L_jbic
!
                     g_ljck_p => g_ljci
                     g_jlkc_p => g_jlic
                     L_jbkc_p => L_jbic
!
                  else
!
                     call wf%omega_cc3_ov_vv_reader(batch_k,batch_i,g_lkci,g_klic,L_kbic)
                     g_lkci_p => g_lkci
                     g_klic_p => g_klic
                     L_kbic_p => L_kbic
!
                     g_lick_p => g_lkci
                     g_ilkc_p => g_klic
                     L_ibkc_p => L_kbic
!
                     g_lkcj_p => g_licj
                     g_kljc_p => g_iljc
                     L_kbjc_p => L_ibjc
!
                     g_ljck_p => g_ljci
                     g_jlkc_p => g_jlic
                     L_jbkc_p => L_jbic
!
                  endif
!
               else if (k_batch .eq. j_batch) then
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
                  L_kbic_p => L_jbic
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
!
                  call wf%omega_cc3_ov_vv_reader(batch_k,batch_j,g_lkcj,g_kljc,L_kbjc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_kbjc
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
                        call wf%omega_cc3_W_calc(i, j, k, t_abc, u_abc, t_abji, &
                                                 g_bdci_p(:,:,:,i_rel), &
                                                 g_bdcj_p(:,:,:,j_rel), &
                                                 g_bdck_p(:,:,:,k_rel), &
                                                 g_ljci_p(:,:,j_rel,i_rel), &
                                                 g_lkci_p(:,:,k_rel,i_rel), &
                                                 g_lkcj_p(:,:,k_rel,j_rel), &
                                                 g_licj_p(:,:,i_rel,j_rel), &
                                                 g_lick_p(:,:,i_rel,k_rel), &
                                                 g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_eps(i, j, k, t_abc)
!
                        call wf%omega_cc3_omega1(i, j, k, t_abc, u_abc, omega1, omega2, F_kc, &
                                                 L_jbic_p(:,:,j_rel,i_rel), &
                                                 L_kbic_p(:,:,k_rel,i_rel), &
                                                 L_kbjc_p(:,:,k_rel,j_rel), &
                                                 L_ibjc_p(:,:,i_rel,j_rel), &
                                                 L_ibkc_p(:,:,i_rel,k_rel), &
                                                 L_jbkc_p(:,:,j_rel,k_rel))
!
!
                        call wf%omega_cc3_omega2(i, j, k, t_abc, u_abc, v_abc, omega2, &
                                                 g_dbic_p(:,:,:,i_rel), &
                                                 g_dbjc_p(:,:,:,j_rel), &
                                                 g_dbkc_p(:,:,:,k_rel), &
                                                 g_jlic_p(:,:,j_rel,i_rel), &
                                                 g_klic_p(:,:,k_rel,i_rel), &
                                                 g_kljc_p(:,:,k_rel,j_rel), &
                                                 g_iljc_p(:,:,i_rel,j_rel), &
                                                 g_ilkc_p(:,:,i_rel,k_rel), &
                                                 g_jlkc_p(:,:,j_rel,k_rel))
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!     Close files
!
      call disk%close_file(wf%g_bdck_t)
      call disk%close_file(wf%g_ljck_t)
      call disk%close_file(wf%g_dbkc_t)
      call disk%close_file(wf%g_jlkc_t)
      call disk%close_file(wf%L_jbkc_t)
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o) 
         call mem%dealloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o) 
         call mem%dealloc(g_ljci,wf%n_o,wf%n_v,wf%n_o,wf%n_o) 
         call mem%dealloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o) 
         call mem%dealloc(L_jbic,wf%n_v,wf%n_v,wf%n_o,wf%n_o) 
!
      else
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%dealloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%dealloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
!
         call mem%dealloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%dealloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%dealloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
!
         call mem%dealloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
!
         call mem%dealloc(g_jlic,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%dealloc(g_klic,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%dealloc(g_kljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%dealloc(g_iljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%dealloc(g_ilkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%dealloc(g_jlkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
!
         call mem%dealloc(L_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(L_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(L_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(L_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(L_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%dealloc(L_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
!
      endif 
!
      call mem%dealloc(F_kc,wf%n_v,wf%n_o)
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%dealloc(u_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%dealloc(v_abc,wf%n_v,wf%n_v,wf%n_v)
!
      call mem%dealloc(t_abji,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Omega and store on disk
!!    (bd|ck) ordered as dbc,k
!!    (db|kc) ordered as bcd,k
!!    (lj|ck) ordered as lc,jk
!!    (jl|kc) ordered as cl,jk
!!    (jb|kc) stored as L_jbkc = 2g_jbkc - g_jckb ordered as bc,jk
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs !Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs !Array for sorted integrals
      real(dp), dimension(:,:), allocatable     :: v2_help !Help array for constructing L_jbkc
!
      integer(i15) :: k, j, record 
      type(batching_index) :: batch_k
!
      integer(i15) :: req_0, req_k
      integer(i15) :: current_k_batch
!
      integer(i15) :: ioerror=-1
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
   end subroutine omega_cc3_integrals_cc3
!
!
   module subroutine omega_cc3_vvv_reader_cc3(wf,batch_x,g_bdcx,g_dbxc)
!!
!!    Read the bdck and dbkc integrals needed in the current batch
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_bdcx
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_dbxc
!
      integer(i15) :: ioerror
      integer(i15) :: x, x_abs
!
      do x = 1,batch_x%length
!
         x_abs = batch_x%first + x - 1
!
         read(wf%g_bdck_t%unit,rec=x_abs, iostat=ioerror) g_bdcx(:,:,:,x)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a)') 'Failed to read bdck file'
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
            write(output%unit,'(t3,a)') 'Failed to read dbkc file'
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!
   end subroutine omega_cc3_vvv_reader_cc3
!
!
   module subroutine omega_cc3_ov_vv_reader_cc3(wf,batch_y,batch_x,g_lycx,g_ylxc,L_ybxc)
!!
!!    Read the ljck, jlkc and jbkc integrals needed in the current batches
!!
!!    Rolf H. Myhre, January 2019
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
!
      integer(i15) :: ioerror, record
      integer(i15) :: x, y, x_abs, y_abs
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
               write(output%unit,'(t3,a)') 'Failed to read ljck file'
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
!
   end subroutine omega_cc3_ov_vv_reader_cc3
!
!
   module subroutine omega_cc3_W_calc_cc3(wf, i, j, k, t_abc, u_abc, t_abji, &
                                          g_bdci, g_bdcj, g_bdck, &
                                          g_ljci, g_lkci, g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Calculate the the contributions to the t_3 amplitudes 
!!    for occupied indices i,j,k
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abji
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci 
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj 
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck 
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljci 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkci 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkcj 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_licj 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lick 
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljck 
!
!
!     Contributions to W
!     W^abc_ijk = P^abc_ijk(\sum_d t^ad_ij(bd|ck) - \sum_l t^ab_il(lj|ck))
!
!
!     t^ad_ij*(bd|ck)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 t_abji(:,:,j,i), &
                 wf%n_v, &
                 g_bdck, &
                 wf%n_v, &
                 zero, &
                 t_abc, &
                 wf%n_v)
!
!     t^ab_il*(lj|ck)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 t_abji(:,:,:,i), &
                 wf%n_v**2, &
                 g_ljck, &
                 wf%n_o, &
                 one, &
                 t_abc, &
                 wf%n_v**2)
!
!
!
!     t^bd_ji*(ad|ck)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 t_abji(:,:,i,j), &
                 wf%n_v, &
                 g_bdck, &
                 wf%n_v, &
                 zero, &
                 u_abc, &
                 wf%n_v)
!
!     t^ba_jl*(li|ck)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 t_abji(:,:,:,j), &
                 wf%n_v**2, &
                 g_lick, &
                 wf%n_o, &
                 one, &
                 u_abc, &
                 wf%n_v**2)
!
      call sort_123_to_213_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
!     t^ad_ik*(cd|bj)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 t_abji(:,:,k,i), &
                 wf%n_v, &
                 g_bdcj, &
                 wf%n_v, &
                 zero, &
                 u_abc, &
                 wf%n_v)
!
!     t^ac_il*(lk|bj)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 t_abji(:,:,:,i), &
                 wf%n_v**2, &
                 g_lkcj, &
                 wf%n_o, &
                 one, &
                 u_abc, &
                 wf%n_v**2)
!
      call sort_123_to_132_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
!     t^cd_ki*(ad|bj)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 t_abji(:,:,i,k), &
                 wf%n_v, &
                 g_bdcj, &
                 wf%n_v, &
                 zero, &
                 u_abc, &
                 wf%n_v)
!
!     t^ca_kl*(li|bj)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 t_abji(:,:,:,k), &
                 wf%n_v**2, &
                 g_licj, &
                 wf%n_o, &
                 one, &
                 u_abc, &
                 wf%n_v**2)
!
      call sort_123_to_231_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
!     t^bd_jk*(cd|ai)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 t_abji(:,:,k,j), &
                 wf%n_v, &
                 g_bdci, &
                 wf%n_v, &
                 zero, &
                 u_abc, &
                 wf%n_v)
!
!     t^bc_jl*(lk|ai)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 t_abji(:,:,:,j), &
                 wf%n_v**2, &
                 g_lkci, &
                 wf%n_o, &
                 one, &
                 u_abc, &
                 wf%n_v**2)
!
      call sort_123_to_312_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
!     t^cd_kj*(bd|ai)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 t_abji(:,:,j,k), &
                 wf%n_v, &
                 g_bdci, &
                 wf%n_v, &
                 zero, &
                 u_abc, &
                 wf%n_v)
!
!     t^cb_kl*(lj|ai)
!     ---------------
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 t_abji(:,:,:,k), &
                 wf%n_v**2, &
                 g_ljci, &
                 wf%n_o, &
                 one, &
                 u_abc, &
                 wf%n_v**2)
!
      call sort_123_to_321_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
   end subroutine omega_cc3_W_calc_cc3
!
!
   module subroutine omega_cc3_eps_cc3(wf, i, j, k, t_abc)
!!
!!    Divide W^abc_ijk with -epsilon^abc_ijk to obtain T^abc_ijk
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v), intent(inout) :: t_abc
!
      integer(i15) a, b, c
!
      real(dp) :: epsilon_ijk, epsilon_c, epsilon_cb
!
!
      epsilon_ijk = wf%fock_diagonal(i,1) + wf%fock_diagonal(j,1) + wf%fock_diagonal(k,1)
!
!$omp parallel do schedule(static) private(a)
      do a = 1,wf%n_v
!
         t_abc(a,a,a) = zero
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
               t_abc(a,b,c) = t_abc(a,b,c)*one/(epsilon_cb - wf%fock_diagonal(wf%n_o + a, 1))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!
   end subroutine omega_cc3_eps_cc3
!
!
   module subroutine omega_cc3_omega1_cc3(wf, i, j, k, t_abc, u_abc, omega1, omega2, F_kc, &
                                          L_jbic, L_kbic, L_kbjc, L_ibjc, L_ibkc, L_jbkc)
!!
!!    Calculate the triples contribution to omega1 and
!!    the Fock contribution to omega2
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_kc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!
!
!     Construct u_abc = t_abc - t_cba
!
      call construct_123_minus_321(t_abc, u_abc, wf%n_v)
!
!     omega_ai += sum_bc (t^abc - t^cba)*L_jbkc
!
      call dgemv('N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 one, &
                 u_abc, &
                 wf%n_v, &
                 L_jbkc, &
                 1, &
                 one, &
                 omega1(:,i), &
                 1)
!
!     omega_abij += sum_c (t^abc - t^cba)*F_kc
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 u_abc, &
                 wf%n_v**2, &
                 F_kc(:,k), &
                 1, &
                 one, &
                 omega2(:,:,i,j), &
                 1)
!
!
!     omega_ak += sum_cb (t^cba - t^abc)*L_jbic
!     Same as above if i == k, but this is never true
!
      call dgemv('N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 -one, &
                 u_abc, &
                 wf%n_v, &
                 L_jbic, &
                 1, &
                 one, &
                 omega1(:,k), &
                 1)
!
!     omega_abkj += sum_c (t^cba - t^abc)*F_ic
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 -one, &
                 u_abc, &
                 wf%n_v**2, &
                 F_kc(:,i), &
                 1, &
                 one, &
                 omega2(:,:,k,j), &
                 1)
!
!
!
      if (j .ne. k) then
!
!        Construct u_abc = t_acb - t_cab
!
         call construct_132_minus_312(t_abc, u_abc, wf%n_v)
!
!        omega_ai += sum_cb (t^acb - t^cab)*L_kbjc
!
         call dgemv('N', &
                    wf%n_v, &
                    wf%n_v**2, &
                    one, &
                    u_abc, &
                    wf%n_v, &
                    L_kbjc, &
                    1, &
                    one, &
                    omega1(:,i), &
                    1)
!
!        omega_abji += sum_c (t^acb - t^cab)*F_jc
!
         call dgemv('N', &
                    wf%n_v**2, &
                    wf%n_v, &
                    one, &
                    u_abc, &
                    wf%n_v**2, &
                    F_kc(:,j), &
                    1, &
                    one, &
                    omega2(:,:,i,k), &
                    1)
!
!
         if (i .ne. j) then
!
!           omega_aj += sum_cb (t^cab - t^acb)*L_kbic
!
            call dgemv('N', &
                       wf%n_v, &
                       wf%n_v**2, &
                       -one, &
                       u_abc, &
                       wf%n_v, &
                       L_kbic, &
                       1, &
                       one, &
                       omega1(:,j), &
                       1)
!
!           omega_abjk += sum_c (t^cab - t^acb)*F_ic
!
            call dgemv('N', &
                       wf%n_v**2, &
                       wf%n_v, &
                       -one, &
                       u_abc, &
                       wf%n_v**2, &
                       F_kc(:,i), &
                       1, &
                       one, &
                       omega2(:,:,j,k), &
                       1)
!
!
         end if
!
      end if
!
!
      if (i .ne. j .and. j .ne. k) then
!
!        Construct u_abc = t_bac - t_bca
!        This is zero if j == k
!
         call construct_213_minus_231(t_abc, u_abc, wf%n_v)
!
!        omega_aj += sum_cb (t^bac - t^bca)*L_ibkc
!
         call dgemv('N', &
                    wf%n_v, &
                    wf%n_v**2, &
                    one, &
                    u_abc, &
                    wf%n_v, &
                    L_ibkc, &
                    1, &
                    one, &
                    omega1(:,j), &
                    1)
!
!        omega_abij += sum_c (t^bac - t^bca)*F_kc
!
         call dgemv('N', &
                    wf%n_v**2, &
                    wf%n_v, &
                    one, &
                    u_abc, &
                    wf%n_v**2, &
                    F_kc(:,k), &
                    1, &
                    one, &
                    omega2(:,:,j,i), &
                    1)
!
!
!        omega_ak += sum_cb (t^bca - t^bac)*L_ibjc
!
         call dgemv('N', &
                    wf%n_v, &
                    wf%n_v**2, &
                    -one, &
                    u_abc, &
                    wf%n_v, &
                    L_ibjc, &
                    1, &
                    one, &
                    omega1(:,k), &
                    1)
!
!        omega_abki += sum_c (t^bca - t^bac)*F_jc
!
         call dgemv('N', &
                    wf%n_v**2, &
                    wf%n_v, &
                    -one, &
                    u_abc, &
                    wf%n_v**2, &
                    F_kc(:,j), &
                    1, &
                    one, &
                    omega2(:,:,k,i), &
                    1)
!
!
      end if
!
   end subroutine omega_cc3_omega1_cc3
!
!
   module subroutine omega_cc3_omega2_cc3(wf, i, j, k, t_abc, u_abc, v_abc, omega2, &
                                          g_dbic, g_dbjc, g_dbkc, &
                                          g_jlic, g_klic, g_kljc, g_iljc, g_ilkc, g_jlkc)
!!
!!    Calculate the triples contribution to omega1 and
!!    the Fock contribution to omega2
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
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
!     omega_adij += sum_bc (2*t_abc - t_acb - t_cba)*g_dbkc
!
      call dgemm('N','N', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v**2, &
                 alpha, &
                 u_abc, &
                 wf%n_v, &
                 g_dbkc, &
                 wf%n_v**2, &
                 one, &
                 omega2(:,:,i,j), &
                 wf%n_v)
!
      if (i .ne. j) then
!
!        omega_ablj += \sum_c (2*t_abc - t_acb - t_cba)*g_ilkc
!
         call dgemm('N','N', &
                    wf%n_v**2, &
                    wf%n_o, &
                    wf%n_v, &
                    -one, &
                    u_abc, &
                    wf%n_v**2, &
                    g_ilkc, &
                    wf%n_v, &
                    one, &
                    omega2(:,:,:,j), &
                    wf%n_v**2)
!
      end if
!
!
!     resort to u_abc = 2*t_bac - t_bca - t_cab
!
      call sort_123_to_213(u_abc,v_abc,wf%n_v,wf%n_v,wf%n_v)
!
!     omega_adji += sum_bc (2*t_bac - t_bca - t_cab)*g_dbkc
!
      call dgemm('N','N', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v**2, &
                 alpha, &
                 v_abc, &
                 wf%n_v, &
                 g_dbkc, &
                 wf%n_v**2, &
                 one, &
                 omega2(:,:,j,i), &
                 wf%n_v)
!
!     omega_abli += \sum_c (2*t_bac - t_bca - t_cab)*g_jlkc
!
      call dgemm('N','N', &
                 wf%n_v**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 v_abc, &
                 wf%n_v**2, &
                 g_jlkc, &
                 wf%n_v, &
                 one, &
                 omega2(:,:,:,i), &
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
!        omega_adik += sum_bc (2*t_acb - t_abc - t_cab)*g_dbjc
!
         call dgemm('N','N', &
                    wf%n_v, &
                    wf%n_v, &
                    wf%n_v**2, &
                    alpha, &
                    u_abc, &
                    wf%n_v, &
                    g_dbjc, &
                    wf%n_v**2, &
                    one, &
                    omega2(:,:,i,k), &
                    wf%n_v)
!
         if (i .ne. k .and. i .ne. j) then
!
!           omega_ablk += \sum_c (2*t_acb - t_abc - t_cab)*g_iljc
!
            call dgemm('N','N', &
                       wf%n_v**2, &
                       wf%n_o, &
                       wf%n_v, &
                       -one, &
                       u_abc, &
                       wf%n_v**2, &
                       g_iljc, &
                       wf%n_v, &
                       one, &
                       omega2(:,:,:,k), &
                       wf%n_v**2)
!
         end if
!
!        resort to u_abc = 2*t_bca - t_bac - t_cba
!
         call sort_123_to_213(u_abc,v_abc,wf%n_v,wf%n_v,wf%n_v)
!
!        omega_adki += sum_bc (2*t_bca - t_bac - t_cba)*g_dbjc
!
         call dgemm('N','N', &
                    wf%n_v, &
                    wf%n_v, &
                    wf%n_v**2, &
                    alpha, &
                    v_abc, &
                    wf%n_v, &
                    g_dbjc, &
                    wf%n_v**2, &
                    one, &
                    omega2(:,:,k,i), &
                    wf%n_v)
!
!        omega_abli += \sum_c (2*t_bca - t_bac - t_cba)*g_kljc
!
         call dgemm('N','N', &
                    wf%n_v**2, &
                    wf%n_o, &
                    wf%n_v, &
                    -one, &
                    v_abc, &
                    wf%n_v**2, &
                    g_kljc, &
                    wf%n_v, &
                    one, &
                    omega2(:,:,:,i), &
                    wf%n_v**2)
!
      end if
!
!
      if(i .ne. k .and. (j .ne. k .or. i .ne. j)) then
!
!        construct u_abc = 2*t_acb - t_abc - t_cab
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
!           omega_adkj += sum_bc (2*t_acb - t_abc - t_cab)*g_dbic
!
            call dgemm('N','N', &
                       wf%n_v, &
                       wf%n_v, &
                       wf%n_v**2, &
                       alpha, &
                       u_abc, &
                       wf%n_v, &
                       g_dbic, &
                       wf%n_v**2, &
                       one, &
                       omega2(:,:,k,j), &
                       wf%n_v)
!
!           omega_ablj += \sum_c (2*t_acb - t_abc - t_cab)*g_kljic
!
            call dgemm('N','N', &
                       wf%n_v**2, &
                       wf%n_o, &
                       wf%n_v, &
                       -one, &
                       u_abc, &
                       wf%n_v**2, &
                       g_klic, &
                       wf%n_v, &
                       one, &
                       omega2(:,:,:,j), &
                       wf%n_v**2)
!
         end if
!
!        resort to u_abc = 2*t_cab - t_acb - t_bac
!
         call sort_123_to_213(u_abc,v_abc,wf%n_v,wf%n_v,wf%n_v)
!
         if (i .ne. j) then
!
!           omega_adjk += sum_bc (2*t_cab - t_acb - t_bac)*g_dbic
!
            call dgemm('N','N', &
                       wf%n_v, &
                       wf%n_v, &
                       wf%n_v**2, &
                       alpha, &
                       v_abc, &
                       wf%n_v, &
                       g_dbic, &
                       wf%n_v**2, &
                       one, &
                       omega2(:,:,j,k), &
                       wf%n_v)
!
         end if
!
         if (j .ne. k) then
!
!           omega_ablk += sum_c (2*t_cab - t_acb - t_bac)*g_jlic
!
            call dgemm('N','N', &
                       wf%n_v**2, &
                       wf%n_o, &
                       wf%n_v, &
                       -one, &
                       v_abc, &
                       wf%n_v**2, &
                       g_jlic, &
                       wf%n_v, &
                       one, &
                       omega2(:,:,:,k), &
                       wf%n_v**2)
!
         end if
!
      end if
!
!
   end subroutine omega_cc3_omega2_cc3
!
!
end submodule
