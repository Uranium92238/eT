submodule (cc3_class) omega_cc3
!
!!
!!    Omega submodule (cc3)
!!    Alex C. Paul and Rolf H. Myhre 2018
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
      real(dp), dimension(:,:), allocatable :: omega1
      real(dp), dimension(:), allocatable   :: omega2
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
      call wf%omega_cc3_a(omega1,omega2)
!
      call dcopy(wf%n_t1, omega1, 1, omega, 1)
      call dcopy(wf%n_t2, omega2, 1, omega(wf%n_t1+1, 1), 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
      call mem%dealloc(omega2, wf%n_t2)
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
      real(dp), dimension(wf%n_t2), intent(inout) :: omega2
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc 
      real(dp), dimension(:,:,:), allocatable :: u_abc 
      real(dp), dimension(:,:,:), allocatable :: v_abc 
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), pointer              :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), pointer              :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), pointer              :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc
      real(dp), dimension(:,:,:,:), pointer              :: g_jlic_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_klic_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_kljc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_iljc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_ilkc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_jlkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), pointer              :: g_jbic_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_kbic_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_kbjc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), pointer              :: g_jbkc_p => null()
!
      integer(i15) :: i, j, k
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
      if (batch_i%num_batches .eq. 1) then
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o) 
         g_bdci_p => g_bdci
         g_bdcj_p => g_bdci
         g_bdck_p => g_bdci
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o) 
         g_dbic_p => g_dbic
         g_dbjc_p => g_dbjc
         g_dbkc_p => g_dbkc
!
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,wf%n_o,wf%n_o) 
         g_ljci_p => g_ljci
         g_lkci_p => g_ljci
         g_lkcj_p => g_ljci
         g_licj_p => g_ljci
         g_lick_p => g_ljci
         g_ljck_p => g_ljci
!
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o) 
         g_jlic_p => g_jlic
         g_klic_p => g_jlic
         g_kljc_p => g_jlic
         g_iljc_p => g_jlic
         g_ilkc_p => g_jlic
         g_jlkc_p => g_jlic
!
         call mem%alloc(g_jbic,wf%n_v,wf%n_v,wf%n_o,wf%n_o) 
         g_jbic_p => g_jbic
         g_kbic_p => g_jbic
         g_kbjc_p => g_jbic
         g_ibjc_p => g_jbic
         g_ibkc_p => g_jbic
         g_jbkc_p => g_jbic
!
      else
         call batch_i%determine_limits(1)
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         g_bdci_p => g_bdci
         g_bdcj_p => g_bdcj
         g_bdck_p => g_bdck
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         call mem%alloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%length) 
         g_dbic_p => g_bdci
         g_dbjc_p => g_bdcj
         g_dbkc_p => g_bdck
!
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length) 
         g_ljci_p => g_ljci
         g_lkci_p => g_lkci
         g_lkcj_p => g_lkcj
         g_licj_p => g_licj
         g_lick_p => g_lick
         g_ljck_p => g_ljck
!
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_klic,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_kljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_iljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_ilkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         call mem%alloc(g_jlkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length) 
         g_jlic_p => g_jlic
         g_klic_p => g_klic
         g_kljc_p => g_kljc
         g_iljc_p => g_iljc
         g_ilkc_p => g_ilkc
         g_jlkc_p => g_jlkc
!
         call mem%alloc(g_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         call mem%alloc(g_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length) 
         g_jbic_p => g_jbic
         g_kbic_p => g_kbic
         g_kbjc_p => g_kbjc
         g_ibjc_p => g_ibjc
         g_ibkc_p => g_ibkc
         g_jbkc_p => g_jbkc
!
      endif 
!
      g_bdci(1,1,1,1) = 9.3
      write(output%unit,*) "g_bdci_p =", g_bdci_p(1,1,1,1)
      write(output%unit,*) "g_bdcj_p =", g_bdcj_p(1,1,1,1)
!
      do i_batch = 1,batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         do j_batch = 1,batch_j%num_batches
!
            call batch_j%determine_limits(j_batch)
!
            do k_batch = 1,batch_k%num_batches
!
               call batch_k%determine_limits(k_batch)
!
               do i = 1,batch_i%length
                  do j = 1,batch_j%length
                     do k = 1,batch_k%length
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%dealloc(u_abc,wf%n_v,wf%n_v,wf%n_v)
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
!!    (jb|kc) ordered as bc,jk
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs !Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs !Array for sorted integrals
!
      integer(i15) :: k, j, record 
      type(batching_index) :: batch_k
!
      integer(i15) :: req_0, req_k
      integer(i15) :: current_k_batch
!
      integer(i15) :: ioerror=-1
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
               record  = (batch_k%first + k - 1)*wf%n_o + j
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
               record  = (batch_k%first + k - 1)*wf%n_o + j
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
      call wf%g_jbkc_t%init('g_jbkc_t','direct','unformatted',dp*wf%n_v**2)
      call disk%open_file(wf%g_jbkc_t,'write')
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
               record  = (batch_k%first + k - 1)*wf%n_o + j
               write(wf%g_jbkc_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,j,k)
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
      call disk%close_file(wf%g_jbkc_t,'keep')
!
!
   end subroutine omega_cc3_integrals_cc3
!
!
end submodule
