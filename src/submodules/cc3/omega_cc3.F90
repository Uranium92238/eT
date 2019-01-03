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
      real(dp), dimension(:,:), allocatable :: omega1, omega2
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call mem%alloc(omega2, wf%n_t2, 1)
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
      write(output%unit,*) "In construct_cc3_omega"
!
      call wf%omega_cc3_a(omega1,omega2)
!
      call dcopy(wf%n_t1, omega1, 1, omega, 1)
      call dcopy(wf%n_t2, omega2, 1, omega(wf%n_t1+1, 1), 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
      call mem%dealloc(omega2, wf%n_t2, 1)
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
      real(dp), dimension(wf%n_t2, 1), intent(inout) :: omega2
!
!
!     Set up required integrals on disk
      call wf%omega_cc3_integrals() 
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
      real(dp), dimension(:,:), allocatable :: g_pqrs !Array for constructed integrals
      real(dp), dimension(:,:), allocatable :: h_pqrs !Array for sorted integrals
!
      type(batching_index) :: batch_k
!
      integer(i15) :: req_0, req_k
      integer(i15) :: current_k_batch
!
      write(output%unit,*) "In omega_cc3_integrals"
!
      req_0 = wf%integrals%n_J*wf%n_v**2
      req_k = 2*wf%n_v**3 + wf%integrals%n_J*wf%n_v
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
!
!     (bd|ck)
      do current_k_batch = 1,batch_k%num_batches
!
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(g_pqrs,wf%n_v**2,wf%n_v*batch_k%length)
         call mem%alloc(h_pqrs,wf%n_v**2,wf%n_v*batch_k%length)
!
         call wf%get_vvvo(g_pqrs, &
                           1,wf%n_v, &
                           1,wf%n_v, &
                           1,wf%n_v, &
                           batch_k%first,batch_k%last)
!
         call sort_1234_to_2134(g_pqrs,h_pqrs,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
!
      enddo
!
   end subroutine omega_cc3_integrals_cc3
!
!
end submodule
