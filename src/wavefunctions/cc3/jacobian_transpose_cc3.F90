!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
submodule (cc3_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (cc3)
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * c_i,
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
   module subroutine prepare_for_jacobian_transpose_cc3(wf)
!!
!!    Prepare for Jacobian transpose (CC3)
!!    Write some integrals and intermediates to disk
!!    Written by Rolf H. Myhre and Alexander Paul, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      call wf%prepare_cc3_jacobian_transpose_integrals
      call wf%prepare_cc3_jacobian_transpose_intermediates
!
   end subroutine prepare_for_jacobian_transpose_cc3
!
!
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, c)
!!
!!    Effective Jacobian transpose transformation (CC3)
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
!!    Directs the transformation by the transpose of the  CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >,
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine. On exit, the vector c is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj, sigma_abij
!
      integer :: i, j, a, b, ai, bj, aibj, b_end ! Index
!
      type(timings) :: cc3_timer
      type(timings) :: ccsd_timer
!
      call cc3_timer%init('CC3 contribution)')
      call ccsd_timer%init('CCSD contribution)')
!
!     Allocate and zero the transformed singles vector
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      sigma_ai = zero
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai)
!
         enddo
      enddo
!
!     :: CCS contributions to the transformed singles vector ::
!
      call ccsd_timer%start()
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, c_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, c_ai)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_transpose_ccsd_a1(sigma_ai, c_ai)
      call wf%jacobian_transpose_ccsd_b1(sigma_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj, b_end)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, i
!
               if (i .ne. j) then
                  b_end = wf%n_v
               else
                  b_end = a
               endif
!
               do b = 1, b_end
!
                  bj = wf%n_v*(j - 1) + b
!
                  aibj = ai*(ai-3)/2 + ai + bj
!
                  c_aibj(a,i,b,j) = c(wf%n_t1 + aibj)
                  c_aibj(b,j,a,i) = c(wf%n_t1 + aibj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_transpose_ccsd_c1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_d1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_e1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_f1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_g1(sigma_ai, c_aibj)
!
!     :: CCSD contributions to the transformed doubles vector ::
!     Allocate unpacked transformed vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      sigma_aibj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_transpose_ccsd_a2(sigma_aibj, c_ai)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_transpose_ccsd_b2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_c2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_d2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_e2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_f2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_g2(sigma_aibj, c_aibj)
!
!     Compute CC3 contributions to sigma_ai and sigma_aibj and symmetrise sigma_aibj
!     CCSD H2 and I2 are already symmetric and will be computed afterwards
!
      call ccsd_timer%freeze()
!
      call mem%alloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(sigma_aibj, sigma_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call cc3_timer%start()
      call wf%jacobian_transpose_cc3_A(omega, c_ai, c_abij, sigma_ai, sigma_abij)
      call cc3_timer%freeze()
      call cc3_timer%switch_off()
!
!     Done with singles vector c; Overwrite the incoming singles c vector for exit
!
      call ccsd_timer%start()
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, sigma_ai, 1, c, 1)
!
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     Last two CCSD-terms (H2, I2) are already symmetric.
!     Perform the symmetrization sigma_aibj = P_ij^ab sigma_aibj
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(sigma_abij, sigma_aibj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call symmetric_sum(sigma_aibj, (wf%n_v)*(wf%n_o))
!
      call sort_1234_to_1324(sigma_aibj, sigma_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Compute CCSD H2 and I2 contributions
!
      call wf%jacobian_transpose_ccsd_h2(sigma_abij, c_abij)
      call wf%jacobian_transpose_ccsd_i2(sigma_abij, c_abij)
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call ccsd_timer%freeze()
      call ccsd_timer%switch_off()
!
!     overwrite the incoming, packed doubles c vector for exit
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj, b_end)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, i
!
               if (j .ne. i) then
                  b_end = wf%n_v
               else
                  b_end = a
               endif
!
               do b = 1, b_end
!
                  bj = wf%n_v*(j - 1) + b
!
                  aibj = ai*(ai-3)/2 + ai + bj
!
                  c(wf%n_t1 + aibj) = sigma_abij(a,b,i,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_A_cc3(wf, omega, c_ai, c_abij, sigma_ai, sigma_abij)
!!
!!    Terms of the transpose of the  CC3 Jacobi matrix
!!
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
   end subroutine jacobian_transpose_cc3_A_cc3
!
!
   module subroutine prepare_cc3_jacobian_transpose_integrals_cc3(wf)
!!
!!    Construct integrals needed in CC3 jacobian transpose and store on disk
!!    (ab|cd) ordered as abc,d
!!    (mi|lk) ordered as lm,ik
!!    (ib|kd) ordered as bd,ik
!!    (le|ck) ordered as lce,k
!!    (cd|mk) ordered as dcm,k
!!
!!    written by Rolf H. Myhre and Alexander Paul, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs !Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs !Array for sorted integrals
!
      integer :: i, k, d, record
      type(batching_index) :: batch_d, batch_k
!
      integer :: req_0, req_d, req_k
      integer :: d_batch, k_batch
!
      integer :: ioerror=-1
!
!     (be|cd)
!
      call wf%get_g_pqrs_required(req_0,req_d,wf%n_v,wf%n_v,wf%n_v,1)
      req_d = req_d + wf%n_v**3
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d,req_0,req_d)
      call batch_d%determine_limits(1)
!
      call mem%alloc(g_pqrs,wf%n_v,wf%n_v,wf%n_v,batch_d%length)
!
      call wf%g_becd_t%init('g_becd_t','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%g_becd_t,'write')
!
      do d_batch = 1,batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
         call wf%get_vvvv(g_pqrs, &
                          1,wf%n_v, &
                          1,wf%n_v, &
                          1,wf%n_v, &
                          batch_d%first,batch_d%last)
!
         do d = 1,batch_d%length
!
            record = batch_d%first + d - 1
!
            write(wf%g_becd_t%unit,rec=record,iostat=ioerror) g_pqrs(:,:,:,d)
!
            if(ioerror .ne. 0) then
               call output%error_msg('Failed to write becd_t file')
            endif
!
         enddo
!
      enddo
!
      call disk%close_file(wf%g_becd_t,'keep')
!
      call batch_d%determine_limits(1)
      call mem%dealloc(g_pqrs,wf%n_v,wf%n_v,wf%n_v,batch_d%length)
!
!
!     (mi|lk) 
!
      call wf%get_g_pqrs_required(req_0,req_k,wf%n_o,wf%n_o,wf%n_o,1)
      req_k = req_k + 2*wf%n_o**3
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
      call batch_k%determine_limits(1)
!
      call mem%alloc(g_pqrs , wf%n_o , wf%n_o , wf%n_o , batch_k%length)
      call mem%alloc(h_pqrs , wf%n_o , wf%n_o , wf%n_o , batch_k%length)
!
      call wf%g_milk_t%init('g_milk_t','direct','unformatted',dp*wf%n_o**2)
      call disk%open_file(wf%g_milk_t,'write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call wf%get_oooo(g_pqrs, &
                          1,wf%n_o, &
                          1,wf%n_o, &
                          1,wf%n_o, &
                          batch_k%first,batch_k%last)
!
         call sort_1234_to_3124(g_pqrs , h_pqrs , wf%n_o , wf%n_o , wf%n_o , batch_k%length)
!
         do k = 1,batch_k%length
            do i = 1, wf%n_o
!
               record = (batch_k%first + k - 2)*wf%n_o + i
               write(wf%g_milk_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,i,k)
!
               if(ioerror .ne. 0) then
                  call output%error_msg('Failed to write milk_t file')
               endif
!
            enddo
         enddo
!
      enddo
!
      call disk%close_file(wf%g_milk_t,'keep')
!
      call batch_k%determine_limits(1)
      call mem%dealloc(g_pqrs,wf%n_o,wf%n_o,wf%n_o,batch_k%length)
      call mem%dealloc(h_pqrs,wf%n_o,wf%n_o,wf%n_o,batch_k%length)
!
!
!     (ib|kd) 
!
      call wf%get_g_pqrs_required(req_0,req_k,wf%n_v,wf%n_o,1,wf%n_o)
      req_k = req_k + 2*wf%n_v**2*wf%n_o
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%g_ibkd_t%init('g_ibkd_t','direct','unformatted',dp*wf%n_v**2)
      call disk%open_file(wf%g_ibkd_t,'write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o , wf%n_v , batch_k%length , wf%n_v)
         call mem%alloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
         call wf%get_ovov(g_pqrs, &
                          1,wf%n_o, &
                          1,wf%n_v, &
                          batch_k%first,batch_k%last, &
                          1,wf%n_v)
!
         call sort_1234_to_2413(g_pqrs , h_pqrs , wf%n_o , wf%n_v , batch_k%length , wf%n_v)
!
         do k = 1,batch_k%length
            do i = 1, wf%n_o
!
               record = (batch_k%first + k - 2)*wf%n_o + i
               write(wf%g_ibkd_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,i,k)
!
               if(ioerror .ne. 0) then
                  call output%error_msg('Failed to write ibkd_t file')
               endif
!
            enddo
         enddo
!
         call mem%dealloc(g_pqrs, wf%n_o , wf%n_o , batch_k%length , wf%n_o)
         call mem%dealloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_ibkd_t,'keep')
!
!
!     (le|ck) 
!
      call wf%get_g_pqrs_required(req_0, req_k, wf%n_o, wf%n_v, wf%n_v, 1)
      req_k = req_k + 2*wf%n_v**2*wf%n_o
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
      call batch_k%determine_limits(1)
!
      call mem%alloc(g_pqrs, wf%n_o, wf%n_v, wf%n_v, batch_k%length)
      call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_v, batch_k%length)
!
      call wf%g_leck_t%init('g_leck_t','direct','unformatted',dp*wf%n_o*wf%n_v**2)
      call disk%open_file(wf%g_leck_t,'write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call wf%get_ovov(g_pqrs, &
                          1,wf%n_o, &
                          1,wf%n_v, &
                          1,wf%n_v, &
                          batch_k%first,batch_k%last)
!
         call sort_1234_to_1324(g_pqrs , h_pqrs , wf%n_o , wf%n_v , wf%n_v , batch_k%length)
!
         do k = 1,batch_k%length
!
            record = batch_k%first + k - 1
            write(wf%g_leck_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,:,k)
!
            if(ioerror .ne. 0) then
               call output%error_msg('Failed to write leck_t file')
            endif
!
         enddo
!
      enddo
!
      call disk%close_file(wf%g_leck_t,'keep')
!
      call batch_k%determine_limits(1)
      call mem%dealloc(g_pqrs, wf%n_o, wf%n_v, wf%n_v, batch_k%length)
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
!
!
!     (cd|mk) 
!
      call wf%get_g_pqrs_required(req_0, req_k, wf%n_v, wf%n_v, wf%n_o, 1)
      req_k = req_k + 2*wf%n_v**2*wf%n_o
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
      call batch_k%determine_limits(1)
!
      call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
      call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
!
      call wf%g_cdmk_t%init('g_cdmk_t','direct','unformatted',dp*wf%n_o*wf%n_v**2)
      call disk%open_file(wf%g_cdmk_t,'write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call wf%get_vvoo(g_pqrs, &
                          1,wf%n_v, &
                          1,wf%n_v, &
                          1,wf%n_o, &
                          batch_k%first,batch_k%last)
!
         call sort_1234_to_2134(g_pqrs , h_pqrs , wf%n_o , wf%n_v , wf%n_v , batch_k%length)
!
         do k = 1,batch_k%length
!
            record = batch_k%first + k - 1
            write(wf%g_cdmk_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,:,k)
!
            if(ioerror .ne. 0) then
               call output%error_msg('Failed to write cdmk_t file')
            endif
!
         enddo
!
      enddo
!
      call disk%close_file(wf%g_cdmk_t,'keep')
!
      call batch_k%determine_limits(1)
      call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%length)
!
!
   end subroutine prepare_cc3_jacobian_transpose_integrals_cc3
!
!
   module subroutine prepare_cc3_jacobian_transpose_intermediates_cc3(wf)
!!
!!    Construct some intermediates needed in CC3 jacobian transpose and store on disk
!!
!!    written by Rolf H. Myhre and Alexander Paul, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abji
!
!     Arrays for intermediates
!     cannot hold the whole X_acdi array
      real(dp), dimension(:,:,:,:), allocatable, target :: X_acdi
      real(dp), dimension(:,:,:,:), allocatable, target :: X_acdj
      real(dp), dimension(:,:,:,:), allocatable, target :: X_acdk
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_acdi_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_acdj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_acdk_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_aikl
      real(dp), dimension(:,:,:,:), allocatable :: Y_akil
!
!     Integrals and Pointers
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch ! used for the current batch
      integer :: req_0, req_1, req_2, req_3
      real(dp) :: batch_buff = 0.0
!
      integer :: ioerror = -1
      integer :: l
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(t_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1342(wf%t2, t_abji, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      req_0 = 0
      req_1 = wf%n_v**3
      req_2 = wf%n_o*wf%n_v + 2*wf%n_v**2
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
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%alloc(g_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%alloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%alloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
!
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
!
         call mem%alloc(g_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
!
      endif
!
      call mem%alloc(Y_aikl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call disk%open_file(wf%g_bdck_t,'read')
      call disk%open_file(wf%g_ljck_t,'read')
      call disk%open_file(wf%g_ibkd_t,'read')
!
      call wf%X_acdi%init('X_acdi','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%X_acdi,'readwrite')
!
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%single_batch_reader(batch_i, wf%g_bdck_t, g_bdci)
         g_bdci_p => g_bdci
!
!           cannot hold X_acdi - read in previous X, add contributions, write to disk again
!
            if (i_batch .gt. 1) then
               call wf%single_batch_reader(batch_i, wf%X_acdi, X_acdi)
               X_acdi_p => X_acdi
            end if
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%double_batch_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci, wf%g_ibkd_t, g_jbic)
            g_ljci_p => g_ljci
            g_jbic_p => g_jbic
!
            if (j_batch .ne. i_batch) then ! read for switched i - j
!
               call wf%single_batch_reader(batch_j, wf%g_bdck_t, g_bdcj)
               g_bdcj_p => g_bdcj
!
!              Don't read X in the first iteration - X_acdi file empty
               if (i_batch .gt. 1 .or. j_batch .gt. 1) then
                  call wf%single_batch_reader(batch_j, wf%X_acdi, X_acdj)
                  X_acdj_p => X_acdj
               end if
!
               call wf%double_batch_reader(batch_i, batch_j, wf%g_ljck_t, g_licj, wf%g_ibkd_t, g_ibjc)
               !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_i, batch_j, g_licj, g_ibjc, L_ibjc)
               g_licj_p => g_licj
               g_ibjc_p => g_ibjc
!
            else
!
               g_bdcj_p => g_bdci
!
               X_acdj_p => X_acdi
!
               g_licj_p => g_ljci
               g_ibjc_p => g_jbic
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. i_batch .and. k_batch .ne. j_batch) then
!
                  call wf%single_batch_reader(batch_k, wf%g_bdck_t, g_bdck)
                  g_bdck_p => g_bdck
!
!                 Don't read X in the first iteration - X_acdi file empty
                  if (i_batch .gt. 1 .or. j_batch .gt. 1 .or. k_batch .gt. 1) then
                     call wf%single_batch_reader(batch_k, wf%X_acdi, X_acdk)
                     X_acdk_p => X_acdk
                  endif
!
                  call wf%double_batch_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, wf%g_ibkd_t, g_kbic)
                  !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_k, batch_i, g_lkci, g_kbic, L_kbic)
                  g_lkci_p => g_lkci
                  g_kbic_p => g_kbic
!
                  call wf%double_batch_reader(batch_i, batch_k, wf%g_ljck_t, g_lick, wf%g_ibkd_t, g_ibkc)
                  !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_i, batch_k, g_lick, g_ibkc, L_ibkc)
                  g_lick_p => g_lick
                  g_ibkc_p => g_ibkc
!
                  call wf%double_batch_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_ibkd_t, g_kbjc)
                  !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_k, batch_j, g_lkcj, g_kbjc, L_kbjc)
                  g_lkcj_p => g_lkcj
                  g_kbjc_p => g_kbjc
!
                  call wf%double_batch_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck, wf%g_ibkd_t, g_jbkc)
                  !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_j, batch_k, g_ljck, g_jbkc, L_jbkc)
                  g_ljck_p => g_ljck
                  g_jbkc_p => g_jbkc
!
               else if (k_batch .eq. i_batch) then
!
                  g_bdck_p => g_bdci
!
                  X_acdk_p => X_acdi
!
                  if (j_batch .eq. i_batch) then
!
                     g_lkci_p => g_ljci
                     g_kbic_p => g_jbic
!
                     g_lick_p => g_ljci
                     g_ibkc_p => g_jbic
!
                     g_lkcj_p => g_ljci
                     g_kbjc_p => g_jbic
!
                     g_ljck_p => g_ljci
                     g_jbkc_p => g_jbic
!
                  else
!
                     call wf%double_batch_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, wf%g_ibkd_t, g_kbic)
                     !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_k, batch_i, g_lkci, g_kbic, L_kbic)
                     g_lkci_p => g_lkci
                     g_kbic_p => g_kbic
!
                     g_lick_p => g_lkci
                     g_ibkc_p => g_kbic
!
                     g_lkcj_p => g_licj
                     g_kbjc_p => g_ibjc
!
                     g_ljck_p => g_ljci
                     g_jbkc_p => g_jbic
!
                  endif
!
               else if (k_batch .eq. j_batch) then
!
                  g_bdck_p => g_bdcj
!
                  X_acdk_p => X_acdj
!
                  g_lkci_p => g_ljci
                  g_kbic_p => g_jbic
!
                  g_lick_p => g_ljci
                  g_ibkc_p => g_jbic
!
                  call wf%double_batch_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_ibkd_t, g_kbjc)
                  !call wf%jacobian_transpose_cc3_ov_vv_reader(batch_k, batch_j, g_lkcj, g_kbjc, L_kbjc)
                  g_lkcj_p => g_lkcj
                  g_kbjc_p => g_kbjc
!
                  g_ljck_p => g_lkcj
                  g_jbkc_p => g_kbjc
!
               endif
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct t^{abc}_{ijk} for given i, j, k
!
                        call wf%omega_cc3_W_calc(i, j, k, t_abc, u_abc, t_abji,  &
                                                g_bdci_p(:,:,:,i_rel),           &
                                                g_bdcj_p(:,:,:,j_rel),           &
                                                g_bdck_p(:,:,:,k_rel),           &
                                                g_ljci_p(:,:,j_rel,i_rel),       &
                                                g_lkci_p(:,:,k_rel,i_rel),       &
                                                g_lkcj_p(:,:,k_rel,j_rel),       &
                                                g_licj_p(:,:,i_rel,j_rel),       &
                                                g_lick_p(:,:,i_rel,k_rel),       &
                                                g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_eps(i, j, k, t_abc)
!
                        call wf%construct_X_and_Y(i, j, k, t_abc, u_abc,      &
                                                   X_acdi_p(:,:,:,i_rel),     &
                                                   X_acdj_p(:,:,:,j_rel),     &
                                                   X_acdk_p(:,:,:,k_rel),     &
                                                   Y_aikl,                    &
                                                   g_jbic_p(:,:,j_rel,i_rel), &
                                                   g_kbic_p(:,:,k_rel,i_rel), &
                                                   g_kbjc_p(:,:,k_rel,j_rel), &
                                                   g_ibjc_p(:,:,i_rel,j_rel), &
                                                   g_ibkc_p(:,:,i_rel,k_rel), &
                                                   g_jbkc_p(:,:,j_rel,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
               call wf%jacobian_transpose_cc3_write_X(batch_k, X_acdk)
!
            enddo ! batch_k
!
            call wf%jacobian_transpose_cc3_write_X(batch_j, X_acdj)
!
         enddo ! batch_j
!
         call wf%jacobian_transpose_cc3_write_X(batch_i, X_acdi)
!
      enddo ! batch_i
!
!     Close files: 
!
      call disk%close_file(wf%g_bdck_t)
      call disk%close_file(wf%g_ljck_t)
      call disk%close_file(wf%g_ibkd_t)
      !call disk%close_file(wf%L_jbkc_t)
!
!     Allocate integral arrays and assign pointers.
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
         !call mem%dealloc(L_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%dealloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%dealloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
!
         call mem%dealloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
!
         call mem%dealloc(g_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
!
      endif
!
!     Resort X_acdi to X_caid for the final contraction with C^ac_il to sigma_dl
!
      call wf%sort_X_to_caid_and_write()
!
!     sort Y_aikl to akil and write to disk 
!
      call mem%alloc(Y_akil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(Y_aikl, Y_akil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(Y_aikl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%Y_akil%init('Y_akil','direct','unformatted', dp*(wf%n_v)*(wf%n_o)**2)
      call disk%open_file(wf%Y_akil,'write')
!
      do l = 1, wf%n_o
!
         write(wf%Y_akil%unit, rec=l, iostat=ioerror) Y_akil(:,:,:,l)
!
      enddo
!
      if(ioerror .ne. 0) then
         call output%error_msg('Failed to write Y_akil file')
      endif
!
      call mem%dealloc(Y_akil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call disk%close_file(wf%Y_akil)
!
   end subroutine prepare_cc3_jacobian_transpose_intermediates_cc3
!
!
   module subroutine construct_X_and_Y_cc3(wf, i, j, k, t_abc, u_abc, X_acdi, X_acdj, X_acdk, Y_aikl, &
                                             g_jbic, g_kbic, g_kbjc, g_ibjc, g_ibkc, g_jbkc)
!!
!!    Constructs the intermediates X_acdi and Y_akil used to compute the contributions to sigma_ai
!!
!!    X_acdi = sum_bjk (t^bac_ijk * g_jbkd + t^abc_ijk * g_jdkb - 2 * t^abc_ijk * g_jbkd)
!!    Y_akil = sum_bjc (t^bac_ijk * g_jblc + t^abc_ijk * g_jclb - 2 * t^abc_ijk * g_jblc)
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_acdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_acdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_acdk
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout)   :: Y_aikl
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jbkc
!                       
   end subroutine construct_X_and_Y_cc3
!
!
   module subroutine jacobian_transpose_cc3_write_X_cc3(wf, batch_x, X_acdx)
!!
!!    Write the contributions to the X_acdi intermediate to file in the respective batches
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_x%length), intent(in) :: X_acdx
!
      integer :: ioerror
      integer :: x, record
!
      do x = 1, batch_x%length
!
         record = batch_x%first + x -1
         write(wf%X_acdi%unit, rec=record, iostat=ioerror) X_acdx(:,:,:,x)
!
      enddo
!
   end subroutine jacobian_transpose_cc3_write_X_cc3
!
!
   module subroutine sort_X_to_caid_and_write_cc3(wf)
!!
!!    Read in intermediate X_acdi from file, resort to X_caid and write to file again
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: X_acdi
      real(dp), dimension(:,:,:,:), allocatable :: X_caid
!
      type(batching_index) :: batch_i
      integer :: record
      integer :: i_batch, i, i_abs, d
      integer :: req_0, req_i
!
      integer :: ioerror
      character(len=100) :: iom
!
      req_0 = 0
      req_i = 2*wf%n_v**3
!
      call batch_i%init(wf%n_o)
!
      call mem%batch_setup(batch_i, req_0, req_i)
!
      call wf%X_caid%init('X_caid','direct','unformatted', dp*(wf%n_v)**2)
      call disk%open_file(wf%X_caid,'write')
!
      do i_batch = 1, batch_i%num_batches
!
         call mem%alloc(X_acdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
!        Read from file
!
         do i = 1, batch_i%length
!
            i_abs = batch_i%first + i - 1
!
            read(wf%X_acdi%unit, rec=i_abs, iostat=ioerror, iomsg=iom) X_acdi(:,:,:,i)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read X_acdi file in sort_X_to_caid_and_write_cc3'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
!        Sort X_acdi to X_caid
!
         call mem%alloc(X_caid, wf%n_v, wf%n_v, batch_i%length, wf%n_v)
!
         call sort_1234_to_2143(X_acdi, X_caid, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(X_acdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
!        Write to file X_caid
!
         do d = 1, wf%n_v
            do i = 1, batch_i%length
!
               record  = (d - 1)*wf%n_o + batch_i%first + i - 1
               write(wf%X_caid%unit, rec=record, iostat=ioerror) X_caid(:,:,i,d)
!
            enddo
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write X_caid file')
         endif
!
         call mem%dealloc(X_caid, wf%n_v, wf%n_v, batch_i%length, wf%n_v)
!
      enddo ! batch_i
!
      call disk%close_file(wf%X_caid)
!
   end subroutine sort_X_to_caid_and_write_cc3
!
!
end submodule jacobian_transpose
