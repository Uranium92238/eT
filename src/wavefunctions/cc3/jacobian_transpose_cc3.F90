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
!!
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
!!
      class(cc3) :: wf
!
   end subroutine prepare_cc3_jacobian_transpose_intermediates_cc3
!
!
end submodule jacobian_transpose
