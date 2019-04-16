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
!     Compute CC3 contributions to sigma_ai and sigma_aibj and symmetrize sigma_aibj
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
!
!     CC3-Contributions from the T3-amplitudes
      call wf%jacobian_transpose_cc3_sigma1_T3_A1(c_abij, sigma_ai)
      call wf%jacobian_transpose_cc3_sigma1_T3_B1(c_abij, sigma_ai)
!
!     CC3-Contributions from the C3-amplitudes
      call wf%jacobian_transpose_cc3_C3_terms(omega, c_ai, c_abij, sigma_ai, sigma_abij)
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
   module subroutine jacobian_transpose_cc3_sigma1_t3_A1_cc3(wf, c_abij, sigma_ai)
!!
!!    Computes first contribution of the T3 amplitudes to sigma_1
!!
!!    Reads in the intermediates X_abid and Y_akil prepared in prepare_jacobian_transpose
!!    contracts with c_abij and adds to sigma_ai
!!
!!    sigma_dl =  sum_abi X_abid * C_abil + sum_aik C_daki * Y_akil
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
      real(dp), dimension(:,:,:,:), allocatable :: Y_akil
!
      type(batching_index) :: batch_d
      integer :: d_batch
      integer :: req_0, req_d
!
!     :: X_abid term ::
!
      call batch_d%init(wf%n_v)
!
      req_0 = 0
      req_d = wf%n_o * wf%n_v**2
!
      call mem%batch_setup(batch_d, req_0, req_d)
!
      call disk%open_file(wf%X_abid,'read')
!
      call batch_d%determine_limits(1)
      call mem%alloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      do d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
!        Read in X_abid written with compound index "id" as record
!
         call compound_record_reader(wf%n_o, batch_d, wf%X_abid, X_abid)
!
         call dgemm('T','N',                    & ! X is transposed
                     batch_d%length,            &
                     wf%n_o,                    &
                     wf%n_o * wf%n_v**2,        &
                     one,                       &
                     X_abid,                    & ! X_d_abi
                     wf%n_o * wf%n_v**2,        &
                     c_abij,                    & ! c_abi_l
                     wf%n_o * wf%n_v**2,        &
                     one,                       &
                     sigma_ai(batch_d%first,1), & ! sigma_dl
                     wf%n_v)
!
      enddo
!
      call batch_d%determine_limits(1)
      call mem%dealloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      call disk%close_file(wf%X_abid)
!
!     :: Y_akil term ::
!
      call disk%open_file(wf%Y_akil,'read')
!
      call mem%alloc(Y_akil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call single_record_reader(wf%n_o, wf%Y_akil, Y_akil)
!
      call disk%close_file(wf%Y_akil)
!
      call dgemm('N','N',              &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v * wf%n_o**2,  &
                  one,                 &
                  c_abij,              & ! C_d_aik
                  wf%n_v,              &
                  Y_akil,              & ! Y_aik_l
                  wf%n_v * wf%n_o**2,  &
                  one,                 &
                  sigma_ai,            & ! sigma_dl
                  wf%n_v)
!
      call mem%dealloc(Y_akil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call disk%close_file(wf%Y_akil)
!
   end subroutine jacobian_transpose_cc3_sigma1_t3_A1_cc3
!
!
   module subroutine jacobian_transpose_cc3_sigma1_t3_B1_cc3(wf, c_abij, sigma_ai)
!!
!!    Computes first contribution of the T3 amplitudes to sigma_1
!!
!!    Constructs t^abc_ijk for fixed ijk and contracts with c_abij
!!    The intermediate X_ck is then contracted with L_kcld
!!
!!    sigma_dl =  sum_abcijk C^ab_ij (t^abc_ijk - t^acb_ijk) L_kcld
!!             =  sum_ck X_ck * L_kcld
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abji
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
      real(dp), dimension(:,:,:,:), allocatable :: L_kcld
      real(dp), dimension(:,:,:,:), allocatable :: L_jbkc
!
!     Intermediate
      real(dp), dimension(:,:), allocatable :: X_ck
!
      type(batching_index) :: batch_i, batch_j, batch_k, batch_l
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: i_batch, j_batch, k_batch, l_batch ! used for the current batch
      integer :: req_0, req_1, req_2, req_3
      real(dp) :: batch_buff = 0.0
!
!     :: Construct intermediate X_ck ::
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(t_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1342(wf%t2, t_abji, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Array for the whole intermediate X_ck
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
      X_ck = zero
!
!     Setup and Batching loops
!
      req_0 = 0
      req_1 = wf%n_v**3
      req_2 = wf%n_o * wf%n_v
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
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
      call disk%open_file(wf%g_bdck_t,'read')
      call disk%open_file(wf%g_ljck_t,'read')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci)
         g_bdci_p => g_bdci
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci)
            g_ljci_p => g_ljci
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj)
               g_bdcj_p => g_bdcj
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj)
               g_licj_p => g_licj
!
            else
!
               g_bdcj_p => g_bdci
!
               g_licj_p => g_ljci
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. i_batch .and. k_batch .ne. j_batch) then
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck)
                  g_bdck_p => g_bdck
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci)
                  g_lkci_p => g_lkci
!

                  call compound_record_reader(batch_i, batch_k, wf%g_ljck_t, g_lick)
                  g_lick_p => g_lick
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
                  g_lkcj_p => g_lkcj
!

                  call compound_record_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck)
                  g_ljck_p => g_ljck
!
!
               else if (k_batch .eq. i_batch) then
!
                  g_bdck_p => g_bdci
!
                  if (j_batch .eq. i_batch) then
!
                     g_lkci_p => g_ljci
!
                     g_lick_p => g_ljci
!
                     g_lkcj_p => g_ljci
!
                     g_ljck_p => g_ljci
!
                  else
!
                     call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci)
!
                     g_lkci_p => g_lkci
!
                     g_lick_p => g_lkci
!
                     g_lkcj_p => g_licj
!
                     g_ljck_p => g_ljci
!
                  endif
!
               else if (k_batch .eq. j_batch) then
!
                  g_lkci_p => g_ljci
!
                  g_lick_p => g_licj
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
!
                  g_lkcj_p => g_lkcj
!
                  g_ljck_p => g_lkcj
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
                                                   g_bdci_p(:,:,:,i_rel),        &
                                                   g_bdcj_p(:,:,:,j_rel),        &
                                                   g_bdck_p(:,:,:,k_rel),        &
                                                   g_ljci_p(:,:,j_rel,i_rel),    &
                                                   g_lkci_p(:,:,k_rel,i_rel),    &
                                                   g_lkcj_p(:,:,k_rel,j_rel),    &
                                                   g_licj_p(:,:,i_rel,j_rel),    &
                                                   g_lick_p(:,:,i_rel,k_rel),    &
                                                   g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_eps(i, j, k, t_abc)
!
                        call wf%jacobian_transpose_cc3_X_ck_calc(i, j, k, t_abc, u_abc,   &
                                                                  X_ck, c_abij)
!
                     enddo ! loop over k
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
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!
!     :: sigma_dl = sum_ck X_ck * L_kcld ::
!
!
      req_0 = 0
      req_1 = wf%n_v**2 * wf%n_o
!
      call batch_l%init(wf%n_o)
      call mem%batch_setup(batch_l, req_0, req_1)
      call batch_l%determine_limits(1)
!
      call mem%alloc(L_jbkc, wf%n_v, wf%n_v, wf%n_o, batch_l%length) ! ordered bcjk
      call mem%alloc(L_kcld, wf%n_v, wf%n_o, wf%n_v, batch_l%length) ! ordered ckdl
!
      do l_batch = 1, batch_l%num_batches
!
         call compound_record_reader(wf%n_o, batch_l, wf%L_jbkc_t, L_jbkc)
!
         call sort_1234_to_1324(L_jbkc, L_kcld, wf%n_v, wf%n_v, wf%n_o, batch_l%length)
!
!        sigma_dl += sum_ck X_ck * L_kcld
!
         call dgemm('T','N',                    &
                     wf%n_v * batch_l%length,   &
                     1,                         &
                     wf%n_v * wf%n_o,           &
                     one,                       &
                     L_kcld,                    & ! L_dl_ck
                     wf%n_v * wf%n_o,           &
                     X_ck,                      & ! X_ck
                     wf%n_v * wf%n_o,           &
                     one,                       &
                     sigma_ai(1,batch_l%first), & ! sigma_dl
                     wf%n_v * wf%n_o)
!
      enddo
!
      call batch_l%determine_limits(1)
      call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, wf%n_o, batch_l%length)
      call mem%dealloc(L_kcld, wf%n_v, wf%n_o, wf%n_v, batch_l%length)
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_cc3_sigma1_t3_B1_cc3
!
!
   module subroutine jacobian_transpose_cc3_X_ck_calc_cc3(wf, i, j, k, t_abc, u_abc, X_ck, c_abij)
!!
!!    Constructs the intermediate X_ck
!!
!!    X_ck =  sum_abcijk C^ab_ij (t^abc_ijk - t^acb_ijk)
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: X_ck
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: c_abij
!
      if (j .ne. k) then ! t_abc - t_acb = 0 if j==k
!
!     Construct u_abc = t_abc - t_acb
!
      call construct_123_minus_132(t_abc, u_abc, wf%n_v)
!
!     X_ck += sum_ab (t^abc - t^acb) * C_abij
!
         call dgemv('T',               &
                     wf%n_v**2,        & ! dim of c
                     wf%n_v,           & ! dim of X
                     one,              &
                     u_abc,            & ! u_c_ab
                     wf%n_v**2,        &
                     c_abij(:,:,i,j),  & ! c_ab_ij
                     1,                &
                     one,              &
                     X_ck(:,k),        & ! X_ck
                     1)
!
!        X_cj += sum_ab C_abik * (t^acb - t^abc)
!
         call dgemv('T',               &
                     wf%n_v**2,        & ! dim of c
                     wf%n_v,           & ! dim of X
                     -one,             &
                     u_abc,            & ! u_c_ab
                     wf%n_v**2,        &
                     c_abij(:,:,i,k),  & ! c_ab_ik
                     1,                &
                     one,              &
                     X_ck(:,j),        & ! X_cj
                     1)
!
      end if
!
      if (i .ne. j) then
!
!        Construct u_abc = t_bac - t_cab
!
         call construct_213_minus_312(t_abc, u_abc, wf%n_v)
!
!        X_ck += sum_ab C_abji * (t^bac - t^cab)
!
         call dgemv('T',               &
                     wf%n_v**2,        & ! dim of c
                     wf%n_v,           & ! dim of X
                     one,              &
                     u_abc,            & ! u_c_ab
                     wf%n_v**2,        &
                     c_abij(:,:,j,i),  & ! c_ab_ji
                     1,                &
                     one,              &
                     X_ck(:,k),        & ! X_ck
                     1)
!
         if(j .ne. k) then
!
!           X_ci += sum_ab C_abjk * (t^cab - t^bac)
!
            call dgemv('T',               &
                        wf%n_v**2,        & ! dim of c
                        wf%n_v,           & ! dim of X
                        -one,             &
                        u_abc,            & ! u_c_ab
                        wf%n_v**2,        &
                        c_abij(:,:,j,k),  & ! c_ab_jk
                        1,                &
                        one,              &
                        X_ck(:,i),        & ! X_ci
                        1)
!
         end if
!
!        Construct u_abc = t_cba - t_bca
!
         call construct_321_minus_231(t_abc, u_abc, wf%n_v) ! t_cba - t_bca = 0 if i==j
!
!        X_ci += sum_ab C_abkj * (t_cba - t_bca)
!
         call dgemv('T',               &
                     wf%n_v**2,        & ! dim of c
                     wf%n_v,           & ! dim of X
                     one,              &
                     u_abc,            & ! u_c_ab
                     wf%n_v**2,        &
                     c_abij(:,:,k,j),  & ! c_ab_kj
                     1,                &
                     one,              &
                     X_ck(:,i),        & ! X_ci
                     1)
!
         if(j .ne. k) then
!
!           X_cj += sum_ab C_abki * (t^bca - t^cba)
!
            call dgemv('T',               &
                        wf%n_v**2,        & ! dim of c
                        wf%n_v,           & ! dim of X
                        one,              &
                        u_abc,            & ! u_c_ab
                        wf%n_v**2,        &
                        c_abij(:,:,k,i),  & ! c_ab_ki
                        1,                &
                        one,              &
                        X_ck(:,j),        & ! X_cj
                        1)
!
         end if
      end if
!
   end subroutine jacobian_transpose_cc3_X_ck_calc_cc3


   module subroutine jacobian_transpose_cc3_C3_terms_cc3(wf, omega, c_ai, c_abij, sigma_ai, sigma_abij)
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The construction of C3 is split into contributions 
!!    from outer products and matrix multiplications
!!
!!    1 array for each Permutation of C_abc will be used 
!!    to reduce the amount of N^7-contractions and sorting
!!
!!    c_μ3 = (ω - ε^abc_ijk)^-1 (c_μ1 < μ1 | [H,tau_ν3] | R > + c_μ2 < μ2 | [H,tau_ν3] | R >
!!
!!    σ1 += c_μ3 < μ3 | [[H,T_2],tau_ν1] | R >
!!    σ2 += c_μ3 < μ3 | [H,tau_ ν2] | R >
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
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
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: c_abc
      real(dp), dimension(:,:,:), allocatable :: c_bac
      real(dp), dimension(:,:,:), allocatable :: c_cba
      real(dp), dimension(:,:,:), allocatable :: c_acb
      real(dp), dimension(:,:,:), allocatable :: c_cab
      real(dp), dimension(:,:,:), allocatable :: c_bca
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:), allocatable :: F_kc
!
!
!     Arrays for intermediates
!     cannot hold the whole X_bcek array
      real(dp), dimension(:,:,:,:), allocatable, target :: X_bcei
      real(dp), dimension(:,:,:,:), allocatable, target :: X_bcej
      real(dp), dimension(:,:,:,:), allocatable, target :: X_bcek
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_bcei_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_bcej_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_bcek_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_cmjk
!
!     Integrals and Pointers
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
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_1, req_2, req_3
      real(dp)             :: batch_buff = 0.0
!
!     Set up arrays for amplitudes
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Arrays for the triples amplitudes and intermediates
!
      call mem%alloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(c_bac, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(c_cba, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(c_acb, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(c_cab, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(c_bca, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      c_abc = zero
      c_bac = zero
      c_cba = zero
      c_acb = zero
      c_cab = zero
      c_bca = zero
!
!     Fock matrix subblock: Resorting for easier contractions later
!
      call mem%alloc(F_kc, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_kc, wf%n_o, wf%n_v)
!
!     vooo Intermediate
      call mem%alloc(Y_cmjk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     Setup and Batching loops
!
      req_0 = 0
      req_1 = 2*(wf%n_v)**3
      req_2 = 2*(wf%n_o)*(wf%n_v) + (wf%n_v)**2
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
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(L_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(X_bcei, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         X_bcei = zero
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
!        Ordered such that batching indices are at the end
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
         call mem%alloc(X_bcei, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(X_bcej, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(X_bcek, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         X_bcei = zero
         X_bcej = zero
         X_bcek = zero
!
      endif
!
      call disk%open_file(wf%g_bdck_t,'read')
      call disk%open_file(wf%g_dbkc_t,'read')
      call disk%open_file(wf%g_ljck_t,'read')
      call disk%open_file(wf%g_jlkc_t,'read')
      call disk%open_file(wf%L_jbkc_t,'read')
!
      call wf%X_bcek%init('X_bcek','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%X_bcek,'readwrite')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci, wf%g_dbkc_t, g_dbic)
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
!        cannot hold X_bcei - read in previous X, add contributions, write to disk again
         if (i_batch .gt. 1) then
            call single_record_reader(batch_i, wf%X_bcek, X_bcei)
            X_bcei_p => X_bcei
         end if
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci, wf%g_jlkc_t,  &
                                       g_jlic, wf%L_jbkc_t, L_jbic)
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_jbic_p => L_jbic
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj, wf%g_dbkc_t, g_dbjc)
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
!              Don't read X in the first iteration - X_bcei file empty
               if (i_batch .gt. 1 .or. j_batch .gt. 1) then
                  call single_record_reader(batch_j, wf%X_bcek, X_bcej)
                  X_bcej_p => X_bcej
               end if
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj, wf%g_jlkc_t,  &
                                          g_iljc, wf%L_jbkc_t, L_ibjc)
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
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. i_batch .and. k_batch .ne. j_batch) then
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck, wf%g_dbkc_t, g_dbkc)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
!
!                 Don't read X in the first iteration - X_bcei file empty
                  if (i_batch .gt. 1 .or. j_batch .gt. 1 .or. k_batch .gt. 1) then
                     call single_record_reader(batch_k, wf%X_bcek, X_bcek)
                     X_bcek_p => X_bcek
                  endif
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, wf%g_jlkc_t,  &
                                             g_klic, wf%L_jbkc_t, L_kbic)
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
                  L_kbic_p => L_kbic
!
                  call compound_record_reader(batch_i, batch_k, wf%g_ljck_t, g_lick, wf%g_jlkc_t,  &
                                             g_ilkc, wf%L_jbkc_t, L_ibkc)
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_jlkc_t,  &
                                             g_kljc, wf%L_jbkc_t, L_kbjc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  call compound_record_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck, wf%g_jlkc_t,  &
                                             g_jlkc, wf%L_jbkc_t, L_jbkc)
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
                     call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, wf%g_jlkc_t,  &
                                                g_klic, wf%L_jbkc_t, L_kbic)
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
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_jlkc_t,  &
                                             g_kljc, wf%L_jbkc_t, L_kbjc)
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
!                       Construct C^abc_ijk for given i, j, k
!                       and construct the intermediates X_bcek, Y_cmjk 
!                       and calculate contributions to sigma2
!
                     !   call wf%jacobian_transpose__cc3_calc_outer()
!
                        call wf%jacobian_transpose__cc3_calc_c3_matmul(i, j, k, c_abij, c_abc, c_bac,    & 
                                                                     c_cba, c_acb, c_cab, c_bca, u_abc,  &
                                                                     g_dbic, g_dbjc, g_dbkc,             &
                                                                     g_jlic, g_klic, g_kljc,             &
                                                                     g_iljc, g_ilkc, g_jlkc)
!
                     !   call wf%jacobian_transpose__cc3_collect_c3(omega, c_abc, c_bac,  &
                                                                     !c_cba, c_acb, c_cab, c_bca)
!
                     !   call wf%jacobian_transpose_cc3_sigma2()
!
                     !   call wf%construct_intermediates_c3()
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
               call single_record_writer(batch_k, wf%X_bcek, X_bcek)
!
            enddo ! batch_k
!
            call single_record_writer(batch_j, wf%X_bcek, X_bcej)
!
         enddo ! batch_j
!
         call single_record_writer(batch_i, wf%X_bcek, X_bcei)
!
      enddo ! batch_i
!
      call disk%close_file(wf%g_bdck_t)
      call disk%close_file(wf%g_dbkc_t)
      call disk%close_file(wf%g_ljck_t)
      call disk%close_file(wf%g_jlkc_t)
      call disk%close_file(wf%L_jbkc_t)
!
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
         call mem%dealloc(X_bcei, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
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
         call mem%dealloc(X_bcei, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(X_bcej, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(X_bcek, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
      endif
!
!     Deallocate amplitudes arrays and Fock matrix
!
      call mem%dealloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(c_bac, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(c_cba, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(c_acb, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(c_cab, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(c_bca, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(F_kc, wf%n_v, wf%n_o)
!
      call mem%dealloc(Y_cmjk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call disk%close_file(wf%X_bcek)
!
   end subroutine jacobian_transpose_cc3_C3_terms_cc3
!
!
   module subroutine jacobian_transpose__cc3_calc_c3_matmul_cc3(wf, i, j, k, c_abij, c_abc, c_bac, & 
                                                               c_cba, c_acb, c_cab, c_bca, u_abc,  &
                                                               g_dbic, g_dbjc, g_dbkc,             &
                                                               g_jlic, g_klic, g_kljc,             &
                                                               g_iljc, g_ilkc, g_jlkc)
!!
!!    Calculate the contributions from matrix multiplications 
!!    to the  C3 amplitudes for fixed indices i,j,k
!!
!!    C^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc
!!    + sum_l (C_ablk g_iljc - C_abil L_jlkc) - sum_d (C_adjk g_ibdc - C_adij L_dbkc)
!!
!!    Contibutions:
!!    sum_l (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_bac
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_cba
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_acb
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_cab
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_bca
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
!     g_dbkc ordered bcd,k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbkc
!
!     g_jlkc ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlkc
!
!     :: Contribution 1 ::
!
!     The same contraction contributes to 3 permutations of the indices in c_abc
!     c_acb <- u_abc = sum_l c_aclj g_ilkb
!     c_cba <- u_abc = sum_l c_cblj g_ilka
!     c_abc <- u_abc = sum_l c_ablj g_ilkc
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  c_abij(:,:,:,j),  & ! c_aclj
                  wf%n_v**2,        &
                  g_ilkc,           & ! g_ilkb ordered bl,ik
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     The same contraction contributes to 3 permutations of the indices in c_abc
!     c_acb <- u_abc = - sum_d c_adij g_dckb
!     c_cba <- u_abc = - sum_d c_cdij g_dbka
!     c_abc <- u_abc = - sum_d c_adij g_dbkc
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  -one,             &
                  c_abij(:,:,i,j),  & ! c_adij
                  wf%n_v,           &
                  g_dbkc,           & ! g_dckb ordered cbd,k
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
      c_acb = c_acb + u_abc
      c_cba = c_cba + u_abc
      c_abc = c_abc - two*u_abc
!
!     :: Contribution 2 ::
!
!     c_bca <- u_abc = sum_l c_bcli g_jlka
!     c_cab <- u_abc = sum_l c_cali g_jlkb
!     c_bac <- u_abc = sum_l c_bali g_jlkc
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  c_abij(:,:,:,i),  & ! c_bcli
                  wf%n_v**2,        &
                  g_jlkc,           & ! g_jlka ordered al,jk
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_bca <- u_abc = - sum_d c_bdji g_dcka
!     c_cab <- u_abc = - sum_d c_cdji g_dakb
!     c_bac <- u_abc = - sum_d c_bdji g_dakc
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  -one,             &
                  c_abij(:,:,j,i),  & ! c_bdji
                  wf%n_v,           &
                  g_dbkc,           & ! g_dcka ordered cad,k
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
      c_bca = c_bca + u_abc
      c_cab = c_cab + u_abc
      c_bac = c_bac - two*u_abc
!
!     :: Contribution 3 ::
!
!     c_cab <- u_abc = sum_l c_calj g_klib
!     c_abc <- u_abc = sum_l c_ablj g_klic
!     c_cba <- u_abc = sum_l c_cblj g_klia
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  c_abij(:,:,:,j),  & ! c_calj
                  wf%n_v**2,        &
                  g_klic,           & ! g_klib ordered bl,ki
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_cab <- u_abc = - sum_d c_cdkj g_daib
!     c_abc <- u_abc = - sum_d c_adkj g_dbic
!     c_cba <- u_abc = - sum_d c_cdkj g_dbia
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  -one,             &
                  c_abij(:,:,k,j),  & ! c_cdkj
                  wf%n_v,           &
                  g_dbic,           & ! g_daib ordered abd,i
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
      c_cab = c_cab + u_abc
      c_abc = c_abc + u_abc
      c_cba = c_cba - two*u_abc
!
!     :: Contribution 4 ::
!
!     c_abc <- u_abc = sum_l c_ablk g_iljc
!     c_bca <- u_abc = sum_l c_bclk g_ilja
!     c_acb <- u_abc = sum_l c_aclk g_iljb
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  c_abij(:,:,:,k),  & ! c_ablk
                  wf%n_v**2,        &
                  g_iljc,           & ! g_iljc ordered cl,ij
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_abc <- u_abc = - sum_d c_adik g_dbjc
!     c_bca <- u_abc = - sum_d c_bdik g_dcja
!     c_acb <- u_abc = - sum_d c_adik g_dcjb
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  -one,             &
                  c_abij(:,:,i,k),  & ! c_adik
                  wf%n_v,           &
                  g_dbjc,           & ! g_dbjc ordered bcd,j
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
      c_abc = c_abc + u_abc
      c_bca = c_bca + u_abc
      c_acb = c_acb - two*u_abc
!
!     :: Contribution 5 ::
!
!     c_cba <- u_abc = sum_l c_cbli g_klja
!     c_bac <- u_abc = sum_l c_bali g_kljc
!     c_cab <- u_abc = sum_l c_cali g_kljb
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  c_abij(:,:,:,i),  & ! c_cbli
                  wf%n_v**2,        &
                  g_kljc,           & ! g_klja ordered al,kj
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_cba <- u_abc = - sum_d c_cdki g_dbja
!     c_bac <- u_abc = - sum_d c_bdki g_dajc
!     c_cab <- u_abc = - sum_d c_cdki g_dajb
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  -one,             &
                  c_abij(:,:,k,i),  & ! c_cdki
                  wf%n_v,           &
                  g_dbjc,           & ! g_dbja ordered bad,j
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
      c_cba = c_cba + u_abc
      c_bac = c_bac + u_abc
      c_cab = c_cab - two*u_abc
!
!     :: Contribution 6 ::
!
!     c_bac <- u_abc = sum_l c_balk g_jlic
!     c_acb <- u_abc = sum_l c_aclk g_jlib
!     c_bca <- u_abc = sum_l c_bclk g_jlia
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  c_abij(:,:,:,k),  & ! c_balk
                  wf%n_v**2,        &
                  g_jlic,           & ! g_jlic ordered cl,ji
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_bac <- u_abc = - sum_d c_bdjk g_daic
!     c_acb <- u_abc = - sum_d c_adjk g_dcib
!     c_bca <- u_abc = - sum_d c_bdjk g_dcia
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  -one,             &
                  c_abij(:,:,j,k),  & ! c_bdjk
                  wf%n_v,           &
                  g_dbic,           & ! g_daic ordered acd,i
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
      c_bac = c_bac + u_abc
      c_acb = c_acb + u_abc
      c_bca = c_bca - two*u_abc
!
   end subroutine jacobian_transpose__cc3_calc_c3_matmul_cc3
!
!
end submodule jacobian_transpose
