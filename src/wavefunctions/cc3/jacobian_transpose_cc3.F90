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
      call disk%open_file(wf%L_kcld_t,'read')
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
      call disk%close_file(wf%L_kcld_t)
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
   end subroutine jacobian_transpose_cc3_C3_terms_cc3
!
!
end submodule jacobian_transpose
