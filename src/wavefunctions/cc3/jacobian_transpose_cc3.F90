!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
!!    Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    sigma_i = A^T * c_i,
!!
!!    where
!!
!!    A_mu,nu = < mu|exp(-T) [H, tau_nu] exp(T)|R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, b, sigma, &
                                                                     cvs, rm_core)
!!
!!    Effective Jacobian transpose transformation (CC3)
!!    Alexander C. Paul and Rolf H. Myhre, March 2019
!!    Adapted to use a covariant representation of L2
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    Directs the transformation by the transpose of the  CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T)|R >,
!!
!!    The transformation is performed as sigma^T = b^T A, where b is the vector
!!    sent to the routine.
!!
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      use array_utilities, only: copy_and_scale, zero_array
      use reordering, only: construct_covariant_1324, symmetrize_add_to_packed
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: b
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: sigma
!
!     Same routines used for tbar3 and L3
      logical, intent(in) :: cvs, rm_core
!
      real(dp), dimension(:,:,:,:), allocatable :: b_abij
      real(dp), dimension(:,:,:,:), allocatable :: sigma_abij
!
      type(timings) :: timer
!
      timer = timings('Jacobian transpose CC3', pl='normal')
      call timer%turn_on()
!
!     Zero the transformed vector
!
      call zero_array(sigma, wf%n_t1 + wf%n_t2)
!
      call wf%ccsd%jacobian_transpose_transformation(b, sigma)
!
!     Compute CC3 contributions to sigma
!
      call mem%alloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(sigma_abij, (wf%n_v*wf%n_o)**2)
!
!     Construct covariant _b_abij = 1/3 (2 b^ab_ij + b^ba_ij)
!
      call mem%alloc(b_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call construct_covariant_1324(b(wf%n_t1+1:), b_abij, wf%n_v, wf%n_o)
!
!     CC3-Contributions from the C3-amplitudes
      call wf%jacobian_transpose_cc3_c3_a(omega, b(1:wf%n_t1), b_abij, &
                                          sigma(1:wf%n_t1), sigma_abij, cvs, rm_core)
!
!     Done with the doubles part pack sigma_abij into sigma
!
      call symmetrize_add_to_packed(sigma_abij, sigma(wf%n_t1 + 1 : wf%n_es_amplitudes), &
                                    wf%n_v, wf%n_o)

!
      call mem%dealloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     CC3-Contributions from the T3-amplitudes
      call wf%jacobian_transpose_cc3_t3_a1(b_abij, sigma(1:wf%n_t1))
      call wf%jacobian_transpose_cc3_t3_b1(b_abij, sigma(1:wf%n_t1), cvs, rm_core)
!
      call mem%dealloc(b_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_t3_a1_cc3(wf, c_abij, sigma_ai)
!!
!!    Jacobian transpose T3 A1 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!    Adapted to use a covariant representation of L2
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    Computes the first contribution of the T3 amplitudes to sigma_1
!!
!!    Reads in the intermediates X_abid and X_ajil prepared in prepare_jacobian_transpose
!!    contracts with a covariant c_abij (_c_abij) and adds to sigma_ai
!!
!!    sigma_dl +=  sum_abi X_abid * _C_abil + sum_aik _C_daji * X_ajil
!!
!!    where:  _c^ab_ij = 1/3 (2c^ab_ij + c^ba_ij)
!!
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
      real(dp), dimension(:,:,:,:), allocatable :: X_ajil
!
      type(batching_index) :: batch_d
      integer :: d_batch
      integer :: req_0, req_d
!
      type(timings) :: timer
!
      timer = timings('Time in CC3 T3 a1', pl='verbose')
      call timer%turn_on()
!
!     :: X_abid term ::
!
      batch_d = batching_index(wf%n_v)
!
      req_0 = 0
      req_d = wf%n_o * wf%n_v**2
!
      call mem%batch_setup(batch_d, req_0, req_d, 'jacobian_transpose_cc3_t3_a1')
!
      call wf%X_abid%open_('read')
!
      call mem%alloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%max_length)
!
      do d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
!        Read in X_abid written with compound index "id" as record
!
         call wf%X_abid%read_compound_full_batch(X_abid, wf%n_o, batch_d)
!
         call dgemm('T','N',                    &
                     batch_d%length,            &
                     wf%n_o,                    &
                     wf%n_o*wf%n_v**2,          &
                     one,                       &
                     X_abid,                    & ! X_d_abi
                     wf%n_o*wf%n_v**2,          &
                     c_abij,                    & ! c_abi_l
                     wf%n_o*wf%n_v**2,          &
                     one,                       &
                     sigma_ai(batch_d%first,1), & ! sigma_dl
                     wf%n_v)
!
      enddo
!
      call mem%dealloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%max_length)
!
      call wf%X_abid%close_()
!
      call mem%batch_finalize()
!
!     :: X_ajil term ::
!
      call wf%X_ajil%open_('read')
!
      call mem%alloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%X_ajil%read_(X_ajil, 1, wf%n_o)
!
      call wf%X_ajil%close_()
!
      call dgemm('N','N',              &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v * wf%n_o**2,  &
                  one,                 &
                  c_abij,              & ! C_d_aji
                  wf%n_v,              &
                  X_ajil,              & ! X_aji_l
                  wf%n_v * wf%n_o**2,  &
                  one,                 &
                  sigma_ai,            & ! sigma_dl
                  wf%n_v)
!
      call mem%dealloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_cc3_t3_a1_cc3
!
!
   module subroutine jacobian_transpose_cc3_t3_b1_cc3(wf, c_abij, sigma_ai, cvs, rm_core)
!!
!!    Jacobian transpose T3 B1 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!    Adapted to use a covariant representation of L2 and a contravariant t3
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    Constructs t^abc_ijk for fixed ijk and contracts with c_abij
!!    The intermediate X_ai is then contracted with L_iald
!!
!!    sigma_dl +=  sum_abcijk C^bc_jk (t^abc_ijk - t^bac_ijk) L_iald
!!             +=  sum_abcijk _C^bc_jk u^abc_ijk L_iald
!!             +=  sum_ai X_ai * L_iald
!!
!!    where: _c^ab_ij  = 1/3 (2 c^ab_ij + c^ba_ij)
!!           u^abc_ijk = 4t^abc_ijk + t_bca_ijk + t_cab_ijk
!!                     - 2t^acb_ijk - 2t_cba_ijk - 2t_bac_ijk
!!
      use array_utilities, only: zero_array
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      logical, intent(in) :: cvs, rm_core
!
      real(dp), dimension(:,:,:), allocatable :: t_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
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
      real(dp), dimension(:,:), allocatable :: X_ai
!
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: i_batch, j_batch, k_batch ! used for the current batch
      integer :: req_0, req_i, req_1, req_2, req_3, req_1_eri
      integer :: req_single_batch
!
      type(timings) :: timer
!
      timer = timings('Time in CC3 T3 b1', pl='verbose')
      call timer%turn_on()
!
!     :: Construct intermediate X_ai ::
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Array for the whole intermediate X_ai
      call mem%alloc(X_ai, wf%n_v, wf%n_o)
      call zero_array(X_ai, wf%n_t1)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + wf%n_v**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + wf%n_v**3*wf%n_o &
                       + wf%n_v*wf%n_o**3
!
      req_1 = wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 2*wf%n_o*wf%n_v
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,       &
                           req_0, req_i, req_1, req_1,      &
                           req_2, req_2, req_2, req_3,      &
                           'jacobian_transpose_cc3_t3_b1',  &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(g_bdci, g_bdci_p, sorting, batch_i)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(g_bdcj, g_bdcj_p, sorting, batch_j)
!
               call wf%setup_oovo(g_licj, g_licj_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%setup_vvvo(g_bdck, g_bdck_p, sorting, batch_k)
!
                  call wf%setup_oovo(g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
!
               else if (k_batch .eq. i_batch) then !k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
!
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
               endif
!
!              sigma_dl +=  sum_abcijk C^bc_jk (t^abc_ijk - t^bac_ijk) L_iald
!              CVS: in principle check j,k and l but due to the symmetry in L_iald
!                   we can also check i,j,k
!
               do i = batch_i%first, batch_i%get_last()
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%get_last(), i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%get_last(), j)
!
!                       Check for core orbitals:
!                       cvs: i,j,k cannot all correspond to valence orbitals
!                       rm_core: i,j,k may not contain any core orbital
!
!                       Here t3 is contracted with L2 and can, thus, be restricted as well
                        if (wf%ijk_amplitudes_are_zero(i, j, k, cvs, rm_core)) cycle
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct t^{abc}_{ijk} for given i, j, k
!
                        call wf%construct_W(i, j, k, sorting, t_abc, t_abij, &
                                            g_bdci_p(:,:,:,i_rel),         &
                                            g_bdcj_p(:,:,:,j_rel),         &
                                            g_bdck_p(:,:,:,k_rel),         &
                                            g_ljci_p(:,:,j_rel,i_rel),     &
                                            g_lkci_p(:,:,k_rel,i_rel),     &
                                            g_lkcj_p(:,:,k_rel,j_rel),     &
                                            g_licj_p(:,:,i_rel,j_rel),     &
                                            g_lick_p(:,:,i_rel,k_rel),     &
                                            g_ljck_p(:,:,j_rel,k_rel))
!
                        call construct_contravariant_t3(t_abc, sorting, wf%n_v)
                        call wf%divide_by_orbital_differences(i, j, k, t_abc)
!
                        call wf%construct_x_ai_intermediate(i, j, k, t_abc, sorting, c_abij, X_ai)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif
!
      call mem%batch_finalize()
!
      call wf%jacobian_transpose_cc3_b1_x_ai(X_ai, sigma_ai)
!
      call mem%dealloc(X_ai, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_cc3_t3_b1_cc3
!
!
   module subroutine construct_x_ai_intermediate_cc3(wf, i, j, k, u_abc, v_abc, L2, X_ai)
!!
!!    Construct X_ai intermediate
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    X_ai = sum_bcjk (t^abc_ijk - t^bac_ijk) L^bc_jk
!!
!!    This can be reformulated to:
!!    X_ai = sum_bcjk u_abc L^bc_jk
!!
!!    where:
!!    u_abc = (4t_abc - 2t_acb - 2t_cba - 2t_bac + t_bca + t_cab)
!!    L^ab_ij = 1/3 (2 L'^ab_ij + L'^ba_ij)
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      use reordering, only: sort_123_to_312, sort_123_to_213
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)      :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)      :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                :: X_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: L2
!
      real(dp) :: factor_ij, factor_jk
      logical  :: skip
!
      if (i .ne. j .and. j .ne. k) then
         factor_ij = one
         factor_jk = one
         skip = .false.
      else if (j .eq. k )  then
         factor_ij = one
         factor_jk = half
         skip = .true.
      else ! i == j
         factor_ij = half
         factor_jk = one
         skip = .true.
      end if
!
!     X_ai += sum_bc u_abc L^bc_jk
      call wf%construct_X_vo_permutation(i, u_abc, L2(:,:,j,k), X_ai, factor_jk)
!
!     v_abc = 4t_bca - 2t_cba - 2t_bac - 2t_acb + t_abc + t_cab
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     X_ak += sum_cb u_abc L^bc_ij
      call wf%construct_X_vo_permutation(k, v_abc, L2(:,:,i,j), X_ai, factor_ij)
!
      if (.not. skip) then
!
!        v_abc = 4t_bac - 2t_bca - 2t_abc - 2t_cab + t_cba + t_acb
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        X_aj += sum_cb v_abc L^bc_ik
         call wf%construct_X_vo_permutation(j, v_abc, L2(:,:,i,k), X_ai, one)
!
      end if
!
   end subroutine construct_x_ai_intermediate_cc3
!
!
   module subroutine jacobian_transpose_cc3_b1_x_ai_cc3(wf, X_ai, sigma_ai)
!!
!!    Jacobian transpose CC3 B1 X_ai
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Handles the final contractions of L_iald with the intermediate
!!       X_ai = _C^bc_jk u^abc_ijk
!!
!!    sigma_dl +=  sum_ai X_ai * L_iald
!!
!!    where: _c^ab_ij  = 1/3 (2 c^ab_ij + c^ba_ij)
!!           u^abc_ijk = 4t^abc_ijk + t_bca_ijk + t_cab_ijk
!!                     - 2t^acb_ijk - 2t_cba_ijk - 2t_bac_ijk
!!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: X_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(:,:,:), allocatable :: temp, L_J_oo, L_J_ov, L_J_ov_ai
!
!     :: sigma_dl += sum_ai X_ai * L_iald ::
!
!     sigma_dl += sum_ai 2X_ai L^J_ia L^J_ld - sum_ai X_ai L^J_la L^J_id
!
      call mem%alloc(L_J_ov, wf%eri_t1%n_J, wf%n_o, wf%n_v)
      call mem%alloc(L_J_ov_ai, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call wf%L_t1%get(L_J_ov, 1, wf%n_o, wf%n_o+1, wf%n_mo)
!
      call sort_123_to_132(L_J_ov, L_J_ov_ai, wf%eri_t1%n_J, wf%n_o, wf%n_v)
!
      call mem%alloc(temp, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
!     First term
!
      call dgemv('N',            &
                 wf%eri_t1%n_J,  &
                 wf%n_t1,        &
                 two,            &
                 L_J_ov_ai,      & ! L_J_ai
                 wf%eri_t1%n_J,  &
                 X_ai, 1,        & ! X_ai
                 zero,           &
                 temp, 1)      ! L_J
!
      call dgemv('T',            &
                 wf%eri_t1%n_J,  &
                 wf%n_t1,        &
                 one,            &
                 L_J_ov_ai,      & ! L_dl_J
                 wf%eri_t1%n_J,  &
                 temp, 1,        & ! L_J
                 one,            &
                 sigma_ai, 1)

      call mem%dealloc(L_J_ov_ai, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
!     Second term
!
      call dgemm('N', 'N',             &
                 wf%eri_t1%n_J*wf%n_o, &
                 wf%n_o,               &
                 wf%n_v,               &
                 one,                  &
                 L_J_ov,               & ! L_J_la
                 wf%eri_t1%n_J*wf%n_o, &
                 X_ai,                 & ! X_ai
                 wf%n_v,               &
                 zero,                 &
                 temp,                 & ! L_J_li
                 wf%eri_t1%n_J*wf%n_o)
!
      call mem%alloc(L_J_oo, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
!     L_J_li -> L_J_il
      call sort_123_to_132(temp, L_J_oo, wf%eri_t1%n_J, wf%n_o, wf%n_o)

      call mem%dealloc(temp, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                 wf%n_v,               &
                 wf%n_o,               &
                 wf%eri_t1%n_J*wf%n_o, &
                 -one,                 &
                 L_J_ov,               & ! L_Ji_d
                 wf%eri_t1%n_J*wf%n_o, &
                 L_J_oo,               & ! L_Ji_l
                 wf%eri_t1%n_J*wf%n_o, &
                 one,                  &
                 sigma_ai,             &
                 wf%n_v)
!
      call mem%dealloc(L_J_ov, wf%eri_t1%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(L_J_oo, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_cc3_b1_x_ai_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a_cc3(wf, omega,            &
                                                     c_ai, c_abij,         &
                                                     sigma_ai, sigma_abij, &
                                                     cvs, rm_core)
!!
!!    Jacobian transpose CC3 c3 A
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!    Adapted to first construct a covariant intermediate for L3 from
!!    which the full L3 is obtained by a linear combination
!!    by Alexander C. Paul and Rolf H. Myhre, Sep 2020
!!
!!    For that c^ab_ij = 1/3 (2 L^ab_ij + L^ba_ij)
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
!!    c_mu3 = (omega - epsilon^abc_ijk)^-1 (c_mu1 < mu1| [H,tau_nu3] |R >
!!                                        + c_mu2 < mu2| [H,tau_nu3] |R >)
!!
!!    sigma_1 += c_mu3 < mu3| [[H,T_2],tau_nu1] |R >
!!    sigma_2 += c_mu3 < mu3| [H,tau_ nu2] |R >
!!
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: construct_contravariant_t3
      use array_utilities, only: zero_array
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
      logical, intent(in)  :: cvs, rm_core
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Arrays for intermediates cannot hold the whole Y_ebck array
      real(dp), dimension(:,:,:,:), allocatable, target :: Y_ebci
      real(dp), dimension(:,:,:,:), allocatable, target :: Y_ebcj
      real(dp), dimension(:,:,:,:), allocatable, target :: Y_ebck
      real(dp), dimension(:,:,:,:), contiguous, pointer :: Y_ebci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: Y_ebcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: Y_ebck_p => null()
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer              :: req_single_batch
!
      type(timings) :: timer
!
      timer = timings('Time in CC3 C3', pl='verbose')
      call timer%turn_on()
!
!     Set up arrays for amplitudes
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + wf%n_v**3 + wf%n_v*wf%n_o**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + 3*wf%n_v**3*wf%n_o &
                       + 2*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 3*wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 4*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,       &
                           req_0, req_i, req_1, req_1,      &
                           req_2, req_2, req_2, req_3,      &
                           'jacobian_transpose_cc3_c3_a',   &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(Y_cmjk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(Y_cmjk, wf%n_v*wf%n_o**3)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(Y_ebci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(Y_ebci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(Y_ebcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(Y_ebck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      wf%Y_ebck = direct_stream_file('Y_ebck', wf%n_v**3)
      call wf%Y_ebck%open_()
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(g_bdci, g_bdci_p, sorting, batch_i, left=.true.)
!
         call wf%setup_vvov(g_dbic, g_dbic_p, sorting, batch_i, left=.true.)
!
         call zero_array(Y_ebci, batch_i%length*wf%n_v**3)
         Y_ebci_p => Y_ebci
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            call wf%setup_ooov(g_jlic, g_jlic_p, sorting, batch_j, batch_i)
!
            call wf%setup_ovov(g_ibjc, g_ibjc_p, sorting, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(g_bdcj, g_bdcj_p, sorting, batch_j, left=.true.)
!
               call wf%setup_vvov(g_dbjc, g_dbjc_p, sorting, batch_j, left=.true.)
!
               call wf%Y_ebck%read_range(Y_ebcj, batch_j)
               Y_ebcj_p => Y_ebcj
!
               call wf%setup_oovo(g_licj, g_licj_p, sorting, batch_i, batch_j)
!
               call wf%setup_ooov(g_iljc, g_iljc_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
!
               call wf%point_vvvo(g_dbjc_p, g_dbic, batch_j%length)
!
               Y_ebcj_p => Y_ebci
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
!
               call wf%point_vooo(g_iljc_p, g_jlic, batch_i%length, batch_j%length)
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%setup_vvvo(g_bdck, g_bdck_p, sorting, batch_k, left=.true.)
!
                  call wf%setup_vvov(g_dbkc, g_dbkc_p, sorting, batch_k, left=.true.)
!
                  call wf%Y_ebck%read_range(Y_ebck, batch_k)
                  Y_ebck_p => Y_ebck
!
                  call wf%setup_oovo(g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ooov(g_ilkc, g_ilkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ooov(g_jlkc, g_jlkc_p, sorting, batch_j, batch_k)
                  call wf%setup_ooov(g_klic, g_klic_p, sorting, batch_k, batch_i)
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ovov(g_ibkc, g_ibkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ovov(g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbic, batch_k%length)
!
                  Y_ebck_p => Y_ebci
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_ilkc_p, g_jlic, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_jlic, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_kljc_p, g_jlic, batch_k%length, batch_j%length)
!
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
                  call wf%point_vvoo(g_jbkc_p, g_ibjc, batch_j%length, batch_k%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbjc, batch_k%length)
!
                  Y_ebck_p => Y_ebcj
!
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_ilkc_p, g_iljc, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_kljc, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
!
                  call wf%setup_ovov(g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
!
               endif
!
               do i = batch_i%first, batch_i%get_last()
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%get_last(), i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%get_last(), j)
!
!                       Check for core orbitals
                        if (wf%ijk_amplitudes_are_zero(i, j, k, cvs, rm_core)) cycle
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct covariant V_abc for given i,j,k
!                       L_abc is obtained by the linear combination
!                       4V_abc - 2V_bac - 2V_cba - 2V_acb + V_bca + V_cab
!                       and divide by omega - eps^abc_ijk
!
                        call wf%construct_W(i, j, k, sorting, u_abc, c_abij, &
                                            g_dbic_p(:,:,:,i_rel),     &
                                            g_dbjc_p(:,:,:,j_rel),     &
                                            g_dbkc_p(:,:,:,k_rel),     &
                                            g_jlic_p(:,:,j_rel,i_rel), &
                                            g_klic_p(:,:,k_rel,i_rel), &
                                            g_kljc_p(:,:,k_rel,j_rel), &
                                            g_iljc_p(:,:,i_rel,j_rel), &
                                            g_ilkc_p(:,:,i_rel,k_rel), &
                                            g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%outer_product_terms_l3(i, j, k, c_ai, c_abij,     &
                                                       u_abc, wf%fock_ia,         &
                                                       g_ibjc_p(:,:,i_rel,j_rel), &
                                                       g_ibkc_p(:,:,i_rel,k_rel), &
                                                       g_jbkc_p(:,:,j_rel,k_rel))
!
                        call construct_contravariant_t3(u_abc, sorting, wf%n_v)
!
                        call wf%divide_by_orbital_differences(i, j, k, u_abc, omega)
!
                        call wf%jacobian_transpose_cc3_contractions(i, j, k, t_abij,           &
                                                                    u_abc, sorting,            &
                                                                    sigma_abij, Y_cmjk,        &
                                                                    Y_ebci_p(:,:,:,i_rel),     &
                                                                    Y_ebcj_p(:,:,:,j_rel),     &
                                                                    Y_ebck_p(:,:,:,k_rel),     &
                                                                    g_bdci_p(:,:,:,i_rel),     &
                                                                    g_bdcj_p(:,:,:,j_rel),     &
                                                                    g_bdck_p(:,:,:,k_rel),     &
                                                                    g_ljci_p(:,:,j_rel,i_rel), &
                                                                    g_lkci_p(:,:,k_rel,i_rel), &
                                                                    g_lkcj_p(:,:,k_rel,j_rel), &
                                                                    g_licj_p(:,:,i_rel,j_rel), &
                                                                    g_lick_p(:,:,i_rel,k_rel), &
                                                                    g_ljck_p(:,:,j_rel,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
!              write the intermediate Y_ebck to file.
!              Will be read in after the loops for the contractions to sigma_ai
!
               if (k_batch .ne. j_batch) then !k_batch != j_batch, k_batch != i_batch
                  call wf%Y_ebck%write_range(Y_ebck, batch_k)
               endif
!
            enddo ! batch_k
!
            if (j_batch .ne. i_batch) then
               call wf%Y_ebck%write_range(Y_ebcj, batch_j)
            endif
!
         enddo ! batch_j
!
         call wf%Y_ebck%write_range(Y_ebci, batch_i)
!
      enddo ! batch_i
!
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(Y_ebci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(Y_ebci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(Y_ebcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(Y_ebck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%batch_finalize()
!
!     Contribution of the Y_cmjk to sigma1
!
      call wf%jacobian_transpose_cc3_c3_a1_y_o(sigma_ai, Y_cmjk)
!
      call mem%dealloc(Y_cmjk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     Contribution of the Y_ebck to sigma1
!
      call wf%jacobian_transpose_cc3_c3_a1_y_v(sigma_ai)
!
      call wf%Y_ebck%close_()
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_cc3_c3_a_cc3
!
!
   module subroutine outer_product_terms_L3_cc3(wf, i, j, k, L1, L2, L_abc, &
                                                F_ov, g_ibjc, g_ibkc, g_jbkc)
!!
!!    Direct product terms L3
!!    Written by Alexander C. Paul and Rolf H. Myhre, Sep 2020
!!
!!    Contributions from outer products to the covariant L3:
!!
!!    P^abc_ijk (L^a_i g_jbkc + L^b_j g_iakc + L^c_k g_iajb
!!              + L^ab_ij F_kc + L^ac_ik F_jb + L^bc_jk F_ia)
!!
!!    computed in unrolled loops which is faster than 6 dger.
!!
!!    NB: The covariant L2 = 1/3 (2L^ab_ij + L^ba_ij) is needed.
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: L1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: L2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: L_abc
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: F_ov
!
!     g_ibjc ordered bc,ij
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_jbkc
!
      integer :: a, b, c
!
!$omp parallel do private(a,b,c) collapse(2)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               L_abc(a,b,c) = L_abc(a,b,c)          &
                            + L1(a,i)*g_jbkc(b,c)   &
                            + L1(b,j)*g_ibkc(a,c)   &
                            + L1(c,k)*g_ibjc(a,b)   &
                            + L2(a,b,i,j)*F_ov(k,c) &
                            + L2(a,c,i,k)*F_ov(j,b) &
                            + L2(b,c,j,k)*F_ov(i,a)
!
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine outer_product_terms_L3_cc3
!
!
   module subroutine jacobian_transpose_cc3_contractions_cc3(wf, i, j, k, t2, u_abc, v_abc, &
                                                             sigma2,  Y_cmjk,               &
                                                             Y_ebci, Y_ebcj, Y_ebck,        &
                                                             g_bdci, g_bdcj, g_bdck,        &
                                                             g_ljci, g_lkci, g_lkcj,        &
                                                             g_licj, g_lick, g_ljck)
!!
!!    Jacobian transpose CC3 contractions
!!    Written by Alexander C. Paul and Rolf H. Myhre, Okt 2020
!!
!!    sigma_adij =   sum_ckd c^abc_ijk g_bdck
!!    sigma_abil = - sum_cki c^bac_ijk g_lick
!!    Y_ebck = sum_aij c^abc_ijk * t^ae_ij
!!    Y_cmjk = sum_abj c^abc_ijk * t^ab_mj
!!
!!    All permutations for i,j,k have to be considered
!!    due to the restrictions in the i,j,k loops
!!
      use reordering, only: sort_123_to_312, sort_123_to_213
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)       :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)       :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)  :: t2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out) :: sigma2
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out) :: Y_cmjk
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)         :: Y_ebci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)         :: Y_ebcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)         :: Y_ebck
!
!     g_bdck ordered as bcd,k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)          :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)          :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)          :: g_bdck
!
!     g_ljck ordered as cl,jk
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                  :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                  :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                  :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                  :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                  :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                  :: g_ljck
!
!     sigma_adki += sum_bc,j u_bca g_bdcj
!     sigma_ablk -= sum_c u_cab g_ljci
!
      call wf%omega2_cc3_permutation(k, i, u_abc, g_bdcj, g_ljci, sigma2)
!
!     Y_ebck += sum_a,ij u_abc t^ae_ij
      call wf%construct_Y_vvvo_permutation(u_abc, t2(:,:,j,i), Y_ebck, one)
!
!     Y_cmik += sum_ab,j u_abc t^ab_mj
      call wf%construct_Y_vooo_permutation(i, k, u_abc, t2(:,:,:,j), Y_cmjk, one)
!
!     cab -> bca
!     abc -> cab
!     bca -> abc
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     sigma_adjk += sum_bc,i v_bca g_bdci
!     sigma_ablj -= sum_c v_cab g_lick
!
      call wf%omega2_cc3_permutation(j, k, v_abc, g_bdci, g_lick, sigma2)
!
!     Y_ebcj += sum_a,ki v^bca t^ae_ki
      call wf%construct_Y_vvvo_permutation(v_abc, t2(:,:,i,k), Y_ebcj, one)
!
!     Y_cmkj += sum_ab,i v^abc t^ab_mi
      call wf%construct_Y_vooo_permutation(k, j, v_abc, t2(:,:,:,i), Y_cmjk, one)
!
!     abc -> bca
!     bca -> cab
!     cab -> abc
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     sigma_adij += sum_bc,k u_bca g_bdck
!     sigma_abli -= sum_c u_cab g_lkcj
!
      call wf%omega2_cc3_permutation(i, j, u_abc, g_bdck, g_lkcj, sigma2)
!
!     Y_ebci += sum_a,jk u_abc t^ea_kj
      call wf%construct_Y_vvvo_permutation(u_abc, t2(:,:,k,j), Y_ebci, one)
!
!     Y_cmji += sum_ab,k u_abc t^ab_mk
      call wf%construct_Y_vooo_permutation(j, i, u_abc, t2(:,:,:,k), Y_cmjk, one)
!
      if (i .ne. j .and. j .ne. k) then
!
!        acb -> bca
!        bac -> cab
!        cba -> abc
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        sigma_adik += sum_bc v_bca g_bdcj
!        sigma_abli -= sum_c v_cab g_ljck
!
         call wf%omega2_cc3_permutation(i, k, v_abc, g_bdcj, g_ljck, sigma2)
!
!        Y_ebci = sum_a,kj v_abc t^ea_jk
         call wf%construct_Y_vvvo_permutation(v_abc, t2(:,:,j,k), Y_ebci, one)
!
!        Y_cmki = sum_ab,j v_abc t^ab_mj
         call wf%construct_Y_vooo_permutation(k, i, v_abc, t2(:,:,:,j), Y_cmjk, one)
!
!        bac -> bca
!        cba -> cab
!        acb -> abc
         call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        sigma_adji += sum_bc,k u_bca g_bdck
!        sigma_ablj += - sum_c u^cab g_lkci
!
         call wf%omega2_cc3_permutation(j, i, u_abc, g_bdck, g_lkci, sigma2)
!
!        Y_ebcj = sum_a,ik u_abc t^ea_ki
         call wf%construct_Y_vvvo_permutation(u_abc, t2(:,:,k,i), Y_ebcj, one)
!
!        Y_cmij = sum_ab,k u_abc t^ab_mk
         call wf%construct_Y_vooo_permutation(i, j, u_abc, t2(:,:,:,k), Y_cmjk, one)
!
!        cba -> bca
!        acb -> cab
!        bac -> abc
         call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        sigma_adkj += sum_bc v_bca g_bdci
!        sigma_ablk += - sum_c v_cab g_licj
!
         call wf%omega2_cc3_permutation(k, j, v_abc, g_bdci, g_licj, sigma2)
!
!        Y_ebck = sum_a,ji v_abc t^ea_ij
         call wf%construct_Y_vvvo_permutation(v_abc, t2(:,:,i,j), Y_ebck, one)
!
!        Y_cmjk = sum_ab,i v_abc t^ab_mi
         call wf%construct_Y_vooo_permutation(j, k, v_abc, t2(:,:,:,i), Y_cmjk, one)
!
      end if
!
   end subroutine jacobian_transpose_cc3_contractions_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a1_y_o_cc3(wf, sigma_ai, Y_cmjk)
!!
!!    Jacobian tranpose contribution of Y_vooo
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_1 += sum_mjk Y_cmjk * g_mjlk
!!    sigma_1 += sum_cmj Y_cmjk * g_mjcd
!!    sigma_1 += sum_cmk Y_cmjk * g_mdck
!!
      use reordering, only: sort_1234_to_2314, sort_123_to_132
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: Y_cmjk
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_mjck
!
      real(dp), dimension(:,:,:), allocatable :: X_J_ck
      real(dp), dimension(:,:,:), allocatable :: X_Jk_c
!
      real(dp), dimension(:,:,:), allocatable :: X_J_lk
      real(dp), dimension(:,:,:), allocatable :: X_J_kl
!
      real(dp), dimension(:,:,:), allocatable :: L_J_vv
!
      type(batching_index) :: batch_k, batch_v
      integer :: k_batch, v_batch
      integer :: req_0, req_k, req_v, req_2
!
      batch_k = batching_index(wf%n_o)
      batch_v = batching_index(wf%n_v)
!
      req_0 = 0
      req_k = wf%n_v*wf%n_o**2 + (2*wf%n_v + wf%n_o)*wf%eri_t1%n_J
      req_v = max(wf%n_v,wf%n_o)*wf%eri_t1%n_J
      req_2 = 0
!
      call mem%alloc(X_J_lk, wf%eri_t1%n_J, wf%n_o, wf%n_o)
      call mem%alloc(X_J_kl, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
      call mem%batch_setup(batch_k, batch_v, req_0, req_k, req_v, req_2, &
                           'jacobian_transpose_cc3_c3_a1_y_o')
!
      call mem%alloc(Y_mjck, wf%n_v, batch_k%max_length, wf%n_o, wf%n_o)
!
      call mem%alloc(X_J_ck, wf%eri_t1%n_J, wf%n_v, batch_k%max_length)
      call mem%alloc(X_Jk_c, wf%eri_t1%n_J, batch_k%max_length, wf%n_v)
!
      call mem%alloc(L_J_vv, wf%eri_t1%n_J, max(wf%n_v, wf%n_o), batch_v%max_length)
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call sort_1234_to_2314(Y_cmjk(:,:,:,batch_k%first:), Y_mjck, &
                                wf%n_v, wf%n_o, wf%n_o, batch_k%length)
!
         call wf%L_t1%get(X_J_lk, 1, wf%n_o, 1, wf%n_o)
!
         call dgemm('N', 'N',              &
                    wf%eri_t1%n_J,         &
                    wf%n_v*batch_k%length, &
                    wf%n_o**2,             &
                    one,                   &
                    X_J_lk,                & ! L_J_mj
                    wf%eri_t1%n_J,         &
                    Y_mjck,                & ! Y_mj_c#k
                    wf%n_o**2,             &
                    zero,                  &
                    X_J_ck,                & ! X_J_c#k
                    wf%eri_t1%n_J)
!
!        :: Term 1: sigma_dk -= sum_cmj Y_cmjk * g_cdmj  ::
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
            call wf%L_t1%get(L_J_vv, wf%n_o+1, wf%n_o+wf%n_v, &
                                        wf%n_o+batch_v%first, wf%n_o+batch_v%get_last())
!
            call dgemm('T','N',                                 &
                        batch_v%length,                         &
                        batch_k%length,                         &
                        wf%n_v*wf%eri_t1%n_J,                   &
                        -one,                                   &
                        L_J_vv,                                 & ! L_Jc_#d
                        wf%n_v*wf%eri_t1%n_J,                   &
                        X_J_ck,                                 & ! X_Jc_#k
                        wf%n_v*wf%eri_t1%n_J,                   &
                        one,                                    &
                        sigma_ai(batch_v%first:,batch_k%first), & ! sigma_#d_#k
                        wf%n_v)
!
         enddo
!
!        :: Term 2: sigma_cl += sum_mjk Y_cmjk * g_mjlk ::
!
         call wf%L_t1%get(X_J_lk, 1, wf%n_o, batch_k%first, batch_k%get_last())
!
         call sort_123_to_132(X_J_lk, X_J_kl, wf%eri_t1%n_J, wf%n_o, batch_k%length)
         call sort_123_to_132(X_J_ck, X_Jk_c, wf%eri_t1%n_J, wf%n_v, batch_k%length)
!
         call dgemm('T', 'N',                      &
                    wf%n_v,                        &
                    wf%n_o,                        &
                    wf%eri_t1%n_J*batch_k%length,  &
                    one,                           &
                    X_Jk_c,                        & ! X_J#k_c
                    wf%eri_t1%n_J*batch_k%length,  &
                    X_J_kl,                        & ! L_J#k_l
                    wf%eri_t1%n_J*batch_k%length,  &
                    one,                           &
                    sigma_ai,                      & ! sigma_c_l
                    wf%n_v)
!
!        :: Term 3: sigma_dj -= sum_cmj Y_cmjk * g_ckmd ::
!
         call wf%L_t1%get(X_J_ck, wf%n_o+1, wf%n_o+wf%n_v, &
                                     batch_k%first, batch_k%get_last())
!
         call dgemm('N', 'T',              &
                    wf%eri_t1%n_J,         &
                    wf%n_o**2,             &
                    wf%n_v*batch_k%length, &
                    one,                   &
                    X_J_ck,                & ! L_J_c#k
                    wf%eri_t1%n_J,         &
                    Y_mjck,                & ! Y_mj_c#k
                    wf%n_o**2,             &
                    zero,                  &
                    X_J_lk,                & ! X_J_mj
                    wf%eri_t1%n_J)
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
            call wf%L_t1%get(L_J_vv, 1, wf%n_o, &
                                        wf%n_o+batch_v%first, wf%n_o+batch_v%get_last())
!
            call dgemm('T','N',                     &
                        batch_v%length,             &
                        wf%n_o,                     &
                        wf%n_o*wf%eri_t1%n_J,       &
                        -one,                       &
                        L_J_vv,                     & ! L_Jm_#d
                        wf%n_o*wf%eri_t1%n_J,       &
                        X_J_lk,                     & ! X_Jm_j
                        wf%n_o*wf%eri_t1%n_J,       &
                        one,                        &
                        sigma_ai(batch_v%first:,1), & ! sigma_#dj
                        wf%n_v)
!
         enddo
      enddo
!
      call mem%dealloc(L_J_vv, wf%eri_t1%n_J, max(wf%n_v,wf%n_o), batch_v%max_length)
!
      call mem%dealloc(X_Jk_c, wf%eri_t1%n_J, batch_k%max_length, wf%n_v)
      call mem%dealloc(X_J_ck, wf%eri_t1%n_J, wf%n_v, batch_k%max_length)
!
      call mem%dealloc(Y_mjck, wf%n_v, batch_k%max_length, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_J_kl, wf%eri_t1%n_J, wf%n_o, wf%n_o)
      call mem%dealloc(X_J_lk, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
      call mem%batch_finalize()
!
   end subroutine jacobian_transpose_cc3_c3_a1_y_o_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a1_y_v_cc3(wf, sigma_ai)
!!
!!    Jacobian transpose contribution Y_vvvo
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_dk += sum_bec Y_ebck * g_becd
!!    sigma_cl += sum_bek Y_ebck * g_lkbe
!!    sigma_bl += sum_cek Y_ebck * g_leck
!!
      use array_utilities, only: zero_array
      use reordering, only: sort_1234_to_2134, sort_123_to_132, sort_123_to_312
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_ebck
      real(dp), dimension(:,:,:,:), allocatable :: Y_beck
!
      real(dp), dimension(:,:,:), allocatable :: L_J_lk
      real(dp), dimension(:,:,:), allocatable :: L_Jk_l
!
      real(dp), dimension(:,:,:), allocatable :: L_J_le
      real(dp), dimension(:,:,:), allocatable :: L_eJ_l
!
      real(dp), dimension(:,:,:), allocatable :: X_J_ck
      real(dp), dimension(:,:,:), allocatable :: X_J_kc
!
      real(dp), dimension(:,:,:), allocatable :: X_J_vv
!
      type(batching_index) :: batch_k, batch_v
      integer :: k_batch, v_batch
      integer :: req_0, req_k, req_v, req_2
!
      batch_k = batching_index(wf%n_o)
      batch_v = batching_index(wf%n_v)
!
      req_0 = 0
      req_k = 2*wf%n_v**3 + 2*wf%n_v*wf%eri_t1%n_J + 2*wf%n_o*wf%eri_t1%n_J
      req_v = wf%n_v*wf%eri_t1%n_J + 2*wf%n_o*wf%eri_t1%n_J
      req_2 = 0
!
      call mem%batch_setup(batch_k, batch_v, req_0, req_k, req_v, req_2, &
                           'jacobian_transpose_cc3_c3_a1_y_v')
!
      call mem%alloc(Y_ebck, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
      call mem%alloc(Y_beck, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
!
      call mem%alloc(L_Jk_l, wf%eri_t1%n_J, batch_k%max_length, wf%n_o)
      call mem%alloc(L_J_lk, wf%eri_t1%n_J, wf%n_o, batch_k%max_length)
!
      call mem%alloc(L_J_le, wf%eri_t1%n_J, wf%n_o, batch_v%max_length)
      call mem%alloc(L_eJ_l, batch_v%max_length, wf%eri_t1%n_J, wf%n_o)
!
      call mem%alloc(X_J_ck, wf%eri_t1%n_J, wf%n_v, batch_k%max_length)
      call mem%alloc(X_J_kc, wf%eri_t1%n_J, batch_k%max_length, wf%n_v)
!
      call mem%alloc(X_J_vv, wf%eri_t1%n_J, wf%n_v, batch_v%max_length)
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call wf%Y_ebck%read_range(Y_ebck, batch_k)
         call sort_1234_to_2134(Y_ebck, Y_beck, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call zero_array(X_J_ck, wf%eri_t1%n_J*wf%n_v*batch_k%length)
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
            call wf%L_t1%get(X_J_vv, wf%n_o+1, wf%n_o+wf%n_v, &
                                        wf%n_o+batch_v%first, wf%n_o+batch_v%get_last())
!
            call dgemm('N','N',                      &
                        wf%eri_t1%n_J,               &
                        wf%n_v*batch_k%length,       &
                        wf%n_v*batch_v%length,       &
                        one,                         &
                        X_J_vv,                      & ! L_J_b#e
                        wf%eri_t1%n_J,               &
                        Y_beck(:,batch_v%first,1,1), & ! Y_b#e_c#k
                        wf%n_v**2,                   &
                        one,                         &
                        X_J_ck,                      & ! X_J_c#k
                        wf%eri_t1%n_J)
!
         enddo ! v_batch
!
!        :: Term 1: sigma_dk += sum_bec Y_ebck * g_becd ::
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
            call wf%L_t1%get(X_J_vv, wf%n_o+1, wf%n_o+wf%n_v, &
                                        wf%n_o+batch_v%first, wf%n_o+batch_v%get_last())
!
            call dgemm('T','N',                                 &
                        batch_v%length,                         &
                        batch_k%length,                         &
                        wf%n_v*wf%eri_t1%n_J,                   &
                        one,                                    &
                        X_J_vv,                                 & ! L_Jc_#d
                        wf%n_v*wf%eri_t1%n_J,                   &
                        X_J_ck,                                 & ! X_Jc_#k
                        wf%n_v*wf%eri_t1%n_J,                   &
                        one,                                    &
                        sigma_ai(batch_v%first:,batch_k%first), & ! sigma_#d#k
                        wf%n_v)
!
         enddo ! v_batch
!
!        :: Term 2: sigma_cl -= sum_bek Y_ebck * g_belk ::
!
         call sort_123_to_132(X_J_ck, X_J_kc, wf%eri_t1%n_J, wf%n_v, batch_k%length)
!
         call wf%L_t1%get(L_J_lk, 1, wf%n_o, batch_k%first, batch_k%get_last())
         call sort_123_to_132(L_J_lk, L_Jk_l, wf%eri_t1%n_J, wf%n_o, batch_k%length)
!
         call dgemm('T', 'N',                      &
                    wf%n_v, wf%n_o,                &
                    wf%eri_t1%n_J*batch_k%length,  &
                    -one,                          &
                    X_J_kc,                        & !X_J#k_c
                    wf%eri_t1%n_J*batch_k%length,  &
                    L_Jk_l,                        & !L_J#k_l
                    wf%eri_t1%n_J*batch_k%length,  &
                    one,                           &
                    sigma_ai,                      & !sigma_cl
                    wf%n_v)
!
!
!        :: Term 3: sigma_bl += - sum_cek Y_ebck * g_ckle ::
!
         call wf%L_t1%get(X_J_ck, 1+wf%n_o, wf%n_v+wf%n_o, &
                                     batch_k%first, batch_k%get_last())
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
            call dgemm('N','T',                      &
                        wf%n_v*batch_v%length,       &
                        wf%eri_t1%n_J,               &
                        wf%n_v*batch_k%length,       &
                        one,                         &
                        Y_beck(:,batch_v%first,1,1), & ! Y_b#e_c#k
                        wf%n_v**2,                   &
                        X_J_ck,                      & ! L_J_c#k
                        wf%eri_t1%n_J,               &
                        zero,                        &
                        X_J_vv,                      & ! X_b#e_J
                        wf%n_v*batch_v%length)
!
            call wf%L_t1%get(L_J_le, 1, wf%n_o, &
                                        batch_v%first+wf%n_o, batch_v%get_last()+wf%n_o)
!
            call sort_123_to_312(L_J_le, L_eJ_l, wf%eri_t1%n_J, wf%n_o, batch_v%length)
!
            call dgemm('N','N',                       &
                        wf%n_v, wf%n_o,               &
                        wf%eri_t1%n_J*batch_v%length, &
                        -one,                         &
                        X_J_vv,                       & ! X_b_#eJ
                        wf%n_v,                       &
                        L_eJ_l,                       & ! L_#eJ_l
                        wf%eri_t1%n_J*batch_v%length, &
                        one,                          &
                        sigma_ai,                     & ! sigma_bl
                        wf%n_v)
!
         enddo
!
      enddo ! k_batch
!
      call mem%dealloc(X_J_vv, wf%eri_t1%n_J, wf%n_v, batch_v%max_length)
!
      call mem%dealloc(X_J_kc, wf%eri_t1%n_J, batch_k%max_length, wf%n_v)
      call mem%dealloc(X_J_ck, wf%eri_t1%n_J, wf%n_v, batch_k%max_length)
!
      call mem%dealloc(L_eJ_l, wf%eri_t1%n_J, batch_v%max_length, wf%n_o)
      call mem%dealloc(L_J_le, wf%eri_t1%n_J, wf%n_o, batch_v%max_length)
!
      call mem%dealloc(L_J_lk, wf%eri_t1%n_J, wf%n_o, batch_k%max_length)
      call mem%dealloc(L_Jk_l, wf%eri_t1%n_J, batch_k%max_length, wf%n_o)
!
      call mem%dealloc(Y_beck, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
      call mem%dealloc(Y_ebck, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
!
      call mem%batch_finalize()
!
   end subroutine jacobian_transpose_cc3_c3_a1_y_v_cc3
!
!
   module subroutine construct_Y_vvvo_permutation_cc3(wf, L3, t2, Y_vvvo, factor)
!!
!!    Construct Y_vvvo permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to the Y_vvvo intermediates of the type:
!!    Y_ebck += sum_a,ij u_abc t^ae_ij
!!    X_abid -= sum_c,jk u_abc g_kcjd
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: L3
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: t2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: Y_vvvo
!
      real(dp) :: factor
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_v**2,   &
                  wf%n_v,      &
                  factor,      &
                  t2(:,:),     &
                  wf%n_v,      &
                  L3,          &
                  wf%n_v,      &
                  one,         &
                  Y_vvvo,      &
                  wf%n_v)
!
   end subroutine construct_Y_vvvo_permutation_cc3
!
!
   module subroutine construct_Y_vooo_permutation_cc3(wf, o1, o2, L3, t2, &
                                                      Y_vooo, factor)
!!
!!    Construct Y_vooo permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to the Y_vooo intermediates of the type:
!!    Y_cmik += sum_ab,j u_abc t^ab_mj
!!    X_alji -= sum_bc,k u_abc g_lbkc
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: o1, o2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: L3
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in) :: t2
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out) :: Y_vooo
!
      real(dp) :: factor
!
      call dgemm('T','N',            &
                  wf%n_v,            &
                  wf%n_o,            &
                  wf%n_v**2,         &
                  factor,            &
                  L3,                &
                  wf%n_v**2,         &
                  t2(:,:,:),         &
                  wf%n_v**2,         &
                  one,               &
                  Y_vooo(:,:,o1,o2), &
                  wf%n_v)
!
   end subroutine construct_Y_vooo_permutation_cc3
!
!
   module subroutine construct_X_vo_permutation_cc3(wf, o1, t3, L2, X_vo, factor)
!!
!!    Construct X_vo permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to CC3 omega1 vector, e.g.:
!!    X_ai += sum_bc u_abc L^bc_jk
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: o1
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: t3
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: X_vo
!
      real(dp), intent(in) :: factor
!
      call dgemv('N',         &
                 wf%n_v,      &
                 wf%n_v**2,   &
                 factor,      &
                 t3,          & ! u_a_bc
                 wf%n_v,      &
                 L2, 1,       & ! L_bc_oo
                 one,         &
                 X_vo(:,o1), 1)
!
   end subroutine construct_X_vo_permutation_cc3
!
!
end submodule jacobian_transpose
