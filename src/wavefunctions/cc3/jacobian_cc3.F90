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
submodule (cc3_class) jacobian
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    rho_i = A * c_i,
!!
!!    where
!!
!!    A_mu,nu = < mu|exp(-T) [H, tau_nu] exp(T) |R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine effective_jacobian_transformation_cc3(wf, omega, c, rho)
!!
!!    Effective Jacobian transformation (CC3)
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!    Adapted to construct a contravariant representation of R2 for CC3 which
!!    is in the end transformed back to its covariant form.
!!    by Rolf H. Myhre and Alexander C. Paul, Sep 2020
!!
!!    Directs the transformation by the CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu|exp(-T) [H, tau_nu] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu
!!       = sum_ck A_mu,ck c_ck + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl)
!!
!!    On exit, c is overwritten by rho. That is, c(ai) = rho_a_i,
!!    and c(aibj) = rho_aibj.
!!
      use array_utilities, only: scale_diagonal, zero_array
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: symmetrize_add_contra_to_packed
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: c
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: rho
!
      real(dp), dimension(:,:,:,:), allocatable :: c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_abij
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CC3 transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(rho, wf%n_t1 + wf%n_t2)
!
      call wf%ccsd%jacobian_transformation(c, rho)
!
!     :: Compute CC3 contributions to rho ::
!
      call wf%construct_c1_cholesky(c(1:wf%n_t1), wf%L_t1, wf%L_c1)
!
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(rho_abij, (wf%n_v*wf%n_o)**2)
!
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(c(wf%n_t1+1:), c_abij, wf%n_v, wf%n_o, &
                                          wf%n_v, wf%n_o)
      call scale_diagonal(two, C_abij, wf%n_v, wf%n_o)
!
!     rho_abij is obtained in contravariant form (2 rho^ab_ij - rho^ba_ij)
!     and has to be converted back to the covariant form when packing in.
!
!     CC3-Contributions from the T3-amplitudes
      call wf%jacobian_cc3_t3_a2(c(1:wf%n_t1), rho_abij)
      call wf%jacobian_cc3_t3_b2(rho_abij, wf%cvs, wf%rm_core)
!
!     CC3-Contributions from the C3-amplitudes
      call wf%jacobian_cc3_c3_a(omega, c_abij, rho, rho_abij, &
                                wf%cvs, wf%rm_core)
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call scale_diagonal(half, rho_abij, wf%n_v, wf%n_o)
      call symmetrize_add_contra_to_packed(rho_abij, rho(wf%n_t1+1:wf%n_es_amplitudes), &
                                           wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transformation_cc3
!
!
   module subroutine jacobian_cc3_t3_a2_cc3(wf, c_ai, rho_abij)
!!
!!    Jacobian CC3 T3 A2
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!    Adapted to give a contravariant representation of rho2 due to the
!!    way the X-intermediates are constructed
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    Reads in the intermediates X_abid and X_ajil prepared in
!!    prepare_jacobian_transform contracts with c_ai and adds to rho_abij
!!
!!    ~rho_abil += sum_abi X_abid * C_dl
!!    ~rho_daji += sum_aik C_dl * X_ajil
!!
!!    where: ~rho^ab_ij = 2 rho^ab_ij - rho^ba_ij
!!           X_abid = - sum_jck u_abc * g_kcjd
!!           X_ajil = - sum_bck u_abc * g_lbkc
!!           u^abc_ijk = 4t^abc_ijk + t_bca_ijk + t_cab_ijk
!!                     - 2t^acb_ijk - 2t_cba_ijk - 2t_bac_ijk
!!
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out) :: rho_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
      real(dp), dimension(:,:,:,:), allocatable :: X_ajil
!
      type(batching_index) :: batch_d
      integer :: d_batch
      integer :: req_0, req_d
!
      type(timings), allocatable :: cc3_timer_t3_a2
!
      cc3_timer_t3_a2 = timings('Time in CC3 T3 a2', pl='verbose')
      call cc3_timer_t3_a2%turn_on()
!
!     :: X_abid term ::
!
      batch_d = batching_index(wf%n_v)
!
      req_0 = 0
      req_d = wf%n_o * wf%n_v**2
!
      call mem%batch_setup(batch_d, req_0, req_d, 'jacobian_cc3_t3_a2')
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
         call dgemm('N','N',                 &
                     wf%n_o*wf%n_v**2,       &
                     wf%n_o,                 &
                     batch_d%length,         &
                     one,                    &
                     X_abid,                 & ! X_abi_d
                     wf%n_o*wf%n_v**2,       &
                     c_ai(batch_d%first, 1), & ! c_d_l
                     wf%n_v,                 &
                     one,                    &
                     rho_abij,               & ! rho_abil
                     wf%n_o*wf%n_v**2)
!
      enddo
!
      call mem%dealloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%max_length)
!
      call mem%batch_finalize()
!
      call wf%X_abid%close_()
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
      call dgemm('N','T',           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  wf%n_o,           &
                  one,              &
                  c_ai,             & ! C_d_l
                  wf%n_v,           &
                  X_ajil,           & ! X_aji_l
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  rho_abij,         & ! rho_daji
                  wf%n_v)
!
      call mem%dealloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call cc3_timer_t3_a2%turn_off()
!
   end subroutine jacobian_cc3_t3_a2_cc3
!
!
   module subroutine jacobian_cc3_t3_b2_cc3(wf, rho_abij, cvs, rm_core)
!!
!!    Jacobian CC3 T3 B2
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!    Adapted to give a contravariant representation of rho2 due to the
!!    use of a contravariant representation of t3
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    ~rho_abij +=  sum_ckdl C^d_l L_kcld (t^abc_ijk - t^bac_ijk)
!!              +=  sum_ck F_kc_c1 * (t^abc_ijk - t^bac_ijk)
!!
!!    where: ~rho^ab_ij = 2 rho^ab_ij - rho^ba_ij
!!
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3) :: wf
!
      logical, intent(in) :: cvs, rm_core
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
      real(dp), dimension(:,:,:), allocatable :: t_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck_c1
!
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
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_i, req_1, req_2, req_3, req_1_eri
      integer              :: req_single_batch
!
      type(timings), allocatable :: cc3_timer_t3_b2
!
      cc3_timer_t3_b2 = timings('Time in CC3 T3 b2', pl='verbose')
      call cc3_timer_t3_b2%turn_on()
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     C1 transformed Fock matrix
!
      call mem%alloc(F_ov_ck_c1, wf%n_v, wf%n_o)
      call wf%construct_c1_fock(F_ov_ck_c1)
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
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           'jacobian_cc3_t3_b2',       &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
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
!     Loop over the batches in i,j,k
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(wf%eri_t1, g_bdci, g_bdci_p, sorting, batch_i)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(wf%eri_t1, g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(wf%eri_t1, g_bdcj, g_bdcj_p, sorting, batch_j)
!
               call wf%setup_oovo(wf%eri_t1, g_licj, g_licj_p, sorting, batch_i, batch_j)
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
                  call wf%setup_vvvo(wf%eri_t1, g_bdck, g_bdck_p, sorting, batch_k)
!
                  call wf%setup_oovo(wf%eri_t1, g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
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
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
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
!                       Check for core orbitals:
!                       cvs: i,j,k cannot all correspond to valence orbitals
!                       rm_core: i,j,k may not contain any core orbital
!
!                       Here t3 contributes to rho2 and can, thus, be restricted as well
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
                        call wf%jacobian_cc3_b2_fock(i, j, k, t_abc, sorting, &
                                                     rho_abij, F_ov_ck_c1)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
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
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(F_ov_ck_c1, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call cc3_timer_t3_b2%turn_off()
!
   end subroutine jacobian_cc3_t3_b2_cc3
!
!
   module subroutine jacobian_cc3_c3_a_cc3(wf, omega, c_abij, rho_ai, rho_abij, &
                                           cvs, rm_core)
!!
!!    Jacobian CC3 C3 A
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!    Adapted to give a contravariant representation of rho2 due to the
!!    use of a contravariant representation of R3
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The triples amplitudes are expressed in terms of doubles amplitudes:
!!    C_3 = (omega - epsilon_mu3)^-1 (< mu3| [H,C_2] | HF >
!!                                  + < mu3| [[H,C_1],T_2] |HF >)
!!
!!    c^abc = (omega - epsilon^abc_ijk)^-1 * P^abc_ijk
!!             (sum_d c^ad_ij g_ckbd - sum_l c^ab_il g_cklj
!!            + sum_d t^ad_ij g'_bdck - sum_l t^ab_il g'_cklj
!!
!!     rho1 += < mu1| [H,C_3] |R >
!!    ~rho2 += < mu2| [H,C_3] |R >
!!
!!    where: ~rho^ab_ij = 2 rho^ab_ij - rho^ba_ij
!!
!!    Based on omega_cc3_a_cc3 written by Rolf H. Myhre
!!
      use reordering, only: squareup_and_sort_1234_to_1324, sort_12_to_21
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      logical, intent(in) :: cvs, rm_core
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
      real(dp), dimension(:,:,:), allocatable :: c_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck
!
!     T1 transformed integrals
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
!     C1 transformed integrals
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_c1_p => null()
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
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer              :: req_single_batch
!
      type(timings), allocatable :: cc3_timer_c3
!
      cc3_timer_c3 = timings('Time in CC3 C3', pl='verbose')
      call cc3_timer_c3%turn_on()
!
      call mem%alloc(F_ov_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_ck, wf%n_o, wf%n_v)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_c1_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + wf%n_v**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + 3*wf%n_v**3*wf%n_o &
                       + 3*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 3*wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 6*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           'jacobian_cc3_c3_a',        &
                           req_single_batch=req_single_batch)
!
!
      call mem%alloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci_c1, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
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
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif
!
!     Loop over the batches in i,j,k
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(wf%eri_t1, g_bdci, g_bdci_p, sorting, batch_i)
         call wf%setup_vvvo(wf%eri_c1, g_bdci_c1, g_bdci_c1_p, sorting, batch_i)
!
         call wf%setup_vvov(g_dbic, g_dbic_p, sorting, batch_i)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(wf%eri_t1, g_ljci, g_ljci_p, sorting, batch_j, batch_i)
            call wf%setup_oovo(wf%eri_c1, g_ljci_c1, g_ljci_c1_p, sorting, batch_j, batch_i)
!
            call wf%setup_ooov(g_jlic, g_jlic_p, sorting, batch_j, batch_i)
!
            call wf%setup_ovov(wf%eri_t1, g_ibjc, g_ibjc_p, sorting, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(wf%eri_t1, g_bdcj, g_bdcj_p, sorting, batch_j)
               call wf%setup_vvvo(wf%eri_c1, g_bdcj_c1, g_bdcj_c1_p, sorting, batch_j)
!
               call wf%setup_vvov(g_dbjc, g_dbjc_p, sorting, batch_j)
!
               call wf%setup_oovo(wf%eri_t1, g_licj, g_licj_p, sorting, batch_i, batch_j)
               call wf%setup_oovo(wf%eri_c1, g_licj_c1, g_licj_c1_p, sorting, batch_i, batch_j)
!
               call wf%setup_ooov(g_iljc, g_iljc_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
               call wf%point_vvvo(g_bdcj_c1_p, g_bdci_c1, batch_j%length)
!
               call wf%point_vvvo(g_dbjc_p, g_dbic, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
               call wf%point_vooo(g_licj_c1_p, g_ljci_c1, batch_i%length, batch_j%length)
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
                  call wf%setup_vvvo(wf%eri_t1, g_bdck, g_bdck_p, sorting, batch_k)
                  call wf%setup_vvvo(wf%eri_c1, g_bdck_c1, g_bdck_c1_p, sorting, batch_k)
!
                  call wf%setup_vvov(g_dbkc, g_dbkc_p, sorting, batch_k)
!
                  call wf%setup_oovo(wf%eri_t1, g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
                  call wf%setup_oovo(wf%eri_c1, g_lick_c1, g_lick_c1_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(wf%eri_c1, g_ljck_c1, g_ljck_c1_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(wf%eri_c1, g_lkci_c1, g_lkci_c1_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(wf%eri_c1, g_lkcj_c1, g_lkcj_c1_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ooov(g_ilkc, g_ilkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ooov(g_jlkc, g_jlkc_p, sorting, batch_j, batch_k)
                  call wf%setup_ooov(g_klic, g_klic_p, sorting, batch_k, batch_i)
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ovov(wf%eri_t1, g_ibkc, g_ibkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ovov(wf%eri_t1, g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
                  call wf%point_vvvo(g_bdck_c1_p, g_bdci_c1, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbic, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_lick_c1_p, g_ljci_c1, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_c1_p, g_ljci_c1, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_c1_p, g_ljci_c1, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_c1_p, g_ljci_c1, batch_k%length, batch_j%length)
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
                  call wf%point_vvvo(g_bdck_c1_p, g_bdcj_c1, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbjc, batch_k%length)
!
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
                  call wf%setup_oovo(wf%eri_c1, g_lkcj_c1, g_lkcj_c1_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_c1_p, g_licj_c1, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_c1_p, g_lkcj_c1, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_c1_p, g_ljci_c1, batch_k%length, batch_i%length)
!
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_ilkc_p, g_iljc, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_kljc, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
!
                  call wf%setup_ovov(wf%eri_t1, g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
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
!                       Check for core orbitals:
!                       cvs: i,j,k cannot all correspond to valence orbitals
!                       rm_core: i,j,k may not contain any core orbital
                        if (wf%ijk_amplitudes_are_zero(i, j, k, cvs, rm_core)) cycle
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct C^{abc}_{ijk} for given i, j, k
!                       and calculate contributions to rho1 and rho2
!                       Using c1-transformed integrals the terms have the same form
!                       as the omega terms (where t_abc = c_abc)
!
                        call wf%construct_V(i, j, k, sorting, c_abc,      &
                                            t_abij, c_abij,               &
                                            g_bdci_p(:,:,:,i_rel),        &
                                            g_bdcj_p(:,:,:,j_rel),        &
                                            g_bdck_p(:,:,:,k_rel),        &
                                            g_bdci_c1_p(:,:,:,i_rel),     &
                                            g_bdcj_c1_p(:,:,:,j_rel),     &
                                            g_bdck_c1_p(:,:,:,k_rel),     &
                                            g_ljci_p(:,:,j_rel,i_rel),    &
                                            g_lkci_p(:,:,k_rel,i_rel),    &
                                            g_lkcj_p(:,:,k_rel,j_rel),    &
                                            g_licj_p(:,:,i_rel,j_rel),    &
                                            g_lick_p(:,:,i_rel,k_rel),    &
                                            g_ljck_p(:,:,j_rel,k_rel),    &
                                            g_ljci_c1_p(:,:,j_rel,i_rel), &
                                            g_lkci_c1_p(:,:,k_rel,i_rel), &
                                            g_lkcj_c1_p(:,:,k_rel,j_rel), &
                                            g_licj_c1_p(:,:,i_rel,j_rel), &
                                            g_lick_c1_p(:,:,i_rel,k_rel), &
                                            g_ljck_c1_p(:,:,j_rel,k_rel))
!
                        call construct_contravariant_t3(c_abc, sorting, wf%n_v)
                        call wf%divide_by_orbital_differences(i, j, k, c_abc, omega)
!
                        call wf%omega_cc3_contractions(i, j, k, c_abc, sorting,   &
                                                       rho_ai, rho_abij, F_ov_ck, &
                                                       g_dbic_p(:,:,:,i_rel),     &
                                                       g_dbjc_p(:,:,:,j_rel),     &
                                                       g_dbkc_p(:,:,:,k_rel),     &
                                                       g_jlic_p(:,:,j_rel,i_rel), &
                                                       g_klic_p(:,:,k_rel,i_rel), &
                                                       g_kljc_p(:,:,k_rel,j_rel), &
                                                       g_iljc_p(:,:,i_rel,j_rel), &
                                                       g_ilkc_p(:,:,i_rel,k_rel), &
                                                       g_jlkc_p(:,:,j_rel,k_rel), &
                                                       g_ibjc_p(:,:,i_rel,j_rel), &
                                                       g_ibkc_p(:,:,i_rel,k_rel), &
                                                       g_jbkc_p(:,:,j_rel,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      call mem%dealloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci_c1, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
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
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
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
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      call mem%dealloc(F_ov_ck, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%batch_finalize()
!
      call cc3_timer_c3%turn_off()
!
   end subroutine jacobian_cc3_c3_a_cc3
!
!
   module subroutine construct_c1_fock_cc3(wf, F_ia_c1)
!!
!!    Construct C1-transformed fock matrix
!!    Written by Alexander C. Paul, Feb 2019
!!
!!    Calculates C1-transformed occupied-virtual elements of the Fock matrix
!!    required for the CC3 jacobian and returns it ordered as n_v, n_o
!!
!!    F_ia_c1 = sum_j L_iajj' = sum_j 2 g_iajj' - g_ij'ja
!!
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: F_ia_c1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajk
!
      integer :: i, a, j
!
      call zero_array(F_ia_c1, wf%n_v*wf%n_o)
!
!     Construct the g_ovoo_c1
!
      call mem%alloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%eri_c1%get('ovoo', g_iajk)
!
!     Add contributions and resort to F_ia_c1(a,i)
!
!$omp parallel do private(a,i,j)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
               F_ia_c1(a,i) = F_ia_c1(a,i) + two*g_iajk(i,a,j,j) - g_iajk(j,a,i,j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine construct_c1_fock_cc3
!
!
   module subroutine jacobian_cc3_b2_fock_cc3(wf, i, j, k, u_abc, v_abc, rho_ablj, F_ov_ck)
!!
!!    Jacobian CC3 B2 fock
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!    Adapted to give a contravariant representation of rho2 due to the
!!    use of a contravariant representation of t3
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    rho_2 =+ P^{ab}_{ij} sum_kc (t^abc_ijk - t^cba_ijk) F_kc
!!
!!    This can be reformulated to use the "contravariant" t3
!!
!!    ~rho^ab_ij = 2rho^ab_ij - rho^ba_ij =+ P^{ab}_{ij} sum_kc u_abc F_kc
!!
!!    where:
!!    u_abc = (4t_abc - 2t_acb - 2t_cba - 2t_bac + t_bca + t_cab)
!!
!!    Note that rho_abij has to be symmetrized outside of this routine
!!    This routine is also used for Z_bcjk = tbar^abc_ijk R^a_i
!!
      use reordering, only: sort_123_to_312, sort_123_to_213
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_ablj
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: F_ov_ck
!
      real(dp) :: factor_ij, factor_jk
!
      if (i .ne. j .and. j .ne. k) then
         factor_ij = one
         factor_jk = one
      else if (j .eq. k )  then
         factor_ij = one
         factor_jk = half
      else ! i == j
         factor_ij = half
         factor_jk = one
      end if
!
!     rho_abjk += sum_c u_abc*F_ic
      call wf%omega2_fock_cc3_permutation(u_abc, F_ov_ck(:,i), rho_ablj(:,:,j,k), factor_jk)
!
!     abc -> cab
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     rho^ab_ij += sum_c v_abc F^C1_kc
      call wf%omega2_fock_cc3_permutation(v_abc, F_ov_ck(:,k), rho_ablj(:,:,i,j), factor_ij)
!
      if (k .ne. j .and. j .ne. i) then
!
!        acb -> cab
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        rho^ab_ik += sum_c v_abc*F_jc
         call wf%omega2_fock_cc3_permutation(v_abc, F_ov_ck(:,j), rho_ablj(:,:,i,k), one)
!
      end if
!
   end subroutine jacobian_cc3_b2_fock_cc3
!
!
   module subroutine rho2_fock_cc3_permutation_cc3(wf, o1, o2, R3, F_ov_c1, rho2, factor)
!!
!!    Rho doubles Fock CC3 permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to CC3 omega2 vector, e.g.:
!!    ~rho_abjk += sum_c u_abc F_ic
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: o1, o2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: R3
!
      real(dp), dimension(wf%n_v), intent(in) :: F_ov_c1
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho2
!
      real(dp), intent(in) :: factor
!
      call dgemv('N',            &
                 wf%n_v**2,      &
                 wf%n_v,         &
                 factor,         &
                 R3,             &
                 wf%n_v**2,      &
                 F_ov_c1, 1,     &
                 one,            &
                 rho2(:,:,o1,o2), 1)
!
   end subroutine rho2_fock_cc3_permutation_cc3
!
!
end submodule jacobian
