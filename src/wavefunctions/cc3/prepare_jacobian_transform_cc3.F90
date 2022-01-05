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
submodule (cc3_class) prepare_jacobian_transform
!
!!
!!    Prepare jacobian transformation
!!
!!    Routines setting up the files containing intermediates for the linear
!!    transform of trial vectors by the Jacobian matrix and its transpose.
!!
!!    X_dbai = - sum_cjk (2 * t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!    X_ajil = - sum_bck (2 * t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_cc3(wf)
!!
!!    Prepare for jacobian
!!    Written by Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      type(timings) :: prep_timer
!
      call output%printf('v', 'Preparing for (a0) right excited state equations', &
                         chars=[trim(wf%name_)], fs='(/t3,a)')
!
      call wf%ccsd%prepare_for_jacobian()
!
      prep_timer = timings("Time preparing for CC3 Jacobian", pl='normal')
      call prep_timer%turn_on()
!
      call wf%prepare_cc3_jacobian_intermediates()
!
      call prep_timer%turn_off()
!
   end subroutine prepare_for_jacobian_cc3
!
!
   module subroutine prepare_for_jacobian_transpose_cc3(wf)
!!
!!    Prepare for jacobian transpose transformation
!!    Written by Rolf H. Myhre, April 2019
!!
!!    Modified by Tor S. Haugland
!!
!!    Also prepares CCSD jacobian_transpose.
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      type(timings) :: prep_timer
!
      call output%printf('v', 'Preparing for (a0) left excited state equations' &
                         &, chars=[trim(wf%name_)], fs='(/t3,a)')
!
      call wf%ccsd%prepare_for_jacobian_transpose()
!
      prep_timer = timings("Time preparing for CC3 Jacobian transpose", pl='normal')
      call prep_timer%turn_on()
!
      call wf%prepare_cc3_jacobian_intermediates()
!
      call prep_timer%turn_off()
!
   end subroutine prepare_for_jacobian_transpose_cc3
!
!
   module subroutine prepare_cc3_jacobian_intermediates_cc3(wf)
!!
!!
!!    Prepare intermediates for jacobian CC3 transformations
!!    written by Rolf H. Myhre and Alexander C. Paul, April 2019
!!    Adapted to a contravariant representation of t3
!!    by Rolf H. Myhre and Alexander C. Paul, Sep 2020
!!
!!    Construct X_abdi and X_ajil needed in CC3 jacobian transpose and store on disk
!!    For that: construct t^abc_ijk in single batches of ijk
!!    and contract with the respective integrals
!!
!!    t^abc = - (eps^abc)^-1 P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il(lj|ck))
!!
!!    X_abid -= sum_jck u^abc_ijk g_kcjd
!!    X_ajil -= sum_bck u^abc_ijk g_lbkc
!!    where
!!       u^abc_ijk = 4t^abc_ijk + t_bca_ijk + t_cab_ijk
!!                 - 2t^acb_ijk - 2t_cba_ijk - 2t_bac_ijk
!!
      use reordering, only: squareup_and_sort_1234_to_1324, sort_1234_to_1342
      use reordering, only: construct_contravariant_t3
      use array_utilities, only: zero_array
!
      implicit none
!
      class(cc3) :: wf
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Arrays for intermediates cannot hold the whole X_dbai array
      real(dp), dimension(:,:,:,:), allocatable, target :: X_dbai
      real(dp), dimension(:,:,:,:), allocatable, target :: X_dbaj
      real(dp), dimension(:,:,:,:), allocatable, target :: X_dbak
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_dbai_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_dbaj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_dbak_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable :: X_alji
      real(dp), dimension(:,:,:,:), allocatable :: X_ajil
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lbkc_p => null()
!
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k, batch_full
      integer :: i_batch, j_batch, k_batch ! used for the current batch
      integer :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer :: req_single_batch
!
!     Alloc and squareup the t2 amplitudes
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Need a batching index of full dimension to get the ovov integral
      batch_full = batching_index(wf%n_o)
      call batch_full%do_single_batch
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
      req_single_batch = req_0 + req_1_eri*wf%n_o + 2*wf%n_v**3*wf%n_o &
                       + wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 2*wf%n_v**3 + wf%n_v**2*wf%n_o
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 2*wf%n_o*wf%n_v
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,            &
                           req_0, req_i, req_1, req_1,           &
                           req_2, req_2, req_2, req_3,           &
                           'prepare_cc3_jacobian_intermediates', &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(X_alji, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(X_alji, wf%n_v*wf%n_o**3)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%alloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(X_dbai, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
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
         call mem%alloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, batch_i%max_length)
         call mem%alloc(g_lbjc, wf%n_v, wf%n_v, wf%n_o, batch_i%max_length)
         call mem%alloc(g_lbkc, wf%n_v, wf%n_v, wf%n_o, batch_i%max_length)
!
         call mem%alloc(X_dbai, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(X_dbaj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(X_dbak, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      wf%X_dbai = direct_stream_file('X_dbai', wf%n_v**3)
      call wf%X_dbai%open_()
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(g_bdci, g_bdci_p, sorting, batch_i)
!
         call wf%setup_ovov(g_lbic, g_lbic_p, sorting, batch_full, batch_i)
!
         call zero_array(X_dbai, batch_i%length*wf%n_v**3)
         X_dbai_p => X_dbai
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            if (j_batch .ne. i_batch) then ! read for switched i - j
!
               call wf%setup_vvvo(g_bdcj, g_bdcj_p, sorting, batch_j)
!
               call wf%setup_ovov(g_lbjc, g_lbjc_p, sorting, batch_full, batch_j)
!
               call wf%X_dbai%read_range(X_dbaj, batch_j)
               X_dbaj_p => X_dbaj
!
               call wf%setup_oovo(g_licj, g_licj_p, sorting, batch_i, batch_j)
!
            else ! j_batch == i_batch
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
!
               call wf%point_vvoo(g_lbjc_p, g_lbic, wf%n_o, batch_j%length)
!
               X_dbaj_p => X_dbai
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
                  call wf%setup_ovov(g_lbkc, g_lbkc_p, sorting, batch_full, batch_k)
!
                  call wf%X_dbai%read_range(X_dbak, batch_k)
                  X_dbak_p => X_dbak
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
                  call wf%point_vvoo(g_lbkc_p, g_lbic, wf%n_o, batch_k%length)
!
                  X_dbak_p => X_dbai
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
                  call wf%point_vvoo(g_lbkc_p, g_lbjc, wf%n_o, batch_k%length)
!
                  X_dbak_p => X_dbaj
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
                        if (k .eq. i) then ! k == j == i
                           cycle
                        end if
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
                        call wf%construct_x_intermediates(i, j, k, t_abc, sorting, &
                                                          X_alji,                 &
                                                          X_dbai_p(:,:,:,i_rel),  &
                                                          X_dbaj_p(:,:,:,j_rel),  &
                                                          X_dbak_p(:,:,:,k_rel),  &
                                                          g_lbic_p(:,:,:,i_rel),  &
                                                          g_lbjc_p(:,:,:,j_rel),  &
                                                          g_lbkc_p(:,:,:,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
                  call wf%X_dbai%write_range(X_dbak, batch_k)
               endif
!
            enddo ! batch_k
!
            if (j_batch .ne. i_batch) then
               call wf%X_dbai%write_range(X_dbaj, batch_j)
            endif
!
         enddo ! batch_j
!
         call wf%X_dbai%write_range(X_dbai, batch_i)
!
      enddo ! batch_i
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(X_dbai, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
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
         call mem%dealloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, batch_i%max_length)
         call mem%dealloc(g_lbjc, wf%n_v, wf%n_v, wf%n_o, batch_i%max_length)
         call mem%dealloc(g_lbkc, wf%n_v, wf%n_v, wf%n_o, batch_i%max_length)
!
         call mem%dealloc(X_dbai, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(X_dbaj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(X_dbak, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
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
!     Resort X_dbai to X_abid for the final contraction with C^ac_il to sigma_dl
!
      call wf%sort_x_to_abid_and_write()
!
      call wf%X_dbai%close_('delete')
!
!     sort X_alji to ajil and write to disk
!
      call mem%alloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1342(X_alji, X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_alji, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      wf%X_ajil = direct_stream_file('X_ajil',wf%n_v*wf%n_o**2)
      call wf%X_ajil%open_('write')
!
      call wf%X_ajil%write_(X_ajil, 1, wf%n_o)
!
      call mem%dealloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%X_ajil%close_()
!
   end subroutine prepare_cc3_jacobian_intermediates_cc3
!
!
   module subroutine construct_x_intermediates_cc3(wf, i, j, k, u_abc, v_abc, X_alji, &
                                                   X_dbai, X_dbaj, X_dbak, &
                                                   g_lbic, g_lbjc, g_lbkc)
!!
!!    Construct X intermediates
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!    Adapted to a contravariant representation of t3
!!    by Rolf H. Myhre and Alexander C. Paul, Sep 2020
!!
!!    Constructs the intermediates X_vvvo and X_vooo
!!    used to compute the contributions to sigma_ai
!!
!!    X_dbai = sum_cjk (2t^bac + 2t^cba + 2t^acb - t_bca - t_cab - 4 * t_abc) * g_kcjd
!!    X_alji = sum_cjk (2t^bac + 2t^cba + 2t^acb - t_bca - t_cab - 4 * t_abc) * g_lbkc
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
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)       :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)       :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out) :: X_alji
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)         :: X_dbai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)         :: X_dbaj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)         :: X_dbak
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)          :: g_lbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)          :: g_lbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)          :: g_lbkc
!
!     X_dbak -= sum_c,ij u_cba g_jdic
      call wf%construct_Y_vvvo_permutation(u_abc, g_lbic(:,:,j), X_dbak, -one)
!
!     X_alik -= sum_bc,j u_bca g_lbjc
      call wf%construct_Y_vooo_permutation(i, k, u_abc, g_lbjc, X_alji, -one)
!
!     bac -> cba
!     cab -> bca
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     X_dbaj -= sum_c,ik u_cab g_idkc
      call wf%construct_Y_vvvo_permutation(v_abc, g_lbkc(:,:,i), X_dbaj, -one)
!
!     X_alkj -= sum_bc,i u_cab g_lbic
      call wf%construct_Y_vooo_permutation(k, j, v_abc, g_lbic, X_alji, -one)
!
!     acb -> bac -> cba
!     abc -> cab -> bca
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     X_dbai -= sum_c,jk u_acb g_kdjc
      call wf%construct_Y_vvvo_permutation(u_abc, g_lbjc(:,:,k), X_dbai, -one)
!
!     X_alji -= sum_bc,k u_abc g_lbkc
      call wf%construct_Y_vooo_permutation(j, i, u_abc, g_lbkc, X_alji, -one)
!
      if (i .ne. j .and. j .ne. k) then
!
!        (b <-> c)
!        abc -> acb -> bac -> cba
!        acb -> abc -> cab -> bca
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        X_dbai -= sum_c,jk u_acb g_jdkc
         call wf%construct_Y_vvvo_permutation(v_abc, g_lbkc(:,:,j), X_dbai, -one)
!
!        X_alki -= sum_bc,j u_acb g_lbjc
         call wf%construct_Y_vooo_permutation(k, i, v_abc, g_lbjc, X_alji, -one)
!
!        cab -> abc -> acb -> bac -> cba
!        bac -> acb -> abc -> cab -> bca
         call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        X_dbaj -= sum_c,ik u_cab g_kdic
         call wf%construct_Y_vvvo_permutation(u_abc, g_lbic(:,:,k), X_dbaj, -one)
!
!        X_alij = sum_ab,k u_abc g_lbkc
         call wf%construct_Y_vooo_permutation(i, j, u_abc, g_lbkc, X_alji, -one)
!
!        bca -> cab -> abc -> acb -> bac -> cba
!        cba -> bac -> acb -> abc -> cab -> bca
         call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        X_dbak -= sum_c,ij u_bca g_idjc
         call wf%construct_Y_vvvo_permutation(v_abc, g_lbjc(:,:,i), X_dbak, -one)
!
!        X_aljk = sum_ab,i u_cba g_lbic
         call wf%construct_Y_vooo_permutation(j, k, v_abc, g_lbic, X_alji, -one)
!
      end if
!
   end subroutine construct_x_intermediates_cc3
!
!
   module subroutine sort_x_to_abid_and_write_cc3(wf)
!!
!!    Read in intermediate X_dbai from file, resort to X_abid and write to file again
!!
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      use reordering, only: sort_1234_to_3241
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: X_dbai
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
!
      type(batching_index) :: batch_i
      integer :: i_batch
      integer :: req_0, req_i
!
      req_0 = 0
      req_i = 2*wf%n_v**3
!
      batch_i = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, req_0, req_i, 'sort_x_to_abid_and_write_cc3')
!
      wf%X_abid = direct_stream_file('X_abid',wf%n_v**2)
      call wf%X_abid%open_('write')
!
      call mem%alloc(X_dbai, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
      call mem%alloc(X_abid, wf%n_v, wf%n_v, batch_i%max_length, wf%n_v)
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
!        Read from file
!
         call wf%X_dbai%read_range(X_dbai, batch_i)
!
!        Sort X_dbai to X_abid
!
         call sort_1234_to_3241(X_dbai, X_abid, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
!        Write to file
!
         call wf%X_abid%write_compound_batch_full(X_abid, batch_i, wf%n_v)
!
      enddo ! batch_i
!
      call mem%dealloc(X_abid, wf%n_v, wf%n_v, batch_i%max_length, wf%n_v)
      call mem%dealloc(X_dbai, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
      call mem%batch_finalize()
!
      call wf%X_abid%close_()
!
   end subroutine sort_x_to_abid_and_write_cc3
!
!
end submodule prepare_jacobian_transform
