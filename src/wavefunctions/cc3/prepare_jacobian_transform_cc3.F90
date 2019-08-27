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
submodule (cc3_class) prepare_jacobian_transform
!
!!
!!    Prepare jacobian transformation (cc3)
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
!!    Routines setting up the files containing intermediates for the linear 
!!    transform of trial vectors by the Jacobian matrix and its transpose.
!!
!!    X_abdi = - sum_cjk (2 * t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!    X_ajil = - sum_bck (2 * t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!
!!
!!    Also contains the routine setting up the integral files 
!!    for the transpose transformation:
!!    g_vvvv ordered 1324
!!    g_oooo ordered 1243
!!    g_voov ordered 1324 
!!    g_vvoo ordered 1342
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prep_cc3_jacobian_trans_integrals_cc3(wf)
!!
!!    Construct integrals needed in CC3 jacobian transpose and store on disk
!!    (be|cd) ordered as bce,d
!!    (mj|lk) ordered as mjk,l
!!    (ck|ld) ordered as cl,kd
!!    (cd|lk) ordered as cl,kd
!!
!!    written by Rolf H. Myhre and Alexander Paul, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      type(batching_index) :: batch_d, batch_l
!
      integer :: req_0, req_d, req_l
      integer :: d_batch, l_batch
!
!     (be|cd) stored as bce#d
!
      call wf%get_g_pqrs_required(req_0, req_d, wf%n_v, wf%n_v, wf%n_v, 1)
      req_d = req_d + 2*wf%n_v**3
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d,req_0,req_d)
      call batch_d%determine_limits(1)
!
      call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_d%length)
      call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_d%length)
!
      wf%g_becd_t = direct_file('g_becd_t',wf%n_v**3)
      call wf%g_becd_t%open_('write')
!
      do d_batch = 1,batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
         call wf%get_vvvv(g_pqrs,   &
                          1,wf%n_v, &
                          1,wf%n_v, &
                          1,wf%n_v, &
                          batch_d%first,batch_d%last)
!
         call sort_1234_to_1324(g_pqrs, h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_d%length)
!
         call single_record_writer(batch_d, wf%g_becd_t, h_pqrs)
!
      enddo
!
      call wf%g_becd_t%close_()
!
      call batch_d%determine_limits(1)
      call mem%dealloc(g_pqrs,wf%n_v,wf%n_v,wf%n_v,batch_d%length)
      call mem%dealloc(h_pqrs,wf%n_v,wf%n_v,wf%n_v,batch_d%length)
!
!
!     (mj|lk) stored as mjk#l
!
      call wf%get_g_pqrs_required(req_0,req_l,wf%n_o,wf%n_o,1,wf%n_o)
      req_l = req_l + 2*wf%n_o**3
!
      call batch_l%init(wf%n_o)
      call mem%batch_setup(batch_l,req_0,req_l)
      call batch_l%determine_limits(1)
!
      call mem%alloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_l%length)
!
      wf%g_mjlk_t = direct_file('g_mjlk_t',wf%n_o**3)
      call wf%g_mjlk_t%open_('write')
!
      do l_batch = 1,batch_l%num_batches
!
         call batch_l%determine_limits(l_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, batch_l%length, wf%n_o)
!
         call wf%get_oooo(g_pqrs,   &
                          1,wf%n_o, &
                          1,wf%n_o, &
                          batch_l%first,batch_l%last, &
                          1,wf%n_o)
!
         call sort_1234_to_1243(g_pqrs, h_pqrs, wf%n_o, wf%n_o, batch_l%length, wf%n_o)
!
         call single_record_writer(batch_l, wf%g_mjlk_t, h_pqrs)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, batch_l%length, wf%n_o)
!
      enddo
!
      call wf%g_mjlk_t%close_()
!
      call batch_l%determine_limits(1)
      call mem%dealloc(h_pqrs,wf%n_o,wf%n_o,wf%n_o,batch_l%length)
!
!
!     (ck|ld) !stored as cl#k#d
!
      call wf%get_g_pqrs_required(req_0, req_d, wf%n_v, wf%n_o, wf%n_o, 1)
      req_d = req_d + wf%n_v*wf%n_o**2
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d,req_0,req_d)
      call batch_d%determine_limits(1)
!
      call mem%alloc(g_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length) 
      call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length) 
!
      wf%g_ckld_t = direct_file('g_ckld_t',wf%n_v*wf%n_o)
      call wf%g_ckld_t%open_('write')
!
      do d_batch = 1,batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
         call wf%get_voov(g_pqrs,   &
                          1,wf%n_v, &
                          1,wf%n_o, &
                          1,wf%n_o, &
                          batch_d%first,batch_d%last)
!
         call sort_1234_to_1324(g_pqrs, h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length)
         call compound_record_writer(wf%n_o, batch_d, wf%g_ckld_t, h_pqrs) !store as cl#k#d
!
      enddo
!
      call wf%g_ckld_t%close_()
!
      call batch_d%determine_limits(1)
      call mem%dealloc(g_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length) 
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length) 
!
!
!     (cd|lk) !stored as cl#k#d
!
      call wf%get_g_pqrs_required(req_0, req_d, wf%n_v, 1, wf%n_o, wf%n_o)
      req_d = req_d + 2*wf%n_v*wf%n_o**2
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d,req_0,req_d)
      call batch_d%determine_limits(1)
!
      call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length)
!
      wf%g_cdlk_t = direct_file('g_cdlk_t',wf%n_v*wf%n_o)
      call wf%g_cdlk_t%open_('write')
!
!     Going to batch over both d and k, so both are used as records, store as cl#k#d
!
      do d_batch = 1,batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
         call mem%alloc(g_pqrs, wf%n_v, batch_d%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_pqrs,                     &
                          1,wf%n_v,                   &
                          batch_d%first,batch_d%last, &
                          1,wf%n_o,                   &
                          1,wf%n_o)
!
         call sort_1234_to_1342(g_pqrs, h_pqrs, wf%n_v, batch_d%length, wf%n_o, wf%n_o) !sort to clkd
!
         call compound_record_writer(wf%n_o, batch_d, wf%g_cdlk_t, h_pqrs)
!
         call mem%dealloc(g_pqrs, wf%n_v, batch_d%length, wf%n_o, wf%n_o)
!
      enddo
!
      call wf%g_cdlk_t%close_()
!
      call batch_d%determine_limits(1)
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_d%length)
!
   end subroutine prep_cc3_jacobian_trans_integrals_cc3
!
!
   module subroutine prep_cc3_g_lbkc_t_file_cc3(wf)
!!
!!    Construct ovov-integral needed only in the construction of the intermediates 
!!    for the CC3 jacobian transformations and store on disk
!!
!!    (lb|kc) ordered as bcl,k
!!
!!    written by Rolf H. Myhre and Alexander Paul, Mai 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      type(batching_index) :: batch_k
!
      integer :: req_0, req_k
      integer :: k_batch
!
!     (lb|kc)  ! stored as bcl#k
!
      call wf%get_g_pqrs_required(req_0,req_k,wf%n_o,wf%n_v,1,wf%n_v)
      req_k = req_k + 2*wf%n_v**2*wf%n_o
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
      call batch_k%determine_limits(1)
!
      call mem%alloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
      wf%g_lbkc_t = direct_file('g_lbkc_t',wf%n_o*wf%n_v**2)
      call wf%g_lbkc_t%open_('write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o , wf%n_v , batch_k%length , wf%n_v)
!
         call wf%get_ovov(g_pqrs,   &
                          1,wf%n_o, &
                          1,wf%n_v, &
                          batch_k%first,batch_k%last, &
                          1,wf%n_v)
!
         call sort_1234_to_2413(g_pqrs , h_pqrs , wf%n_o , wf%n_v , batch_k%length , wf%n_v) ! sort to bclk
!
         call single_record_writer(batch_k, wf%g_lbkc_t, h_pqrs)
!
         call mem%dealloc(g_pqrs, wf%n_o , wf%n_v , batch_k%length , wf%n_v)
!
      enddo
!
      call batch_k%determine_limits(1)
      call mem%dealloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
      call wf%g_lbkc_t%close_()

!
   end subroutine prep_cc3_g_lbkc_t_file_cc3
!
!
   module subroutine prep_cc3_jacobian_intermediates_cc3(wf)
!!
!!    Construct X_abdi and X_ajil needed in CC3 jacobian transpose and store on disk
!!    For that: construct t^abc_ijk in single batches of ijk 
!!    and contract with the respective integrals
!!
!!    t^abc_ijk = - (ε^abc_ijk)^-1 P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il(lj|ck))
!!
!!    X_abid = - sum_jck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!    X_ajil = - sum_bck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
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
      real(dp), dimension(:,:,:), allocatable :: v_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Arrays for intermediates
!     cannot hold the whole X_abdi array
      real(dp), dimension(:,:,:,:), allocatable, target :: X_abdi
      real(dp), dimension(:,:,:,:), allocatable, target :: X_abdj
      real(dp), dimension(:,:,:,:), allocatable, target :: X_abdk
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_abdi_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_abdj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_abdk_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable :: X_alij
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
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch ! used for the current batch
      integer :: req_0, req_1, req_2, req_3
!
!     Construct the g_lbkc_t file (only needed for the intermediates)
!     g_ljck_t and g_bdck_t are already on disk from the ground state calculation
!
      call wf%prep_cc3_g_lbkc_t_file()
!
!     Alloc and squareup the t2 amplitudes
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Setup and Batching loops
!
      req_0 = 3*wf%n_v**3 + wf%n_v*wf%n_o**3
      req_1 = 2*wf%n_v**3 + wf%n_o*wf%n_v**2
      req_2 = wf%n_o*wf%n_v
      req_3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                                 req_0, req_1, req_2, req_3, zero)
!
!     Allocate integral arrays
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then ! no batching
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
!     Arrays for the triples amplitudes
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Remaining integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, batch_i%length)
         call mem%alloc(g_lbjc, wf%n_v, wf%n_v, wf%n_o, batch_i%length)
         call mem%alloc(g_lbkc, wf%n_v, wf%n_v, wf%n_o, batch_i%length)
!
         call mem%alloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(X_abdj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(X_abdk, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
      end if
!
!     Array for the whole intermediate X_alij
      call mem%alloc(X_alij, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(X_alij, wf%n_v*wf%n_o**3)
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_lbkc_t%open_('read')
!
      wf%X_abdi = direct_file('X_abdi',wf%n_v**3)
      call wf%X_abdi%open_()
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci, wf%g_lbkc_t, g_lbic)
         g_bdci_p => g_bdci
         g_lbic_p => g_lbic
!
         call zero_array(X_abdi, wf%n_o*wf%n_v**3)
         X_abdi_p => X_abdi
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci)
            g_ljci_p => g_ljci
!
            if (j_batch .ne. i_batch) then ! read for switched i - j
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj, wf%g_lbkc_t, g_lbjc)
               g_bdcj_p => g_bdcj
               g_lbjc_p => g_lbjc
!
               call single_record_reader(batch_j, wf%X_abdi, X_abdj)
               X_abdj_p => X_abdj
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj)
               g_licj_p => g_licj
!
            else ! j_batch == i_batch
!
               g_bdcj_p => g_bdci
               g_lbjc_p => g_lbic
!
               X_abdj_p => X_abdi
!
               g_licj_p => g_ljci
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck, wf%g_lbkc_t, g_lbkc)
                  g_bdck_p => g_bdck
                  g_lbkc_p => g_lbkc
!
                  call single_record_reader(batch_k, wf%X_abdi, X_abdk)
                  X_abdk_p => X_abdk
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
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
                  g_lbkc_p => g_lbic
!
                  X_abdk_p => X_abdi
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_ljci
!
                  g_lkcj_p => g_ljci
                  g_ljck_p => g_ljci
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
                  g_lbkc_p => g_lbjc
!
                  X_abdk_p => X_abdj
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_licj
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
                  g_lkcj_p => g_lkcj
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
                        if (k .eq. i) then ! k == j == i
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct t^{abc}_{ijk} for given i, j, k
!
                        call wf%omega_cc3_W_calc(i, j, k, t_abc, u_abc, t_abij,  &
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
                        call wf%construct_x_intermediates(i, j, k, t_abc, u_abc, v_abc,  &
                                                          X_alij,                        &
                                                          X_abdi_p(:,:,:,i_rel),         &
                                                          X_abdj_p(:,:,:,j_rel),         &
                                                          X_abdk_p(:,:,:,k_rel),         &
                                                          g_lbic_p(:,:,:,i_rel),         &
                                                          g_lbjc_p(:,:,:,j_rel),         &
                                                          g_lbkc_p(:,:,:,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
                  call single_record_writer(batch_k, wf%X_abdi, X_abdk)
               endif
!
            enddo ! batch_k
!
            if (j_batch .ne. i_batch) then
               call single_record_writer(batch_j, wf%X_abdi, X_abdj)
            endif
!
         enddo ! batch_j
!
         call single_record_writer(batch_i, wf%X_abdi, X_abdi)
!
      enddo ! batch_i
!
!     Close files: 
!
      call wf%g_bdck_t%close_()
      call wf%g_ljck_t%close_()
!
!     g_lbkc_t file only needed for the construction of the intermediates
      call wf%g_lbkc_t%close_('delete')
!
!     Deallocate integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
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
         call mem%dealloc(g_lbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lbjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(X_abdj, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(X_abdk, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
      endif
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Resort X_abdi to X_abid for the final contraction with C^ac_il to sigma_dl
!
      call wf%sort_x_to_abid_and_write()
!
      call wf%X_abdi%close_('delete')
!
!     sort X_alij (ordered as alij) to ajil and write to disk 
!
      call mem%alloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1432(X_alij, X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_alij, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      wf%X_ajil = direct_file('X_ajil',wf%n_v*wf%n_o**2)
      call wf%X_ajil%open_('write')
!
      call single_record_writer(wf%n_o, wf%X_ajil, X_ajil)
!
      call mem%dealloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%X_ajil%close_()
!
   end subroutine prep_cc3_jacobian_intermediates_cc3
!
!
   module subroutine construct_x_intermediates_cc3(wf, i, j, k, t_abc, u_abc, v_abc, x_alij,      &
                                                   X_abdi, X_abdj, X_abdk, g_lbic, g_lbjc, g_lbkc)
!!
!!    Constructs the intermediates X_abdi and X_ajil used to compute the contributions to sigma_ai
!!
!!    X_abdi = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_kcjd
!!    X_ajil = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_lbkc
!!
!!    g_lbic, g_lbjc, g_lbkc can be used for g_pcqd as well: 
!!    The p(i,j,k) can be set in dgemm and q(i,j,k) is defined by the array used
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out)  :: X_alij ! ordered alij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: X_abdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: X_abdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: X_abdk
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)           :: g_lbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)           :: g_lbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o), intent(in)           :: g_lbkc
!
!     construct u_abc = 2*t_abc - t_acb - t_cba
!
      call construct_123_min_132_min_321(t_abc, u_abc, wf%n_v)
!
!     X_abdi += - sum_c,kj (2*t_abc - t_acb - t_cba)*g_kcjd
!
      call dgemm('N','N',           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_v,           &
                  -one,             &
                  u_abc,            & ! u_ab_c
                  wf%n_v**2,        &
                  g_lbjc(:,:,k),    & ! g_c_d,kj
                  wf%n_v,           &
                  one,              &
                  X_abdi,           & ! X_ab_d,i
                  wf%n_v**2)
!
!     X_alij += - sum_bc,k (2*t_abc - t_acb - t_cba)*g_lbkc
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  -one,             &
                  u_abc,            & ! u_a_bc
                  wf%n_v,           &
                  g_lbkc,           & ! g_bc_l,k
                  wf%n_v**2,        &
                  one,              &
                  X_alij(:,:,i,j),  & ! X_a_l,ij
                  wf%n_v)
!
      if (j .ne. i) then
!
!        resort to v_abc = 2*t_bac - t_bca - t_cab
!
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        X_abdj += - sum_c,ki (2*t_bac - t_cab - t_bca)*g_kcid
!
         call dgemm('N','N',           &
                     wf%n_v**2,        &
                     wf%n_v,           &
                     wf%n_v,           &
                     -one,             &
                     v_abc,            & ! v_ab_c
                     wf%n_v**2,        &
                     g_lbic(:,:,k),    & ! g_c_d,ki
                     wf%n_v,           &
                     one,              &
                     X_abdj,           & ! X_ab_d,j
                     wf%n_v**2)
!
!        X_ajil += - sum_bc,k (2*t_bac - t_cab - t_bca)*g_lbkc
!
         call dgemm('N','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     -one,             &
                     v_abc,            & ! v_a_bc
                     wf%n_v,           &
                     g_lbkc,           & ! g_bc_l,k
                     wf%n_v**2,        &
                     one,              &
                     X_alij(:,:,j,i),  & ! X_a_l,ji
                     wf%n_v)
!
      end if ! j .ne. i
!
!     construct u_cba = 2*t_cba - t_bca - t_abc
!
      call construct_321_min_231_min_123(t_abc, u_abc, wf%n_v)
!
!     X_abdj += - sum_c,ij (2*t_cba - t_bca - t_abc)*g_icjd
!
      call dgemm('N','N',           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_v,           &
                  -one,             &
                  u_abc,            & ! u_ab_c
                  wf%n_v**2,        &
                  g_lbjc(:,:,i),    & ! g_c_d,ij
                  wf%n_v,           &
                  one,              &
                  X_abdk,           & ! X_ab_d,k
                  wf%n_v**2)
!
!     X_akjl += - sum_bc,i (2*t_cba - t_bca - t_abc)*g_lbic
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  -one,             &
                  u_abc,            & ! u_a_bc
                  wf%n_v,           &
                  g_lbic,           & ! g_bc_l,i
                  wf%n_v**2,        &
                  one,              &
                  X_alij(:,:,k,j),  & ! X_a_l,kj
                  wf%n_v)
!
      if (k .ne. j) then
!
!        resort to u_cab = 2*t_cab - t_acb - t_bac
!
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        X_abdj += - sum_c,ik (2*t_cab - t_acb - t_bac)*(g_ickd)
!
         call dgemm('N','N',           & 
                     wf%n_v**2,        &
                     wf%n_v,           &
                     wf%n_v,           &
                     -one,             &
                     v_abc,            & ! v_ab_c
                     wf%n_v**2,        &
                     g_lbkc(:,:,i),    & ! g_c_d,ik
                     wf%n_v,           &
                     one,              &
                     X_abdj,           & ! X_ab_d,j
                     wf%n_v**2)
!
!        X_ajkl += - sum_bc,i (2*t_cab - t_acb - t_bac)*g_lbic
!
         call dgemm('N','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     -one,             &
                     v_abc,            & ! v_a_bc
                     wf%n_v,           &
                     g_lbic,           & ! g_bc_l,i
                     wf%n_v**2,        &
                     one,              &
                     X_alij(:,:,j,k),  & ! X_a_l,jk
                     wf%n_v)
!
         if (j .ne. i) then
!
!           construct u_acb = 2*t_acb - t_abc - t_cab
!
            call construct_132_min_123_min_312(t_abc, u_abc, wf%n_v)
!
!           X_abdi += - sum_c,jk (2*t_acb - t_abc - t_cab)*(g_jckd)
!
            call dgemm('N','N',           & 
                        wf%n_v**2,        &
                        wf%n_v,           &
                        wf%n_v,           &
                        -one,             &
                        u_abc,            & ! u_ab_c
                        wf%n_v**2,        &
                        g_lbkc(:,:,j),    & ! g_c_d,jk
                        wf%n_v,           &
                        one,              &
                        X_abdi,           & ! X_ab_d,i
                        wf%n_v**2)
!
!           X_aikl += - sum_bc,j (2*t_acb - t_abc - t_cab)*g_lbjc
!
            call dgemm('N','N',           &
                        wf%n_v,           &
                        wf%n_o,           &
                        wf%n_v**2,        &
                        -one,             &
                        u_abc,            & ! u_a_bc
                        wf%n_v,           &
                        g_lbjc,           & ! g_bc_l,j
                        wf%n_v**2,        &
                        one,              &
                        X_alij(:,:,i,k),  & ! X_a_l,ik
                        wf%n_v)
!
!           resort to u_bca = 2*t_bca - t_bac - t_cba
!
            call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!           X_abdk += - sum_c,ji (2*t_bca - t_bac - t_cba)*(g_jcid)
!
            call dgemm('N','N',           & 
                        wf%n_v**2,        &
                        wf%n_v,           &
                        wf%n_v,           &
                        -one,             &
                        v_abc,            & ! v_ab_c
                        wf%n_v**2,        &
                        g_lbic(:,:,j),    & ! g_c_d,ji
                        wf%n_v,           &
                        one,              &
                        X_abdk,           & ! X_ab_d,k
                        wf%n_v**2)
!
!           X_akil += - sum_bc,j (2*t_bca - t_bac - t_cba)*g_lbjc
!
            call dgemm('N','N',           &
                        wf%n_v,           &
                        wf%n_o,           &
                        wf%n_v**2,        &
                        -one,             &
                        v_abc,            & ! v_a_bc
                        wf%n_v,           &
                        g_lbjc,           & ! g_bc_l,j
                        wf%n_v**2,        &
                        one,              &
                        X_alij(:,:,k,i),  & ! X_a_l,ki
                        wf%n_v)
!
         end if ! j .ne. i
      end if ! k .ne. j
!                
   end subroutine construct_x_intermediates_cc3
!
!
   module subroutine sort_x_to_abid_and_write_cc3(wf)
!!
!!    Read in intermediate X_abdi from file, resort to X_abid and write to file again
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abdi
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
!
      type(batching_index) :: batch_i
      integer :: i_batch
      integer :: req_0, req_i
!
      req_0 = 0
      req_i = 2*wf%n_v**3
!
      call batch_i%init(wf%n_o)
!
      call mem%batch_setup(batch_i, req_0, req_i)
!
      wf%X_abid = direct_file('X_abid',wf%n_v**2)
      call wf%X_abid%open_('write')
!
      call batch_i%determine_limits(1)
      call mem%alloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
      call mem%alloc(X_abid, wf%n_v, wf%n_v, batch_i%length, wf%n_v)
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
!        Read from file
!
         call single_record_reader(batch_i, wf%X_abdi, X_abdi)
!
!        Sort X_abdi to X_abid
!
         call sort_1234_to_1243(X_abdi, X_abid, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
!        Write to file
!
         call compound_record_writer(wf%n_v, batch_i, wf%X_abid, X_abid, .true.)
!
      enddo ! batch_i
!
      call batch_i%determine_limits(1)
      call mem%dealloc(X_abid, wf%n_v, wf%n_v, batch_i%length, wf%n_v)
      call mem%dealloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
      call wf%X_abid%close_()
!
   end subroutine sort_x_to_abid_and_write_cc3
!
!
end submodule prepare_jacobian_transform
