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
submodule (cc3_class) prepare_jacobian_transpose
!
!!
!!    Prepare jacobian transpose submodule (cc3)
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
!!    Routines setting up the files containing intermediates and integrals 
!!    for the linear transform of trial vectors by the transpose of the Jacobian matrix
!!
!!    Sets up the integrals
!!    g_vvvv, g_ovov ordered 2413, g_voov ordered 2413, g_oooo ordered 1324
!!
!!    And the intermediates
!!    X_abdi = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_kcjd
!!    Y_akil = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_lbkc
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
!!    Called from solve - Write some integrals and intermediates to disk
!!
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
   module subroutine prepare_cc3_jacobian_transpose_integrals_cc3(wf)
!!
!!    Construct integrals needed in CC3 jacobian transpose and store on disk
!!    (ab|cd) ordered as abc,d
!!    (mi|lk) ordered as lm,ik
!!    (lb|kc) ordered as bcl,k
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
!     (lb|kc) 
!
      call wf%get_g_pqrs_required(req_0,req_k,wf%n_v,wf%n_o,1,wf%n_o)
      req_k = req_k + 2*wf%n_v**2*wf%n_o
!
      call batch_k%init(wf%n_o)
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call wf%g_lbkc_t%init('g_lbkc_t','direct','unformatted',dp*wf%n_o*wf%n_v**2)
      call disk%open_file(wf%g_lbkc_t,'write')
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
         call sort_1234_to_2413(g_pqrs , h_pqrs , wf%n_o , wf%n_v , batch_k%length , wf%n_v) ! sort to bclk
!
         do k = 1,batch_k%length
!
               record = batch_k%first + k - 1
               write(wf%g_lbkc_t%unit,rec=record,iostat=ioerror) h_pqrs(:,:,:,k)
!
               if(ioerror .ne. 0) then
                  call output%error_msg('Failed to write lbkc_t file')
               endif
!
         enddo
!
         call mem%dealloc(g_pqrs, wf%n_o , wf%n_v , batch_k%length , wf%n_v)
         call mem%dealloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
      enddo
!
      call disk%close_file(wf%g_lbkc_t,'keep')
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
         call sort_1234_to_1324(g_pqrs , h_pqrs , wf%n_o , wf%n_v , wf%n_v , batch_k%length) ! lcek
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
      call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_v, batch_k%length)
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
!!    Construct X_abdi and Y_akil needed in CC3 jacobian transpose and store on disk
!!    For that: construct t^abc_ijk in single batches of ijk 
!!    and contract with the respective integrals
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
!     cannot hold the whole X_abdi array
      real(dp), dimension(:,:,:,:), allocatable, target :: X_abdi
      real(dp), dimension(:,:,:,:), allocatable, target :: X_abdj
      real(dp), dimension(:,:,:,:), allocatable, target :: X_abdk
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_abdi_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_abdj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: X_abdk_p => null()
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
      req_1 = wf%n_v**3 + (wf%n_o)*(wf%n_v)**2
      req_2 = wf%n_v**2
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
         call mem%alloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
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
         call mem%alloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, batch_i%length)
         call mem%alloc(g_lbjc, wf%n_v, wf%n_v, wf%n_o, batch_i%length)
         call mem%alloc(g_lbkc, wf%n_v, wf%n_v, wf%n_o, batch_i%length)
!
      endif
!
      call mem%alloc(Y_aikl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call disk%open_file(wf%g_bdck_t,'read')
      call disk%open_file(wf%g_ljck_t,'read')
      call disk%open_file(wf%g_lbkc_t,'read')
!
      call wf%X_abdi%init('X_abdi','direct','unformatted',dp*wf%n_v**3)
      call disk%open_file(wf%X_abdi,'readwrite')
!
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%single_batch_reader(batch_i, wf%g_bdck_t, g_bdci, wf%g_lbkc_t, g_lbic)
         g_bdci_p => g_bdci
         g_lbic_p => g_lbic
!
!           cannot hold X_abdi - read in previous X, add contributions, write to disk again
!
            if (i_batch .gt. 1) then
               call wf%single_batch_reader(batch_i, wf%X_abdi, X_abdi)
               X_abdi_p => X_abdi
            end if
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%double_batch_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci)
            g_ljci_p => g_ljci
!
            if (j_batch .ne. i_batch) then ! read for switched i - j
!
               call wf%single_batch_reader(batch_j, wf%g_bdck_t, g_bdcj, wf%g_lbkc_t, g_lbjc)
               g_bdcj_p => g_bdcj
               g_lbjc_p => g_lbjc
!
!              Don't read X in the first iteration - X_abdi file empty
               if (i_batch .gt. 1 .or. j_batch .gt. 1) then
                  call wf%single_batch_reader(batch_j, wf%X_abdi, X_abdj)
                  X_abdj_p => X_abdj
               end if
!
               call wf%double_batch_reader(batch_i, batch_j, wf%g_ljck_t, g_licj)
               g_licj_p => g_licj
!
            else
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
               if (k_batch .ne. i_batch .and. k_batch .ne. j_batch) then
!
                  call wf%single_batch_reader(batch_k, wf%g_bdck_t, g_bdck, wf%g_lbkc_t, g_lbkc)
                  g_bdck_p => g_bdck
                  g_lbkc_p => g_lbkc
!
!                 Don't read X in the first iteration - X_abdi file empty
                  if (i_batch .gt. 1 .or. j_batch .gt. 1 .or. k_batch .gt. 1) then
                     call wf%single_batch_reader(batch_k, wf%X_abdi, X_abdk)
                     X_abdk_p => X_abdk
                  endif
!
                  call wf%double_batch_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci)
                  g_lkci_p => g_lkci
!
                  call wf%double_batch_reader(batch_i, batch_k, wf%g_ljck_t, g_lick)
                  g_lick_p => g_lick
!
                  call wf%double_batch_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
                  g_lkcj_p => g_lkcj
!
                  call wf%double_batch_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck)
                  g_ljck_p => g_ljck
!
               else if (k_batch .eq. i_batch) then
!
                  g_bdck_p => g_bdci
                  g_lbkc_p => g_lbic
!
                  X_abdk_p => X_abdi
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
                     call wf%double_batch_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci)
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
                  g_bdck_p => g_bdcj
                  g_lbkc_p => g_lbjc
!
                  X_abdk_p => X_abdj
!
                  g_lkci_p => g_ljci
!
                  g_lick_p => g_ljci
!
                  call wf%double_batch_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
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
                        call wf%construct_X_and_Y(i, j, k, t_abc, u_abc,   &
                                                   X_abdi_p(:,:,:,i_rel),  &
                                                   X_abdj_p(:,:,:,j_rel),  &
                                                   X_abdk_p(:,:,:,k_rel),  &
                                                   Y_aikl,                 &
                                                   g_lbic_p(:,:,:,i_rel),  &
                                                   g_lbjc_p(:,:,:,j_rel),  &
                                                   g_lbkc_p(:,:,:,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
               call wf%jacobian_transpose_cc3_write_X(batch_k, X_abdk)
!
            enddo ! batch_k
!
            call wf%jacobian_transpose_cc3_write_X(batch_j, X_abdj)
!
         enddo ! batch_j
!
         call wf%jacobian_transpose_cc3_write_X(batch_i, X_abdi)
!
      enddo ! batch_i
!
!     Close files: 
!
      call disk%close_file(wf%g_bdck_t)
      call disk%close_file(wf%g_ljck_t)
      call disk%close_file(wf%g_lbkc_t)
!
!     Allocate integral arrays and assign pointers.
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_lbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
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
      endif
!
!     Resort X_abdi to X_baid for the final contraction with C^ac_il to sigma_dl
!
      call wf%sort_X_to_baid_and_write()
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
   module subroutine construct_X_and_Y_cc3(wf, i, j, k, t_abc, u_abc, X_abdi, X_abdj, X_abdk,   &
                                          Y_aikl, g_lbic, g_lbjc, g_lbkc)
!!
!!    Constructs the intermediates X_abdi and Y_akil used to compute the contributions to sigma_ai
!!
!!    X_abdi = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_kcjd
!!    Y_akil = sum_cjk (t^cba_ijk + t^acb_ijk - 2 * t^abc_ijk) * g_lbkc
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the loops
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
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdi
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_abdk
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout)   :: Y_aikl
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_lbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_lbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_lbkc
!                       
   end subroutine construct_X_and_Y_cc3
!
!
   module subroutine jacobian_transpose_cc3_write_X_cc3(wf, batch_x, X_abdx)
!!
!!    Write the contributions to the X_abdi intermediate to file in the respective batches
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
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_x%length), intent(in) :: X_abdx
!
      integer :: ioerror
      integer :: x, record
!
      do x = 1, batch_x%length
!
         record = batch_x%first + x -1
         write(wf%X_abdi%unit, rec=record, iostat=ioerror) X_abdx(:,:,:,x)
!
      enddo
!
   end subroutine jacobian_transpose_cc3_write_X_cc3
!
!
   module subroutine sort_X_to_baid_and_write_cc3(wf)
!!
!!    Read in intermediate X_abdi from file, resort to X_baid and write to file again
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abdi
      real(dp), dimension(:,:,:,:), allocatable :: X_baid
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
      call wf%X_baid%init('X_baid','direct','unformatted', dp*(wf%n_v)**2)
      call disk%open_file(wf%X_baid,'write')
!
      do i_batch = 1, batch_i%num_batches
!
         call mem%alloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
!        Read from file
!
         do i = 1, batch_i%length
!
            i_abs = batch_i%first + i - 1
!
            read(wf%X_abdi%unit, rec=i_abs, iostat=ioerror, iomsg=iom) X_abdi(:,:,:,i)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a)') 'Failed to read X_abdi file in sort_X_to_baid_and_write_cc3'
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
!
         enddo
!
!        Sort X_abdi to X_baid
!
         call mem%alloc(X_baid, wf%n_v, wf%n_v, batch_i%length, wf%n_v)
!
         call sort_1234_to_2143(X_abdi, X_baid, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(X_abdi, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
!        Write to file
!
         do d = 1, wf%n_v
            do i = 1, batch_i%length
!
               record  = (d - 1)*wf%n_o + batch_i%first + i - 1
               write(wf%X_baid%unit, rec=record, iostat=ioerror) X_baid(:,:,i,d)
!
            enddo
         enddo
!
         if(ioerror .ne. 0) then
            call output%error_msg('Failed to write X_baid file')
         endif
!
         call mem%dealloc(X_baid, wf%n_v, wf%n_v, batch_i%length, wf%n_v)
!
      enddo ! batch_i
!
      call disk%close_file(wf%X_baid)
!
   end subroutine sort_X_to_baid_and_write_cc3
!
!
end submodule prepare_jacobian_transpose