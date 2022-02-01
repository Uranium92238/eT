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
submodule (cc3_class) omega_cc3
!
!!
!!    Omega submodule
!!
!!    Routines to construct
!!
!!    Omega =  < mu| exp(-T) H exp(T) |R >
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
!!    Written by Rolf H. Myhre, January 2019
!!    Adapted to construct a contravariant representation of omega2 for CC3 which
!!    is in the end transformed back to its covariant form.
!!    by Rolf H. Myhre and Alexander C. Paul, Sep 2020
!!
!!    Directs the construction of the projection vector < mu| exp(-T) H exp(T) |R >
!!    for the current amplitudes of the object wfn
!!
      use array_utilities, only: scale_diagonal, zero_array
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: symmetrize_add_contra_to_packed
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: omega_abij
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct CC3 Omega', pl='normal')
!
      call timer%turn_on()
!
      call zero_array(omega, wf%n_t1 + wf%n_t2)
!
      call wf%ccsd%construct_omega(omega)
!
      call mem%alloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(omega_abij, wf%n_t1**2)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%omega_cc3_a(omega(1:wf%n_t1), omega_abij, t_abij)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call scale_diagonal(half, omega_abij, wf%n_v, wf%n_o)
      call symmetrize_add_contra_to_packed(omega_abij, omega(wf%n_t1+1:wf%n_gs_amplitudes), &
                                           wf%n_v, wf%n_o)
!
      call mem%dealloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2, t_abij)
!!
!!    CC3 Omega A
!!    Written by Rolf H. Myhre, January 2019
!!    Adapted to give a contravariant representation of omega2 due to the
!!    use of a contravariant representation of t3
!!    by Rolf H. Myhre and Alexander C. Paul, Sep 2020
!!
!!    t_3 = -< mu3|{U,T2}|HF > (epsilon_mu3)^-1
!!
!!    omega_1 += < mu1|[H,T3]|HF >
!!
!!    omega_2 += < mu2|[H,T3]|HF >
!!
!!    The doubles part is returned as
!!    ~Omega^ab_ij = 2 Omega^ab_ij - Omega^ba_ij
!!
      use reordering, only: sort_12_to_21
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)    :: t_abij
!
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck ! Transpose the fock matrix ov block
!
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
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer :: req_single_batch
!
!     Resorting of the Fock-Matrix for easier contractions later
      call mem%alloc(F_ov_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_ck, wf%n_o, wf%n_v)
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
      req_single_batch = req_0 + req_1_eri*wf%n_o + 2*wf%n_v**3*wf%n_o &
                       + 2*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 2*wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 4*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           'omega_cc3_a',              &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%alloc(g_ljci,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         call mem%alloc(g_ibjc,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         else
            call mem%alloc(sorting,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%alloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%alloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
!
         call mem%alloc(g_ljci,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_lkci,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_lkcj,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_licj,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_lick,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_ljck,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%alloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%alloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
!
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_klic,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_kljc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_iljc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_ilkc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_jlkc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
!
         call mem%alloc(g_ibjc,wf%n_v,wf%n_v,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_ibkc,wf%n_v,wf%n_v,batch_i%max_length,batch_i%max_length)
         call mem%alloc(g_jbkc,wf%n_v,wf%n_v,batch_i%max_length,batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting,wf%n_v,wf%n_o,wf%n_o,batch_i%max_length)
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
         call wf%setup_vvov(g_dbic, g_dbic_p, sorting, batch_i)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(wf%eri_t1, g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            call wf%setup_ooov(g_jlic, g_jlic_p, sorting, batch_j, batch_i)
!
            call wf%setup_ovov(wf%eri_t1, g_ibjc, g_ibjc_p, sorting, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(wf%eri_t1, g_bdcj, g_bdcj_p, sorting, batch_j)
!
               call wf%setup_vvov(g_dbjc, g_dbjc_p, sorting, batch_j)
!
               call wf%setup_oovo(wf%eri_t1, g_licj, g_licj_p, sorting, batch_i, batch_j)
!
               call wf%setup_ooov(g_iljc, g_iljc_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
!
               call wf%point_vvvo(g_dbjc_p, g_dbic, batch_j%length)
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
                  call wf%setup_vvvo(wf%eri_t1, g_bdck, g_bdck_p, sorting, batch_k)
!
                  call wf%setup_vvov(g_dbkc, g_dbkc_p, sorting, batch_k)
!
                  call wf%setup_oovo(wf%eri_t1, g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ooov(g_ilkc, g_ilkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ooov(g_jlkc, g_jlkc_p, sorting, batch_j, batch_k)
                  call wf%setup_ooov(g_klic, g_klic_p, sorting, batch_k, batch_i)
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ovov(wf%eri_t1, g_ibkc, g_ibkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ovov(wf%eri_t1, g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
!
               else if (k_batch .eq. i_batch) then ! k_batch = j_batch = i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbic, batch_k%length)
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
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
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
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
                        call wf%construct_W(i, j, k, sorting, u_abc, t_abij, &
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
                        call construct_contravariant_t3(u_abc, sorting, wf%n_v)
!
                        call wf%divide_by_orbital_differences(i, j, k, u_abc)
!
                        call wf%omega_cc3_contractions(i, j, k, u_abc, sorting,   &
                                                       omega1, omega2, F_ov_ck,   &
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
                     enddo ! k
                  enddo ! j
               enddo ! i
            enddo ! k_batch
         enddo ! j_batch
      enddo ! i_batch
!
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%dealloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%dealloc(g_ljci,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         call mem%dealloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         call mem%dealloc(g_ibjc,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         else
            call mem%dealloc(sorting,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         end if
!
      else
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%dealloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%dealloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
!
         call mem%dealloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%dealloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         call mem%dealloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
!
         call mem%dealloc(g_ljci,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_lkci,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_lkcj,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_licj,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_lick,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_ljck,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
!
         call mem%dealloc(g_jlic,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_klic,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_kljc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_iljc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_ilkc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_jlkc,wf%n_v,wf%n_o,batch_i%max_length,batch_i%max_length)
!
         call mem%dealloc(g_ibjc,wf%n_v,wf%n_v,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_ibkc,wf%n_v,wf%n_v,batch_i%max_length,batch_i%max_length)
         call mem%dealloc(g_jbkc,wf%n_v,wf%n_v,batch_i%max_length,batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting,wf%n_v,wf%n_v,wf%n_v,batch_i%max_length)
         else
            call mem%dealloc(sorting,wf%n_v,wf%n_o,wf%n_o,batch_i%max_length)
         end if
!
      endif
!
      call mem%dealloc(F_ov_ck,wf%n_v,wf%n_o)
!
      call mem%batch_finalize()
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_contractions_cc3(wf, i, j, k, u_abc, v_abc, &
                                                omega1, omega2, F_ov_ck,   &
                                                g_dbic, g_dbjc, g_dbkc,    &
                                                g_jlic, g_klic, g_kljc,    &
                                                g_iljc, g_ilkc, g_jlkc,    &
                                                g_ibjc, g_ibkc, g_jbkc)
!!
!!    omega cc3 contractions
!!
!!    Calculate the triples contribution to omega1 and omega2
!!
!!    Written by Rolf H. Myhre and Alexander C. Paul, August 2020
!!    Adapted to give a contravariant representation of rho2 due to the
!!    use of a contravariant representation of t3
!!    by Rolf H. Myhre and Alexander C. Paul, Okt 2020
!!
!!    u_abc = (1 - 1/2 delta_ij - 1/2 delta_jk)
!             (4t_abc - 2t_bac - 2t_cba - 2t_acb + t_bca + t_cab)
!!
!!    omega^a_i += sum_bcjk u^abc_ijk g_jbkc
!!
!!    ~omega^ab_ij += sum_ck  u^abc_ijk F_kc
!!    ~omega^ab_lj += sum_cik u^abc_ijk g_ilkc
!!    ~omega^ad_ij += sum_bck u^abc_ijk g_dbkc
!!
!!    with
!!    ~omega^ab_ij = 2 omega^ab_ij - omega^ba_ij
!!
      use reordering, only: sort_123_to_312, sort_123_to_213
!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)         :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)         :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                 :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: F_ov_ck
!
!     (db|kc) ordered bcd,k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)            :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)            :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)            :: g_dbkc
!
!     (jl|kc) ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                    :: g_jlkc
!
!     (jb|kc) ordered bc,jk
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                    :: g_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                    :: g_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                    :: g_jbkc
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
!     ~omega_adki += sum_bc u_abc g_dbjc
!     ~omega_ablk -= sum_c u_abc g_jlic
!
      call wf%omega2_cc3_permutation(k, i, u_abc, g_dbjc, g_jlic, omega2)
!
!     omega_ak += sum_cb u_abc g_ibjc
      call wf%omega1_cc3_permutation(k, u_abc, g_ibjc, omega1, factor_ij)
!
!     ~omega_abjk += sum_c u_abc F_ic
      call wf%omega2_fock_cc3_permutation(u_abc, F_ov_ck(:,i), omega2(:,:,j,k), factor_jk)
!
!     abc -> cab
!     cab -> bca
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     ~omega_adjk += sum_bc v_abc g_dbic
!     ~omega_ablj -= sum_c u_abc g_ilkc
!
      call wf%omega2_cc3_permutation(j, k, v_abc, g_dbic, g_ilkc, omega2)
!
!     ~omega_abij += sum_c u_abc F_kc
      call wf%omega2_fock_cc3_permutation(v_abc, F_ov_ck(:,k), omega2(:,:,i,j), factor_ij)
!
!     bca -> cab
!     abc -> bca
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     ~omega_adij += sum_bc u_abc g_dbkc
!     ~omega_abli -= sum_c u_abc g_kljc
!
      call wf%omega2_cc3_permutation(i, j, u_abc, g_dbkc, g_kljc, omega2)
!
!     omega_ai += sum_bc u_abc g_jbkc
      call wf%omega1_cc3_permutation(i, u_abc, g_jbkc, omega1, factor_jk)
!
      if (i .ne. j .and. j .ne. k) then
!
!        bac -> cab
!        acb -> bca
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        ~omega_adik += sum_bc v_abc g_dbjc
!        ~omega_abli -= sum_c v_abc g_jlkc
!
         call wf%omega2_cc3_permutation(i, k, v_abc, g_dbjc, g_jlkc, omega2)
!
!        cba -> cab
!        bac -> bca
         call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        ~omega_adji += sum_bc u_abc g_dbkc
!        ~omega_ablj -= sum_c u_abc g_klic
!
         call wf%omega2_cc3_permutation(j, i, u_abc, g_dbkc, g_klic, omega2)
!
!        omega_aj += sum_cb u_abc g_ibkc
         call wf%omega1_cc3_permutation(j, u_abc, g_ibkc, omega1, one)
!
         call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        ~omega_adkj += sum_bc v_abc g_dbic
!        ~omega_ablk -= sum_c v_abc g_iljc
!
         call wf%omega2_cc3_permutation(k, j, v_abc, g_dbic, g_iljc, omega2)
!
!        ~omega_abik += sum_c v_abc F_jc
         call wf%omega2_fock_cc3_permutation(v_abc, F_ov_ck(:,j), omega2(:,:,i,k), one)
!
      end if
!
   end subroutine omega_cc3_contractions_cc3
!
!
   module subroutine omega2_cc3_permutation_cc3(wf, o1, o2, t3, g_vvv_o, g_vo_oo, omega2)
!!
!!    Omega doubles CC3 permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to CC3 omega2 vector:
!!    These contributions are:
!!    ~omega^ab_lj += sum_cik u^abc_ijk g_ilkc
!!    ~omega^ad_ij += sum_bck u^abc_ijk g_dbkc
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: o1, o2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: t3
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_vvv_o
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_vo_oo
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
!
      call dgemm('T','N',           &
                 wf%n_v,            &
                 wf%n_v,            &
                 wf%n_v**2,         &
                 one,               &
                 t3,                &
                 wf%n_v**2,         &
                 g_vvv_o,           &
                 wf%n_v**2,         &
                 one,               &
                 omega2(:,:,o1,o2), &
                 wf%n_v)
!
      call dgemm('T','N',          &
                 wf%n_v**2,        &
                 wf%n_o,           &
                 wf%n_v,           &
                 -one,             &
                 t3,               &
                 wf%n_v,           &
                 g_vo_oo,          &
                 wf%n_v,           &
                 one,              &
                 omega2(:,:,:,o1), &
                 wf%n_v**2)
!
   end subroutine omega2_cc3_permutation_cc3
!
!
   module subroutine omega2_fock_cc3_permutation_cc3(wf, t3, F_ov, omega2, factor)
!!
!!    Omega doubles Fock CC3 permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to CC3 omega2 vector, e.g.:
!!    ~omega_abjk += sum_c u_abc F_ic
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: t3
!
      real(dp), dimension(wf%n_v), intent(in) :: F_ov
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: omega2
!
      real(dp), intent(in) :: factor
!
      call dgemv('T',       &
                 wf%n_v,    &
                 wf%n_v**2, &
                 factor,    &
                 t3,        & ! u_c_ab
                 wf%n_v,    &
                 F_ov, 1,   & ! F_c_i
                 one,       &
                 omega2, 1) ! omega_ab,jk
!
   end subroutine omega2_fock_cc3_permutation_cc3
!
!
   module subroutine omega1_cc3_permutation_cc3(wf, o1, t3, g_vv_oo, omega1, factor)
!!
!!    Omega Singles CC3 permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Constructs one permutation to CC3 omega1 vector, e.g.:
!!    omega_ak += sum_cb u_abc g_ibjc
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: o1
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: t3
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_vv_oo
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
!
      real(dp), intent(in) :: factor
!
      call dgemv('T',          &
                 wf%n_v**2,    &
                 wf%n_v,       &
                 factor,       &
                 t3,           & ! u_bc_a
                 wf%n_v**2,    &
                 g_vv_oo, 1,   & ! g_bc_ij
                 one,          &
                 omega1(:,o1), 1)
!
   end subroutine omega1_cc3_permutation_cc3
!
!
end submodule
