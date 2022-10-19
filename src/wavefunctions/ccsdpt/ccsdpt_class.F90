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
module ccsdpt_class
!
!!
!!    CCSD(T) class module
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2018-2019
!!
!!    Computes the CCSD(T) energy correction to the CCSD ground state energy
!!
!!       E(CCSD(T)) = < tbar(CCSD)| [U, T_3] |CC >
!!
!
   use triples_class, only: triples
!
   use parameters
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use direct_stream_file_class, only : direct_stream_file
!
   implicit none
!
   type, extends(triples) :: ccsdpt
!
!     Integral files
!
      type(direct_stream_file) :: g_bdck
      type(direct_stream_file) :: g_ljck
      type(direct_stream_file) :: g_jbkc
!
      real(dp) :: ccsdpt_energy_correction
!
   contains
!
      procedure :: calculate_energy_correction => calculate_energy_correction_ccsdpt
      procedure :: t_dot_x                     => t_dot_x_ccsdpt
      procedure :: construct_X                 => construct_X_ccsdpt
!
      procedure :: estimate_mem_integral_setup => estimate_mem_integral_setup_ccsdpt
!
      procedure :: print_gs_summary            => print_gs_summary_ccsdpt
!
      procedure :: initialize                  => initialize_ccsdpt
!
   end type ccsdpt
!
!
contains
!
!
   subroutine initialize_ccsdpt(wf, template_wf)
!!
!!    Initialize
!!    Written by Rolf H. Myhre, 2018
!!
      use wavefunction_class, only : wavefunction
!
      implicit none
!
      class(ccsdpt), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'ccsd(t)'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%n_t1 = wf%n_o*wf%n_v
      wf%n_t2 = wf%n_o*wf%n_v*(wf%n_o*wf%n_v + 1)/2
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
      wf%need_g_abcd     = .true.
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
   end subroutine initialize_ccsdpt
!
!
   subroutine calculate_energy_correction_ccsdpt(wf)
!!
!!    Calculate Energy correction CCSD(T)
!!    Written by Alexander C. Paul, Rolf H. Myhre, Nov 2019
!!
!!    The energy correction for singles i,j,k is:
!!
!!    E = 1/3 (4W^abc + W^bca + W^cab)(X^abc - X^cba)/epsilon^abc
!!      = 1/3 (4t^abc + t^bca + t^cab)(X^abc - X^cba)
!!
!!    Where 2 intermediates are used:
!!       W_abc = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il (lj|ck))
!!
!!       X^abc = W^abc + t_ai(jb|kc) + t_bj(ia|kc) + t_ck(ia|jb)
!!
      use batching_index_class, only: batching_index
      use reordering, only: squareup_and_sort_1234_to_1324
      use eri_tool_class, only: eri_tool
      use eri_adapter_class, only: eri_adapter
!
      implicit none
!
      class(ccsdpt) :: wf
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: x_abc
!
      real(dp), dimension(:,:,:,:), allocatable :: sorted
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
      type(batching_index), allocatable :: batch_i, batch_j, batch_k
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_i, req_1, req_2, req_3, req_1_eri
      integer :: req_single_batch
!
      type(eri_adapter), allocatable :: eri
!
      type(timings) :: ccsdpt_timer
!
      ccsdpt_timer = timings('CCSD(T) correction', pl='m')
      call ccsdpt_timer%turn_on()
!
      eri = eri_adapter(eri_tool(wf%L_mo), wf%n_o, wf%n_v)
!
      wf%ccsdpt_energy_correction = zero
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2,t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + 2*wf%n_v**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + wf%n_v**3*wf%n_o &
                       + wf%n_o**3*wf%n_v + wf%n_v**2*wf%n_o**2
!
      req_1 = wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 2*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,             &
                           req_0, req_i, req_1, req_1,            &
                           req_2, req_2, req_2, req_3,            &
                           'calculate_energy_correction_ccsdpt',  &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(x_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorted, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorted, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
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
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorted, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%alloc(sorted, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
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
         call wf%setup_vvvo(eri, g_bdci, g_bdci_p, sorted, batch_i)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(eri, g_ljci, g_ljci_p, sorted, batch_j, batch_i)
!
            call wf%setup_ovov(eri, g_ibjc, g_ibjc_p, sorted, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(eri, g_bdcj, g_bdcj_p, sorted, batch_j)
!
               call wf%setup_oovo(eri, g_licj, g_licj_p, sorted, batch_i, batch_j)
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
                  call wf%setup_vvvo(eri, g_bdck, g_bdck_p, sorted, batch_k)
!
                  call wf%setup_oovo(eri, g_lick, g_lick_p, sorted, batch_i, batch_k)
                  call wf%setup_oovo(eri, g_ljck, g_ljck_p, sorted, batch_j, batch_k)
                  call wf%setup_oovo(eri, g_lkci, g_lkci_p, sorted, batch_k, batch_i)
                  call wf%setup_oovo(eri, g_lkcj, g_lkcj_p, sorted, batch_k, batch_j)
!
                  call wf%setup_ovov(eri, g_ibkc, g_ibkc_p, sorted, batch_i, batch_k)
                  call wf%setup_ovov(eri, g_jbkc, g_jbkc_p, sorted, batch_j, batch_k)
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
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
                  call wf%point_vvoo(g_jbkc_p, g_ibjc, batch_j%length, batch_k%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
!
                  call wf%setup_oovo(eri, g_lkcj, g_lkcj_p, sorted, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
                  call wf%setup_ovov(eri, g_jbkc, g_jbkc_p, sorted, batch_j, batch_k)
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
                        if (k .eq. i) then ! k == j == i
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct w^abc_ijk for given i, j, k
!
                        call wf%construct_W(i, j, k, sorted, t_abc, t_abij, &
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
                        call wf%construct_X(i, j, k, wf%t1, t_abc, x_abc, &
                                            g_ibjc_p(:,:,i_rel,j_rel),    &
                                            g_ibkc_p(:,:,i_rel,k_rel),    &
                                            g_jbkc_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, t_abc)
!
                        call wf%t_dot_x(t_abc, x_abc, sorted)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      wf%ccsdpt_energy_correction = third*wf%ccsdpt_energy_correction
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(x_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorted, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorted, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
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
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorted, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%dealloc(sorted, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      call mem%batch_finalize()
!
      call ccsdpt_timer%turn_off()
!
   end subroutine calculate_energy_correction_ccsdpt
!
!
   subroutine construct_X_ccsdpt(wf, i, j, k, t1, w_abc, x_abc, g_ibjc, g_ibkc, g_jbkc)
!!
!!    Construct X
!!    written by Alexander C. Paul and Rolf H. Myhre, Nov 2019
!!
!!    Constructs the following intermediate for single i,j,k
!!
!!       X_abc - X_cba
!!
!!    where:
!!       X_abc = W_abc + t_ai(jb|kc) + t_bj(ia|kc) + t_ck(ia|jb)
!!    and
!!       W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il (lj|ck))
!!
      implicit none
!
      class(ccsdpt), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)    :: w_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: x_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: g_jbkc
!
      integer, intent(in) :: i, j , k
!
      real(dp) :: abc_factor, cba_factor, acb_factor, bac_factor
!
      integer :: a,b,c
!
!     Have to consider all permutations. Depending on i,j,k some of the
!     permutations are equal leading to the following prefactors
!
!
      if(i .eq. j) then
         abc_factor = two
         acb_factor = one
         bac_factor = zero
         cba_factor = one
      else if(j .eq. k) then
         abc_factor = two
         acb_factor = zero
         bac_factor = one
         cba_factor = one
      else
         abc_factor = six
         acb_factor = two
         bac_factor = two
         cba_factor = two
      end if
!
!     Contribution of W to (x_abc - x_cba)
!     ------------------------------------
!
!$omp parallel do schedule(static) private(a,b,c)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               x_abc(a,b,c) = abc_factor*W_abc(a,b,c) - cba_factor*W_abc(c,b,a)
!
               if(i .ne. j) then
!
                  x_abc(a,b,c) = x_abc(a,b,c) - bac_factor*W_abc(b,a,c)
!
               end if
!
               if(j .ne. k) then
!
                  x_abc(a,b,c) = x_abc(a,b,c) - acb_factor*W_abc(a,c,b)
!
               end if
!
            end do
         end do
      end do
!$omp end parallel do
!
!     Contribution of the outer products
!     ----------------------------------
!
!$omp parallel do schedule(static) private(a,b,c)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               x_abc(a,b,c) = x_abc(a,b,c)                    &
                            + abc_factor*(t1(a,i)*g_jbkc(b,c) &
                                        + t1(b,j)*g_ibkc(a,c) &
                                        + t1(c,k)*g_ibjc(a,b))&
!
                              - cba_factor*(t1(c,i)*g_jbkc(b,a) &
                                          + t1(b,j)*g_ibkc(c,a) &
                                          + t1(a,k)*g_ibjc(c,b))
!
               if(i .ne. j) then
!
                  x_abc(a,b,c) = x_abc(a,b,c)                    &
                               - bac_factor*(t1(b,i)*g_jbkc(a,c) &
                                           + t1(a,j)*g_ibkc(b,c) &
                                           + t1(c,k)*g_ibjc(b,a))
!
               end if
!
               if(j .ne. k) then
!
                  x_abc(a,b,c) = x_abc(a,b,c)                    &
                               - acb_factor*(t1(a,i)*g_jbkc(c,b) &
                                           + t1(c,j)*g_ibkc(a,b) &
                                           + t1(b,k)*g_ibjc(a,c))
!
               end if
!
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine construct_x_ccsdpt
!
!
   subroutine t_dot_x_ccsdpt(wf, t_abc, x_abc, u_abc)
!!
!!    T dot X
!!    written by Alexander C. Paul and Rolf H. Myhre, Nov 2019
!!
!!    Construct linear combination of t and dot with
!!    linear combination of V (constructed in construct_V)
!!
!!       (4t^abc + t^bca + t^cab)*(X_abc - X_cba)
!!
!!    The result has to be scaled by 1/3 to obtain the final energy correction
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccsdpt), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)    :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)    :: x_abc
!
      real(dp) :: ddot
!
      integer :: a,b,c
!
      call zero_array(u_abc, wf%n_v**3)
!
!$omp parallel do schedule(static) private(a,b,c)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               u_abc(a,b,c) = u_abc(a,b,c) + four*t_abc(a,b,c) &
                            + t_abc(b,c,a) + t_abc(c,a,b)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%ccsdpt_energy_correction = wf%ccsdpt_energy_correction &
                                  + ddot(wf%n_v**3, u_abc, 1, x_abc, 1)
!
   end subroutine t_dot_x_ccsdpt
!
!
   subroutine print_gs_summary_ccsdpt(wf)
!!
!!    Print ground state summary
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(ccsdpt), intent(inout) :: wf
!
      call wf%calculate_energy_correction()
!
      call output%printf('m', '- Ground state summary:', fs='(/t3,a)')
!
      call output%printf('m', 'Final ground state CCSD energy (a.u.):    (f18.12)', &
                         reals=[wf%energy], fs='(/t6,a)')
!
      call output%printf('m', 'Correlation energy CCSD (a.u.):           (f18.12)', &
                         reals=[wf%correlation_energy], fs='(/t6,a)')
!
      wf%energy = wf%energy + wf%ccsdpt_energy_correction
!
      call output%printf('m', 'Final ground state CCSD(T) energy (a.u.): (f18.12)', &
                         reals=[wf%energy], fs='(/t6,a)')
!
      call output%printf('m', 'CCSD(T) energy correction (a.u.):         (f18.12)', &
                         reals=[wf%ccsdpt_energy_correction], fs='(/t6,a)')
!
      call wf%print_dominant_amplitudes()
!
   end subroutine print_gs_summary_ccsdpt
!
!
   subroutine estimate_mem_integral_setup_ccsdpt(wf, req0, req1)
!!
!!    Estimate memory integral setup
!!    Written by Alexander C. Paul, Dec 2020
!!
!!    Estimate maximum memory needed for ccsd(t) integral setup
!!
!!    get_eri_mo_mem returns the memory needed to construct the requested integral
!!    The dimensions sent in specify if an index is batched (1) or of
!!    full dimension (n_o/n_v)
!!    The memory estimate for the first and second pair of indices
!!    is added to the integers req*.
!!
      implicit none
!
      class(ccsdpt), intent(in) :: wf
!
      integer, intent(out) :: req0, req1
!
      integer, dimension(2) :: req_vvvo, req_ovov, req_oovo
!
      req_vvvo = wf%eri_t1%get_memory_estimate('vvvo', wf%n_v, wf%n_v, wf%n_v, 1)
      req_ovov = wf%eri_t1%get_memory_estimate('ovov', 1, wf%n_v, 1, wf%n_v)
      req_oovo = wf%eri_t1%get_memory_estimate('oovo', wf%n_o, 1, wf%n_v, 1)
!
      req0 = req_vvvo(1)
      req1 = max(req_vvvo(2), sum(req_ovov), sum(req_oovo))
!
   end subroutine estimate_mem_integral_setup_ccsdpt
!
!
end module ccsdpt_class
