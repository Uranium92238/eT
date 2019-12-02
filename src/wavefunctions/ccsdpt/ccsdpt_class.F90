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
   use memory_manager_class, only: mem
   use batching_index_class, only : batching_index
   use global_out, only: output
   use timings_class, only: timings
   use direct_file_class, only : direct_file
   use io_utilities, only : single_record_reader, compound_record_reader 
   use io_utilities, only : single_record_writer, compound_record_writer
   use array_utilities, only: entrywise_product, zero_array
   use reordering
!
   implicit none
!
   type, extends(triples) :: ccsdpt
!
!     Integral files
!
      type(direct_file) :: g_bdck
      type(direct_file) :: g_ljck
      type(direct_file) :: g_jbkc
!
      real(dp) :: ccsdpt_energy_correction
!
   contains
!
      procedure :: prepare_integrals           => prepare_integrals_ccsdpt
!
      procedure :: calculate_energy_correction => calculate_energy_correction_ccsdpt
      procedure :: construct_V                 => construct_V_ccsdpt
      procedure :: t_dot_v                     => t_dot_v_ccsdpt
!
      procedure :: print_gs_summary            => print_gs_summary_ccsdpt
!
      procedure :: cleanup                     => cleanup_ccsdpt
!
   end type ccsdpt
!
   interface ccsdpt
!
      procedure :: new_ccsdpt
!
   end interface ccsdpt
!
contains
!
!
   function new_ccsdpt(system) result(wf)
!!
!!    New CCSD(T)
!!    Written by Rolf H. Myhre, 2018
!!
      use molecular_system_class, only: molecular_system
!
      implicit none
!
      type(ccsdpt) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%name_ = 'ccsd(t)'
!
      call wf%general_cc_preparations(system)
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
      wf%n_t2 = (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
      wf%need_g_abcd     = .true.
!
      call wf%initialize_fock()
!
   end function new_ccsdpt
!   
!
   subroutine cleanup_ccsdpt(wf)
!!
!!    Cleanup
!!    written by Alexander C. Paul and Rolf H. Myhre, Nov 2019
!!
      implicit none
!
      class(ccsdpt) :: wf
!
      call wf%ccs%cleanup()
!
      call wf%g_jbkc%delete_
      call wf%g_bdck%delete_
      call wf%g_ljck%delete_
!
   end subroutine cleanup_ccsdpt
!
!
   subroutine prepare_integrals_ccsdpt(wf)
!!
!!    Prepate integrals for CCSD(T)
!!    written by Alexander C. Paul and Rolf H. Myhre, Nov 2019
!!
!!    Construct and reorder integrals files for CCSD(T)
!!    (using t1-transformed integrals even though it's usually not done)
!!
!!    (bd|ck) ordered dbc,k
!!    (lj|ck) ordered lc,jk
!!    (jb|kc) ordered bc,jk
!!
!!    Algorithm adapted from: Rendell, Lee, Komornicki, Chem. Phys. Lett., 1991, 178
!!
      implicit none
!
      class(ccsdpt) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      type(batching_index), allocatable :: batch_k
!
      integer :: req_0, req_k
      integer :: k_batch
!
!     (jb|kc) ! stored as bcj#k
!
      req_0 = wf%integrals%n_J*wf%n_o*wf%n_v
      req_k = max(wf%integrals%n_J*wf%n_v + wf%n_o*wf%n_v**2, 2*wf%n_o*wf%n_v**2)
!
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_k, req_0, req_k)
      call batch_k%determine_limits(1)
!
      call mem%alloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
      wf%g_jbkc = direct_file('g_jbkc', wf%n_v**2)
      call wf%g_jbkc%open_('write')
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o , wf%n_v , batch_k%length , wf%n_v)
!
         call wf%integrals%construct_g_pqrs_mo(g_pqrs,    &
                                               1, wf%n_o, &
                                               wf%n_o+1, wf%n_o+wf%n_v,     &
                                               batch_k%first, batch_k%last, &
                                               wf%n_o+1, wf%n_o+wf%n_v)
!
         call sort_1234_to_2413(g_pqrs , h_pqrs , wf%n_o , wf%n_v , batch_k%length , wf%n_v) ! sort to bcjk
!
         call compound_record_writer(wf%n_o, batch_k, wf%g_jbkc, h_pqrs)
!
         call mem%dealloc(g_pqrs, wf%n_o , wf%n_v , batch_k%length , wf%n_v)
!
      enddo
!
      call batch_k%determine_limits(1)
      call mem%dealloc(h_pqrs, wf%n_v , wf%n_v , wf%n_o , batch_k%length)
!
      call wf%g_jbkc%close_()
!
!     (bd|ck) ordered dbc,k
!
      req_0 = wf%integrals%n_J*wf%n_v**2
      req_k = 2*wf%n_v**3 + wf%integrals%n_J*wf%n_v
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      wf%g_bdck = direct_file('g_bdck',wf%n_v**3)
      call wf%g_bdck%open_('write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call wf%integrals%construct_g_pqrs_mo(g_pqrs,    &
                                               wf%n_o+1, wf%n_o+wf%n_v,     &
                                               wf%n_o+1, wf%n_o+wf%n_v,     &
                                               wf%n_o+1, wf%n_o+wf%n_v,     &
                                               batch_k%first, batch_k%last)
!
         call sort_1234_to_2134(g_pqrs,h_pqrs,wf%n_v,wf%n_v,wf%n_v,batch_k%length)
!
         call single_record_writer(batch_k, wf%g_bdck, h_pqrs)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo
!
      call wf%g_bdck%close_()
!
!     (lj|ck) ordered lcjk
!
      req_0 = wf%integrals%n_J*wf%n_o**2
      req_k = 2*wf%n_o**2*wf%n_v + wf%integrals%n_J*wf%n_v
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      wf%g_ljck = direct_file('g_ljck',wf%n_v*wf%n_o)
      call wf%g_ljck%open_('write')
!
      do k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
         call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o ,batch_k%length)
!
         call wf%integrals%construct_g_pqrs_mo(g_pqrs,    &
                                               1, wf%n_o, &
                                               1, wf%n_o, &
                                               wf%n_o+1, wf%n_v+wf%n_o, &
                                               batch_k%first, batch_k%last)
!
         call sort_1234_to_1324(g_pqrs,h_pqrs,wf%n_o,wf%n_o,wf%n_v,batch_k%length)
!
         call compound_record_writer(wf%n_o, batch_k, wf%g_ljck, h_pqrs)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_k%length)
!
      enddo
!
      call wf%g_ljck%close_()
!
   end subroutine prepare_integrals_ccsdpt
!
!
   subroutine calculate_energy_correction_ccsdpt(wf)
!!
!!    Calculate Energy correction CCSD(T)
!!    Written by Alexander C. Paul, Rolf H. Myhre, Nov 2019
!!
!!    The energy correction for singles i,j,k is:
!!
!!    E = 1/3 (4W^abc + W^bca + W^cab)(V^abc - V^cba)/epsilon^abc
!!      = 1/3 (4t^abc + t^bca + t^cab)(V^abc - V^cba)
!!
!!    Where 2 intermediates are used:
!!       W_abc = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il (lj|ck))
!!
!!       V^abc = W^abc + t_ai(jb|kc) + t_bj(ia|kc) + t_ck(ia|jb)
!!
      implicit none
!
      class(ccsdpt) :: wf
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
      real(dp), dimension(:,:,:), allocatable :: v_abc
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
      integer :: req_0, req_1, req_2, req_3
!
      type(timings) :: ccsdpt_timer
!
      ccsdpt_timer = timings('CCSD(T) correction', pl='m')
      call ccsdpt_timer%turn_on()
!
      wf%ccsdpt_energy_correction = zero
!
!     prepare integral files
!
      call wf%prepare_integrals
!
      call mem%alloc(t_abij,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2,t_abij,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
      req_0 = 3*wf%n_v**3
      req_1 = wf%n_v**3
      req_2 = wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, buffer_size = zero)
!
!     Allocate integral arrays and assign pointers.
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
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
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Arrays for the triples amplitudes
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call wf%g_bdck%open_('read')
      call wf%g_ljck%open_('read')
      call wf%g_jbkc%open_('read')
!
!     Loop over the batches in i,j,k
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck, g_bdci)
         g_bdci_p => g_bdci
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck, g_ljci)
            call compound_record_reader(batch_i, batch_j, wf%g_jbkc, g_ibjc)
!
            g_ljci_p => g_ljci
            g_ibjc_p => g_ibjc
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck, g_bdcj)
               g_bdcj_p => g_bdcj
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck, g_licj)
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
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call single_record_reader(batch_k, wf%g_bdck, g_bdck)
                  g_bdck_p => g_bdck
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck, g_lkci)
                  g_lkci_p => g_lkci
!
                  call compound_record_reader(batch_i, batch_k, wf%g_ljck, g_lick, &
                                                                wf%g_jbkc, g_ibkc)
                  g_lick_p => g_lick
                  g_ibkc_p => g_ibkc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck, g_lkcj)
                  g_lkcj_p => g_lkcj
!
                  call compound_record_reader(batch_j, batch_k, wf%g_ljck, g_ljck, &
                                                                wf%g_jbkc, g_jbkc)
                  g_ljck_p => g_ljck
                  g_jbkc_p => g_jbkc
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_ljci
                  g_ibkc_p => g_ibjc
!
                  g_lkcj_p => g_ljci
                  g_ljck_p => g_ljci
                  g_jbkc_p => g_ibjc
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_licj
                  g_ibkc_p => g_ibjc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck, g_lkcj, &
                                                                wf%g_jbkc, g_jbkc) ! g_jbkc = g_kbjc
                  g_lkcj_p => g_lkcj
                  g_ljck_p => g_lkcj
                  g_jbkc_p => g_jbkc
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
!                       Construct w^abc_ijk for given i, j, k
!
                        call wf%construct_W(i, j, k, t_abc, u_abc, t_abij, &
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
                        call wf%construct_V(i, j, k, wf%t1, t_abc, v_abc, &
                                            g_ibjc_p(:,:,i_rel,j_rel),    &
                                            g_ibkc_p(:,:,i_rel,k_rel),    &
                                            g_jbkc_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, t_abc)
!
                        call wf%t_dot_v(t_abc, v_abc, u_abc)
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
!     Close files
!
      call wf%g_bdck%close_()
      call wf%g_ljck%close_()
      call wf%g_jbkc%close_()
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
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
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
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
      call ccsdpt_timer%turn_off()
!
   end subroutine calculate_energy_correction_ccsdpt
!
!
   subroutine construct_V_ccsdpt(wf, i, j, k, t1, W_abc, V_abc, g_ibjc, g_ibkc, g_jbkc)
!!
!!    Construct V intermediate
!!    written by Alexander C. Paul and Rolf H. Myhre, Nov 2019
!!
!!    Constructs the intermediate the following intermediate for single i,j,k
!!
!!       V_abc - V_cba
!!
!!    where:
!!       V_abc = W_abc + t_ai(jb|kc) + t_bj(ia|kc) + t_ck(ia|jb)
!!    and
!!       W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il (lj|ck))
!!
      implicit none
!
      class(ccsdpt), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)    :: W_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: V_abc
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
!     Contribution of W to (V_abc - V_cba)
!     ------------------------------------
!
!$omp parallel do schedule(static) private(a,b,c)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!              
               V_abc(a,b,c) = abc_factor*W_abc(a,b,c) - cba_factor*W_abc(c,b,a) 
!              
               if(i .ne. j) then
!              
                  V_abc(a,b,c) = V_abc(a,b,c) - bac_factor*W_abc(b,a,c) 
!
               end if
!              
               if(j .ne. k) then
!              
                  V_abc(a,b,c) = V_abc(a,b,c) - acb_factor*W_abc(a,c,b) 
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
               V_abc(a,b,c) = V_abc(a,b,c)                    &
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
                  V_abc(a,b,c) = V_abc(a,b,c)                  &
                             - bac_factor*(t1(b,i)*g_jbkc(a,c) &
                                         + t1(a,j)*g_ibkc(b,c) &
                                         + t1(c,k)*g_ibjc(b,a))
!
               end if
!              
               if(j .ne. k) then
!              
                  V_abc(a,b,c) = V_abc(a,b,c)                  &
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
   end subroutine construct_V_ccsdpt
!
!
   subroutine t_dot_v_ccsdpt(wf, t_abc, v_abc, u_abc)
!!
!!    T dot V
!!    written by Alexander C. Paul and Rolf H. Myhre, Nov 2019
!!
!!    Construct linear combination of t and dot with 
!!    linear combination of V (constructed in construct_V)
!!
!!       (4t^abc + t^bca + t^cab)*(V_abc - V_cba)
!!
!!    The result has to be scaled by 1/3 to obtain the final energy correction
!!
      implicit none
!
      class(ccsdpt), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)    :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)    :: v_abc
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
               u_abc(a,b,c) = u_abc(a,b,c) + four*t_abc(a,b,c) + t_abc(b,c,a) + t_abc(c,a,b)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%ccsdpt_energy_correction = wf%ccsdpt_energy_correction &
                                  + ddot(wf%n_v**3, u_abc, 1, v_abc, 1)
!
   end subroutine t_dot_v_ccsdpt
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
      real(dp) :: t1_diagnostic
!
      call wf%calculate_energy_correction()
!
      call output%printf('- Ground state summary:', fs='(/t3,a)', pl='m')
!
      call output%printf('Final ground state CCSD energy (a.u.):    (f18.12)', &
                         reals=[wf%energy], fs='(/t6,a)', pl='m')
!
      call output%printf('Correlation energy CCSD (a.u.):           (f18.12)', &
                         reals=[wf%correlation_energy], fs='(/t6,a)', pl='m')
!
      wf%energy = wf%energy + wf%ccsdpt_energy_correction
!
      call output%printf('Final ground state CCSD(T) energy (a.u.): (f18.12)', &
                         reals=[wf%energy], fs='(/t6,a)', pl='m')
!
      call output%printf('CCSD(T) energy correction (a.u.):         (f18.12)', &
                          pl='m', reals=[wf%ccsdpt_energy_correction], fs='(/t6,a)')
!
      call wf%print_dominant_amplitudes()
!
      t1_diagnostic = wf%get_t1_diagnostic() 
      call output%printf('T1 diagnostic (|T1|/sqrt(N_e)): (f14.12)', &
                         reals=[t1_diagnostic], fs='(/t6,a)', pl='m')
!
   end subroutine print_gs_summary_ccsdpt
!
!
end module ccsdpt_class
