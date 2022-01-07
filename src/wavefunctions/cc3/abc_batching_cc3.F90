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
submodule (cc3_class) batching_abc
!
!!
!!    Batching abc submodule
!!
!!    Routines that construct the triples T-amplitudes and
!!    components of the excitation vectors in batches
!!    of the virtual indices a,b,c
!!
!
   implicit none
!
!
contains
!
!
   module subroutine setup_vvvo_abc_cc3(wf, g_vvvo, point, unordered_g_vvvo, &
                                        batch_p, batch_r, eri_c1)
!!
!!    Setup vvvo abc
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Setups the vvvo integral and the pointer for a triples calculation.
!!    The integral is given in batches over two virtual indices
!!
!!    The integral is returned in 2413 order
!!
!!    g_vvvo:           final sorted integral array
!!    point:            pointer to g_vvvo integral
!!    unordered_g_vvvo: help array for reordering
!!    c_ai:             use the c1-transformed integral
!!
      use reordering, only: sort_1234_to_2413
!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_p, batch_r
!
      real(dp), dimension(wf%n_v, wf%n_o, batch_p%length, batch_r%length), &
                                                      target, intent(out)  :: g_vvvo
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_vvvo
      type(eri_adapter), intent(inout), optional  :: eri_c1
!
      if (.not. present(eri_c1)) then
!
         call wf%eri_t1%get('vvvo', unordered_g_vvvo, &
                                 first_p=batch_p%first, last_p=batch_p%get_last(), &
                                 first_r=batch_r%first, last_r=batch_r%get_last())
!
      else
!
         call eri_c1%get('vvvo', unordered_g_vvvo, &
                                 first_p=batch_p%first, last_p=batch_p%get_last(), &
                                 first_r=batch_r%first, last_r=batch_r%get_last())
!
      end if
!
       call sort_1234_to_2413(unordered_g_vvvo, g_vvvo, batch_p%length, wf%n_v, &
                                                        batch_r%length, wf%n_o)
!
      point(1:wf%n_v, 1:wf%n_o, 1:batch_p%length, 1:batch_r%length) => g_vvvo
!
   end subroutine setup_vvvo_abc_cc3
!
!
   module subroutine point_vvvo_abc_cc3(wf, point, g_vvvo, length1, length2)
!!
!!    Point vvvo abc
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Sets a pointer to a vvvo integral
!!    where the two virtual indices can be batched over
!!    NB: The batching indices need to be the last index
!!
!!    point:   pointer to g_vvvo integral
!!    g_vvvo:  Integral array in vovv order
!!    length1: length of the batch for the second to last index
!!    length2: length of the batch for the last index
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: length1, length2
!
      real(dp), dimension(wf%n_v, wf%n_o, length1, length2), target, intent(in)  :: g_vvvo
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_v, 1:wf%n_o, 1:length1, 1:length2) => g_vvvo
!
   end subroutine point_vvvo_abc_cc3
!
!
   module subroutine setup_vvov_abc_cc3(wf, g_vvov, point, unordered_g_vvov, &
                                        batch_q, batch_s)
!!
!!    Setup vvov abc
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Setups the vvov integral and the pointer for a triples calculation.
!!    The integral is given in batches over two virtual indices
!!
!!    The integral is returned in 1324 order
!!
!!    g_vvov:           final sorted integral array
!!    point:            pointer to g_vvov integral
!!    unordered_g_vvov: help array for reordering
!!
      use reordering, only: sort_1234_to_1324
!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_q, batch_s
!
      real(dp), dimension(wf%n_o, wf%n_v, batch_q%length, batch_s%length), &
                                                      target, intent(out)  :: g_vvov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_vvov
!
      call wf%eri_t1%get('vvov', unordered_g_vvov, &
                              first_q=batch_q%first, last_q=batch_q%get_last(), &
                              first_s=batch_s%first, last_s=batch_s%get_last())
!
      call sort_1234_to_1324(unordered_g_vvov, g_vvov, wf%n_v, batch_q%length, &
                                                       wf%n_o, batch_s%length)
!
      point(1:wf%n_v, 1:wf%n_o, 1:batch_q%length, 1:batch_s%length) => g_vvov
!
   end subroutine setup_vvov_abc_cc3
!
!
   module subroutine point_vvov_abc_cc3(wf, point, g_vvov, length1, length2)
!!
!!    Point vvov abc
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Sets a pointer to a vvov integral
!!    where the two virtual indices can be batched over
!!    NB: The batching indices need to be the last index
!!
!!    point:   pointer to g_vvov integral
!!    g_vvov:  Integral array in vovv order
!!    length1: length of the batch for the second to last index
!!    length2: length of the batch for the last index
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: length1, length2
!
      real(dp), dimension(wf%n_o, wf%n_v, length1, length2), target, intent(in)  :: g_vvov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_v, 1:wf%n_o, 1:length1, 1:length2) => g_vvov
!
   end subroutine point_vvov_abc_cc3
!
!
   module subroutine setup_ovov_abc_cc3(wf, g_ovov, point, unordered_g_ovov, &
                                        batch_q, batch_s)
!!
!!    Setup ovov abc
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Setups the ovov integral and the pointer for a triples calculation.
!!    The integral is given in batches over two virtual indices
!!
!!    The integral is returned in 1324 order
!!
!!    g_ovov:           final sorted integral array
!!    point:            pointer to g_vvov integral
!!    unordered_g_ovov: help array for reordering
!!
      use reordering, only: sort_1234_to_1324
!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_q, batch_s
!
      real(dp), dimension(wf%n_o, wf%n_o, batch_q%length, batch_s%length), &
                                                      target, intent(out)  :: g_ovov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_ovov
!
      call wf%eri_t1%get('ovov', unordered_g_ovov, &
                              first_q=batch_q%first, last_q=batch_q%get_last(), &
                              first_s=batch_s%first, last_s=batch_s%get_last())
!
       call sort_1234_to_1324(unordered_g_ovov, g_ovov, wf%n_o, batch_q%length, &
                                                        wf%n_o, batch_s%length)
!
      point(1:wf%n_o, 1:wf%n_o, 1:batch_q%length, 1:batch_s%length) => g_ovov
!
   end subroutine setup_ovov_abc_cc3
!
!
   module subroutine point_ovov_abc_cc3(wf, point, g_ovov, length1, length2)
!!
!!    Point ovov abc
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Sets a pointer to a ovov integral
!!    where the two virtual indices can be batched over
!!    NB: The batching indices need to be the last index
!!
!!    point:   pointer to g_ovov integral
!!    g_ovov:  Integral array in vovv order
!!    length1: length of the batch for the second to last index
!!    length2: length of the batch for the last index
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: length1, length2
!
      real(dp), dimension(wf%n_o, wf%n_o, length1, length2), target, intent(in)  :: g_ovov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_o, 1:wf%n_o, 1:length1, 1:length2) => g_ovov
!
   end subroutine point_ovov_abc_cc3
!
!
   module subroutine setup_oovo_abc_cc3(wf, g_oovo, point, unordered_g_oovo, &
                                        batch_r, eri_c1)
!!
!!    Setup oovo abc
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Setups the oovo integral and the pointer for a triples calculation
!!    The integral is given in batches over the virtual index
!!
!!    The integral is returned in 1243 order:
!!
!!    g_oovo:           final sorted integral array
!!    point:            pointer to g_oovo integral
!!    unordered_g_oovo: help array for reordering
!!    c_ai:             use the c1-transformed integral
!!
      use reordering, only: sort_1234_to_1243
!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_r
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o, batch_r%length), &
                                                target, intent(out)  :: g_oovo
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_oovo
      type(eri_adapter), intent(inout), optional  :: eri_c1
!
      if (.not. present(eri_c1)) then
!
         call wf%eri_t1%get('oovo', unordered_g_oovo, &
                           first_r=batch_r%first,     &
                           last_r=batch_r%get_last())
!
      else
!
         call eri_c1%get('oovo', unordered_g_oovo, &
                         first_r=batch_r%first,    &
                         last_r=batch_r%get_last())
!
      end if
!
      call sort_1234_to_1243(unordered_g_oovo, g_oovo, wf%n_o, wf%n_o, &
                                                       batch_r%length, wf%n_o)
!
      point(1:wf%n_o, 1:wf%n_o, 1:wf%n_o, 1:batch_r%length) => g_oovo
!
   end subroutine setup_oovo_abc_cc3
!
!
   module subroutine point_oovo_abc_cc3(wf, point, g_oovo, length)
!!
!!    Point oovo abc
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Sets a pointer to a oovo integral
!!    where the virtual index can be batched over
!!    NB: The virtual index needs to be the last index
!!
!!    point:  pointer to g_oovo integral
!!    g_oovo: Integral array in ooov order
!!    length: length of the batch
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: length
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o, length), target, intent(in)  :: g_oovo
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_o, 1:wf%n_o, 1:wf%n_o, 1:length) => g_oovo
!
   end subroutine point_oovo_abc_cc3
!
!
   module subroutine setup_ooov_abc_cc3(wf, g_ooov, point, unordered_g_ooov, batch_s)
!!
!!    Setup ooov abc
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Setups the ooov integral and the pointer for a triples calculation
!!    The integral is given in batches over the virtual index
!!
!!    The integral is returned in 1324 order:
!!
!!    g_ooov:           final sorted integral array
!!    point:            pointer to g_ooov integral
!!    unordered_g_ooov: help array for reordering
!!
      use reordering, only: sort_1234_to_2134
!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_s
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o, batch_s%length), &
                                                target, intent(out)  :: g_ooov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_ooov
!
      call wf%eri_t1%get('ooov', unordered_g_ooov, &
                        first_s=batch_s%first,     &
                        last_s=batch_s%get_last())
!
      call sort_1234_to_2134(unordered_g_ooov, g_ooov, &
                             wf%n_o, wf%n_o, wf%n_o, batch_s%length)
!
      point(1:wf%n_o, 1:wf%n_o, 1:wf%n_o, 1:batch_s%length) => g_ooov
!
   end subroutine setup_ooov_abc_cc3
!
!
   module subroutine point_ooov_abc_cc3(wf, point, g_ooov, length)
!!
!!    Point ooov abc
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Sets a pointer to a ooov integral
!!    where the virtual index can be batched over
!!    NB: The virtual index needs to be the last index
!!
!!    point:  pointer to g_ooov integral
!!    g_ooov: Integral array in ooov order
!!    length: length of the batch
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: length
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o, length), target, intent(in)  :: g_ooov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_o, 1:wf%n_o, 1:wf%n_o, 1:length) => g_ooov
!
   end subroutine point_ooov_abc_cc3
!
!
   module subroutine estimate_mem_integral_setup_abc_cc3(wf, req0, req1)
!!
!!    Estimate memory integrals setup abc
!!    Written by Alexander C. Paul, Dec 2020
!!
!!    Estimate maximum memory needed for cc3 integral setup for abc batching
!!
!!    get_eri_t1_mem returns the memory needed to construct the requested integral
!!    The dimensions sent in specify if an index is batched (1) or of
!!    full dimension (n_o/n_v)
!!    The memory estimate for the first and second pair of indices
!!    is added to the integers req*.
!!
!!    NB: The memory needed to get vvov and vvvo is identical
!!        The memory needed to get oovo and ooov is identical
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(out) :: req0, req1
!
      integer, dimension(2) :: req_vvvo, req_ovov, req_oovo
!
      req_vvvo = wf%eri_t1%get_memory_estimate('vvvo', 1, wf%n_v, 1, wf%n_o)
      req_ovov = wf%eri_t1%get_memory_estimate('ovov', wf%n_o, 1, wf%n_o, 1)
      req_oovo = wf%eri_t1%get_memory_estimate('oovo', wf%n_o, wf%n_o, 1, wf%n_o)
!
      req0 = req_oovo(1)
      req1 = max(req_vvvo(1) + req_vvvo(2), req_ovov(1) + req_ovov(2), req_oovo(2))
!
   end subroutine estimate_mem_integral_setup_abc_cc3
!
!
   module subroutine estimate_mem_c1_integral_setup_abc_cc3(wf, req0, req1, eri_c1)
!!
!!    Estimate memory C1 transformed integrals setup abc
!!    Written by Alexander C. Paul, Dec 2020
!!
!!    Estimate maximum memory needed for cc3 integral setup
!!    for the C1-transformed integrals with abc batching
!!
!!   _c1 get_eri_t1_mem returns the memory needed to construct the requested
!!    c1-transformed integral
!!    The dimensions sent in specify if an index is batched (1) or of
!!    full dimension (n_o/n_v)
!!
!!    6 memory estimates are returned:
!!    1 for each index added to the first 4 req* variables sent in
!!    1 for the first and 1 for the second pair of indices (last 2 req* variables)
!!
!!    NB: The memory requirement is overestimated by the routines.
!!
!
      use eri_adapter_class, only: eri_adapter
!
      implicit none
      class(cc3), intent(in) :: wf
!
      integer, intent(out) :: req0, req1
      type(eri_adapter), intent(in) :: eri_c1
!
      integer, dimension(2) :: req_vvvo, req_oovo
!
      req_vvvo = eri_c1%get_memory_estimate('vvvo', 1, wf%n_v, 1, wf%n_o)
      req_oovo = eri_c1%get_memory_estimate('oovo', wf%n_o, wf%n_o, 1, wf%n_o)
!
      req0 = req_oovo(1)
      req1 = max(req_vvvo(1) + req_vvvo(2), req_oovo(2))
!
   end subroutine estimate_mem_c1_integral_setup_abc_cc3
!
!
   module subroutine construct_W_abc_cc3(wf, a, b, c,            &
                                         t_ijk, u_ijk, t_ijab,   &
                                         g_ljak, g_ljbk, g_ljck, &
                                         g_bdak, g_cdak, g_cdbk, &
                                         g_adbk, g_adck, g_bdck, &
                                         overwrite)
!!
!!    Construct W abc batching
!!    Written by Rolf H. Myhre and Alexander C. Paul July 2019
!!
!!    Contributions to W
!!       W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il(lj|ck))
!!
!!    based on omega_cc3_W_calc written by Rolf H. Myhre
!!
      use reordering, only: add_132_to_123, add_213_to_123, add_312_to_123
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)      :: t_ijk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)      :: u_ijk
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: t_ijab
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)         :: g_ljak
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)         :: g_ljbk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)         :: g_ljck
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: g_bdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: g_cdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: g_cdbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: g_adbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: g_adck
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: g_bdck
!
!     If present and false, t_abc is not overwritten by first dgemm
      logical, optional, intent(in) :: overwrite
!
      real(dp) :: alpha
!
      if (.not. present(overwrite)) then
         alpha = zero
      elseif (overwrite) then
         alpha = zero
      else
         alpha = one
      endif
!
!
!     u_ijk terms
!     -----------
!
!     t^db_ij*(ad|ck)
!     ---------------
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,b),  & ! t_ij_d,b
                  wf%n_o**2,        &
                  g_adck,           & ! g_dk,ac
                  wf%n_v,           &
                  alpha,            &
                  t_ijk,            &
                  wf%n_o**2)
!
!     -t^ab_il*(lj|ck)
!     ----------------
!
      call dgemm('N', 'N',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  t_ijab(:,:,a,b),  & ! t_i_l,ab
                  wf%n_o,           &
                  g_ljck,           & ! g_l_jk,c
                  wf%n_o,           &
                  one,              &
                  t_ijk,            &
                  wf%n_o)
!
!     t^dc_jk*(bd|ai)
!     ---------------
!
      call dgemm('T', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_v,           &
                  one,              &
                  g_bdak,           & ! g_d_i,ba
                  wf%n_v,           &
                  t_ijab(:,:,:,c),  & ! t_jk_d,c
                  wf%n_o**2,        &
                  one,              &
                  t_ijk,            &
                  wf%n_o)
!
!     -t^ca_kl*(li|bj)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljbk,           & ! g_l_ij,b
                  wf%n_o,           &
                  t_ijab(:,:,c,a),  & ! t_k_l,ca
                  wf%n_o,           &
                  one,              &
                  t_ijk,            &
                  wf%n_o**2)
!
!
!     u_ikj terms
!     -----------
!
!     t^dc_ik*(ad|bj)
!     ---------------
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,c),  & ! t_ik_d,c
                  wf%n_o**2,        &
                  g_adbk,           & ! g_dj,ab
                  wf%n_v,           &
                  zero,             &
                  u_ijk,            &
                  wf%n_o**2)
!
!     -t^ac_il*(lk|bj)
!     ----------------
!
      call dgemm('N', 'N',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  t_ijab(:,:,a,c),  & ! t_i_l,ac
                  wf%n_o,           &
                  g_ljbk,           & ! g_l_kj,b
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o)
!
!     t^db_kj*(cd|ai)
!     ---------------
!
      call dgemm('T', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_v,           &
                  one,              &
                  g_cdak,           & ! g_d_i,ca
                  wf%n_v,           &
                  t_ijab(:,:,:,b),  & ! t_kj_d,b
                  wf%n_o**2,        &
                  one,              &
                  u_ijk,            &
                  wf%n_o)
!
!     -t^ba_jl*(li|ck)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljck,           & ! g_l_ik,c
                  wf%n_o,           &
                  t_ijab(:,:,b,a),  & ! t_j_l,ba
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
!
      call add_132_to_123(one, u_ijk, t_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!
!     u_jik terms
!     -----------
!
!     t^da_ji*(bd|ck)
!     ---------------
!
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,a),  & ! t_ji_d,a
                  wf%n_o**2,        &
                  g_bdck,           & ! g_dk,bc
                  wf%n_v,           &
                  zero,             &
                  u_ijk,            &
                  wf%n_o**2)
!
!     -t^cb_kl*(lj|ai)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljak,           & ! g_l_ji,a
                  wf%n_o,           &
                  t_ijab(:,:,c,b),  & ! t_k_l,cb
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
      call add_213_to_123(one, u_ijk, t_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!
!     u_kij terms
!     -----------
!
!     t^da_ki*(cd|bj)
!     ---------------
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,a),  & ! t_ki_d,a
                  wf%n_o**2,        &
                  g_cdbk,           & ! g_dj,cb
                  wf%n_v,           &
                  zero,             &
                  u_ijk,            &
                  wf%n_o**2)
!
!     -t^bc_jl*(lk|ai)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljak,           & ! g_l_ki,a
                  wf%n_o,           &
                  t_ijab(:,:,b,c),  & ! t_j_l,bc
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
      call add_312_to_123(one, u_ijk, t_ijk, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine construct_W_abc_cc3
!
!
   module subroutine divide_by_orbital_differences_abc_cc3(wf, a, b, c, t_ijk, &
                                                           omega, cvs, rm_core)
!!
!!    Divide by orbital energy differences (abc batching)
!!    Written by Alexander C. Paul, July 2019
!!
!!    Divide an array of triples amplitudes (single a,b,c)
!!    by the respective orbital energy differences eps^abc_ijk
!!
!!    omega: If present divide by eps^abc_ijk - omega instead of eps^abc_ijk
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(inout) :: t_ijk
!
      real(dp), optional, intent(in) :: omega
      logical, optional, intent(in)  :: cvs, rm_core
!
      integer  :: i, j, k
      logical  :: cvs_, rm_core_
      real(dp) :: epsilon_abc, epsilon_k, epsilon_kj
!
      cvs_ = .false.
      rm_core_ = .false.
      if (present(cvs)) cvs_ = cvs
      if (present(rm_core)) rm_core_ = rm_core
!
      epsilon_abc =  - wf%orbital_energies(wf%n_o + a) &
                     - wf%orbital_energies(wf%n_o + b) &
                     - wf%orbital_energies(wf%n_o + c)
!
      if (present(omega)) epsilon_abc =  epsilon_abc + omega
!
!$omp parallel do schedule(static) private(i)
      do i = 1, wf%n_o
!
         t_ijk(i,i,i) = zero
!
      enddo
!$omp end parallel do
!
!$omp parallel do schedule(static) private(i,j,k,epsilon_k,epsilon_kj)
      do k = 1, wf%n_o
!
         epsilon_k = epsilon_abc + wf%orbital_energies(k)
!
         do j = 1, wf%n_o
!
            epsilon_kj = epsilon_k + wf%orbital_energies(j)
!
            do i = 1, wf%n_o
!
!              Check for core orbitals (used for excited states):
!              cvs: i,j,k cannot all correspond to valence orbitals
!              rm_core: i,j,k may not contain any core orbital
!
               if (wf%ijk_amplitudes_are_zero(i, j, k, cvs_, rm_core_)) then
                  t_ijk(i,j,k) = zero
                  cycle
               end if
!
               t_ijk(i,j,k) = t_ijk(i,j,k)*one/(epsilon_kj + wf%orbital_energies(i))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine divide_by_orbital_differences_abc_cc3
!
!
   module subroutine outer_product_terms_L3_abc_cc3(wf, a, b, c, L1, L2, L_ijk, &
                                                    F_ov, g_jakb, g_jakc, g_jbkc)
!!
!!    Outer product terms L3 abc-batching
!!    Written by Alexander C. Paul, Okt 2020
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
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: L1
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: L2
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(inout) :: L_ijk
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: F_ov
!
!     (ja|kb) ordered jk,ab
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: g_jakb
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: g_jakc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: g_jbkc
!
      integer :: i, j, k
!
!$omp parallel do private(i,j,k) collapse(2)
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               L_ijk(i,j,k) = L_ijk(i,j,k)          &
                            + L1(i,a)*g_jbkc(j,k)   &
                            + L1(j,b)*g_jakc(i,k)   &
                            + L1(k,c)*g_jakb(i,j)   &
                            + L2(i,j,a,b)*F_ov(k,c) &
                            + L2(i,k,a,c)*F_ov(j,b) &
                            + L2(j,k,b,c)*F_ov(i,a)
!
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine  outer_product_terms_L3_abc_cc3
!
!
end submodule batching_abc
