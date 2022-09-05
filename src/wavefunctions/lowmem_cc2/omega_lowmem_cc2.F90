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
submodule (lowmem_cc2_class) omega_lowmem_cc2
!
!!
!!    Omega submodule
!!
!!    Routines to construct
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
   module subroutine construct_omega_lowmem_cc2(wf, omega)
!!
!!    Construct omega
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Direqts the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(out) :: omega
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct omega lowmem-cc2', pl='normal')
      call timer%turn_on()
!
      call zero_array(omega, wf%n_t1)
!
      call wf%omega_ccs_a1(omega)
!
      call wf%omega_cc2(omega)
!
      call timer%turn_off()
!
   end subroutine construct_omega_lowmem_cc2
!
!
   module subroutine setup_L_Jvo_lowmem_cc2(wf, L_J_vo, point, batch_o)
!!
!!    Setup L_Jvo
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Sets up the vo Cholesky vector and the pointer for a CC2 calculation
!!    The Cholesky vector is given in batches over the occupied index
!!
!!    L_J_vo: final sorted Cholesky vector
!!    point:  pointer to L_J_vo vector
!!
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      type(batching_index), intent(in) :: batch_o
!
      real(dp), dimension(wf%eri_t1%n_J, wf%n_v, batch_o%length), target, intent(out) :: L_J_vo
      real(dp), dimension(:,:,:), pointer, contiguous, intent(out) :: point
!
      call wf%L_t1%get(L_J_vo, wf%n_o+1, wf%n_mo, batch_o%first, batch_o%get_last())
!
      point(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_o%length) => L_J_vo
!
   end subroutine setup_L_Jvo_lowmem_cc2
!
!
   module subroutine setup_L_Jov_lowmem_cc2(wf, L_J_ov_reordered, point, &
                                            L_Jov, batch_o)
!!
!!    Setup L_Jov
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Sets up the ov Cholesky vector and the pointer for a CC2 calculation
!!    The Cholesky vector is given in batches over the occupied index
!!
!!    The Cholesky vector is returned in 132 order:
!!
!!    L_J_ov:          final sorted Cholesky vector
!!    point:           pointer to L_J_ov vector
!!    L_Jov: help array for reordering
!!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      type(batching_index), intent(in) :: batch_o
!
      real(dp), dimension(wf%eri_t1%n_J, wf%n_v, batch_o%length), &
                                          target, intent(out) :: L_J_ov_reordered
      real(dp), dimension(:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(wf%eri_t1%n_J, batch_o%length, wf%n_v), intent(out) :: L_Jov
!
      call wf%L_t1%get(L_Jov, batch_o%first, batch_o%get_last(), &
                                  wf%n_o+1, wf%n_mo)
!
      call sort_123_to_132(L_Jov, L_J_ov_reordered, &
                           wf%eri_t1%n_J, batch_o%length, wf%n_v)
!
      point(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_o%length) => L_J_ov_reordered
!
   end subroutine setup_L_Jov_lowmem_cc2
!
!
   module subroutine setup_L_Joo_lowmem_cc2(wf, L_J_oo_reordered, point, &
                                            L_Joo, batch_o)
!!
!!    Setup L_Joo
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Sets up the oo Cholesky vector and the pointer for a CC2 calculation
!!    The Cholesky vector is given in batches over the occupied index
!!
!!    The Cholesky vector is returned in 132 order:
!!
!!    L_J_oo:          final sorted Cholesky vector
!!    point:           pointer to L_J_oo vector
!!    L_Joo: help array for reordering
!!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      type(batching_index), intent(in) :: batch_o
!
      real(dp), dimension(wf%eri_t1%n_J, wf%n_o, batch_o%length), &
                                          target, intent(out) :: L_J_oo_reordered
      real(dp), dimension(:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(wf%eri_t1%n_J, batch_o%length, wf%n_o), intent(out) :: L_Joo
!
      call wf%L_t1%get(L_Joo, batch_o%first, batch_o%get_last(), 1, wf%n_o)
!
      call sort_123_to_132(L_Joo, L_J_oo_reordered, &
                           wf%eri_t1%n_J, batch_o%length, wf%n_o)
!
      point(1:wf%eri_t1%n_J, 1:wf%n_o, 1:batch_o%length) => L_J_oo_reordered
!
   end subroutine setup_L_Joo_lowmem_cc2
!
!
   module subroutine omega_cc2_lowmem_cc2(wf, omega1)
!!
!!    Omega CC2
!!    Written by Alexander C. Paul, December 2021
!!
!!    Compute the CC2 contributions to omega
!!    without storing the full t2 amplitudes
!!
!!    t^ab_ij = -(eps^ab_ij)^(-1) sum_J L_J_ai L_J_bj
!!    u^ab_ij = 2 t^ab_ij - t^ba_ij
!!
!!    omega^a_i = sum_bj  u^ab_ij F_jb
!!    omega^a_k = sum_bij u^ab_ij L_J_jb L_J_ik
!!    omega^c_i = sum_abj u^ab_ij L_J_jb L_J_ca
!!
      use array_initialization, only: zero_array
      use reordering, only: sort_12_to_21, sort_123_to_132
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
!
      real(dp), dimension(:,:), allocatable :: u_ab
!
      real(dp), dimension(:,:), allocatable :: F_ov_bj
!
      real(dp), dimension(:,:,:), allocatable, target :: L_J_vi
      real(dp), dimension(:,:,:), allocatable, target :: L_J_vj
      real(dp), dimension(:,:,:), contiguous, pointer :: L_J_vi_p => null()
      real(dp), dimension(:,:,:), contiguous, pointer :: L_J_vj_p => null()
!
      real(dp), dimension(:,:,:), allocatable, target :: L_J_ib
      real(dp), dimension(:,:,:), allocatable, target :: L_J_jb
      real(dp), dimension(:,:,:), contiguous, pointer :: L_J_ib_p => null()
      real(dp), dimension(:,:,:), contiguous, pointer :: L_J_jb_p => null()
!
      real(dp), dimension(:,:,:), allocatable, target :: X_J_ai
      real(dp), dimension(:,:,:), allocatable, target :: X_J_aj
      real(dp), dimension(:,:,:), contiguous, pointer :: X_J_ai_p => null()
      real(dp), dimension(:,:,:), contiguous, pointer :: X_J_aj_p => null()
!
      real(dp), dimension(:,:,:), allocatable, target :: L_J_ik
      real(dp), dimension(:,:,:), allocatable, target :: L_J_jk
      real(dp), dimension(:,:,:), contiguous, pointer :: L_J_ik_p => null()
      real(dp), dimension(:,:,:), contiguous, pointer :: L_J_jk_p => null()
!
      real(dp), dimension(:), allocatable  :: array_for_reordering
      real(dp), dimension(:,:), allocatable :: temp_J_v
!
!
      type(direct_stream_file) :: X_J_vo_file
!
      type(batching_index) :: batch_i, batch_j
!
      integer :: i_batch, j_batch
      integer :: i, j, i_rel, j_rel
      integer :: req_0, req_1, req_2, req_1_sort, req_single_batch
!
      X_J_vo_file = direct_stream_file('X_J_vo_cc2', wf%eri_t1%n_J*wf%n_v)
      call X_J_vo_file%open_
!
      call mem%alloc(F_ov_bj, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_bj, wf%n_o, wf%n_v)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
!
      req_1_sort = max(wf%eri_t1%n_J, wf%n_v) * max(wf%n_v, wf%n_o)
!
      req_0 = wf%n_v**2 + wf%eri_t1%n_J*wf%n_v
      req_1 = 3*wf%eri_t1%n_J*wf%n_v + wf%eri_t1%n_J*wf%n_o
      req_2 = 0
!
      req_single_batch = req_0 + (3*wf%eri_t1%n_J*wf%n_v + wf%eri_t1%n_J*wf%n_o + req_1_sort)*wf%n_o
!
      call mem%batch_setup(batch_i, batch_j, req_0, req_1 + req_1_sort, req_1, req_2, &
                           'omega_lowmem_cc2', req_single_batch=req_single_batch)
!
      call mem%alloc(u_ab, wf%n_v, wf%n_v)
      call mem%alloc(temp_J_v, wf%eri_t1%n_J, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(L_J_vi, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
         call mem%alloc(L_J_ib, wf%eri_t1%n_J, wf%n_v, wf%n_o)
         call mem%alloc(L_J_ik, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
         call mem%alloc(X_J_ai, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
         call mem%alloc(array_for_reordering, req_1_sort*wf%n_o)
!
      else ! batching
!
         call mem%alloc(L_J_vi, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
         call mem%alloc(L_J_vj, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
         call mem%alloc(L_J_ib, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
         call mem%alloc(L_J_jb, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
         call mem%alloc(L_J_ik, wf%eri_t1%n_J, wf%n_o, batch_i%max_length)
         call mem%alloc(L_J_jk, wf%eri_t1%n_J, wf%n_o, batch_i%max_length)
!
         call mem%alloc(X_J_ai, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
         call mem%alloc(X_J_aj, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
         call mem%alloc(array_for_reordering, req_1_sort*batch_i%max_length)
!
      endif
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_L_Jvo(L_J_vi, L_J_vi_p, batch_i)
         call wf%setup_L_Jov(L_J_ib, L_J_ib_p, array_for_reordering, batch_i)
         call wf%setup_L_Joo(L_J_ik, L_J_ik_p, array_for_reordering, batch_i)
!
         call zero_array(X_J_ai, wf%eri_t1%n_J*wf%n_v*batch_i%length)
         X_J_ai_p(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_i%length) => X_J_ai
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_L_Jvo(L_J_vj, L_J_vj_p, batch_j)
               call wf%setup_L_Jov(L_J_jb, L_J_jb_p, array_for_reordering, batch_j)
               call wf%setup_L_Joo(L_J_jk, L_J_jk_p, array_for_reordering, batch_j)
!
               call X_J_vo_file%read_range(X_J_aj, batch_j)
               X_J_aj_p(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_j%length) => X_J_aj
!
            else
!
               L_J_vj_p(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_j%length) => L_J_vi
               L_J_jb_p(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_j%length) => L_J_ib
               L_J_jk_p(1:wf%eri_t1%n_J, 1:wf%n_o, 1:batch_j%length) => L_J_ik
!
               X_J_aj_p(1:wf%eri_t1%n_J, 1:wf%n_v, 1:batch_j%length) => X_J_ai
!
            end if
!
            do i = batch_i%first, batch_i%get_last()
!
               i_rel = i - batch_i%first + 1
!
               do j = batch_j%first, min(batch_j%get_last(), i)
!
                  j_rel = j - batch_j%first + 1
!
                  call wf%construct_contravariant_t2_single_ij(i, j, array_for_reordering, &
                                                               u_ab, &
                                                               L_J_vi_p(:,:,i_rel), &
                                                               L_J_vj_p(:,:,j_rel))
!
                  call wf%omega_cc2_fock(i, j, u_ab, omega1, F_ov_bj)
!
                  call wf%omega_cc2_v2o2J(i, j, u_ab, temp_J_v, &
                                          omega1,              &
                                          X_J_ai_p(:,:,i_rel), &
                                          X_J_aj_p(:,:,j_rel), &
                                          L_J_ib_p(:,:,i_rel), &
                                          L_J_jb_p(:,:,j_rel), &
                                          L_J_ik_p(:,:,i_rel), &
                                          L_J_jk_p(:,:,j_rel))
!
               end do ! j
            end do ! i
               if (j_batch .ne. i_batch) call X_J_vo_file%write_range(X_J_aj, batch_j)
         end do ! j_batch
         call X_J_vo_file%write_range(X_J_ai, batch_i)
      end do ! i_batch
!
      call mem%dealloc(u_ab, wf%n_v, wf%n_v)
      call mem%dealloc(temp_J_v, wf%eri_t1%n_J, wf%n_v)
      call mem%dealloc(F_ov_bj, wf%n_v, wf%n_o)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(L_J_vi, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
         call mem%dealloc(L_J_ib, wf%eri_t1%n_J, wf%n_v, wf%n_o)
         call mem%dealloc(L_J_ik, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
         call mem%dealloc(X_J_ai, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
         call mem%dealloc(array_for_reordering, req_1_sort*wf%n_o)
!
      else ! batching
!
         call mem%dealloc(L_J_vi, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
         call mem%dealloc(L_J_vj, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(L_J_ib, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
         call mem%dealloc(L_J_jb, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(L_J_ik, wf%eri_t1%n_J, wf%n_o, batch_i%max_length)
         call mem%dealloc(L_J_jk, wf%eri_t1%n_J, wf%n_o, batch_i%max_length)
!
         call mem%dealloc(X_J_ai, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
         call mem%dealloc(X_J_aj, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(array_for_reordering, req_1_sort*batch_i%max_length)
!
      endif
!
      call mem%batch_finalize()
!
      call wf%omega_cc2_Jv2o(X_J_vo_file, omega1)
!
      call X_J_vo_file%close_('delete')
!
   end subroutine omega_cc2_lowmem_cc2
!
!
   module subroutine construct_contravariant_t2_single_ij_lowmem_cc2(wf, i, j, v_ab, u_ab, &
                                                                     L_J_vi, L_J_vj)
!!
!!    Construct contravariant t2 for single i,j
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    t^ab_ij = -(eps^ab_ij)^(-1) sum_J L_J_ai L_J_bj
!!    u^ab_ij = 2 t^ab_ij - t^ba_ij
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      integer, intent(in) :: i, j
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: v_ab
      real(dp), dimension(wf%n_v, wf%n_v), intent(out)   :: u_ab
!
      real(dp), dimension(wf%eri_t1%n_J, wf%n_v), intent(in) :: L_J_vi, L_J_vj
!
      call dgemm('T', 'N',       &
                  wf%n_v,        &
                  wf%n_v,        &
                  wf%eri_t1%n_J, &
                  one,           &
                  L_J_vi,        &
                  wf%eri_t1%n_J, &
                  L_J_vj,        &
                  wf%eri_t1%n_J, &
                  zero,          &
                  v_ab,          &
                  wf%n_v)
!
      call wf%make_contravariant_doubles_single_ij(i, j, v_ab, u_ab)
!
   end subroutine construct_contravariant_t2_single_ij_lowmem_cc2
!
!
   module subroutine make_contravariant_doubles_single_ij_lowmem_cc2(wf, i, j, u_ab, v_ab)
!!
!!    Make contravariant doubles single i,j
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Divide an array of lowmem_cc2 amplitudes (single i,j)
!!    by the respective orbital energy differences eps^ab_ij
!!    And construct the contravariant linear combination
!!       u^ab_ij = 2 t^ab_ij - t^ba_ij
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      integer, intent(in) :: i, j
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)  :: u_ab
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: v_ab
!
      integer a, b
!
      real(dp) :: epsilon_ij, epsilon_b
!
      epsilon_ij = wf%orbital_energies(i) + wf%orbital_energies(j)
!
!$omp parallel do schedule(static) private(a,b,epsilon_b)
      do b = 1, wf%n_v
!
         epsilon_b = epsilon_ij - wf%orbital_energies(wf%n_o + b)
!
         do a = 1, wf%n_v
!
            v_ab(a,b) = (two*u_ab(a,b) - u_ab(b,a))&
                        /(epsilon_b - wf%orbital_energies(wf%n_o + a))
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine make_contravariant_doubles_single_ij_lowmem_cc2
!
!
   module subroutine omega_cc2_fock_lowmem_cc2(wf, i, j, t_ab, omega1, F_ov)
!!
!!    Omega CC2 Fock
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Compute:
!!       omega^a_i = sum_bj  u^ab_ij F_jb
!!    for a given i and j
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      integer, intent(in) :: i, j
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: t_ab
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: F_ov
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
!
      call dgemv('N',           &
                  wf%n_v,       &
                  wf%n_v,       &
                  one,          &
                  t_ab,         &
                  wf%n_v,       &
                  F_ov(:,j), 1, &
                  one,          &
                  omega1(:,i), 1)
!
      if (i .ne. j) then
!
         call dgemv('T',           &
                     wf%n_v,       &
                     wf%n_v,       &
                     one,          &
                     t_ab,         &
                     wf%n_v,       &
                     F_ov(:,i), 1, &
                     one,          &
                     omega1(:,j), 1)
!
      end if
!
   end subroutine omega_cc2_fock_lowmem_cc2
!
!
   module subroutine omega_cc2_v2o2J_lowmem_cc2(wf, i, j, t_ab, temp_J_v, omega1, &
                                                X_J_ai, X_J_aj,           &
                                                L_Jib, L_Jjb, L_Jik, L_Jjk)
!!
!!    Omega CC2 v2o2J
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Compute:
!!       X_J_ai = sum_bj u^ab_ij L_J_jb
!!       omega^a_k = sum_Ji X_J_ai L_J_ik
!!    for a given i and j
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      integer, intent(in) :: i, j
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: t_ab
      real(dp), dimension(wf%eri_t1%n_J, wf%n_v), intent(out) :: temp_J_v
!
      real(dp), dimension(wf%n_v, wf%n_o),     intent(inout) :: omega1
      real(dp), dimension(wf%eri_t1%n_J, wf%n_v), intent(inout) :: X_J_ai, X_J_aj
!
      real(dp), dimension(wf%eri_t1%n_J, wf%n_v), intent(in) :: L_Jib, L_Jjb
      real(dp), dimension(wf%eri_t1%n_J, wf%n_o), intent(in) :: L_Jik, L_Jjk
!
      call dgemm('N', 'T',       &
                  wf%eri_t1%n_J, &
                  wf%n_v,        &
                  wf%n_v,        &
                  one,           &
                  L_Jjb,         & ! L_J_b,j
                  wf%eri_t1%n_J, &
                  t_ab,          & ! t_a_b,ij
                  wf%n_v,        &
                  zero,          &
                  temp_J_v,      & ! X_J_a,i
                  wf%eri_t1%n_J)
!
      call daxpy(wf%eri_t1%n_J*wf%n_v, one, temp_J_v, 1, X_J_ai, 1)
!
      call dgemm('T', 'N',       &
                  wf%n_v,        &
                  wf%n_o,        &
                  wf%eri_t1%n_J, &
                  -one,          &
                  temp_J_v,      & ! X_J_a,i
                  wf%eri_t1%n_J, &
                  L_Jik,         & ! L_J_k,i
                  wf%eri_t1%n_J, &
                  one,           &
                  omega1,        & ! omega_a_k
                  wf%n_v)
!
      if (i .ne. j) then
!
         call dgemm('N', 'N',       &
                     wf%eri_t1%n_J, &
                     wf%n_v,        &
                     wf%n_v,        &
                     one,           &
                     L_Jib,         & ! L_J_b,i
                     wf%eri_t1%n_J, &
                     t_ab,          & ! t_b_a,ij
                     wf%n_v,        &
                     zero,          &
                     temp_J_v,      & ! X_J_a,j
                     wf%eri_t1%n_J)
!
         call daxpy(wf%eri_t1%n_J*wf%n_v, one, temp_J_v, 1, X_J_aj, 1)
!
         call dgemm('T', 'N',       &
                     wf%n_v,        &
                     wf%n_o,        &
                     wf%eri_t1%n_J, &
                     -one,          &
                     temp_J_v,      & ! X_J_a,j
                     wf%eri_t1%n_J, &
                     L_Jjk,         & ! L_J_k,j
                     wf%eri_t1%n_J, &
                     one,           &
                     omega1,        & ! omega_a_k
                     wf%n_v)
!
      end if
!
   end subroutine omega_cc2_v2o2J_lowmem_cc2
!
!
   module subroutine omega_cc2_Jv2o_lowmem_cc2(wf, X_J_vo_file, omega1)
!!
!!    Omega cc2 Jv2o
!!    Written by Alexander C. Paul, Dec 2021
!!
!!    Compute:
!!       omega^c_i = sum_abj X_Jai L_J_ca
!!
!!    where X_J_ai is stored on file (s. omega_cc2_lowmem_cc2)
!!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      type(direct_stream_file), intent(inout) :: X_J_vo_file
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ca, L_J_vv, X_J_vo
      real(dp), dimension(:,:), allocatable :: temp
!
      type(batching_index) :: batch_c, batch_i
!
      integer :: c_batch, i_batch, req_0, req_c, req_i, req_2
      integer :: c, i, c_rel, i_rel
!
      batch_i = batching_index(wf%n_o)
      batch_c = batching_index(wf%n_v)
!
      req_0 = 0
      req_i = wf%eri_t1%n_J*wf%n_v
      req_c = 2*wf%eri_t1%n_J*wf%n_v
      req_2 = 1
!
      call mem%batch_setup(batch_i, batch_c, req_0, req_i, req_c, req_2, 'omega_cc2_Jv2o')
!
      call mem%alloc(X_J_vo, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
         call mem%alloc(L_J_ca, wf%eri_t1%n_J, batch_c%length, wf%n_v)
         call mem%alloc(L_J_vv, wf%eri_t1%n_J, wf%n_v, batch_c%length)
!
         call wf%L_t1%get(L_J_ca, &
                          wf%n_o+batch_c%first, &
                          wf%n_o+batch_c%get_last(), &
                          wf%n_o+1, wf%n_mo)
!
         call sort_123_to_132(L_J_ca, L_J_vv, wf%eri_t1%n_J, batch_c%length, wf%n_v)
         call mem%dealloc(L_J_ca, wf%eri_t1%n_J, wf%n_v, batch_c%length)
!
         do i_batch = 1, batch_i%num_batches
!
            call batch_i%determine_limits(i_batch)
!
            call mem%alloc(temp, batch_c%length, batch_i%length)
!
!           only need to read once if there is only a single batch in i
            if (.not. (batch_i%num_batches == 1 .and. c_batch > 1)) then
               call X_J_vo_file%read_range(X_J_vo, batch_i)
            end if
!
            call dgemm('T', 'N',           &
                        batch_c%length,    &
                        batch_i%length,    &
                        wf%eri_t1%n_J*wf%n_v, &
                        one,               &
                        L_J_vv,            & ! L_Ja_c
                        wf%eri_t1%n_J*wf%n_v, &
                        X_J_vo,            & ! X_Ja_i
                        wf%eri_t1%n_J*wf%n_v, &
                        zero,              &
                        temp,              & ! omega_c_i
                        batch_c%length)
!
            do i = batch_i%first, batch_i%get_last()
               do c = batch_c%first, batch_c%get_last()
!
                  c_rel = c - batch_c%first + 1
                  i_rel = i - batch_i%first + 1
!
                  omega1(c, i) = omega1(c, i) + temp(c_rel, i_rel)
!
               end do
            end do
!
            call mem%dealloc(temp, batch_c%length, batch_i%length)
!
         end do ! i_batch
!
         call mem%dealloc(L_J_vv, wf%eri_t1%n_J, wf%n_v, batch_c%length)
!
      end do ! c_batch
!
      call mem%dealloc(X_J_vo, wf%eri_t1%n_J, wf%n_v, batch_i%max_length)
!
      call mem%batch_finalize
!
   end subroutine omega_cc2_Jv2o_lowmem_cc2
!
!
end submodule omega_lowmem_cc2
