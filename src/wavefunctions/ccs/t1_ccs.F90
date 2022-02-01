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
submodule (ccs_class) t1_ccs
!
!!
!!    t1 submodule
!!
!!    Submodule containing routines that can be used t1-transform arrays.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine t1_transform_ccs(wf, Z_pq)
!!
!!    T1 transform
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Assumes that Z is in the MO basis and performs the T1 transformation,
!!
!!       Z_pq <- sum_rs X_ps Z_sr Y_qr,    i.e.    Z <- X Z Y^T
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    Here, t1 is a full MO matrix whose only non-zero block is the vir-occ
!!    part, where it is equal to t_i^a.
!!
      use array_utilities, only: zero_array, sandwich
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
      real(dp), dimension(:,:), allocatable :: X, Y
!
      integer :: p, i, a
!
!     Construct the X and Y arrays
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      call zero_array(X, (wf%n_mo)**2)
      call zero_array(Y, (wf%n_mo)**2)
!
!$omp parallel do private(p)
      do p = 1, wf%n_mo
!
         X(p, p) = one
         Y(p, p) = one
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            X(wf%n_o + a, i) = -wf%t1(a, i)
            Y(i, wf%n_o + a) = wf%t1(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
!     Construct Z_pq = sum_sr X_ps Z_sr Y_qr
!
      call sandwich(Z_pq, X, Y, wf%n_mo, .false.)
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
!
   end subroutine t1_transform_ccs
!
!
   module subroutine add_t1_terms_ccs(wf, Z_pq)
!!
!!    Add t1 terms
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
!!    Here, Z is assumed to be the density matrix with no T1 contributions on input
!!    - the so-called T1-transformed density matrix.
!!    The routine adds the missing T1 contributions to Z:
!!
!!       Z_pq <- sum_rs X_sp Z_sr Y_rq,    i.e.    Z <- X^T Z Y
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    t1 is a full MO matrix whose only non-zero block is the vir-occ part,
!!    where it is equal to t_i^a.
!!
!!    Based on t1_transform_ccs by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      use array_utilities, only: zero_array, sandwich
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
      real(dp), dimension(:,:), allocatable :: X, Y
!
      integer :: p, i, a
!
!     Construct the X and Y arrays
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      call zero_array(X, (wf%n_mo)**2)
      call zero_array(Y, (wf%n_mo)**2)
!
!$omp parallel do private(p)
      do p = 1, wf%n_mo
!
         X(p, p) = one
         Y(p, p) = one
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            X(wf%n_o + a, i) = -wf%t1(a, i)
            Y(i, wf%n_o + a) = wf%t1(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
!     Construct Z_pq = sum_sr X_sp Z_sr Y_rq
!
      call sandwich(Z_pq, X, Y, wf%n_mo, .true.)
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
!
   end subroutine add_t1_terms_ccs
!
!
   module subroutine add_t1_terms_and_transform_ccs(wf, Z_pq, Z_out)
!!
!!    Add t1 terms and transform
!!    Written by Tor S. Haugland, Nov 2019 (as do_visualization)
!!
!!    Here Z, on input, is assumed to be the density matrix with no T1 contributions
!!    - the so-called T1-transformed density matrix.
!!    The routine adds the missing T1 contributions to Z and transforms it
!!    with the AO coefficients to obtain a density as needed by the visualization tool.
!!
!!       Z_alpha,beta = (sum_pq  Z_pq C_alpha,p C_beta,q)
!!
!!    Renamed and moved here, by Alexander C. Paul, May 2020
!!
      use array_utilities, only: symmetric_sandwich_right_transposition
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)  :: Z_pq
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: Z_out
!
      real(dp), dimension(:,:), allocatable :: Z_mo
!
      call mem%alloc(Z_mo, wf%n_mo, wf%n_mo)
!
      call dcopy(wf%n_mo**2, Z_pq, 1, Z_mo, 1)
!
      call wf%add_t1_terms(Z_mo)
!
!     D_alpha,beta = sum_pq  D_pq C_alpha,p C_beta,q
!
      call symmetric_sandwich_right_transposition(Z_out, Z_mo,             &
                                                  wf%orbital_coefficients, &
                                                  wf%ao%n, wf%n_mo)
!
      call mem%dealloc(Z_mo, wf%n_mo, wf%n_mo)
!
   end subroutine add_t1_terms_and_transform_ccs
!
!
   module subroutine construct_t1_cholesky_ccs(wf, t1, L_mo, L_t1)
!!
!!    Construct T1 Cholesky
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Wrapper for the constructors of the different Cholesky T1 blocks
!!
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) ::t1
      class(abstract_eri_cholesky), intent(inout) :: L_mo
      class(abstract_eri_cholesky), intent(inout) :: L_t1
!
      type(timings) :: timer
      timer = timings('T1 transform cholesky', pl='m')
      call timer%turn_on()
!
!     occupied-occupied block
!
      call wf%construct_cholesky_t1_oo(t1, L_mo, L_t1)
!
!     virtual-occupied block
!
      call wf%construct_cholesky_t1_vo(t1, L_mo, L_t1)
!
!     virtual-virtual block
!
      call wf%construct_cholesky_t1_vv(t1, L_mo, L_t1)
!
      call wf%L_t1%notify_observers()
!
      call timer%turn_off()
!
   end subroutine construct_t1_cholesky_ccs
!
!
   module subroutine construct_cholesky_t1_oo_ccs(wf, t1, L_mo, L_t1)
!!
!!    Construct Cholesky T1 oo
!!    Written by Rolf H. Myhre, Sarai D. Folkestad and Eirik F. Kjønstad 2020-2021
!!
!!    Computes
!!
!!       L_J_ij_T1 = L_J_ij_MO + sum_a L_J_ia_MO*t_aj,
!!
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
      class(abstract_eri_cholesky), intent(inout) :: L_mo
      class(abstract_eri_cholesky), intent(inout) :: L_t1
!
      real(dp), dimension(:), allocatable :: L_J_oo
      real(dp), dimension(:), allocatable :: L_J_ox
!
      type(batching_index) :: batch_o
!
      integer :: o_batch
!
      batch_o = batching_index(wf%n_o)
      call mem%batch_setup(batch_o, 0, L_mo%n_J*(wf%n_v + wf%n_o), &
                           tag='Cholesky t1 oo')
!
      call mem%alloc(L_J_oo, L_mo%n_J*batch_o%max_length*wf%n_o)
      call mem%alloc(L_J_ox, L_mo%n_J*batch_o%max_length*wf%n_v)
!
      do o_batch = 1,batch_o%num_batches
!
         call batch_o%determine_limits(o_batch)
!
!        L_J_ij_t1 = L_J_ij_mo
         call L_mo%get(L_J_oo, batch_o%first, batch_o%get_last(), 1, wf%n_o)
!
!        L_J_ij_t1 += sum_b L_J_ib_mo t_bj
         call L_mo%get(L_J_ox, batch_o%first, batch_o%get_last(), wf%n_o + 1, wf%n_mo)
!
         call dgemm('N', 'N',                                  &
                    L_mo%n_J*batch_o%length, wf%n_o, wf%n_v,   &
                    one,                                       &
                    L_J_ox, L_mo%n_J*batch_o%length,           &
                    t1, wf%n_v,                                &
                    one,                                       &
                    L_J_oo, L_mo%n_J*batch_o%length)
!
         call L_t1%set(L_J_oo, batch_o%first, batch_o%get_last(), 1, wf%n_o)
!
      enddo
!
      call mem%dealloc(L_J_oo, L_mo%n_J*batch_o%max_length*wf%n_o)
      call mem%dealloc(L_J_ox, L_mo%n_J*batch_o%max_length*wf%n_v)
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_t1_oo_ccs
!
!
   module subroutine construct_cholesky_t1_vo_ccs(wf, t1, L_mo, L_t1)
!!
!!    Construct Cholesky T1 vo
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Computes
!!
!!       L_ai_J_T1 = L_J_ai_MO + sum_b L_J_ab_MO*t_bi
!!                             - sum_j t_aj*(L_J_ji_MO + sum_b L_J_jb_MO*t_bi)
!!
!
      use reordering, only: add_132_to_123
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
      class(abstract_eri_cholesky), intent(inout) :: L_mo
      class(abstract_eri_cholesky), intent(inout) :: L_t1
!
      real(dp), dimension(:), allocatable, target :: L_J_ij
      real(dp), dimension(:), allocatable, target :: L_J_ji
      real(dp), dimension(:), allocatable, target :: L_J_ai
      real(dp), dimension(:), allocatable, target :: L_J_ab
!
      real(dp), dimension(:,:,:), pointer :: L_J_ij_p
      real(dp), dimension(:,:,:), pointer :: L_J_ji_p
      real(dp), dimension(:,:,:), pointer :: L_J_ia_p
      real(dp), dimension(:,:,:), pointer :: L_J_ai_p
      real(dp), dimension(:,:,:), pointer :: L_J_ab_p
!
      type(batching_index) :: batch_v, batch_o
!
      integer :: v_batch, o_batch, req
!
      batch_o = batching_index(wf%n_o)
      batch_v = batching_index(wf%n_v)
      req = L_mo%n_J*(wf%n_o + max(wf%n_v,wf%n_o))
      call mem%batch_setup(batch_o, batch_v, 0, req, req, 0, tag='Cholesky t1 vo')
!
      call mem%alloc(L_J_ij, L_mo%n_J*batch_o%max_length*max(wf%n_v, wf%n_o))
      call mem%alloc(L_J_ji, L_mo%n_J*batch_o%max_length*wf%n_o)
!
      call mem%alloc(L_J_ai, L_mo%n_J*batch_v%max_length*wf%n_o)
      call mem%alloc(L_J_ab, L_mo%n_J*batch_v%max_length*max(wf%n_v, wf%n_o))
!
      do v_batch = 1, batch_v%num_batches
!
         call batch_v%determine_limits(v_batch)
!
         L_J_ai_p(1:L_mo%n_J, 1:batch_v%length, 1:wf%n_o) => L_J_ai
         L_J_ab_p(1:L_mo%n_J, 1:batch_v%length, 1:wf%n_v) => L_J_ab
!
!        L_J_ai_t1 = L_J_ai_mo
         call L_mo%get(L_J_ai_p, wf%n_o+batch_v%first, wf%n_o+batch_v%get_last(), 1, wf%n_o)
!
!        L_J_ai_t1 += sum_b L_J_ab_mo t_bi
         call L_mo%get(L_J_ab_p, wf%n_o+batch_v%first, wf%n_o+batch_v%get_last(), wf%n_o+1, wf%n_mo)
!
         call dgemm('N', 'N',                                     &
                    L_mo%n_J*batch_v%length, wf%n_o, wf%n_v,      &
                    one,                                          &
                    L_J_ab_p, L_mo%n_J*batch_v%length,            &
                    t1, wf%n_v,                                   &
                    one,                                          &
                    L_J_ai_p, L_mo%n_J*batch_v%length)
!
         do o_batch = 1, batch_o%num_batches
!
            call batch_o%determine_limits(o_batch)
!
            L_J_ji_p(1:L_mo%n_J, 1:batch_o%length, 1:wf%n_o) => L_J_ji
            L_J_ia_p(1:L_mo%n_J, 1:batch_o%length, 1:wf%n_v) => L_J_ij
!
!           L'_J_ji = sum_b L_J_jb_mo t_bi
            call L_mo%get(L_J_ia_p, batch_o%first, batch_o%get_last(), wf%n_o+1, wf%n_mo)
!
            call dgemm('N', 'N',                                  &
                       L_mo%n_J*batch_o%length, wf%n_o, wf%n_v,   &
                       one,                                       &
                       L_J_ia_p, L_mo%n_J*batch_o%length,         &
                       t1, wf%n_v,                                &
                       zero,                                      &
                       L_J_ji_p, L_mo%n_J*batch_o%length)
!
!           L'_J_ji += L_J_ji_mo
            L_J_ij_p(1:L_mo%n_J, 1:wf%n_o, 1:batch_o%length) => L_J_ij
            call L_mo%get(L_J_ij_p, 1, wf%n_o, batch_o%first, batch_o%get_last())
!
            call add_132_to_123(one, L_J_ji_p, L_J_ij_p, L_mo%n_J, wf%n_o, batch_o%length)
!
!           L_J_ai_t1 -= sum_j t_aj L'_J_ji
            L_J_ia_p(1:L_mo%n_J, 1:wf%n_o, 1:batch_v%length) => L_J_ab
!
            call dgemm('N', 'T',                                        &
                       L_mo%n_J*wf%n_o, batch_v%length, batch_o%length, &
                       -one,                                            &
                       L_J_ij_p, L_mo%n_J*wf%n_o,                       &
                       t1(batch_v%first, batch_o%first), wf%n_v,        &
                       zero,                                            &
                       L_J_ia_p, L_mo%n_J*wf%n_o)
!
            call add_132_to_123(one, L_J_ia_p, L_J_ai_p, L_mo%n_J, batch_v%length, wf%n_o)
!
         enddo
!
         call L_t1%set(L_J_ai_p, wf%n_o+batch_v%first, wf%n_o+batch_v%get_last(), 1, wf%n_o)
!
      enddo
!
      call mem%dealloc(L_J_ab, L_mo%n_J*batch_v%max_length*max(wf%n_v, wf%n_o))
      call mem%dealloc(L_J_ai, L_mo%n_J*batch_v%max_length*wf%n_o)
!
      call mem%dealloc(L_J_ji, L_mo%n_J*batch_o%max_length*wf%n_o)
      call mem%dealloc(L_J_ij, L_mo%n_J*batch_o%max_length*max(wf%n_v, wf%n_o))
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_t1_vo_ccs
!
!
   module subroutine construct_cholesky_t1_vv_ccs(wf, t1, L_mo, L_t1)
!!
!!    Construct Cholesky T1 vv
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Computes
!!
!!       L_J_ab_T1 = L_J_ab_MO - sum_j t_aj*L_J_jb_MO
!!
!
      use reordering, only: add_132_to_123
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
      class(abstract_eri_cholesky), intent(inout) :: L_mo
      class(abstract_eri_cholesky), intent(inout) :: L_t1
!
      real(dp), dimension(:,:,:), allocatable :: L_J_vv
      real(dp), dimension(:,:,:), allocatable, target :: L_J_xv
!
      type(batching_index) :: batch_v
!
      integer :: v_batch
!
      batch_v = batching_index(wf%n_v)
      call mem%batch_setup(batch_v, 0, 2*L_mo%n_J*max(wf%n_v, wf%n_o), tag='Cholesky t1 vv')
!
      call mem%alloc(L_J_xv, L_mo%n_J, batch_v%max_length, max(wf%n_v, wf%n_o))
      call mem%alloc(L_J_vv, L_mo%n_J, batch_v%max_length, max(wf%n_v, wf%n_o))
!
      do v_batch = 1, batch_v%num_batches
!
         call batch_v%determine_limits(v_batch)
!
!        L_J_ab_t1 -= sum_j t_aj L_J_jb_mo
         call L_mo%get(L_J_xv, wf%n_o + batch_v%first, wf%n_o + batch_v%get_last(), 1, wf%n_o)
!
         call dgemm('N', 'T',                                  &
                    L_mo%n_J*batch_v%length, wf%n_v, wf%n_o,   &
                    -one,                                      &
                    L_J_xv, L_mo%n_J*batch_v%length,           &
                    t1, wf%n_v,                                &
                    zero,                                      &
                    L_J_vv, L_mo%n_J*batch_v%length)
!
!        L_J_ab_t1 = L_J_ab_mo
         call L_mo%get(L_J_xv, wf%n_o + 1, wf%n_mo, wf%n_o + batch_v%first, wf%n_o + batch_v%get_last())
!
         call add_132_to_123(one, L_J_vv, L_J_xv, L_mo%n_J, wf%n_v, batch_v%length)
!
         call L_t1%set(L_J_xv, wf%n_o+1, wf%n_mo, wf%n_o + batch_v%first, wf%n_o + batch_v%get_last())
!
      enddo
!
      call mem%dealloc(L_J_xv, L_mo%n_J, batch_v%max_length, max(wf%n_v, wf%n_o))
      call mem%dealloc(L_J_vv, L_mo%n_J, batch_v%max_length, max(wf%n_v, wf%n_o))
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_t1_vv_ccs
!
!
end submodule t1_ccs
