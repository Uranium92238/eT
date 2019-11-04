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
submodule (ccs_class) t1_ccs
!
!!
!!    t1 submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Oct 2019
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
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
      real(dp), dimension(:,:), allocatable :: X, Y
!
      real(dp), dimension(:,:), allocatable :: W ! W_sq = sum_r Z_sr Y_rq^T, intermediate
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
!     Construct intermediate W = Z Y^T and then use it to do transformation
!
      call mem%alloc(W, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'T', &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  Z_pq,    & ! Z_s_r
                  wf%n_mo, &
                  Y,       & ! Y_q_r
                  wf%n_mo, &
                  zero,    &
                  W,       & ! W_sq = sum_r Z_sr Y_rq
                  wf%n_mo)
!
      call dgemm('N', 'N', &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  X,       &
                  wf%n_mo, &
                  W,       &
                  wf%n_mo, &
                  zero,    &
                  Z_pq,    & ! Z_pq = (X W)_pq = sum_s X_ps W_sq = sum_sr X_ps Z_sr Y_rq
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
      call mem%dealloc(W, wf%n_mo, wf%n_mo)
!
   end subroutine t1_transform_ccs
!
!
   module subroutine t1_transform_4_ccs(wf, Z_tuvw, Z_pqrs, t1)
!!
!!    T1 transform 4 index arrays
!!    Written by Andreas Skeidsvoll, Apr 2019
!!
!!    Assumes that Z is in the MO basis and performs the T1 transformation,
!!
!!       Z_pqrs = sum_tuvw X_pt Y_qu X_rm Y_sn Z_tuvw,
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    Here, t1 is a full MO matrix whose only non-zero block is the vir-occ
!!    part, where it is equal to t_i^a.
!!    NB: needs place for an additional 2*wf%n_mo**4 + wf%n_t1 in memory.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo), intent(in) :: Z_tuvw
      real(dp), dimension(wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo), intent(out) :: Z_pqrs
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
!
      real(dp), dimension(:,:), allocatable :: X, Y
      real(dp), dimension(:,:,:,:), allocatable :: Z_tuvs, Z_puvs, Z_vspu, Z_vspq
!
      integer :: p, a, i
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      X = zero
      Y = zero
!
      do p = 1, wf%n_mo
!
         X(p, p) = one
         Y(p, p) = one
!
      enddo
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            X(wf%n_o + a, i) = -t1(a, i)
            Y(i, wf%n_o + a) = t1(a, i)
!
         enddo
      enddo
!
      call mem%alloc(Z_tuvs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'T',   &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 wf%n_mo,    &
                 one,        &
                 Z_tuvw,     &
                 wf%n_mo**3, &
                 Y,          &
                 wf%n_mo,    &
                 zero,       &
                 Z_tuvs,     &
                 wf%n_mo**3)
!
      call mem%alloc(Z_puvs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'N',   &
                 wf%n_mo,    &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 one,        &
                 X,          &
                 wf%n_mo,    &
                 Z_tuvs,     &
                 wf%n_mo,    &
                 zero,       &
                 Z_puvs,     &
                 wf%n_mo)
!
      call mem%dealloc(Z_tuvs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call mem%alloc(Z_vspu, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call sort_1234_to_3412(Z_puvs, Z_vspu, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call mem%dealloc(Z_puvs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call mem%alloc(Z_vspq, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'T',   &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 wf%n_mo,    &
                 one,        &
                 Z_vspu,     &
                 wf%n_mo**3, &
                 Y,          &
                 wf%n_mo,    &
                 zero,       &
                 Z_vspq,     &
                 wf%n_mo**3)
!
      call mem%dealloc(Z_vspu, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
!     Z_rspq = Z_pqrs because of particle permutation symmetry
!
      call dgemm('N', 'N',   &
                 wf%n_mo,    &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 one,        &
                 X,          &
                 wf%n_mo,    &
                 Z_vspq,     &
                 wf%n_mo,    &
                 zero,       &
                 Z_pqrs,     &
                 wf%n_mo)
!
      call mem%dealloc(Z_vspq, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
!
   end subroutine t1_transform_4_ccs
!
!
end submodule t1_ccs