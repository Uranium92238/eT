!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
submodule (ccs_class) t1_ccs_complex
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
   module subroutine t1_transform_ccs_complex(wf, Z_pq)
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
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
      complex(dp), dimension(:,:), allocatable :: X, Y
!
      integer :: p, i, a
!
!     Construct the X and Y arrays
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      call zero_array_complex(X, (wf%n_mo)**2)
      call zero_array_complex(Y, (wf%n_mo)**2)
!
!$omp parallel do private(p)
      do p = 1, wf%n_mo
!
         X(p, p) = one_complex
         Y(p, p) = one_complex
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            X(wf%n_o + a, i) = -wf%t1_complex(a, i)
            Y(i, wf%n_o + a) = wf%t1_complex(a, i)
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
   end subroutine t1_transform_ccs_complex
!
!
   module subroutine t1_transpose_transform_ccs_complex(wf, Z_pq)
!!
!!    T1 transform transpose
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
!!    Assumes that Z is in the T1 basis and performs the transpose T1 transformation,
!!
!!       Z_pq <- sum_rs X_sp Z_sr Y_rq,    i.e.    Z <- X^T Z Y
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    Here, t1 is a full MO matrix whose only non-zero block is the vir-occ
!!    part, where it is equal to t_i^a.
!!
!!    Based on t1_transform_ccs_complex by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
      complex(dp), dimension(:,:), allocatable :: X, Y
!
      integer :: p, i, a
!
!     Construct the X and Y arrays
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      call zero_array_complex(X, (wf%n_mo)**2)
      call zero_array_complex(Y, (wf%n_mo)**2)
!
!$omp parallel do private(p)
      do p = 1, wf%n_mo
!
         X(p, p) = one_complex
         Y(p, p) = one_complex
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            X(wf%n_o + a, i) = -wf%t1_complex(a, i)
            Y(i, wf%n_o + a) = wf%t1_complex(a, i)
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
   end subroutine t1_transpose_transform_ccs_complex
!
!
   module subroutine t1_transform_4_ccs_complex(wf, Z_tuvw, Z_pqrs, t1)
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
      complex(dp), dimension(wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo), intent(in) :: Z_tuvw
      complex(dp), dimension(wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo), intent(out) :: Z_pqrs
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
!
      complex(dp), dimension(:,:), allocatable :: X, Y
      complex(dp), dimension(:,:,:,:), allocatable :: Z_tuvs, Z_puvs, Z_vspu, Z_vspq
!
      integer :: p, a, i
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      X = zero_complex
      Y = zero_complex
!
      do p = 1, wf%n_mo
!
         X(p, p) = one_complex
         Y(p, p) = one_complex
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
      call zgemm('N', 'T',   &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 wf%n_mo,    &
                 one_complex,        &
                 Z_tuvw,     &
                 wf%n_mo**3, &
                 Y,          &
                 wf%n_mo,    &
                 zero_complex,       &
                 Z_tuvs,     &
                 wf%n_mo**3)
!
      call mem%alloc(Z_puvs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call zgemm('N', 'N',   &
                 wf%n_mo,    &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 one_complex,        &
                 X,          &
                 wf%n_mo,    &
                 Z_tuvs,     &
                 wf%n_mo,    &
                 zero_complex,       &
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
      call zgemm('N', 'T',   &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 wf%n_mo,    &
                 one_complex,        &
                 Z_vspu,     &
                 wf%n_mo**3, &
                 Y,          &
                 wf%n_mo,    &
                 zero_complex,       &
                 Z_vspq,     &
                 wf%n_mo**3)
!
      call mem%dealloc(Z_vspu, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
!     Z_rspq = Z_pqrs because of particle permutation symmetry
!
      call zgemm('N', 'N',   &
                 wf%n_mo,    &
                 wf%n_mo**3, &
                 wf%n_mo,    &
                 one_complex,        &
                 X,          &
                 wf%n_mo,    &
                 Z_vspq,     &
                 wf%n_mo,    &
                 zero_complex,       &
                 Z_pqrs,     &
                 wf%n_mo)
!
      call mem%dealloc(Z_vspq, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
!
   end subroutine t1_transform_4_ccs_complex
!
!
end submodule t1_ccs_complex
