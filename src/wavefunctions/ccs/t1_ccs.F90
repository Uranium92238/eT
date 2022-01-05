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
      call symmetric_sandwich_right_transposition(Z_out, Z_mo,              &
                                                  wf%orbital_coefficients, &
                                                  wf%ao%n, wf%n_mo)
!
      call mem%dealloc(Z_mo, wf%n_mo, wf%n_mo)
!
   end subroutine add_t1_terms_and_transform_ccs
!
!
end submodule t1_ccs
