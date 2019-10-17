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
submodule (ccs_class) oei_ccs
!
!!
!!    One-electron integrals submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!!    Submodule containing routines that can be used to construct t1-transformed one-electron integrals.
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
   module subroutine construct_mu_ccs(wf, mu_pqk)
!!
!!    Construct mu
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Constructs 
!!
!!       mu_pqk, k = 1, 2, 3 (x, y, z)
!!
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(out) :: mu_pqk 
!
      call wf%get_mo_mu(mu_pqk)
!
      call wf%t1_transform(mu_pqk(:,:,1))
      call wf%t1_transform(mu_pqk(:,:,2))
      call wf%t1_transform(mu_pqk(:,:,3))
!
   end subroutine construct_mu_ccs
!
!
   module subroutine construct_h_ccs(wf, h_pq)
!!
!!    Construct h
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    Constructs the one-electron Hamiltonian integrals h_pq 
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: h_pq 
!
      call wf%get_mo_h(h_pq)
      call wf%t1_transform(h_pq)
!
   end subroutine construct_h_ccs
!
!
   module subroutine construct_q_ccs(wf, q_pqk)
!!
!!    Construct q
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    Constructs 
!!
!!       q_pqk, k = 1, 2, 3, 4, 5, 6 (xx, xy, xz, yy, yz, and zz)
!!
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 6), intent(out) :: q_pqk 
!
      call wf%get_mo_q(q_pqk)
!
      call wf%t1_transform(q_pqk(:,:,1))
      call wf%t1_transform(q_pqk(:,:,2))
      call wf%t1_transform(q_pqk(:,:,3))
      call wf%t1_transform(q_pqk(:,:,4))
      call wf%t1_transform(q_pqk(:,:,5))
      call wf%t1_transform(q_pqk(:,:,6))
!
   end subroutine construct_q_ccs
!
!
end submodule oei_ccs
