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
submodule (cc3_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (cc3)
!!    Written by Alexander Paul and Rolf H. Myhre March 2019
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | ν >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, c)
!!
!!    Effective Jacobian transpose transformation (CC3)
!!    Written by Alexander Paul and Rolf H. Myhre March 2019
!!
!!    Directs the transformation by the transpose of the  CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >,
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine. On exit, the vector c is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij, c_abji
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj, sigma_abij
!
      integer :: i, j, a, b, ai, bj, aibj, b_end ! Index
!
      type(timings) :: cc3_timer
      type(timings) :: ccsd_timer
!
      call cc3_timer%init('CC3 contribution)')
      call ccsd_timer%init('CCSD contribution)')
!
!     Allocate and zero the transformed singles vector
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      sigma_ai = zero
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai)
!
         enddo
      enddo
!
!     :: CCS contributions to the transformed singles vector ::
!
      call ccsd_timer%start()
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, c_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, c_ai)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_transpose_ccsd_a1(sigma_ai, c_ai)
      call wf%jacobian_transpose_ccsd_b1(sigma_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj, b_end)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, i
!
               if (i .ne. j) then
                  b_end = wf%n_v
               else
                  b_end = a
               endif
!
               do b = 1, b_end
!
                  bj = wf%n_v*(j - 1) + b
!
                  aibj = ai*(ai-3)/2 + ai + bj
!
                  c_aibj(a,i,b,j) = c(wf%n_t1 + aibj)
                  c_aibj(b,j,a,i) = c(wf%n_t1 + aibj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_transpose_ccsd_c1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_d1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_e1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_f1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_g1(sigma_ai, c_aibj)
!
!     :: CCSD contributions to the transformed doubles vector ::
!     Allocate unpacked transformed vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      sigma_aibj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_transpose_ccsd_a2(sigma_aibj, c_ai)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_transpose_ccsd_b2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_c2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_d2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_e2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_f2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_g2(sigma_aibj, c_aibj)
!
!     Compute CC3 contributions to sigma_ai and sigma_aibj and symmetrise sigma_aibj
!     CCSD H2 and I2 are already symmetric and will be computed afterwards
!
      call ccsd_timer%freeze()
!
      call mem%alloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(c_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1342(c_aibj, c_abji, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(sigma_aibj, sigma_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call cc3_timer%start()
      call wf%jacobian_transpose_cc3_A(omega, c_ai, c_abji, sigma_ai, sigma_abij)
      call cc3_timer%freeze()
      call cc3_timer%switch_off()
!
!     Done with singles vector c; Overwrite the incoming singles c vector for exit
!
      call ccsd_timer%start()
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, sigma_ai, 1, c, 1)
!
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     Last two CCSD-terms (H2, I2) are already symmetric.
!     Perform the symmetrization sigma_aibj = P_ij^ab sigma_aibj
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(sigma_abij, sigma_aibj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call symmetric_sum(sigma_aibj, (wf%n_v)*(wf%n_o))
!
      call sort_1234_to_1324(sigma_aibj, sigma_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1243(c_abji, c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(c_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Compute CCSD H2 and I2 contributions
!
      call wf%jacobian_transpose_ccsd_h2(sigma_abij, c_abij)
      call wf%jacobian_transpose_ccsd_i2(sigma_abij, c_abij)
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call ccsd_timer%freeze()
      call ccsd_timer%switch_off()
!
!     overwrite the incoming, packed doubles c vector for exit
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj, b_end)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, i
!
               if (j .ne. i) then
                  b_end = wf%n_v
               else
                  b_end = a
               endif
!
               do b = 1, b_end
!
                  bj = wf%n_v*(j - 1) + b
!
                  aibj = ai*(ai-3)/2 + ai + bj
!
                  c(wf%n_t1 + aibj) = sigma_abij(a,b,i,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_A_cc3(wf, omega, c_ai, c_abji, sigma_ai, sigma_abij)
!!
!!    Terms of the transpose of the  CC3 Jacobi matrix
!!
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abji
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
   end subroutine jacobian_transpose_cc3_A_cc3
!
!
end submodule jacobian_transpose
