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
submodule (cc2_class) triplet_jacobian_cc2
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ |exp(-T) [H, τ_ν] exp(T) | R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_triplet_jacobian_cc2(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      call wf%initialize_t2()
      call wf%construct_t2()
!
   end subroutine prepare_for_triplet_jacobian_cc2
!
!
   module subroutine triplet_jacobian_transformation_cc2(wf, c, rho)
!!
!!    Triplet Jacobian transformation
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
!!
!!       A_μ,ν = < μ |exp(-T) [H, τ_ν] exp(T) | R >,
!!
!!    Where τ_ν is a triplet excitation operator,
!!    and < μ | = 1/2 < R | τ_μ^dagger
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
      use array_initialization, only: zero_array
      use reordering, only: symmetric_sum, antisymmetric_sum, sort_1234_to_1432, squareup
!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_triplet_amplitudes), intent(in) :: c
      real(dp), dimension(wf%n_triplet_amplitudes), intent(out) :: rho
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj_p, c_aibj_m
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj_pm
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj_p, rho_aibj_m, rho_ajbi_p
!
      type(timings), allocatable :: timer
!
      timer = timings('Triplet Jacobian cc2 transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(rho, wf%n_triplet_amplitudes)

      call wf%ccs%triplet_jacobian_transformation(c(1:wf%n_t1), rho(1:wf%n_t1))
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_v*wf%n_o)
      call wf%triplet_jacobian_s_s_b(rho(1:wf%n_t1), c(1:wf%n_t1), t_aibj)
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(c_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%packout_triplet_d(c, c_aibj_p, c_aibj_m)
!
      call mem%alloc(c_aibj_pm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call dcopy(wf%n_o**2*wf%n_v**2, c_aibj_p, 1, c_aibj_pm, 1)
      call daxpy(wf%n_o**2*wf%n_v**2, one, c_aibj_m, 1, c_aibj_pm, 1)
!
      call wf%triplet_jacobian_s_d_a(rho(1:wf%n_t1), c_aibj_pm)
      call wf%triplet_jacobian_s_d_b(rho(1:wf%n_t1), c_aibj_pm)
      call wf%triplet_jacobian_s_d_c(rho(1:wf%n_t1), c_aibj_pm)
!
      call mem%dealloc(c_aibj_pm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
      call mem%alloc(rho_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
      call wf%triplet_jacobian_d_s_a(rho_aibj_p, c(1:wf%n_t1))
      call wf%triplet_jacobian_d_s_a(rho_aibj_m, c(1:wf%n_t1))
!
      call symmetric_sum(rho_aibj_p, wf%n_o*wf%n_v)
      call antisymmetric_sum(rho_aibj_m, wf%n_o*wf%n_v)
!
!     ~P_ij rho_aibj_p = rho_aibj_p - rho_ajbi_p
      call mem%alloc(rho_ajbi_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(rho_aibj_p, rho_ajbi_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy(wf%n_v**2*wf%n_o**2, -one, rho_ajbi_p, 1, rho_aibj_p, 1)
      call mem%dealloc(rho_ajbi_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%triplet_jacobian_d_d_a(rho_aibj_p, c_aibj_p)
      call wf%triplet_jacobian_d_d_a(rho_aibj_m, c_aibj_m)
!
      call mem%dealloc(c_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%packin_triplet_d(rho, rho_aibj_p, rho_aibj_m)
!
      call mem%dealloc(rho_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine triplet_jacobian_transformation_cc2
!
!
   module subroutine triplet_jacobian_d_d_a_cc2(wf, rho_aibj, c_aibj)
!!
!!    Triplet Jacobian  doubles-doubles A
!!    Written by  Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_aibj += c_aibj * (e_a + e_b - e_i - e_j)
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)  :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(:), allocatable :: eps_o, eps_v
!
      integer :: a, i, b, j
!
      call mem%alloc(eps_o, wf%n_o)
      call mem%alloc(eps_v, wf%n_v)
!
      eps_o(1:wf%n_o) = wf%orbital_energies(1:wf%n_o)
      eps_v(1:wf%n_v) = wf%orbital_energies(wf%n_o+1:wf%n_mo)
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  rho_aibj(a, i, b, j) = rho_aibj(a, i, b, j) &
                                       + c_aibj(a, i, b, j)*(eps_v(a) - eps_o(i) &
                                                            + eps_v(b) - eps_o(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(eps_o, wf%n_o)
      call mem%dealloc(eps_v, wf%n_v)
!
   end subroutine triplet_jacobian_d_d_a_cc2
!
!
end submodule triplet_jacobian_cc2
