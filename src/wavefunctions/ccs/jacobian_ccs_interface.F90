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
   module subroutine prepare_for_jacobian_ccs(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_ccs
!
!
   module subroutine jacobian_transformation_ccs(wf, c)
!!
!!    Jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
!!    On exit, c is overwritten by rho.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine jacobian_transformation_ccs
!
!
   module subroutine jacobian_ccs_a1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_b F_ab c_bi - sum_j F_ji c_aj,
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1
!
   end subroutine jacobian_ccs_a1_ccs
!
!
   module subroutine jacobian_ccs_b1_ccs(wf, rho_ai, c_bj)
!!
!!    Jacobian CCS B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 term
!!
!!       A1: sum_bj (2 g_aijb - g_abji) * c_bj
!!
!!    and adds it to rho_ai.
!!
!!    Separate calculation of both terms due to batching.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine jacobian_ccs_b1_ccs
