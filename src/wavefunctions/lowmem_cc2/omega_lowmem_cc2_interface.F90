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
   module subroutine construct_omega_lowmem_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Direqts the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_lowmem_cc2
!
!
   module subroutine omega_cc2_a1_lowmem_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd u_bj_ci * g_abjc,
!!
!!    with 
!!       
!!       u_bj_ci = 2*t_bj_ci - t_bi_cj
!!
!!    and
!!
!!       t_bj_ci = - g_bjci/ε^{bc}_{ji}
!!
!!    and adds it to the projection vector omega
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine omega_cc2_a1_lowmem_cc2
!
!
   module subroutine omega_cc2_b1_lowmem_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_ckl (2g_kb_ji - g_jb_ki) * t_aj_bk,
!!
!!    with
!!
!!       t_aj_bk = - g_ajbk/ε^{ab}_{jk}
!!
!!    and adds it to the projection vector (omega) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine omega_cc2_b1_lowmem_cc2
!
!
   module subroutine omega_cc2_c1_lowmem_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the C1 term,
!!
!!       C1: sum_bj u_ai_bj * F_{jb},
!!
!!    with 
!!       
!!       u_ai_bj = 2*t_ai_bj - t_aj_bi
!!
!!    and
!!
!!       t_ai_bj = - g_aibj/ε^{ab}_{ij}
!!
!!    and adds it to the projection vector (omega) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine omega_cc2_c1_lowmem_cc2
