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
   module subroutine omega_doubles_a1_doubles(wf, omega, u)
!!
!!    Omega doubles A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd u_bicj g_abjc = sum_ckd u_bjc_i * g_a_bjc,
!!
!!    and adds it to the projection vector omega
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
   end subroutine omega_doubles_a1_doubles
!
!
   module subroutine omega_doubles_b1_doubles(wf, omega, u)
!!
!!    Omega doubles B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_ckl g_kb,ji * u_aj,bk,
!!
!!    with
!!
!!      u_aj_bk = 2t_aj,bk - t_ak,bj
!!
!!    and adds it to the projection vector (omega)
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
   end subroutine omega_doubles_b1_doubles
!
!
   module subroutine omega_doubles_c1_doubles(wf, omega, u)
!!
!!    Omega doubles C1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the C1 term,
!!
!!       C1: sum_bj u_ai,bj * F_jb,
!!
!!    with
!!
!!       u_ai_bj = 2*t_ai_bj - t_aj_bi
!!
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
    end subroutine omega_doubles_c1_doubles
