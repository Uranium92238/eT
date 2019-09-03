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
   module subroutine omega_ccsd_a2_ccsd(wf, omega_abij, t_abij)
!!
!!    Omega A2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(inout) :: omega_abij
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine omega_ccsd_b2_ccsd(wf, omega_abij, t_abij)
!!
!!    Omega B2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_abij
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   module subroutine omega_ccsd_c2_ccsd(wf, omega_aibj, t_aibj)
!!
!!    Omega C2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
   end subroutine omega_ccsd_c2_ccsd
!
!
   module subroutine omega_ccsd_d2_ccsd(wf, omega_aibj, t_aibj)
!!
!!    Omega D2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
   end subroutine omega_ccsd_d2_ccsd
!
!
   module subroutine omega_ccsd_e2_ccsd(wf, omega_aibj, t_aibj)
!!
!!    Omega E2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
   end subroutine omega_ccsd_e2_ccsd
!
!
   module subroutine construct_omega_ccsd(wf, omega)
!!
!!    Construct omega (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
         class(ccsd), intent(inout) :: wf
!
         real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_ccsd
!
!
