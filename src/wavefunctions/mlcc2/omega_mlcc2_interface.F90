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
   module subroutine construct_omega_mlcc2(wf, omega)
!!
!!    Construct omega
!!    Written by Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_mlcc2
!
!
   module subroutine omega_cc2_a1_mlcc2(wf, omega, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Omega MLCC2 A1 term
!!    Adapted by Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
!
   end subroutine omega_cc2_a1_mlcc2
!
!
   module subroutine omega_cc2_b1_mlcc2(wf, omega, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Omega MLCC2 B1 term
!!    Adapted by Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
!
   end subroutine omega_cc2_b1_mlcc2
!
!
   module subroutine omega_cc2_c1_mlcc2(wf, omega, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v)
!!
!!    Omega MLCC2 C1 term
!!    Adapted by Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
    end subroutine omega_cc2_c1_mlcc2
!
!
   module subroutine construct_omega_doubles_mlcc2(wf, omega2)
!!
!!    Construct Omega doubles
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_x2), intent(out) :: omega2
!
   end subroutine construct_omega_doubles_mlcc2