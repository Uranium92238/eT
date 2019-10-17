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
submodule (ccs_class) omega_ccs
!
!!
!!    Omega submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_omega_ccs(wf, omega)
!!
!!    Construct Omega (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      call zero_array(omega, wf%n_gs_amplitudes)
      call wf%omega_ccs_a1(omega)
!
   end subroutine construct_omega_ccs
!
!
   module subroutine omega_ccs_a1_ccs(wf, omega)
!!
!!    Omega A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, March 2017
!!
!!    Adds the A1 contribution to omega,
!!
!!       Omega_ai^A1 =+ F_ai_T1.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes) :: omega
!
      type(timings) :: omega_ccs_a1_timer
!
      omega_ccs_a1_timer = timings('omega ccs a1')
      call omega_ccs_a1_timer%turn_on()
!
      call daxpy((wf%n_o)*(wf%n_v), one, wf%fock_ai, 1, omega, 1)
!
      call omega_ccs_a1_timer%turn_off()
!
   end subroutine omega_ccs_a1_ccs
!
!
end submodule omega_ccs
