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
submodule (cc2_class) omega_cc2
!
!!
!!    Omega submodule
!!
!!    Routines to construct
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(out) :: omega
!
      type(timings) :: timer
!
      timer = timings('Construct CC2 omega', pl='normal')
      call timer%turn_on()
!
      call zero_array(omega, wf%n_t1)
!
      call wf%ccs%construct_omega(omega)
!
      call wf%construct_u_aibj()
!
      call wf%omega_doubles_a1(omega, wf%u_aibj)
      call wf%omega_doubles_b1(omega, wf%u_aibj)
      call wf%omega_doubles_c1(omega, wf%u_aibj)
!
      call timer%turn_off()
!
   end subroutine construct_omega_cc2
!
!
end submodule omega_cc2
