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
submodule (cc2_class) omega_cc2
!
!!
!!    Omega submodule (CC2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
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
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      type(timings) :: timer
!
      timer = new_timer('omega CC2')
      call timer%turn_on()
!
      call zero_array(omega, wf%n_gs_amplitudes)
!
      call wf%omega_ccs_a1(omega)
!
      call wf%construct_u()
!
      call wf%omega_doubles_a1(omega, wf%u)
      call wf%omega_doubles_b1(omega, wf%u)
      call wf%omega_doubles_c1(omega, wf%u)
!
      call timer%turn_off()
!
   end subroutine construct_omega_cc2
!
!
   module subroutine construct_omega2_cc2(wf, omega2, t)
!!
!!    Construct Omega2
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the doubles part of omega for CC2
!!
!!    Note that this is not used to solve CC2 equations, 
!!    but used for debug of Jacobian matrix.
!!
!!    omega2 = 1/Δ_aibj (g_aibj + e_aibj t_aibj)
!! 
!!    NOTE: made in the biorthonormal basis
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_t2), intent(out) :: omega2
      real(dp), dimension(wf%n_t2), intent(in) :: t
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer :: a, i, b, j, ai, bj, aibj, aiai
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_vovo(g_aibj)
!
      call packin(omega2, g_aibj, wf%n_v*wf%n_o)
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  ai = wf%n_v*(i - 1) + a
                  bj = wf%n_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai,bj) - 3)/2 + ai + bj
!
                  omega2(aibj) = omega2(aibj) + t(aibj)*&
                        (wf%orbital_energies(a + wf%n_o)&
                        +wf%orbital_energies(b + wf%n_o)&
                        -wf%orbital_energies(i)&
                        -wf%orbital_energies(j))
!
               enddo
            enddo
         enddo
      enddo
!
      do ai = 1, wf%n_v*wf%n_o
!
         aiai = ai*(ai - 3)/2 + ai*2
!
         omega2(aiai) = half*omega2(aiai)
!
      enddo
!
   end subroutine construct_omega2_cc2
!
end submodule omega_cc2
