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
   module subroutine effective_jacobian_transformation_lowmem_cc2(wf, omega, c)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      real(dp), intent(in) :: omega 
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine effective_jacobian_transformation_lowmem_cc2
!
!
   module subroutine jacobian_cc2_A1_lowmem_cc2(wf, rho_ai, c_bj)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: rho_ai
!
    end subroutine jacobian_cc2_A1_lowmem_cc2
!
!
    module subroutine jacobian_cc2_B1_lowmem_cc2(wf, rho_ai, c_bj, eps_o, eps_v)
!!
!!    Jacobian CC2 B1
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_cc2_B1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_a1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Effective Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_a1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_b1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Effective jacobian B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in)  :: eps_o
      real(dp), dimension(wf%n_v), intent(in)  :: eps_v
!
   end subroutine effective_jacobian_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_c1_lowmem_cc2(wf, omega, rho_ai, c_cj, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_cj
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
end subroutine effective_jacobian_cc2_c1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_d1_lowmem_cc2(wf, omega, rho_ai, c_bl, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bl
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_d1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_e1_lowmem_cc2(wf, omega, rho_ai, c_dk, eps_o, eps_v)
!!
!!    Jacobian CC2 F1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_dk
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_e1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_f1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 f1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_f1_lowmem_cc2
!
!
