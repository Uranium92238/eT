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
   module subroutine jacobian_doubles_b1_abstract_doubles(wf, rho_ai, c_aibj)
!!
!!    Jacobian doubles B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_doubles_b1_abstract_doubles
!
!
   module subroutine jacobian_doubles_c1_abstract_doubles(wf, rho_ai, c_aibj)
!!
!!    Jacobian doubles C1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_doubles_c1_abstract_doubles
!
!
   module subroutine jacobian_doubles_d1_abstract_doubles(wf, rho_ai, c_bicj)
!!
!!    Jacobian doubles D1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bicj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_doubles_d1_abstract_doubles
!
!
   module subroutine save_jacobian_a1_intermediates_abstract_doubles(wf)
!!
!!    Save jacobian a1 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
   end subroutine save_jacobian_a1_intermediates_abstract_doubles
!
!
   module subroutine jacobian_doubles_a1_abstract_doubles(wf, rho_ai, c_ai)
!!
!!    Jacobian doubles A1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o)                :: rho_ai
!
   end subroutine jacobian_doubles_a1_abstract_doubles
!
!
   module subroutine jacobian_doubles_a2_abstract_doubles(wf, rho_aibj, c_ai)
!!
!!    Jacobian doubles A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)          :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)      :: rho_aibj
!
   end subroutine jacobian_doubles_a2_abstract_doubles
!
