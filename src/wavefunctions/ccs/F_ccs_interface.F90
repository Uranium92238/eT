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
   module subroutine F_transformation_ccs(wf, c)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019     
!!
!!    Directs the transformation by the F matrix.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine F_transformation_ccs
!
!
   module subroutine F_ccs_a1_0_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,0 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_A1,0_ai = 2 L_iajb c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_a1_0_ccs
!
!
   module subroutine F_ccs_a1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_A1,1_ai = - (F_ib tbar_aj + F_ja tbar_bi) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_a1_1_ccs
!
!
   module subroutine F_ccs_b1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation B1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_B1,1_ai = - (L_ikjb tbar_ak + L_jkia tbar_bk) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_b1_1_ccs
!
!
   module subroutine F_ccs_c1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation C1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_C1,1_ai = (L_cajb tbar_ci + L_cbia tbar_cj) c_bj
!!                = (X_iajb + X_jbia) c_bj
!!
!!    In batches of c
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_c1_1_ccs
