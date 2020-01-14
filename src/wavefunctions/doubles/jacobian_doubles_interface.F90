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
   module subroutine save_jacobian_a1_intermediates_doubles(wf)
!!
!!    Save jacobian a1 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediates 
!!
!!       Y_jl   = t_ckdj L_kcld 
!!       Y_bd   = t_blck L_kcld 
!!
!!    Which are constructed in save_jacobian_a1_intermediates
!!    and stored on files
!!
!!       jacobian_a1_intermediate_oo
!!       jacobian_a1_intermediate_vv
!!
!!    which are wavefunction variables
!!
      implicit none
!
      class(doubles) :: wf
!
   end subroutine save_jacobian_a1_intermediates_doubles
!
!
   module subroutine jacobian_doubles_a1_doubles(wf, rho_ai, c_ai)
!!
!!    Jacobian doubles A1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^A1 = sum_ckdl L_kcld (u_ki^ca c_dl - t_kl^ad c_ci  - t_ki^cd c_al)
!!              = sum_ckdl L_kcld u_ki^ca c_dl - Y_ac c_ci - Y_il c_al)
!!
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o)                :: rho_ai
!
   end subroutine jacobian_doubles_a1_doubles
!
!
 module subroutine jacobian_doubles_b1_doubles(wf, rho_ai, c_aibj)
!!
!!    Jacobian doubles B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    rho_ai^B1 = sum_bj F_jb (2*c_aibj - c_ajbi)
!!              = sum_bj F_jb v_aijb
!!
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_doubles_b1_doubles
!
!
   module subroutine jacobian_doubles_c1_doubles(wf, rho_ai, c_aibj)
!!
!!    Jacobian doubles C1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^C1 = - sum_bjk L_jikb c_ajbk
!!              = - sum_bjk (2*g_jikb - g_kijb) c_ajbk
!!
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_doubles_c1_doubles
!
!
   module subroutine jacobian_doubles_d1_doubles(wf, rho_ai, c_bicj)
!!
!!    Jacobian doubles D1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bicj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_doubles_d1_doubles
!
!
   module subroutine jacobian_doubles_a2_doubles(wf, rho_aibj, c_ai)
!!
!!    Jacobian doubles A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^A2 = sum_c g_aibc c_cj - sum_k g_aikj c_bk
!!
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)          :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)      :: rho_aibj
!
   end subroutine jacobian_doubles_a2_doubles
