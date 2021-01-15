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
   module subroutine prepare_for_jacobian_mlcc2(wf)
!!
!!    Prepare for Jacobian
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
      implicit none 
!
      class(mlcc2), intent(inout) :: wf 
!
   end subroutine prepare_for_jacobian_mlcc2
!
!
   module subroutine jacobian_transformation_mlcc2(wf, c, rho)
!!
!!    Jacobian transformation (mlcc2)
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_ai = rho_ai,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: c
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: rho
!
   end subroutine jacobian_transformation_mlcc2
!
!
   module subroutine jacobian_cc2_a1_mlcc2(wf, rho_ai, c_ai, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v)
!!
!!    Jacobian CC2 A1
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Calculates the A1 term
!!
!!       A1:   sum_bjck L_kcjb u_aick c_bj 
!!           - sum_bjck g_jbkc u_ckbi c_aj 
!!           - sum_bjck g_kcjb u_ckaj c_bi
!!
!!    and adds it to rho_ai.
!!
!!    Index restrictions:
!!
!!       Term 1:
!!       
!!          a, i, c, k : CC2 orbitals
!!
!!          b, j : unrestricted
!!
!!       Term 2:
!!       
!!          b, i, c, k : CC2 orbitals
!!
!!          a, j : unrestricted
!!
!!       Term 3:
!!       
!!          a, j, c, k : CC2 orbitals
!!
!!          b, i : unrestricted
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)   :: rho_ai
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
   end subroutine jacobian_cc2_a1_mlcc2
!
!
   module subroutine jacobian_cc2_a2_mlcc2(wf, rho_aibj, c_ai, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Jacobian CC2 A2
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Calculates the A2 term 
!!
!!       A2: sum_c g_aibc c_cj - sum_k g_aikj c_bk, 
!!
!!    and adds it to rho_aibj.
!!
!!    The first term is calculated in batches over c.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
!!       Term 1: 
!!
!!          c : unrestricted
!!
!!       Term 2: 
!!
!!          k : unrestricted
!!
      implicit none
!
      class(mlcc2) :: wf
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                       :: c_ai
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(out)  :: rho_aibj   
!
   end subroutine jacobian_cc2_a2_mlcc2
!
!
   module subroutine jacobian_cc2_b1_mlcc2(wf, rho_ai, c_aibj, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Jacobian CC2 B1
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Calculates the B1 term
!!
!!       B1: 2 sum_bj F_jb c_aibj - F_jb c_ajbi  
!!           - sum_kjb L_kijb c_akbj + sum_bkc L_abkc c_bick
!!
!!    And adds it to rho_ai.
!!
!!    The fourth term is calculated in batches of index a.
!!    
!!    Index restrictions:
!!
!!       Terms 1 and 2:
!!
!!          a, i, b, j : CC2 orbitals
!!
!!       Term 3:
!!
!!          a, k, b, j : CC2 orbitals
!!
!!          i : unrestricted
!!
!!       Term 4:
!!
!!          c, k, b, i : CC2 orbitals
!!
!!          a : unrestricted
!!
      implicit none
!
      class(mlcc2) :: wf
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(in)  :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                     :: rho_ai   
!
   end subroutine jacobian_cc2_b1_mlcc2
!
!
   module subroutine jacobian_cc2_b2_mlcc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Constructs the B2 term
!!
!!       B2: ε_aibj c_aibj/(1/Δ_aibj) 
!!
!!    and adds it to rho_aibj.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)     :: c_aibj
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)     :: rho_aibj   
!
   end subroutine jacobian_cc2_b2_mlcc2
!
!
   module subroutine save_jacobian_a1_intermediates_mlcc2(wf, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v)
!!
!!    Save jacobian a1 intermediates
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Constructs the intermediates 
!!
!!       Y_ji = - sum_bck g_jbkc u_ckbi
!!       Y_ab = - sum_jck u_ajck g_kcjb
!!      
!!    Which are constructed in save_jacobian_a1_intermediates
!!    and stored on files
!!
!!       jacobian_a1_intermediate_oo
!!       jacobian_a1_intermediate_vv
!!
!!    which are wavefunction variables
!!
!!    Index restrictions:
!!
!!       oo intermediate:
!!    
!!          c, k, b, i : CC2 orbitals
!!
!!          j : unrestricted
!!
!!       vv intermediate:
!!
!!          a, j, k, c : CC2 orbitals
!!
!!          b : unrestricted 
!!
      implicit none
!
      class(mlcc2) :: wf
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
   end subroutine save_jacobian_a1_intermediates_mlcc2
