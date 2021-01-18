!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   module subroutine prepare_for_jacobian_lowmem_cc2(wf)
!!
!!    Prepare for Jacobian
!!    Written by Linda Goletto, Oct 2019
!!
!!    Gets occupied and virtual orbital energies and construcs 
!!    the jacobian_doubles_a1_doubles routine second and 
!!    third intermediates
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_lowmem_cc2
!
!
   module subroutine save_jacobian_a1_2_intermediate_lowmem_cc2(wf, eps_o, eps_v)
!!
!!    Save jacobian a1 second intermediate
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Constructs the intermediate 
!!
!!       X_ji   = L_kcjb t^cb_ki  
!!
!!    and stores it on the file:
!!
!!       jacobian_a1_intermediate_oo
!!
!!    which is a wavefunction variable
!!
!!    Modified by Anders Hutcheson, Oct 2019
!!
!!    Transferred here as separate subroutine in order to only
!!    compute X_ji once at the beginning of the calculation
!!
      implicit none
!
      class(lowmem_cc2) :: wf
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine save_jacobian_a1_2_intermediate_lowmem_cc2
!
!
   module subroutine save_jacobian_a1_3_intermediate_lowmem_cc2(wf, eps_o, eps_v)
!!
!!    Save jacobian a1 term 3 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Constructs the intermediate
!!
!!       X_ab = t_akcj L_kcjb 
!!
!!    and stores it on the file:
!!
!!       jacobian_a1_intermediate_vv
!!
!!    which is a wavefunction variable
!!
!!    Modified by Linda Goletto, Oct 2019
!!
!!    Transferred here as separate subroutine in order to only
!!    compute X_ab once at the beginning of the calculation
!!
      implicit none
!
      class(lowmem_cc2) :: wf
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine save_jacobian_a1_3_intermediate_lowmem_cc2
!
!
   module subroutine effective_jacobian_transformation_lowmem_cc2(wf, omega, c, rho)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Constructs the effective Jacobian transformation
!!    for lowmem CC2 according to
!!    
!!       C. Hättig and F. Weigend, J. Chem. Phys. 113, 5154 (2000).
      implicit none
!
      class(lowmem_cc2) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_t1), intent(in)  :: c
      real(dp), dimension(wf%n_t1), intent(out) :: rho
!
   end subroutine effective_jacobian_transformation_lowmem_cc2
!
!
   module subroutine jacobian_cc2_a1_lowmem_cc2(wf, rho_ai, c_bj, eps_o, eps_v)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 term
!!
!!       A1 = L_kcjb c_bj (2 t^ac_ik - t^ac_ki)
!!           - L_kcjb t^cb_ki c_aj - L_kcjb t^ca_kj c_bi
!!
!!    and adds it to rho_ai
!!
!!    Modified by Linda Goletto and Anders Hutcheson, Oct 2019
!!
!!    Reads two intermediates, which are prepared in prepare_for_jacobian
!!
!!       X_ab  =  L_kcjb t^ac_kj 
!!       X_ji   = L_kcjb t^cb_ki 
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_cc2_a1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_a1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 contribution using an implicit
!!    calculation of the doubles vector
!!
!!       A1 = F_kc * (-eps_ai,ck + w)^-1 * (2 g_aicd c_dk + 2 g_ckad c_di - g_akcd c_di - g_ciad c_dk)
!!          = F_kc * (-eps_ai,ck + w)^-1 * (2 X_aick - X_akci + 2 X_ckai - X_ciak)
!!          = F_kc * (Y_aick + Y_ckai)
!!
!!    and adds it to rho_ai
!!    The term is calculated in batches over the a and c indices.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
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
!!    Linda Goletto, and Alexander C. Paul, Jan 2019
!!
!!
!!    Calculates the B1 contribution using an implicit
!!    calculation of the doubles vector
!!
!!       B1 = - 2 sum_{kcl} F_kc (1/(ε_{aick} + ω)) * (g_ailk c_cl + g_ckli c_al)
!!            + sum_{kcl} F_kc (1/(ε_{akci} + ω)) * (g_akli c_cl + g_cilk c_al)
!!          = 2 sum_{kcl} F_kc (- 2*X_ckai - 2*X_aick + X_ciak + X_akci)
!!
!!    and adds it to rho_ai     
!!
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_o), intent(in)  :: eps_o
      real(dp), dimension(wf%n_v), intent(in)  :: eps_v
!
   end subroutine effective_jacobian_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_c1_lowmem_cc2(wf, omega, rho_ai, c_cj, eps_o, eps_v)
!!
!!    Jacobian CC2 C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the C1 contribution using an implicit
!!    calculation of the doubles vector
!!
!!       C1 = sum_ckbj - L_kijb  (g_akbc * c_cj + g_bjac * c_ck) (omega - ε_akbj)^-1
!!          = sum_kjb - L_kijb  (X_akbj + X_bjak)
!!          = sum_kjb - L_kijb Y_a_kjb
!!
!!    and adds it to the rho_ai vector
!!
!!    Modified by Linda Goletto and Anders Hutcheson, Oct 2019
!!
!!    Integrals, g_akbc and g_bjac, moved out of i batching loop 
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_cj
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_c1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_d1_lowmem_cc2(wf, omega, rho_ai, c_bl, eps_o, eps_v)
!!
!!    Jacobian CC2 D1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the D1 contribution using an implicit
!!    calculation of the doubles vector
!!
!!       D1 = sum_ckbj - L_kijb  (- g_aklj * c_bl - g_bjlk * c_al) (omega - ε_akbj)^-1 
!!          = sum_kjb L_kijb  (X_bjak + X_akbj)
!!          = sum_kjb L_kijb  Y_ajbk
!!
!!    and adds it to the rho_ai vector.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bl
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_d1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_e1_lowmem_cc2(wf, omega, rho_ai, c_dk, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Jan 2019
!!
!!    Calculates the E1 contribution using an implicit
!!    calculation of the doubles vector
!!
!!       E1 = sum_bckd L_abkc  (g_bicd * c_dk + g_ckbd * c_di) (omega - ε_bick)^-1
!!          = sum_bck L_abkc  (X_bick + X_ckbi)
!!
!!    and adds it to the rho_ai vector.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_dk
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_e1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_f1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 effective F1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Jan 2019
!!
!!    Calculates the F1 contribution using an implicit
!!    calculation of the doubles vector
!!
!!       F1 = - L_abkc (1/(ε_{bick} + ω) * (g_lkbi c_cl + g_lick c_bl))
!!
!!    and adds it to the rho_ai vector.
!!
!!    Modified by Linda Goletto and Anders Hutcheson, Oct 2019
!!
!!    Integral, g_abkc, moved out of i batching loop 
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_f1_lowmem_cc2
