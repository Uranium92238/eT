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
   module subroutine construct_left_transition_density_doubles(wf, state)
!!
!!    Construct left one-electron EOM transition density
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^L_pq = < k| E_pq |CC >
!!
!!    where < k| is the left eigenvector of the Jacobian
!!    with amplitudes L_mu
!!
!!          < k| = sum_mu L_{k,mu} < mu| e^-T
!!
      implicit none
!
      class(doubles) :: wf
      integer, intent(in) :: state
!
   end subroutine construct_left_transition_density_doubles
!
!
   module subroutine construct_right_transition_density_doubles(wf, state)
!!
!!    Construct right one-electron EOM transition density
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^R_pq = < Lambda| E_pq |k >
!!
!!    where |k > is the right eigenvector of the Jacobian
!!    with amplitudes R_mu
!!
!!          |k > = sum_mu (tau_mu |CC > R_{k,mu} - tbar_mu |CC > R_{k,mu}) 
!!
      implicit none
!
      class(doubles) :: wf
      integer, intent(in) :: state
!
   end subroutine construct_right_transition_density_doubles
!
!
   module subroutine density_doubles_mu_nu_ov_doubles(wf, density, tbar_aibj, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant ov-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_kc += sum_abij R^a_i tbar^ab_ij (2t^bc_jk - t^bc_kj)
!!                   -sum_abij tbar^ab_ij (R^b_k t^ac_ij + R^c_j t^ab_ik)
!!      
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine density_doubles_mu_nu_ov_doubles
!
!
   module subroutine density_doubles_mu_nu_vo_doubles(wf, density, tbar_aibj, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant vo-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_bj += sum_ai R^a_i tbar^ab_ij
!!      
      implicit none
!
      class(doubles) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine density_doubles_mu_nu_vo_doubles
!
!
   module subroutine construct_eom_etaX_doubles(wf, X, csiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the EOM effective etaX vector, adding the EOM
!!    correction to etaX. 
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_eom_etaX_doubles
!
!
   module subroutine construct_etaX_doubles(wf, X, etaX)
!!
!!    Construct eta^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs left-hand-side vector etaX:
!!
!!       eta^X_mu = < Lambda| [X, tau_mu] |CC >
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_etaX_doubles
!
!
   module subroutine etaX_doubles_a1_doubles(wf, X, etaX_ai)
!!
!!    etaX CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!       A1 = - sum_ckdl (tb_ckal X_id t_ckdl + tb_ckdi X_la t_ckdl)
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
   end subroutine etaX_doubles_a1_doubles
!
!
   module subroutine etaX_doubles_a2_doubles(wf, X, etaX_aibj)
!!
!!    etaX CCSD A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!   
!!    Constructs the A2 term of the etaX vector
!! 
!!       A2 = 2 X_jb tb_ai - X_ib tb_aj
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
   end subroutine etaX_doubles_a2_doubles
!
!
   module subroutine etaX_doubles_b2_doubles(wf, X, etaX_aibj)
!!
!!    etaX CCSD B2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!   
!!    Constructs the B2 term of the etaX vector
!! 
!!       B2 = sum_c tb_aicj X_cb - sum_k tb_aibk X_jk
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
   end subroutine etaX_doubles_b2_doubles
!
!
   module subroutine construct_csiX_doubles(wf, X, csiX)
!!
!!    Construct csiX
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs xi^X_mu :
!!
!!       xi^X_mu = < mu| exp(-T) X exp(T)|R >
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: csiX
!
   end subroutine construct_csiX_doubles
!
!
   module subroutine csiX_doubles_a1_doubles(wf, X, csiX_ai)
!!
!!    csiX CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2010
!!
!!    Constructs the A1 term of csiX
!!
!!       A1 = sum_ck u_aick X_kc,
!!    
!!    where u_aick = 2t_ckai - t_ciak
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: csiX_ai
!
   end subroutine csiX_doubles_a1_doubles
!
!
   module subroutine csiX_doubles_a2_doubles(wf, X, csiX_aibj)
!!
!!    CsiX CCSD A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folekstad
!!
!!    Construct csiX A2 contribution:
!!
!!       A2 = sum_c t_aicj X_bc - sum_k t_aibk X_kj
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: csiX_aibj
!
   end subroutine csiX_doubles_a2_doubles
!
!
   module subroutine etaX_eom_doubles_a1_doubles(wf, X, etaX_ai)
!!
!!    etaX EOM CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs the A1 correction term to eta^X for EOM
!!
!!       A1 = sum_ck tb_aick X_ck + sum_ckdl tb_aick u_ckdl X_ld,
!!
!!    where u_ckdl = 2*t_ckdl - t_cldk
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
   end subroutine etaX_eom_doubles_a1_doubles
!
!
   module subroutine etaX_eom_a_doubles(wf, etaX, csiX)
!!
!!    Get eom contribution
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM contribution to etaX vector
!!
!!       EOM correction:  eta^X,corr_mu += tbar_mu (xi * tbar) 
!!
      implicit none
!
      class(doubles), intent(in) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
   end subroutine etaX_eom_a_doubles
