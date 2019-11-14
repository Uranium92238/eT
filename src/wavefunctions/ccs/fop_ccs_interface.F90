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
   module subroutine construct_right_transition_density_ccs(wf, state)
!!
!!    Construct right one-electron transition density (EOM)
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^R_pq = < Lambda| E_pq |k >
!!
!!    where |k > is the right eigenvector of the Jacobian
!!    with amplitudes R_mu
!!
!!          |k > = sum_mu (tau_mu R_{k,mu} - tbar_mu R_{k,mu}) |CC >
!!
!!    Contributions to the right transition density are split as follows:
!!
!!          D^R_pq = sum_mu D_pq(ref-mu) R_mu + sum_mu tbar_mu D_pq(mu-ref)
!!                   + sum_mu,nu tbar_mu D_pq(mu-nu) R_nu
!!
      implicit none
!
      class(ccs) :: wf
!
      integer, intent(in) :: state
!
   end subroutine construct_right_transition_density_ccs
!
!
   module subroutine density_ccs_mu_nu_oo_ccs(wf, density, tbar_ai, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant oo-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_ij = -sum_a R_ai tbar_aj
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine density_ccs_mu_nu_oo_ccs
!
!
   module subroutine density_ccs_ref_mu_ov_ccs(wf, density, R_ai)
!!
!!    One electron density (EOM) reference/excited-determinant ov-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu < HF| e^(-T) E_pq e^T |nu > Y_mu
!!
!!    explicit term in this routine:
!!          D^R_ia = 2*R^a_i 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine density_ccs_ref_mu_ov_ccs
!
!
   module subroutine density_ccs_mu_nu_vv_ccs(wf, density, tbar_ai, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant vv-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_ab = sum_i R_bi tbar_ai
!!      
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine density_ccs_mu_nu_vv_ccs
!
!
   module subroutine density_mu_mu_oo_ccs(wf, density, scaling_factor)
!!
!!    One electron density (EOM) excited determinant/excited determinant
!!    Written by Alexander C. Paul, June 2019
!!
!!    Should be reusable in the whole CC hierachy
!!
!!    Computes the following term containing the same excited contributions (mu,mu) 
!!    of the left and right vector
!!
!!          D_pq += sum_mu X_mu < mu| tau_mu e^(-T) E_pq e^T |HF > Y_mu
!!               += sum_mu X_mu *Y_mu < HF| e^(-T) E_pq e^T |HF >
!!
!!    This is the Hartree fock density scaled by the overlap X_mu*Y_mu
!!
!!    For the right transition density of any level in CC:
!!
!!       D^R_pq +=  sum_mu R_{k,mu} tbar_mu * < HF| e^-T E_pq e^T |HF >
!!              +=  sum_mu R_{k,mu} tbar_mu * 2 delta_pq delta_p,occ 
!!
!!       scaling_factor: overlap of tbar- and R-vector
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), intent(in) :: scaling_factor
!
   end subroutine density_mu_mu_oo_ccs
!
!
   module subroutine density_mu_ref_ccs(wf, density_out, density_in, scaling_factor)
!!
!!    One electron density (EOM) excited determinant/reference term
!!    Written by Alexander C. Paul, June 2019
!!
!!    Should be reusable in the whole CC hierachy
!!
!!    Computes term that is a density matrix scaled by Y_ref
!!
!!          D_pq += sum_mu X_mu < mu| e^(-T) E_pq e^T |HF > Y_ref
!!
!!    For the right transition density of any level in CC:
!!
!!       D^R_pq -= sum_mu,nu R_{k,mu} tbar_mu* tbar_nu < nu| e^-T E_pq e^T |HF >
!!              -= sum_mu R_{k,mu} tbar_mu D^GS_pq
!!
!!       density_out   : right_transition density
!!       density_in    : ground state density 
!!       scaling_factor: overlap of tbar- and R-vector
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density_out
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)    :: density_in
!
      real(dp), intent(in) :: scaling_factor
!
   end subroutine density_mu_ref_ccs
!
!
   module subroutine construct_left_transition_density_ccs(wf, state)
!!
!!    Construct left one-electron transition density for the state k
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^L_pq = < k| E_pq |CC >
!!
!!    where <k| is the left eigenvector of the Jacobian
!!    with amplitudes L_mu
!!
!!          < k| = sum_mu L_{k,mu} < mu | e^-T
!!
      implicit none
!
      class(ccs) :: wf
!
      integer, intent(in) :: state
!
   end subroutine construct_left_transition_density_ccs
!
!
   module subroutine construct_eom_etaX_ccs(wf, X, csiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the EOM effective etaX vector, adding the EOM
!!    correction to etaX. 
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_eom_etaX_ccs
!
!
   module subroutine construct_etaX_ccs(wf, X, etaX)
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
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_etaX_ccs
!
!
   module subroutine etaX_ccs_a1_ccs(wf, X, etaX_ai)
!!
!!    Construct etaX A1 
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Adds the A1 term of eta_ai^X:
!!
!!       A1 = 2X_ia
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
   end subroutine etaX_ccs_a1_ccs
!
!
   module subroutine etaX_ccs_b1_ccs(wf, X, etaX_ai)
!!
!!    Construct etaX B1
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!    Adds the B1 term of eta_ai^X:
!!
!!       B1 = sum_c tb_ci X_ca - sum_k tb_ak X_ik
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)   :: etaX_ai
!
   end subroutine etaX_ccs_b1_ccs
!
!
   module subroutine construct_csiX_ccs(wf, X, csiX)
!!
!!    Construct csiX
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs x^X_mu :
!!
!!       x^X_mu = < mu| exp(-T) X exp(T)|R >
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: csiX
!
   end subroutine construct_csiX_ccs
!
!
   module subroutine csiX_ccs_a1_ccs(wf, X, csiX_ai)
!!
!!    Construct csiX A1 
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adds the A1 term to csiX:
!! 
!!       x^X_ai += X_ai
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: csiX_ai
!
   end subroutine csiX_ccs_a1_ccs
!
!
   module subroutine etaX_eom_a_ccs(wf, etaX, csiX)
!!
!!    EtaX EOM A
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM A correction to etaX vector:
!!
!!       A:  eta^X,corr_mu += tbar_mu (x * tbar) 
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
   end subroutine etaX_eom_a_ccs
!
!
   module subroutine calculate_transition_strength_ccs(wf, S, etaX, csiX, state, T_l, T_r)
!!
!!    Calculate transition strength
!!    Written by Josefine H. Andersen, February 2019
!!
!!    Given etaX and csiX, this routine calculates the left and right transition 
!!    moments T_l and T_r for the state number "state" and the transition strength 
!!    S = T_l * T_r.
!! 
!!    The left and right states L and R are read from file and made binormal by the routine.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), intent(inout) :: S
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
!
      real(dp), intent(out) :: T_l, T_r
      integer, intent(in)   :: state
!
   end subroutine calculate_transition_strength_ccs
