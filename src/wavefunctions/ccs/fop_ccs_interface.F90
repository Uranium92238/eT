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
   module subroutine construct_right_transition_density_ccs(wf, R_k)
!!
!!    Construct right one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
!!          ρ^R_pq = < Λ | E_pq | k >
!!
!!    where |k> is the right eigenvector of the Jacobian
!!    with amplitudes R_μ
!!
!!          | k > = sum_μ (τ_μ | CC > R_{k,μ} - tbar_μ | CC > R_{k,μ}) 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R_k
!
   end subroutine construct_right_transition_density_ccs
!
!
   module subroutine right_transition_density_ccs_oo_ccs(wf, tbar_ai, R_ai)
!!
!!    Right transition density oo contribution
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_ij = -sum_a R_ai tbar_aj
!!      
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_ccs_oo_ccs
!
!
   module subroutine right_transition_density_ccs_ov_ccs(wf, R_ai)
!!
!!    Right transition density ov
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_pq = 2*R_ai 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_ccs_ov_ccs
!
!
   module subroutine right_transition_density_ccs_vv_ccs(wf, tbar_ai, R_ai)
!!
!!    Right transition density vv contribution
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_ab = sum_i R_bi tbar_ai
!!      
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_ccs_vv_ccs
!
!
   module subroutine right_transition_density_ccs_gs_contr_ccs(wf, tbar_ai, R_ai)
!!
!!    Right transition density, contribution from the ground state density
!!    ρ^R_pq -= sum_μν R_{k,μ}tbar_μ tbar_ν < ν |e^-T E_pq e^T| HF >
!!           -= sum_μ R_{k,μ}tbar_μ (D_GS - D_HF)
!!    CCS:
!!    ρ^R_pq = ρ^R_ai = sum_bj R^b_j tbar^b_j tbar^a_i
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_ccs_gs_contr_ccs
!
!
   module subroutine construct_left_transition_density_ccs(wf, L_k)
!!
!!    Construct left one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
!!          ρ^L_pq = < k | E_pq | CC >
!!
!!    where <k| is the left eigenvector of the Jacobian
!!    with amplitudes L_μ
!!
!!          < k | = sum_μ L_{k,μ} < μ | e^-T
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L_k
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
!!    Construct η^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs left-hand-side vector etaX:
!!
!!       η^X_μ = < Λ | [X, τ_μ] | CC >
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
!!    Adds the A1 term of η_ai^X:
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
!!    Adds the B1 term of η_ai^X:
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
!!    Constructs ξ^X_μ :
!!
!!       ξ^X_μ = < μ | exp(-T) X exp(T)| R >
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
!!       ξ^X_ai =+ X_ai
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
!!       A:  η^X,corr_μ += tbar_μ (ξ * tbar) 
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
