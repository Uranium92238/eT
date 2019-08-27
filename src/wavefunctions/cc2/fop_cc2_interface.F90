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
   module subroutine construct_left_transition_density_cc2(wf, L_k)
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
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L_k
!      
   end subroutine construct_left_transition_density_cc2
!
!
   module subroutine construct_right_transition_density_cc2(wf, R_k)
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
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R_k
!
   end subroutine construct_right_transition_density_cc2
!
!
   module subroutine right_transition_density_cc2_ov_cc2(wf, tbar_aibj, R_ai)
!!
!!    Right transition density ov contribution (cc2)
!!    from the singles part of the excitation vector
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_IA = sum_abij R^a_i tbar^ab_ij (2t^bA_jI - t^bA_Ij)
!!             -sum_abij tbar^ab_ij (R^b_I t^aA_ij + R^A_j t^ab_iI)
!!      
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_cc2_ov_cc2
!
!
   module subroutine right_transition_density_cc2_vo_cc2(wf, tbar_aibj, R_ai)
!!
!!    Right transition density ov contribution (cc2)
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_AI = sum_ai R^a_i tbar^aA_iI
!!      
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_cc2_vo_cc2
!
!
   module subroutine construct_eom_etaX_cc2(wf, X, csiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_eom_etaX_cc2
!
   module subroutine construct_etaX_cc2(wf, X, etaX)
!!
!!    Construct etaX
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folekstad
!!
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_etaX_cc2
!
!
   module subroutine etaX_cc2_a1_cc2(wf, X, etaX_ai)
!!
!!    etaX cc2 A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
   end subroutine etaX_cc2_a1_cc2
!
!
   module subroutine etaX_cc2_a2_cc2(wf, X, etaX_aibj)
!!
!!    etaX cc2 A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
   end subroutine etaX_cc2_a2_cc2
!
!
   module subroutine etaX_cc2_b2_cc2(wf, X, etaX_aibj)
!!
!!    etaX cc2 B2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!! 
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
   end subroutine etaX_cc2_b2_cc2
!
!
   module subroutine construct_csiX_cc2(wf, X, csiX)
!!
!!    Construct csiX (cc2)
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: csiX
!
   end subroutine construct_csiX_cc2
!
!
   module subroutine csiX_cc2_a1_cc2(wf, X, csiX_ai)
!!
!!    csiX cc2 A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: csiX_ai
!
   end subroutine csiX_cc2_a1_cc2
!
!
   module subroutine csiX_cc2_a2_cc2(wf, X, csiX_aibj)
!!
!!    CsiX cc2 A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folekstad
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: csiX_aibj
!
   end subroutine csiX_cc2_a2_cc2
!
!
   module subroutine etaX_eom_cc2_a1_cc2(wf, X, etaX_ai)
!!
!!    etaX EOM cc2 A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
   end subroutine etaX_eom_cc2_a1_cc2
!
!
   module subroutine etaX_eom_a_cc2(wf, etaX, csiX)
!!
!!    Get eom contribution
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM contribution to etaX vector
!!
!!       EOM correction:  η^X,corr_μ += tbar_μ (ξ * tbar) 
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
!
   end subroutine etaX_eom_a_cc2