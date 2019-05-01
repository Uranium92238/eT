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
!
  module subroutine prepare_for_eom_fop_ccs(wf)
!!
!!    Prepare for eom fop
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
  end subroutine prepare_for_eom_fop_ccs
!
!
   module subroutine construct_etaX_ccs(wf, X, etaX)
!!
!!    Construct η^X
!!    Written by Josefine H. Andersen, 2019
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
!!    Construct etaX A1 (CCS)
!!    Written by Josefine H. Andersen, 2019
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
!!    Construct right-hand-side vector csiX 
!!    Written by Josefine H. Andersen, Feb 2019
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
   module subroutine add_etaX_eom_correction_ccs(wf, etaX, csiX, X)
!!
!!    Add EOM conrrection to etaX vector
!!    Written by Josefine H. Andersen, Feb 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
   end subroutine add_etaX_eom_correction_ccs
!
!
   module subroutine etaX_eom_a_ccs(wf, etaX, csiX)
!!
!!    Get eom contribution
!!    Written by Josefine H. Andersen, Feb 2019
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
   module subroutine scale_left_excitation_vector_ccs(wf, L, R)
!!
!!    Make left and right excitation vectors biorthogonal by scaling left vector
!!    Written by Josefine H. Andersen, Feb 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: L
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: R
!
   end subroutine scale_left_excitation_vector_ccs
!
!
   module subroutine calculate_transition_strength_ccs(wf, S, etaX, csiX, state, T_l, T_r)
!!
!!    Calculate transition strength for spectra
!!    Written by Josefine H. Andersen, February 2019
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
