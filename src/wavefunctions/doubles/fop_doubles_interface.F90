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
   module subroutine construct_left_transition_density_doubles(wf, state)
!!
!!    Construct left one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      integer, intent(in) :: state
!      
   end subroutine construct_left_transition_density_doubles
!
!
   module subroutine construct_right_transition_density_doubles(wf, state)
!!
!!    Construct right one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      integer, intent(in) :: state
!
   end subroutine construct_right_transition_density_doubles
!
!
   module subroutine right_transition_density_doubles_ov_doubles(wf, tbar_aibj, R_ai)
!!
!!    Right transition density ov contribution (CCSD)
!!    from the singles part of the excitation vector
!!    Written by Alexander Paul, June 2019
!!    
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_doubles_ov_doubles
!
!
   module subroutine right_transition_density_doubles_vo_doubles(wf, tbar_aibj, R_ai)
!!
!!    Right transition density ov contribution (CCSD)
!!    Written by Alexander Paul, June 2019
!!      
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
   end subroutine right_transition_density_doubles_vo_doubles
!
!
   module subroutine construct_eom_etaX_doubles(wf, X, csiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
   end subroutine construct_eom_etaX_doubles
!
!
   module subroutine construct_etaX_doubles(wf, X, etaX)
!!
!!    Construct η^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
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
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
   end subroutine etaX_eom_a_doubles
!
