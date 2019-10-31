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
   module subroutine initialize_amplitudes_ccs(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_amplitudes_ccs
!
!
   module subroutine destruct_amplitudes_ccs(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_amplitudes_ccs
!
!
   module subroutine initialize_t1_ccs(wf)
!!
!!    Initialize T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_t1_ccs
!
!
   module subroutine destruct_t1_ccs(wf)
!!
!!    Destruct T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_t1_ccs
!
!
   module subroutine initialize_multipliers_ccs(wf)
!!
!!    Initialize multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Allocates the multipliers. This routine must be overwritten in
!!    descendants which have more multipliers.
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_multipliers_ccs
!
!
   module subroutine destruct_multipliers_ccs(wf)
!!
!!    Destruct multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Deallocates the multipliers. This routine must be overwritten in 
!!    descendants which have more multipliers. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
   end subroutine destruct_multipliers_ccs
!
!
   module subroutine initialize_t1bar_ccs(wf)
!!
!!    Initialize T1-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_t1bar_ccs
!
!
   module subroutine destruct_t1bar_ccs(wf)
!!
!!    Destruct T1-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_t1bar_ccs
!
!
   module subroutine initialize_fock_ij_ccs(wf)
!!
!!    Initialize Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_fock_ij_ccs
!
!
   module subroutine destruct_fock_ij_ccs(wf)
!!
!!    Destruct Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_fock_ij_ccs
!
!
   module subroutine initialize_fock_ia_ccs(wf)
!!
!!    Initialize Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_fock_ia_ccs
!
!
   module subroutine destruct_fock_ia_ccs(wf)
!!
!!    Destruct Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_fock_ia_ccs
!
!
   module subroutine initialize_fock_ai_ccs(wf)
!!
!!    Initialize Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_fock_ai_ccs
!
!
   module subroutine destruct_fock_ai_ccs(wf)
!!
!!    Destruct Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_fock_ai_ccs
!
!
   module subroutine initialize_fock_ab_ccs(wf)
!!
!!    Initialize Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_fock_ab_ccs
!
!
   module subroutine destruct_fock_ab_ccs(wf)
!!
!!    Destruct Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_fock_ab_ccs
!
!
   module subroutine initialize_gs_density_ccs(wf)
!!
!!    Initialize density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_gs_density_ccs
!
!
   module subroutine destruct_gs_density_ccs(wf)
!!
!!    Destruct density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_gs_density_ccs
!
!
   module subroutine initialize_transition_densities_ccs(wf)
!!
!!    Initialize left and right transition densities
!!    Written by Alexander C. Paul, June 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_transition_densities_ccs
!
!
   module subroutine destruct_transition_densities_ccs(wf)
!!
!!    Destruct left and right transition densities
!!    Written by Alexander C. Paul, June 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_transition_densities_ccs
!
!
   module subroutine initialize_right_excitation_energies_ccs(wf)
!!
!!    Initialize right excitation energies 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_right_excitation_energies_ccs
!
!
   module subroutine destruct_right_excitation_energies_ccs(wf)
!!
!!    Destruct right excitation energies 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_right_excitation_energies_ccs
!
!
   module subroutine initialize_left_excitation_energies_ccs(wf)
!!
!!    Initialize left excitation energies 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_left_excitation_energies_ccs
!
!
   module subroutine destruct_left_excitation_energies_ccs(wf)
!!
!!    Destruct left excitation energies 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_left_excitation_energies_ccs
!
!
   module subroutine initialize_core_MOs_ccs(wf)
!!
!!    Initialize core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_core_MOs_ccs
!
!
   module subroutine destruct_core_MOs_ccs(wf)
!!
!!    Destruct core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine destruct_core_MOs_ccs
!
!
   module subroutine initialize_fock_ccs(wf)
!!
!!    Initialize Fock
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine initialize_fock_ccs