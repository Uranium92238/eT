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
submodule (ccs_class) initialize_destruct_ccs
!
!!
!!    Initialize destruct submodule
!!
!!    Gathers routines that initialize and destruct the CCS type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_amplitudes_ccs(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Allocates the amplitudes. This routine must be overwritten in
!!    descendants which have more amplitudes.
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%initialize_t1()
!
   end subroutine initialize_amplitudes_ccs
!
!
   module subroutine destruct_amplitudes_ccs(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Deallocates the amplitudes. This routine must be overwritten in
!!    descendants which have more amplitudes.
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%destruct_t1()
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
      if (.not. allocated(wf%t1)) call mem%alloc(wf%t1, wf%n_v, wf%n_o)
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
      if (allocated(wf%t1)) call mem%dealloc(wf%t1, wf%n_v, wf%n_o)
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
      call wf%initialize_t1bar()
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
      call wf%destruct_t1bar()
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
      if (.not. allocated(wf%t1bar)) call mem%alloc(wf%t1bar, wf%n_v, wf%n_o)
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
      if (allocated(wf%t1bar)) call mem%dealloc(wf%t1bar, wf%n_v, wf%n_o)
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
      if (.not. allocated(wf%fock_ij)) call mem%alloc(wf%fock_ij, wf%n_o, wf%n_o)
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
      if (allocated(wf%fock_ij)) call mem%dealloc(wf%fock_ij, wf%n_o, wf%n_o)
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
      if (.not. allocated(wf%fock_ia)) call mem%alloc(wf%fock_ia, wf%n_o, wf%n_v)
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
      if (allocated(wf%fock_ia)) call mem%dealloc(wf%fock_ia, wf%n_o, wf%n_v)
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
      if (.not. allocated(wf%fock_ai)) call mem%alloc(wf%fock_ai, wf%n_v, wf%n_o)
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
      if (allocated(wf%fock_ai)) call mem%dealloc(wf%fock_ai, wf%n_v, wf%n_o)
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
      if (.not. allocated(wf%fock_ab)) call mem%alloc(wf%fock_ab, wf%n_v, wf%n_v)
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
      if (allocated(wf%fock_ab)) call mem%dealloc(wf%fock_ab, wf%n_v, wf%n_v)
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
      call mem%alloc(wf%density, wf%n_mo, wf%n_mo)
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
      if (allocated(wf%density)) call mem%dealloc(wf%density, wf%n_mo, wf%n_mo)
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
      call mem%alloc(wf%left_transition_density, wf%n_mo, wf%n_mo)
!
      call mem%alloc(wf%right_transition_density, wf%n_mo, wf%n_mo)
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
      if (allocated(wf%left_transition_density)) &
         call mem%dealloc(wf%left_transition_density, wf%n_mo, wf%n_mo)
!
      if (allocated(wf%right_transition_density)) &
         call mem%dealloc(wf%right_transition_density, wf%n_mo, wf%n_mo)
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
      if (.not. allocated(wf%right_excitation_energies)) &
                  call mem%alloc(wf%right_excitation_energies, wf%n_singlet_states)
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
      if (allocated(wf%right_excitation_energies)) &
               call mem%dealloc(wf%right_excitation_energies, wf%n_singlet_states)
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
      if (.not. allocated(wf%left_excitation_energies)) &
               call mem%alloc(wf%left_excitation_energies, wf%n_singlet_states)
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
      if (allocated(wf%left_excitation_energies)) &
         call mem%dealloc(wf%left_excitation_energies, wf%n_singlet_states)
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
      if (.not. allocated(wf%core_MOs)) call mem%alloc(wf%core_MOs, wf%n_core_MOs)
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
      if (allocated(wf%core_MOs)) call mem%dealloc(wf%core_MOs, wf%n_core_MOs)
!
   end subroutine destruct_core_MOs_ccs
!
!
   module subroutine initialize_fock_ccs(wf)
!!
!!    Initialize Fock
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Initializes all Fock matrix blocks
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
   end subroutine initialize_fock_ccs
!
!
   module subroutine destruct_fock_ccs(wf)
!!
!!    Destruct Fock
!!    Written by Alexander C. Paul, Dec 2019
!!
!!    Destructs all Fock matrix blocks
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%destruct_fock_ij()
      call wf%destruct_fock_ia()
      call wf%destruct_fock_ai()
      call wf%destruct_fock_ab()
!
   end subroutine destruct_fock_ccs
!
end submodule initialize_destruct_ccs
