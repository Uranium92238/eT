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
submodule (ccs_class) initialize_destruct_ccs_complex
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
   module subroutine initialize_amplitudes_ccs_complex(wf)
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
      call wf%initialize_t1_complex()
!
   end subroutine initialize_amplitudes_ccs_complex
!
!
   module subroutine destruct_amplitudes_ccs_complex(wf)
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
      call wf%destruct_t1_complex()
!
   end subroutine destruct_amplitudes_ccs_complex
!
!
   module subroutine initialize_t1_ccs_complex(wf)
!!
!!    Initialize T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%t1_complex)) call mem%alloc(wf%t1_complex, wf%n_v, wf%n_o)
!
      call zero_array_complex(wf%t1_complex, wf%n_t1)
!
   end subroutine initialize_t1_ccs_complex
!
!
   module subroutine destruct_t1_ccs_complex(wf)
!!
!!    Destruct T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1_complex)) call mem%dealloc(wf%t1_complex, wf%n_v, wf%n_o)
!
   end subroutine destruct_t1_ccs_complex
!
!
   module subroutine initialize_multipliers_ccs_complex(wf)
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
      call wf%initialize_t1bar_complex()
!
   end subroutine initialize_multipliers_ccs_complex
!
!
   module subroutine destruct_multipliers_ccs_complex(wf)
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
      call wf%destruct_t1bar_complex()
!
   end subroutine destruct_multipliers_ccs_complex
!
!
   module subroutine initialize_t1bar_ccs_complex(wf)
!!
!!    Initialize T1-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%t1bar_complex)) call mem%alloc(wf%t1bar_complex, wf%n_v, wf%n_o)
!
   end subroutine initialize_t1bar_ccs_complex
!
!
   module subroutine destruct_t1bar_ccs_complex(wf)
!!
!!    Destruct T1-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1bar_complex)) call mem%dealloc(wf%t1bar_complex, wf%n_v, wf%n_o)
!
   end subroutine destruct_t1bar_ccs_complex
!
!
   module subroutine initialize_fock_ij_ccs_complex(wf)
!!
!!    Initialize Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ij_complex)) call mem%alloc(wf%fock_ij_complex, wf%n_o, wf%n_o)
!
   end subroutine initialize_fock_ij_ccs_complex
!
!
   module subroutine destruct_fock_ij_ccs_complex(wf)
!!
!!    Destruct Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ij_complex)) call mem%dealloc(wf%fock_ij_complex, wf%n_o, wf%n_o)
!
   end subroutine destruct_fock_ij_ccs_complex
!
!
   module subroutine initialize_fock_ia_ccs_complex(wf)
!!
!!    Initialize Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ia_complex)) call mem%alloc(wf%fock_ia_complex, wf%n_o, wf%n_v)
!
   end subroutine initialize_fock_ia_ccs_complex
!
!
   module subroutine destruct_fock_ia_ccs_complex(wf)
!!
!!    Destruct Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ia_complex)) call mem%dealloc(wf%fock_ia_complex, wf%n_o, wf%n_v)
!
   end subroutine destruct_fock_ia_ccs_complex
!
!
   module subroutine initialize_fock_ai_ccs_complex(wf)
!!
!!    Initialize Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ai_complex)) call mem%alloc(wf%fock_ai_complex, wf%n_v, wf%n_o)
!
   end subroutine initialize_fock_ai_ccs_complex
!
!
   module subroutine destruct_fock_ai_ccs_complex(wf)
!!
!!    Destruct Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ai_complex)) call mem%dealloc(wf%fock_ai_complex, wf%n_v, wf%n_o)
!
   end subroutine destruct_fock_ai_ccs_complex
!
!
   module subroutine initialize_fock_ab_ccs_complex(wf)
!!
!!    Initialize Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ab_complex)) call mem%alloc(wf%fock_ab_complex, wf%n_v, wf%n_v)
!
   end subroutine initialize_fock_ab_ccs_complex
!
!
   module subroutine destruct_fock_ab_ccs_complex(wf)
!!
!!    Destruct Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ab_complex)) call mem%dealloc(wf%fock_ab_complex, wf%n_v, wf%n_v)
!
   end subroutine destruct_fock_ab_ccs_complex
!
!
   module subroutine initialize_gs_density_ccs_complex(wf)
!!
!!    Initialize density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      call mem%alloc(wf%density_complex, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_gs_density_ccs_complex
!
!
   module subroutine destruct_gs_density_ccs_complex(wf)
!!
!!    Destruct density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%density_complex)) call mem%dealloc(wf%density_complex, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_gs_density_ccs_complex
!
!
   module subroutine initialize_core_MOs_ccs_complex(wf)
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
   end subroutine initialize_core_MOs_ccs_complex
!
!
   module subroutine destruct_core_MOs_ccs_complex(wf)
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
   end subroutine destruct_core_MOs_ccs_complex
!
!
   module subroutine initialize_fock_ccs_complex(wf)
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
      call wf%initialize_fock_ij_complex()
      call wf%initialize_fock_ia_complex()
      call wf%initialize_fock_ai_complex()
      call wf%initialize_fock_ab_complex()
!
   end subroutine initialize_fock_ccs_complex
!
!
   module subroutine destruct_fock_ccs_complex(wf)
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
      call wf%destruct_fock_ij_complex()
      call wf%destruct_fock_ia_complex()
      call wf%destruct_fock_ai_complex()
      call wf%destruct_fock_ab_complex()
!
   end subroutine destruct_fock_ccs_complex
!
end submodule initialize_destruct_ccs_complex
