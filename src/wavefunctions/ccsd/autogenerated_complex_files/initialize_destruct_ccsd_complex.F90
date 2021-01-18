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
submodule (ccsd_class) initialize_destruct_ccsd_complex
!
!!
!!    Initialize destruct submodule
!!
!!    Gathers routines that initialize and destruct the CCSD type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_amplitudes_ccsd_complex(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Allocates the amplitudes. This routine must be overwritten in
!!    descendants which have more amplitudes.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%initialize_t1_complex()
      call wf%initialize_t2_complex()
      call wf%initialize_u_aibj_complex()
!
   end subroutine initialize_amplitudes_ccsd_complex
!
!
   module subroutine destruct_amplitudes_ccsd_complex(wf)
!!
!!    Destruct amplitudes
!!    Written by Andreas Skeidsvoll, Aug 2019
!!
!!    Deallocates the amplitudes. This routine must be overwritten in
!!    descendants which have more multipliers.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%destruct_t1_complex()
      call wf%destruct_t2_complex()
      call wf%destruct_u_aibj_complex()
!
   end subroutine destruct_amplitudes_ccsd_complex
!
!
   module subroutine initialize_multipliers_ccsd_complex(wf)
!!
!!    Initialize multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Allocates the multipliers. This routine must be overwritten in
!!    descendants which have more multipliers.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%initialize_t1bar_complex()
      call wf%initialize_t2bar_complex()
!
   end subroutine initialize_multipliers_ccsd_complex
!
!
   module subroutine destruct_multipliers_ccsd_complex(wf)
!!
!!    Destruct multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
!!    Deallocates the multipliers. This routine must be overwritten in
!!    descendants which have more multipliers.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%destruct_t1bar_complex()
      call wf%destruct_t2bar_complex()
!
   end subroutine destruct_multipliers_ccsd_complex
!
!
end submodule initialize_destruct_ccsd_complex
