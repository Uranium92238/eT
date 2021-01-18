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
   module subroutine initialize_amplitudes_ccsd(wf)
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
   end subroutine initialize_amplitudes_ccsd
!
!
   module subroutine destruct_amplitudes_ccsd(wf)
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
   end subroutine destruct_amplitudes_ccsd
!
!
   module subroutine initialize_multipliers_ccsd(wf)
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
   end subroutine initialize_multipliers_ccsd
!
!
   module subroutine destruct_multipliers_ccsd(wf)
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
   end subroutine destruct_multipliers_ccsd
