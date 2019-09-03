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
submodule (ccsd_class) initialize_destruct_ccsd
!
!!
!!    Initialize destruct submodule (CCSD)
!!    Set up by Andreas Skeidsvoll, Sep 2019
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
      call wf%initialize_t1()
      call wf%initialize_t2()
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
      call wf%destruct_t1()
      call wf%destruct_t2()
!
   end subroutine destruct_amplitudes_ccsd
!
!
   module subroutine initialize_t2_ccsd(wf)
!!
!!    Initialize t2 amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (.not. allocated(wf%t2)) call mem%alloc(wf%t2, wf%n_t2)
!
   end subroutine initialize_t2_ccsd
!
!
   module subroutine destruct_t2_ccsd(wf)
!!
!!    Destruct t2 amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2)) call mem%dealloc(wf%t2, wf%n_t2)
!
   end subroutine destruct_t2_ccsd
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
      call wf%initialize_t1bar()
      call wf%initialize_t2bar()
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
      call wf%destruct_t1bar()
      call wf%destruct_t2bar()
!
   end subroutine destruct_multipliers_ccsd
!
!
   module subroutine initialize_t2bar_ccsd(wf)
!!
!!    Initialize T2-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (.not. allocated(wf%t2bar)) call mem%alloc(wf%t2bar, wf%n_t2)
!
   end subroutine initialize_t2bar_ccsd
!
!
   module subroutine destruct_t2bar_ccsd(wf)
!!
!!    Destruct T2-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2bar)) call mem%dealloc(wf%t2bar, wf%n_gs_amplitudes)
!
   end subroutine destruct_t2bar_ccsd
!
!
end submodule initialize_destruct_ccsd