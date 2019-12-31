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
submodule (cc2_class) initialize_destruct_cc2
!
!!
!!    Initialize destruct submodule (CC2)
!!    Set up by Andreas Skeidsvoll, Sep 2019
!!
!!    Gathers routines that initialize and destruct the CC2 type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_amplitudes_cc2(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%initialize_t1()
      call wf%initialize_u_aibj()
!
   end subroutine initialize_amplitudes_cc2
!
!
   module subroutine destruct_amplitudes_cc2(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%destruct_t1()
      call wf%destruct_u_aibj()
      call wf%destruct_t2()
!
   end subroutine destruct_amplitudes_cc2
!
!
   module subroutine destruct_multipliers_cc2(wf)
!!
!!    Destruct multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%destruct_t1bar()
      call wf%destruct_t2bar()
!
   end subroutine destruct_multipliers_cc2
!
!
end submodule initialize_destruct_cc2