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
   module subroutine initialize_amplitudes_cc2(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
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
   end subroutine destruct_amplitudes_cc2
!
!
   module subroutine initialize_t2_cc2(wf)
!!
!!    Initialize t2 amplitudes
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(cc2) :: wf
!
   end subroutine initialize_t2_cc2
!
!
   module subroutine destruct_t2_cc2(wf)
!!
!!    Destruct t2 amplitudes
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(cc2) :: wf
!
   end subroutine destruct_t2_cc2
!
!
   module subroutine initialize_t2bar_cc2(wf)
!!
!!    Initialize t2bar amplitudes
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(cc2) :: wf
!
   end subroutine initialize_t2bar_cc2
!
!
   module subroutine destruct_t2bar_cc2(wf)
!!
!!    Destruct t2bar amplitudes
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(cc2) :: wf
!
   end subroutine destruct_t2bar_cc2