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
submodule (doubles_class) initialize_destruct_doubles_complex
!
!!
!!    Initialize and destruct submodule
!!
!!    Gathers routines that get and set wavefunction parameters.
!!
!
   implicit none
!
!
contains
!
   module subroutine initialize_t2_doubles_complex(wf)
!!
!!    Initialize t2 amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!
      implicit none
!
      class(doubles) :: wf
!
      if (.not. allocated(wf%t2_complex)) call mem%alloc(wf%t2_complex, wf%n_t2)
!
   end subroutine initialize_t2_doubles_complex
!
!
   module subroutine destruct_t2_doubles_complex(wf)
!!
!!    Destruct t2 amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!
      implicit none
!
      class(doubles) :: wf
!
      if (allocated(wf%t2_complex)) call mem%dealloc(wf%t2_complex, wf%n_t2)
!
   end subroutine destruct_t2_doubles_complex
!
!
   module subroutine initialize_t2bar_doubles_complex(wf)
!!
!!    Initialize t2bar amplitudes
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(doubles) :: wf
!
      if (.not. allocated(wf%t2bar_complex)) call mem%alloc(wf%t2bar_complex, wf%n_t2)
!
   end subroutine initialize_t2bar_doubles_complex
!
!
   module subroutine destruct_t2bar_doubles_complex(wf)
!!
!!    Destruct t2bar amplitudes
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(doubles) :: wf
!
      if (allocated(wf%t2bar_complex)) call mem%dealloc(wf%t2bar_complex, wf%n_t2)
!
   end subroutine destruct_t2bar_doubles_complex
!
!
   module subroutine initialize_u_aibj_doubles_complex(wf)
!!
!!    Initialize u_aibj
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      if (.not. allocated(wf%u_aibj_complex)) call mem%alloc(wf%u_aibj_complex, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine initialize_u_aibj_doubles_complex
!
!
   module subroutine destruct_u_aibj_doubles_complex(wf)
!!
!!    Destruct u_aibj
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      if (allocated(wf%u_aibj_complex)) call mem%dealloc(wf%u_aibj_complex, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine destruct_u_aibj_doubles_complex
!
end submodule initialize_destruct_doubles_complex
