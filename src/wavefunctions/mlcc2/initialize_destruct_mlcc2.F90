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
submodule (mlcc2_class) initialize_desctruct_mlcc2
!
!!
!!    Initialize destruct submodule (MLCC2)
!!
!!    Gathers routines that initialize and destruct the MLCC2 allocatable variables.
!!
!
   implicit none
!
!
contains
!
   module subroutine initialize_u_mlcc2(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (.not. allocated(wf%u)) call mem%alloc(wf%u, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
   end subroutine initialize_u_mlcc2
!
!
   module subroutine destruct_u_mlcc2(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (allocated(wf%u)) call mem%dealloc(wf%u, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
   end subroutine destruct_u_mlcc2
!
!
   module subroutine initialize_amplitudes_mlcc2(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      call wf%initialize_t1()
      call wf%initialize_u()
      call wf%initialize_x2()
!
   end subroutine initialize_amplitudes_mlcc2
!
!
   module subroutine destruct_amplitudes_mlcc2(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      call wf%destruct_t1()
      call wf%destruct_u()
      call wf%destruct_x2()
!
   end subroutine destruct_amplitudes_mlcc2
!
!
   module subroutine initialize_x2_mlcc2(wf)
!!
!!    Initialize x2 amplitudes
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (.not. allocated(wf%x2)) call mem%alloc(wf%x2, wf%n_x2)
!
   end subroutine initialize_x2_mlcc2
!
!
   module subroutine initialize_t2bar_mlcc2(wf)
!!
!!    Initialize t2bar amplitudes
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (.not. allocated(wf%t2bar)) call mem%alloc(wf%t2bar, wf%n_x2)
!
   end subroutine initialize_t2bar_mlcc2
!
!
   module subroutine destruct_x2_mlcc2(wf)
!!
!!    Destruct x2 amplitudes
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (allocated(wf%x2)) call mem%dealloc(wf%x2, wf%n_x2)
!
   end subroutine destruct_x2_mlcc2
!
!
   module subroutine destruct_t2bar_mlcc2(wf)
!!
!!    Destruct t2bar amplitudes
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (allocated(wf%t2bar)) call mem%dealloc(wf%t2bar, wf%n_x2)
!
   end subroutine destruct_t2bar_mlcc2
!
!
   module subroutine initialize_nto_states_mlcc2(wf)
!!
!!    Initialize nto states
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (.not. allocated(wf%nto_states)) call mem%alloc(wf%nto_states, wf%n_nto_states)
!
   end subroutine initialize_nto_states_mlcc2
!
!
   module subroutine initialize_cnto_states_mlcc2(wf)
!!
!!    Initialize cnto states
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (.not. allocated(wf%cnto_states)) call mem%alloc(wf%cnto_states, wf%n_cnto_states)
!
   end subroutine initialize_cnto_states_mlcc2
!
!
   module subroutine destruct_nto_states_mlcc2(wf)
!!
!!    Destruct nto states
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (allocated(wf%nto_states)) call mem%dealloc(wf%nto_states, wf%n_nto_states)
!
   end subroutine destruct_nto_states_mlcc2
!
!
   module subroutine destruct_cnto_states_mlcc2(wf)
!!
!!    Destruct cnto states
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (allocated(wf%cnto_states)) call mem%dealloc(wf%cnto_states, wf%n_cnto_states)
!
   end subroutine destruct_cnto_states_mlcc2
!
!
end submodule initialize_desctruct_mlcc2