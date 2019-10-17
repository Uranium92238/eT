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
   module subroutine save_excited_state_abstract_doubles(wf, X, n, side)
!!
!!    Save excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(abstract_doubles), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
   end subroutine save_excited_state_abstract_doubles
!
!
   module subroutine save_doubles_vector_abstract_doubles(wf, X, n, file_)
!!
!!    Save doubles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!
      implicit none
!
      class(abstract_doubles), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(sequential_file) :: file_
!
   end subroutine save_doubles_vector_abstract_doubles
!
!
   module subroutine read_excited_state_abstract_doubles(wf, X, n, side)
!!
!!    Read excited state 
!!    Written by Sarai D. Fokestad, Mar 2019 
!!
      implicit none
!
      class(abstract_doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: X
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
   end subroutine read_excited_state_abstract_doubles
!
!
   module subroutine read_doubles_vector_abstract_doubles(wf, X, n, file_)
!!
!!    Read doubles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!
      implicit none 
!
      class(abstract_doubles), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(out) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(sequential_file) :: file_
!
   end subroutine read_doubles_vector_abstract_doubles
!
