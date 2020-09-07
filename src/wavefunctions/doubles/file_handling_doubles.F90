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
submodule (doubles_class) file_handling_doubles
!
!!
!!    File handling submodule
!!
!!    Gathers routines that save wavefunction parameters to file,
!!    and reads them from file, plus other routines related to the 
!!    handling of the files that belong to the wavefunction.
!!
!
   implicit none
!
!
contains
!
   module subroutine read_doubles_vector_doubles(wf, file_, vector)
!!
!!    Read doubles vector
!!    Written by Alexander C. Paul, Oct 2019
!!
      implicit none 
!
      class(doubles), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_t2), intent(out) :: vector
!
      call file_%open_('read')
!
!     The first n_t1 dp numbers are the singles vector,
!     so the first byte of the doubles vector is at n_t1*dp + 1
!
      call file_%read_(vector, wf%n_t2, wf%n_t1*dp + 1)
!
      call file_%close_()
!
   end subroutine read_doubles_vector_doubles
!
!
   module subroutine save_doubles_vector_doubles(wf, file_, vector)
!!
!!    Save doubles vector
!!    Written by Alexander C. Paul, Oct 2019 
!!
      implicit none 
!
      class(doubles), intent(inout) :: wf 
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_t2), intent(in) :: vector 
!
      call file_%open_('write')
!
!     The first n_t1 dp numbers are the singles vector,
!     so the first byte of the doubles vector is written to n_t1*dp+1
!
      call file_%write_(vector, wf%n_t2, wf%n_t1*dp + 1)
!
      call file_%close_()
!
   end subroutine save_doubles_vector_doubles
!
!
   module subroutine read_singles_doubles_vector_doubles(wf, file_, X_singles, X_doubles)
!!
!!    Read singles and doubles vector
!!    Written by Alexander C. Paul, May 2020
!!
!!    Reads singles part, checks if doubles exist and reads them as well
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_t1), intent(out) :: X_singles
      real(dp), dimension(wf%n_t2), intent(out) :: X_doubles
!
      call file_%open_('read', 'rewind')
!
      call file_%read_(X_singles, wf%n_t1)
!
!     Check if doubles exist depending on the size of the stream file
!
      if (file_%get_file_size() == dp*(wf%n_t1 + wf%n_t2)) then
!
         call file_%read_(X_doubles, wf%n_t2, wf%n_t1*dp + 1)
!
      else
!
         call zero_array(X_doubles, wf%n_t2)
!
      end if
!
      call file_%close_()
!
   end subroutine read_singles_doubles_vector_doubles
!
!
end submodule file_handling_doubles 
