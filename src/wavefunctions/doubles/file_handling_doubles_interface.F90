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
   module subroutine read_doubles_vector_doubles(wf, file_, vector)
!!
!!    Read doubles vector
!!    Written by Alexander C. Paul, Oct 2019
!!
      implicit none 
!
      class(doubles), intent(inout) :: wf
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_t2), intent(out) :: vector
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
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_t2), intent(in) :: vector 
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
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_t1), intent(out) :: X_singles
      real(dp), dimension(wf%n_t2), intent(out) :: X_doubles
!
   end subroutine read_singles_doubles_vector_doubles
