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
   module subroutine make_doubles_complex_doubles(wf)
!!
!!    Make doubles complex (doubles)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Allocates complex doubles variables, puts the real variables into the real part of the
!!    complex variables and zero into the imaginary part, and deallocates the real variables.
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
   end subroutine make_doubles_complex_doubles
!
!
   module subroutine cleanup_doubles_complex_doubles(wf)
!!
!!    Cleanup doubles complex (doubles)
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
!!    Dellocates complex doubles variables.
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
   end subroutine cleanup_doubles_complex_doubles
