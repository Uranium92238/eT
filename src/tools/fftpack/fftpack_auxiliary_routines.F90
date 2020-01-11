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
module fftpack_auxiliary_routines
!
!!
!!    FFTPACK auxillaries
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
!!    Original work Copyright (c) 2016 Jon Lo Kim Lin
!!
!!    Redistribution and use of the Software in source and binary forms, with
!!    or without modification, is permitted provided that the following
!!    conditions are met:
!! 
!!    Neither the names of NCAR's Computational and Information Systems
!!    Laboratory, the University Corporation for Atmospheric Research, nor
!!    the names of its sponsors or contributors may be used to endorse or
!!    promote products derived from this Software without specific prior
!!    written permission.
!!    Redistributions of source code must retain the above copyright notices,
!!    this list of conditions, and the disclaimer below.
!!    Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions, and the disclaimer below in the documentation
!!    and/or other materials provided with the distribution.
!!    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!!    OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
!!    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!    IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!!    CLAIM, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!!    OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
!!    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!!    OTHER DEALINGS WITH THE SOFTWARE.
!!
!
   use parameters
!
!
contains
!
!
   pure function get_complex_1d_saved_workspace_length(n) result (return_value)
!!
!!    Get complex 1d saved workspace length
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in) :: n
      integer             :: return_value
!
      associate(lensav => return_value)
!
         lensav = 2*n+int(log(real(n, dp))/log(TWO))+4
!
      end associate
!
   end function get_complex_1d_saved_workspace_length
!
!
   pure function get_complex_1d_workspace_length(n) result (return_value)
!!
!!    Get complex 1d workspace length
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in) :: n
      integer             :: return_value
!
      associate(lenwrk => return_value)
!
         lenwrk = 2*n
!
      end associate
!
   end function get_complex_1d_workspace_length
!
!
end module fftpack_auxiliary_routines