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
module fftpack_complex_initialization_routines
!
!!
!!    FFTPACK complex initialization routines
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!    Changed error handling so that it now is only done through the return variable ierror
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
   use fftpack_auxiliary_routines
!
!
contains
!
!
   subroutine cfft1i(n, wsave, lensav, ierror)
!!
!!    cfft1i: initialization for cfft1b and cfft1f.
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!    Changed error handling so that it now is only done through the return variable ierror
!!
!!    Purpose:
!!
!!    cfft1i initializes array wsave for use in its companion routines
!!    cfft1b and cfft1f. Routine cfft1i must be called before the first
!!    call to cfft1b or cfft1f, and after whenever the value of integer
!!    n changes.
!!
!!    Parameters:
!!
!!    input,
!!    n, the length of the sequence to be
!!    transformed.  the transform is most efficient when n is a product
!!    of small primes.
!!
!!    input
!!    lensav, the dimension of the wsave array.
!!    lensav must be at least 2*n + int(log(real(n))) + 4.
!!
!!    output,
!!    wsave(lensav), containing the prime factors
!!    of n and  also containing certain trigonometric values which will be used
!!    in routines cfft1b or cfft1f.
!!
!!    output
!!    ierror, error_flag.
!!    0, successful exit;
!!    2, input parameter lensav not big enough.
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)   :: n
      integer, intent(in)   :: lensav
      real(dp), intent(out) :: wsave(lensav)
      integer, intent(out)  :: ierror
!
!     Local variables
!
      integer :: iw1, iw2
!
!     Check validity of input arguments
!
      if (size(wsave) < get_complex_1d_saved_workspace_length(n)) then
         ierror = 2
      else
         ierror = 0
      end if
!
      if (ierror .ne. 0) return
!
!     Perform transform
!
      if (n == 1) return
!
!     Set workspace indices
!
      iw1 = 2*n+1
      iw2 = iw1 + 1
!
      call initialize_factors_and_trig_tables(n, wsave, wsave(iw1), wsave(iw2))
!
   end subroutine cfft1i
!
!
   pure subroutine initialize_factors_and_trig_tables(n, wa, fnf, fac)
!!
!!    Initialize factors and trig tables
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
!!    Purpose:
!!
!!    Sets up factors and trig tables
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)   :: n
      real(dp), intent(out) :: wa(*)
      real(dp), intent(out) :: fnf
      real(dp), intent(out) :: fac(*)
!
!     Local variables
!
      integer :: ido, iip, iw, k1, l1, l2, nf
!
!     Get the factorization of n
!
      call compute_factorization(n, nf, fac)
      fnf = real(nf, dp)
      iw = 1
      l1 = 1
!
!     Set up the trigonometric tables
!
      do k1 = 1, nf
         iip = int(fac(k1))
         l2 = l1 * iip
         ido = n / l2
         call compute_trig_table(ido, iip, wa(iw))
         iw = iw + (iip - 1) * (2*ido)
         l1 = l2
      end do
!
   end subroutine initialize_factors_and_trig_tables
!
!
   pure subroutine compute_factorization(n, nf, fac)
!!
!!    Compute factorization
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
!!    Purpose:
!!
!!    Factors of an integer for dp-kind float precision computations.
!!
!!    Parameters:
!!
!!    n, the number for which factorization and other information is needed.
!!
!!    nf, the number of factors.
!!
!!    fac(*), a list of factors of n.
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)   :: n
      integer, intent(out)  :: nf
      real(dp), intent(out) :: fac(*)
!
!     Local variables
!
      integer :: j, nl, nq, nr, ntry
!
      ntry = 0
      nl = n
      nf = 0
      j = 0
!
      do while (1 < nl)
         j = j + 1
         select case (j)
            case (1)
               ntry = 4
            case (2:3)
               ntry = j
            case (4)
               ntry = 5
            case default
               ntry = ntry + 2
         end select
!
         inner_loop: do
            nq = nl / ntry
            nr = nl - ntry * nq

            if ( nr /= 0 ) then
               exit inner_loop
            end if
!
            nf = nf + 1
            fac(nf) = real(ntry, dp)
            nl = nq
!
         end do inner_loop
      end do
!
   end subroutine compute_factorization
!
!
   pure subroutine compute_trig_table(ido, iip, wa)
!!
!!    Compute trig table
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
!!    Purpose:
!!
!!    Computes trigonometric tables for dp-kind float precision arithmetic.
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)   :: ido
      integer, intent(in)   :: iip
      real(dp), intent(out) :: wa(ido,iip-1,2)
!
!     Local variables
!
      integer             :: i, j !! Counters
      real(dp), parameter :: TWO_PI = TWO * acos(-ONE)
      real(dp)            :: argz, arg1, arg2, arg3, arg4
!
      argz = TWO_PI/iip
      arg1 = TWO_PI/( ido * iip)
!
      do j = 2, iip
         arg2 = real(j - 1, dp) * arg1
         do i = 1, ido
            arg3 = real(i - 1, dp) * arg2
            wa(i,j-1,1) = cos(arg3)
            wa(i,j-1,2) = sin(arg3)
         end do
         if (5 < iip) then
            arg4 = real(j - 1, dp) * argz
            wa(1,j-1,1) = cos(arg4)
            wa(1,j-1,2) = sin(arg4)
         end if
      end do
!
   end subroutine compute_trig_table
!
!
end module fftpack_complex_initialization_routines