!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module fftpack_complex_forward_1d
!
!!
!!    FFTPACK complex forward 1d module
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
   subroutine cfft1f(n, inc, c, lenc, wsave, lensav, work, lenwrk, ierror)

      use, intrinsic :: ISO_C_binding, only: c_f_pointer, c_loc
!!
!!    cfft1f: complex forward fast Fourier transform, 1d.
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!    Changed error handling so that it now is only done through the return variable ierror
!!
!!    purpose:
!!
!!    cfft1f computes the one-dimensional Fourier transform of a single
!!    periodic sequence within a complex array. this transform is referred
!!    to as the forward transform or Fourier analysis, transforming the
!!    sequence from physical to spectral space.
!!
!!    this transform is normalized since a call to cfft1f followed
!!    by a call to cfft1b (or vice-versa) reproduces the original
!!    array within roundoff error.
!!
!!    input
!!
!!    integer n, the length of the sequence to be
!!    transformed. the transform is most efficient when
!!    n is a product of small primes.
!!
!!    integer inc, the increment between the locations, in
!!    array c, of two consecutive elements within the sequence to be transformed.
!!
!!    real wsave(lensav).  wsave's contents must be
!!    initialized with a call to cfft1i before the first call to routine cfft1f
!!    or cfft1b for a given transform length n.  wsave's contents may be re-used
!!    for subsequent calls to cfft1f and cfft1b with the same n.
!!
!!    integer lensav, the dimension of the wsave array.
!!    lensav must be at least 2*n + int(log(real(n))) + 4.
!!
!!    integer lenwrk, the dimension of the work array.
!!    lenwrk must be at least 2*n.
!!
!!    input/output
!!    complex c(lenc) containing the sequence to be transformed.
!!
!!    real work(lenwrk), workspace array
!!    integer lenc, the dimension of the c array.
!!    lenc must be at least inc*(n-1) + 1.
!!
!!    output
!!    integer ierror, error_flag.
!!    0, successful exit;
!!    1, input parameter lenc not big enough;
!!    2, input parameter lensav not big enough;
!!    3, input parameter lenwrk not big enough;
!!    20, input error returned by lower level routine.
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)                :: n
      integer, intent(in)                :: inc
      integer, intent(in)                :: lenc
      complex(dp), intent(inout), target :: c(lenc)
      integer, intent(in)                :: lensav
      real(dp),   intent(in)             :: wsave(lensav)
      integer, intent(in)                :: lenwrk
      real(dp),   intent(out)            :: work(lenwrk)
      integer, intent(out)               :: ierror
!
!     Local variables
!
      real(dp), pointer :: real_arg(:) => null()
      integer           :: iw1, iw2
!
!     Check validity of input arguments
!
      call check_calling_arguments(n, inc, c, wsave, work, ierror)
      if (ierror .ne. 0) return
!
!     Perform transform
!
      if (n == 1) return
!
!     Perform a C-language style cast without copying
!
      call c_f_pointer(c_loc(c), real_arg, shape=[2*size(c)])
!
!     Set workspace index pointer
!
      iw1 = (2 * n) + 1
      iw2 = iw1 + 1
!
      call complex_pass_forward(n, inc, real_arg, work, wsave, wsave(iw1), wsave(iw2:))
!
!     Terminate association
!
      nullify(real_arg)
!
   end subroutine cfft1f
!
!
   subroutine check_calling_arguments(n, inc, c, wsave, work, ierror)
!!
!!    Check calling arguments
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!    Changed error handling so that it now is only done through the return variable ierror
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)     :: n
      integer, intent(in)     :: inc
      complex(dp), intent(in) :: c(:)
      real(dp), intent(in)    :: wsave(:)
      real(dp), intent(in)    :: work(:)
      integer, intent(out)    :: ierror
!
!     Check validity of input arguments
!
      if (size(c) < inc*(n-1) + 1) then
         ierror = 1
      else if (size(wsave) < get_complex_1d_saved_workspace_length(n)) then
         ierror = 2
      else if (size(work) < get_complex_1d_workspace_length(n)) then
         ierror = 3
      else
         ierror = 0
      end if
!
   end subroutine check_calling_arguments
!
!
   subroutine complex_pass_forward(n, inc, c, ch, wa, fnf, fac)
!!
!!    Complex pass forward
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)      :: n
      integer, intent(in)      :: inc
      real(dp), intent(inout) :: c(:)
      real(dp), intent(out)    :: ch(:)
      real(dp), intent(in)     :: wa(:)
      real(dp), intent(in)     :: fnf
      real(dp), intent(in)     :: fac(:)
!
!     Local variables
!
      integer  :: ido, inc2, iip, iw
      integer  :: k1, l1, l2, lid
      integer  :: na, nbr, nf
!
      inc2 = 2*inc
      nf = int(fnf)
      na = 0
      l1 = 1
      iw = 1
!
      do k1=1,nf
!
         iip = int(fac(k1))
         l2 = iip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(iip-2,4)
!
         select case (nbr)
            case (1)
               call complex_pass_2_forward(ido,l1,na,c,inc2,ch,2,wa(iw:))
            case (2)
               call complex_pass_2_forward(ido,l1,na,ch,2,c,inc2,wa(iw:))
            case (3)
               call complex_pass_3_forward(ido,l1,na,c,inc2,ch,2,wa(iw:))
            case (4)
               call complex_pass_3_forward(ido,l1,na,ch,2,c,inc2,wa(iw:))
            case (5)
               call complex_pass_4_forward(ido,l1,na,c,inc2,ch,2,wa(iw:))
            case (6)
               call complex_pass_4_forward(ido,l1,na,ch,2,c,inc2,wa(iw:))
            case (7)
               call complex_pass_5_forward(ido,l1,na,c,inc2,ch,2,wa(iw:))
            case (8)
               call complex_pass_5_forward(ido,l1,na,ch,2,c,inc2,wa(iw:))
            case (9)
               call complex_pass_n_forward(ido,iip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw:))
            case (10)
               call complex_pass_n_forward(ido,iip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw:))
         end select
!
         l1 = l2
         iw = iw+(iip-1)*(2*ido)
!
         if (iip <= 5) na = 1-na
!
      end do
!
   end subroutine complex_pass_forward
!
!
   subroutine complex_pass_2_forward(ido, l1, na, cc, in1, ch, in2, wa)
!!
!!    Complex pass 2 forward
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)     :: ido
      integer, intent(in)     :: l1
      integer, intent(in)     :: na
      integer, intent(in)     :: in1
      real(dp), intent(inout) :: cc(in1,l1,ido,2)
      integer, intent(in)     :: in2
      real(dp), intent(inout) :: ch(in2,l1,2,ido)
      real(dp), intent(in)    :: wa(ido,1,2)
!
!     Local variables
!
      integer               :: i !! counter
      real(dp)              :: sn
      real(dp), allocatable :: chold1(:), chold2(:)
      real(dp), allocatable :: ti2(:), tr2(:)
!
      if (1 >= ido) then
         sn = ONE/(2 * l1)
         if (na /= 1) then
!
!           Allocate memory
!
            allocate(chold1(l1))
            allocate(chold2(l1))
!
            chold1 = sn*(cc(1,:,1,1)+cc(1,:,1,2))
            cc(1,:,1,2) = sn*(cc(1,:,1,1)-cc(1,:,1,2))
            cc(1,:,1,1) = chold1
            chold2 = sn*(cc(2,:,1,1)+cc(2,:,1,2))
            cc(2,:,1,2) = sn*(cc(2,:,1,1)-cc(2,:,1,2))
            cc(2,:,1,1) = chold2
!
!           Release memory
!
            deallocate(chold1)
            deallocate(chold2)
         else
            ch(1,:,1,1) = sn*(cc(1,:,1,1)+cc(1,:,1,2))
            ch(1,:,2,1) = sn*(cc(1,:,1,1)-cc(1,:,1,2))
            ch(2,:,1,1) = sn*(cc(2,:,1,1)+cc(2,:,1,2))
            ch(2,:,2,1) = sn*(cc(2,:,1,1)-cc(2,:,1,2))
         end if
      else
         ch(1,:,1,1) = cc(1,:,1,1)+cc(1,:,1,2)
         ch(1,:,2,1) = cc(1,:,1,1)-cc(1,:,1,2)
         ch(2,:,1,1) = cc(2,:,1,1)+cc(2,:,1,2)
         ch(2,:,2,1) = cc(2,:,1,1)-cc(2,:,1,2)
!
!        Allocate memory
!
         allocate(tr2(l1))
         allocate(ti2(l1))

         do i = 2, ido
            ch(1,:,1,i) = cc(1,:,i,1)+cc(1,:,i,2)
            tr2 = cc(1,:,i,1)-cc(1,:,i,2)
            ch(2,:,1,i) = cc(2,:,i,1)+cc(2,:,i,2)
            ti2 = cc(2,:,i,1)-cc(2,:,i,2)
            ch(2,:,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
            ch(1,:,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
         end do
!
!        Release memory
!
         deallocate(tr2)
         deallocate(ti2)
      end if
!
   end subroutine complex_pass_2_forward
!
!
   subroutine complex_pass_3_forward(ido, l1, na, cc, in1, ch, in2, wa)
!!
!!    Complex pass 3 forward
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)     :: ido
      integer, intent(in)     :: l1
      integer, intent(in)     :: na
      integer, intent(in)     :: in1
      real(dp), intent(inout) :: cc(in1,l1,ido,3)
      integer, intent(in)     :: in2
      real(dp), intent(inout) :: ch(in2,l1,3,ido)
      real(dp), intent(in)    :: wa(ido,2,2)
!
!     Local variables
!
      integer               :: i !! Counter
      real(dp), allocatable :: ci2(:), ci3(:)
      real(dp), allocatable :: cr2(:), cr3(:)
      real(dp), allocatable :: di2(:), di3(:)
      real(dp), allocatable :: dr2(:), dr3(:)
      real(dp), allocatable :: ti2(:), tr2(:)
      real(dp), parameter   :: TAUI = -sqrt(THREE)/2 !-0.866025403784439_dp
      real(dp), parameter   :: TAUR = -HALF
      real(dp)              :: sn
!
!     Allocate memory
!
      allocate( ci2(l1), ci3(l1) )
      allocate( cr2(l1), cr3(l1) )
      allocate( ti2(l1), tr2(l1) )
!
      if (1 >= ido) then
         sn = ONE/(3 * l1)
         if (na /= 1) then
            tr2 = cc(1,:,1,2)+cc(1,:,1,3)
            cr2 = cc(1,:,1,1)+TAUR*tr2
            cc(1,:,1,1) = sn*(cc(1,:,1,1)+tr2)
            ti2 = cc(2,:,1,2)+cc(2,:,1,3)
            ci2 = cc(2,:,1,1)+TAUR*ti2
            cc(2,:,1,1) = sn*(cc(2,:,1,1)+ti2)
            cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
            ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
            cc(1,:,1,2) = sn*(cr2-ci3)
            cc(1,:,1,3) = sn*(cr2+ci3)
            cc(2,:,1,2) = sn*(ci2+cr3)
            cc(2,:,1,3) = sn*(ci2-cr3)
         else
            tr2 = cc(1,:,1,2)+cc(1,:,1,3)
            cr2 = cc(1,:,1,1)+TAUR*tr2
            ch(1,:,1,1) = sn*(cc(1,:,1,1)+tr2)
            ti2 = cc(2,:,1,2)+cc(2,:,1,3)
            ci2 = cc(2,:,1,1)+TAUR*ti2
            ch(2,:,1,1) = sn*(cc(2,:,1,1)+ti2)
            cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
            ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
            ch(1,:,2,1) = sn*(cr2-ci3)
            ch(1,:,3,1) = sn*(cr2+ci3)
            ch(2,:,2,1) = sn*(ci2+cr3)
            ch(2,:,3,1) = sn*(ci2-cr3)
         end if
      else
         tr2 = cc(1,:,1,2)+cc(1,:,1,3)
         cr2 = cc(1,:,1,1)+TAUR*tr2
         ch(1,:,1,1) = cc(1,:,1,1)+tr2
         ti2 = cc(2,:,1,2)+cc(2,:,1,3)
         ci2 = cc(2,:,1,1)+TAUR*ti2
         ch(2,:,1,1) = cc(2,:,1,1)+ti2
         cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
         ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
         ch(1,:,2,1) = cr2-ci3
         ch(1,:,3,1) = cr2+ci3
         ch(2,:,2,1) = ci2+cr3
         ch(2,:,3,1) = ci2-cr3
!
!        Allocate memory
!
         allocate( di2(l1), di3(l1) )
         allocate( dr2(l1), dr3(l1) )

         do i=2,ido
            tr2 = cc(1,:,i,2)+cc(1,:,i,3)
            cr2 = cc(1,:,i,1)+TAUR*tr2
            ch(1,:,1,i) = cc(1,:,i,1)+tr2
            ti2 = cc(2,:,i,2)+cc(2,:,i,3)
            ci2 = cc(2,:,i,1)+TAUR*ti2
            ch(2,:,1,i) = cc(2,:,i,1)+ti2
            cr3 = TAUI*(cc(1,:,i,2)-cc(1,:,i,3))
            ci3 = TAUI*(cc(2,:,i,2)-cc(2,:,i,3))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(2,:,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,:,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,:,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,:,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
         end do
!
!        Release memory
!
         deallocate( di2, di3 )
         deallocate( dr2, dr3 )
      end if
!
!     Release memory
!
      deallocate( ci2, ci3 )
      deallocate( cr2, cr3 )
      deallocate( ti2, tr2 )
!
   end subroutine complex_pass_3_forward
!
!
   subroutine complex_pass_4_forward(ido, l1, na, cc, in1, ch, in2, wa)
!!
!!    Complex pass 4 forward
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer, intent(in)     :: ido
      integer, intent(in)     :: l1
      integer, intent(in)     :: na
      integer, intent(in)     :: in1
      real(dp), intent(inout) :: cc(in1,l1,ido,4)
      integer, intent(in)     :: in2
      real(dp), intent(inout) :: ch(in2,l1,4,ido)
      real(dp), intent(in)    :: wa(ido,3,2)
!
!     Local variables
!
      integer  :: i
      real(dp)   :: sn
!
      if (1 >= ido) then
         sn = ONE/(4 * l1)
         if (na /= 1) then
            associate(&
               ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
               ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
               tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
               ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
               tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
               tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
               ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
               tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
               )
               cc(1,:,1,1) = sn*(tr2+tr3)
               cc(1,:,1,3) = sn*(tr2-tr3)
               cc(2,:,1,1) = sn*(ti2+ti3)
               cc(2,:,1,3) = sn*(ti2-ti3)
               cc(1,:,1,2) = sn*(tr1+tr4)
               cc(1,:,1,4) = sn*(tr1-tr4)
               cc(2,:,1,2) = sn*(ti1+ti4)
               cc(2,:,1,4) = sn*(ti1-ti4)
            end associate
         else
            associate( &
               ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
               ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
               tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
               ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
               tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
               tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
               ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
               tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
               )
               ch(1,:,1,1) = sn*(tr2+tr3)
               ch(1,:,3,1) = sn*(tr2-tr3)
               ch(2,:,1,1) = sn*(ti2+ti3)
               ch(2,:,3,1) = sn*(ti2-ti3)
               ch(1,:,2,1) = sn*(tr1+tr4)
               ch(1,:,4,1) = sn*(tr1-tr4)
               ch(2,:,2,1) = sn*(ti1+ti4)
               ch(2,:,4,1) = sn*(ti1-ti4)
            end associate
         end if
      else
         associate( &
            ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
            ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
            tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
            ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
            tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
            tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
            ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
            tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
            )
            ch(1,:,1,1) = tr2+tr3
            ch(1,:,3,1) = tr2-tr3
            ch(2,:,1,1) = ti2+ti3
            ch(2,:,3,1) = ti2-ti3
            ch(1,:,2,1) = tr1+tr4
            ch(1,:,4,1) = tr1-tr4
            ch(2,:,2,1) = ti1+ti4
            ch(2,:,4,1) = ti1-ti4
         end associate
         do i=2,ido
            associate( &
               ti1 => cc(2,:,i,1)-cc(2,:,i,3), &
               ti2 => cc(2,:,i,1)+cc(2,:,i,3), &
               ti3 => cc(2,:,i,2)+cc(2,:,i,4), &
               tr4 => cc(2,:,i,2)-cc(2,:,i,4), &
               tr1 => cc(1,:,i,1)-cc(1,:,i,3), &
               tr2 => cc(1,:,i,1)+cc(1,:,i,3), &
               ti4 => cc(1,:,i,4)-cc(1,:,i,2), &
               tr3 => cc(1,:,i,2)+cc(1,:,i,4) &
               )
               ch(1,:,1,i) = tr2+tr3
               associate( cr3 => tr2-tr3 )
                  ch(2,:,1,i) = ti2+ti3
                  associate( &
                     ci3 => ti2-ti3, &
                     cr2 => tr1+tr4, &
                     cr4 => tr1-tr4, &
                     ci2 => ti1+ti4, &
                     ci4 => ti1-ti4 &
                     )
                     ch(1,:,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
                     ch(2,:,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
                     ch(1,:,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
                     ch(2,:,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
                     ch(1,:,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
                     ch(2,:,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
                  end associate
               end associate
            end associate
         end do
      end if
!
   end subroutine complex_pass_4_forward
!
!
   subroutine complex_pass_5_forward(ido, l1, na, cc, in1, ch, in2, wa)
!!
!!    Complex pass 5 forward
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
!     Dummy arguments
!
      integer  :: ido
      integer  :: in1
      integer  :: in2
      integer  :: l1
      integer  :: na
      real(dp) :: cc(in1,l1,ido,5)
      real(dp) :: ch(in2,l1,5,ido)
      real(dp) :: wa(ido,4,2)
!
!     Local variables
!
      integer             :: i, k
      real(dp)            :: sn
      real(dp)            :: chold1, chold2
      real(dp)            :: ci2, ci3, ci4, ci5
      real(dp)            :: cr2, cr3, cr4, cr5
      real(dp)            :: di2, di3, di4, di5
      real(dp)            :: dr2, dr3, dr4, dr5
      real(dp)            :: ti2, ti3, ti4, ti5
      real(dp)            :: tr2, tr3, tr4, tr5
      real(dp), parameter :: SQRT5 = sqrt(FIVE)
      real(dp), parameter :: SQRT5_PLUS_5 = SQRT5 + FIVE
      real(dp), parameter :: TI11 = -sqrt(SQRT5_PLUS_5/2)/2 !-0.9510565162951536_dp
      real(dp), parameter :: TI12 = -sqrt(FIVE/(TWO*SQRT5_PLUS_5)) !-0.5877852522924731_dp
      real(dp), parameter :: TR11 = (SQRT5 - ONE)/4 ! 0.3090169943749474_dp
      real(dp), parameter :: TR12 = -(ONE + SQRT5)/4 !-0.8090169943749474_dp
!
      if (1 >= ido) then
         sn = ONE/(5 * l1)
         if (na /= 1) then
            do k = 1,l1
               ti5 = cc(2,k,1,2)-cc(2,k,1,5)
               ti2 = cc(2,k,1,2)+cc(2,k,1,5)
               ti4 = cc(2,k,1,3)-cc(2,k,1,4)
               ti3 = cc(2,k,1,3)+cc(2,k,1,4)
               tr5 = cc(1,k,1,2)-cc(1,k,1,5)
               tr2 = cc(1,k,1,2)+cc(1,k,1,5)
               tr4 = cc(1,k,1,3)-cc(1,k,1,4)
               tr3 = cc(1,k,1,3)+cc(1,k,1,4)
               chold1 = sn*(cc(1,k,1,1)+tr2+tr3)
               chold2 = sn*(cc(2,k,1,1)+ti2+ti3)
               cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
               ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
               cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
               ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
               cc(1,k,1,1) = chold1
               cc(2,k,1,1) = chold2
               cr5 = TI11*tr5+TI12*tr4
               ci5 = TI11*ti5+TI12*ti4
               cr4 = TI12*tr5-TI11*tr4
               ci4 = TI12*ti5-TI11*ti4
               cc(1,k,1,2) = sn*(cr2-ci5)
               cc(1,k,1,5) = sn*(cr2+ci5)
               cc(2,k,1,2) = sn*(ci2+cr5)
               cc(2,k,1,3) = sn*(ci3+cr4)
               cc(1,k,1,3) = sn*(cr3-ci4)
               cc(1,k,1,4) = sn*(cr3+ci4)
               cc(2,k,1,4) = sn*(ci3-cr4)
               cc(2,k,1,5) = sn*(ci2-cr5)
            end do
         else
            do k=1,l1
               ti5 = cc(2,k,1,2)-cc(2,k,1,5)
               ti2 = cc(2,k,1,2)+cc(2,k,1,5)
               ti4 = cc(2,k,1,3)-cc(2,k,1,4)
               ti3 = cc(2,k,1,3)+cc(2,k,1,4)
               tr5 = cc(1,k,1,2)-cc(1,k,1,5)
               tr2 = cc(1,k,1,2)+cc(1,k,1,5)
               tr4 = cc(1,k,1,3)-cc(1,k,1,4)
               tr3 = cc(1,k,1,3)+cc(1,k,1,4)
               ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2+tr3)
               ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2+ti3)
               cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
               ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
               cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
               ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
               cr5 = TI11*tr5+TI12*tr4
               ci5 = TI11*ti5+TI12*ti4
               cr4 = TI12*tr5-TI11*tr4
               ci4 = TI12*ti5-TI11*ti4
               ch(1,k,2,1) = sn*(cr2-ci5)
               ch(1,k,5,1) = sn*(cr2+ci5)
               ch(2,k,2,1) = sn*(ci2+cr5)
               ch(2,k,3,1) = sn*(ci3+cr4)
               ch(1,k,3,1) = sn*(cr3-ci4)
               ch(1,k,4,1) = sn*(cr3+ci4)
               ch(2,k,4,1) = sn*(ci3-cr4)
               ch(2,k,5,1) = sn*(ci2-cr5)
            end do
         end if
      else
         do k=1,l1
            ti5 = cc(2,k,1,2)-cc(2,k,1,5)
            ti2 = cc(2,k,1,2)+cc(2,k,1,5)
            ti4 = cc(2,k,1,3)-cc(2,k,1,4)
            ti3 = cc(2,k,1,3)+cc(2,k,1,4)
            tr5 = cc(1,k,1,2)-cc(1,k,1,5)
            tr2 = cc(1,k,1,2)+cc(1,k,1,5)
            tr4 = cc(1,k,1,3)-cc(1,k,1,4)
            tr3 = cc(1,k,1,3)+cc(1,k,1,4)
            ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
            ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
            cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
            ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
            cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
            ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
            cr5 = TI11*tr5+TI12*tr4
            ci5 = TI11*ti5+TI12*ti4
            cr4 = TI12*tr5-TI11*tr4
            ci4 = TI12*ti5-TI11*ti4
            ch(1,k,2,1) = cr2-ci5
            ch(1,k,5,1) = cr2+ci5
            ch(2,k,2,1) = ci2+cr5
            ch(2,k,3,1) = ci3+cr4
            ch(1,k,3,1) = cr3-ci4
            ch(1,k,4,1) = cr3+ci4
            ch(2,k,4,1) = ci3-cr4
            ch(2,k,5,1) = ci2-cr5
         end do
         do i=2,ido
            do k=1,l1
               ti5 = cc(2,k,i,2)-cc(2,k,i,5)
               ti2 = cc(2,k,i,2)+cc(2,k,i,5)
               ti4 = cc(2,k,i,3)-cc(2,k,i,4)
               ti3 = cc(2,k,i,3)+cc(2,k,i,4)
               tr5 = cc(1,k,i,2)-cc(1,k,i,5)
               tr2 = cc(1,k,i,2)+cc(1,k,i,5)
               tr4 = cc(1,k,i,3)-cc(1,k,i,4)
               tr3 = cc(1,k,i,3)+cc(1,k,i,4)
               ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
               ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
               cr2 = cc(1,k,i,1)+TR11*tr2+TR12*tr3
               ci2 = cc(2,k,i,1)+TR11*ti2+TR12*ti3
               cr3 = cc(1,k,i,1)+TR12*tr2+TR11*tr3
               ci3 = cc(2,k,i,1)+TR12*ti2+TR11*ti3
               cr5 = TI11*tr5+TI12*tr4
               ci5 = TI11*ti5+TI12*ti4
               cr4 = TI12*tr5-TI11*tr4
               ci4 = TI12*ti5-TI11*ti4
               dr3 = cr3-ci4
               dr4 = cr3+ci4
               di3 = ci3+cr4
               di4 = ci3-cr4
               dr5 = cr2+ci5
               dr2 = cr2-ci5
               di5 = ci2-cr5
               di2 = ci2+cr5
               ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
               ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
               ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
               ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
               ch(1,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
               ch(2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
               ch(1,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
               ch(2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
            end do
         end do
      end if
!
   end subroutine complex_pass_5_forward
!
!
   subroutine complex_pass_n_forward(ido, iip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa)
!!
!!    Complex pass n forward
!!    Modified from https://github.com/jlokimlin/modern_fftpack by Andreas Skeidsvoll, Feb 2019:
!!    Changed variable types, style of comments, and spacing, to fit the style of eT
!!
      implicit none
!
      integer  :: ido
      integer  :: in1
      integer  :: in2
      integer  :: iip
      integer  :: l1
      integer  :: lid
!
      real(dp) :: cc(in1,l1,iip,ido)
      real(dp) :: cc1(in1,lid,iip)
      real(dp) :: ch(in2,l1,ido,iip)
      real(dp) :: ch1(in2,lid,iip)
      real(dp) :: chold1
      real(dp) :: chold2
      integer  :: i
      integer  :: idlj
      integer  :: iipp2
      integer  :: iipph
      integer  :: j
      integer  :: jc
      integer  :: ki
      integer  :: l
      integer  :: lc
      integer  :: na
      real(dp) :: sn
      real(dp) :: wa(ido,iip-1,2)
      real(dp) :: wai
      real(dp) :: war
!
      iipp2 = iip+2
      iipph = (iip+1)/2
!
      do ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
      end do
!
      do j=2,iipph
         jc = iipp2-j
         do ki=1,lid
            ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
         end do
      end do
!
      do j=2,iipph
         do ki=1,lid
            cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
            cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
         end do
      end do
!
      do l=2,iipph
         lc = iipp2-l
         do ki=1,lid
            cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
            cc1(1,ki,lc) = -wa(1,l-1,2)*ch1(1,ki,iip)
            cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
            cc1(2,ki,lc) = -wa(1,l-1,2)*ch1(2,ki,iip)
         end do
         do j=3,iipph
            jc = iipp2-j
            idlj = mod((l-1)*(j-1),iip)
            war = wa(1,idlj,1)
            wai = -wa(1,idlj,2)
            do ki=1,lid
               cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
               cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
               cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
               cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
            end do
         end do
      end do
!
      if (1 >= ido)then
         sn = ONE/(iip * l1)
         if (na /= 1) then
            do ki=1,lid
               cc1(1,ki,1) = sn*cc1(1,ki,1)
               cc1(2,ki,1) = sn*cc1(2,ki,1)
            end do
            do j=2,iipph
               jc = iipp2-j
               do ki=1,lid
                  chold1 = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
                  chold2 = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
                  cc1(1,ki,j) = chold1
                  cc1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
                  cc1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
                  cc1(1,ki,jc) = chold2
               end do
            end do
         else
            do ki=1,lid
               ch1(1,ki,1) = sn*cc1(1,ki,1)
               ch1(2,ki,1) = sn*cc1(2,ki,1)
            end do
            do j=2,iipph
               jc = iipp2-j
               do ki=1,lid
                  ch1(1,ki,j) = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
                  ch1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
                  ch1(1,ki,jc) = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
                  ch1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
               end do
            end do
         end if
      else
         do ki=1,lid
            ch1(1,ki,1) = cc1(1,ki,1)
            ch1(2,ki,1) = cc1(2,ki,1)
         end do
!
         do j=2,iipph
            jc = iipp2-j
            do ki=1,lid
               ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
               ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
               ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
               ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
            end do
         end do
!
         do i=1,ido
            cc(1,:,1,i) = ch(1,:,i,1)
            cc(2,:,1,i) = ch(2,:,i,1)
         end do
!
         do j=2,iip
            cc(1,:,j,1) = ch(1,:,1,j)
            cc(2,:,j,1) = ch(2,:,1,j)
         end do
!
         do j=2,iip
            do i=2,ido
               cc(1,:,j,i) = wa(i,j-1,1)*ch(1,:,i,j)+wa(i,j-1,2)*ch(2,:,i,j)
               cc(2,:,j,i) = wa(i,j-1,1)*ch(2,:,i,j)-wa(i,j-1,2)*ch(1,:,i,j)
            end do
         end do
!
      end if
!
   end subroutine complex_pass_n_forward
!
!
  end module fftpack_complex_forward_1d