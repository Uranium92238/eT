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
submodule (doubles_class) complex_doubles
!
!!
!!    Complex submodule
!!
!!    Gathers routines that makes the abstract doubles wavefunction complex.
!!
!
   implicit none
!
!
contains
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
      if (allocated(wf%t2)) then
         call wf%initialize_t2_complex()
         wf%t2_complex = cmplx(wf%t2, zero, dp)
         call wf%destruct_t2()
      endif
!
      if (allocated(wf%t2bar)) then
         call wf%initialize_t2bar_complex()
         wf%t2bar_complex = cmplx(wf%t2bar, zero, dp)
         call wf%destruct_t2bar()
      endif
!
      if (allocated(wf%u_aibj)) then
         call wf%initialize_u_aibj_complex()
         wf%u_aibj_complex = cmplx(wf%u_aibj, zero, dp)
         call wf%destruct_u_aibj()
      endif
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
      if (allocated(wf%t2_complex)) then
         call wf%destruct_t2_complex()
      endif
!
      if (allocated(wf%t2bar_complex)) then
         call wf%destruct_t2bar_complex()
      endif
!
      if (allocated(wf%u_aibj_complex)) then
         call wf%destruct_u_aibj_complex()
      endif
!
   end subroutine cleanup_doubles_complex_doubles
!
!
end submodule complex_doubles