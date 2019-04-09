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
submodule (molecular_system_class) ao_integrals
!
!!
!!    AO integrals submodule
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2019
!!
!!    Routines to handle AO integrals. Provides an interface to Libint to eT developers.
!!  
!!    Note: C++ is not fond of 64-bit integers, so ints are explicitly 
!!    translated to 32-bits here before calling Libint routines 
!!
!
   implicit none
!
   include "../../libint/h_wx_cdef.F90"
!
!
contains
!
!
   module subroutine construct_ao_h_wx_molecular_system(molecule, h, s1, s2)
!!
!!    Construct h_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
      integer, intent(in) :: s1, s2 
!
      real(dp), dimension(molecule%shell_limits(s1)%size, molecule%shell_limits(s2)%size), intent(out) :: h 
!
      integer(i6) :: s1_4, s2_4 
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_h_wx_c(h, s1_4, s2_4)
!
   end subroutine construct_ao_h_wx_molecular_system
!
!
   module subroutine construct_ao_s_wx_ao_molecular_system(s, s1, s2)
!!
!!    Construct s_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves parts of the s_αβ integral in the array s. s1 and s2 are the shells
!!    that w and x, respectively belong to.
!!
      implicit none
!
      class(molecular_system), intent(in) :: molecule
!
      real(dp), dimension(molecule%shell_limits(s1),molecule%shell_limits(s1)), intent(inout) :: s
!
      integer, intent(in) :: s1, s2
!
      integer(i6) :: s1_4, s2_4 ! Integers that are passed to libint
!
      s1_4 = int(s1,i6)
      s2_4 = int(s2,i6)
!
      call construct_ao_s_wx_c(s, s1_4, s2_4) 
!
   end subroutine construct_ao_s_wx_ao_molecular_system
!
!
end submodule ao_integrals
