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
!
contains
!
!
   subroutine construct_ao_h_wx_molecular_system(molecule, h, s1, s2)
!!
!!    Construct h_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(molecular_system), intent(in) :: molecule 
!
   
!
   end subroutine construct_ao_h_wx_ao_integral_tool
!
!
end submodule ao_integrals
