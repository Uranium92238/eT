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
module math_utilities
!
!!
!!    Math utilities module
!!
!!    This module contains routines that perform various math operations.
!!
!
   use parameters
!
   implicit none
!
contains
!
!
   integer function double_factorial(angmom)
!!
!!    Double Factorial
!!    Written by Tommaso Giovannini, March 2019
!!
!!    Calculates double factorial n!!
!!
      implicit none 
!
      integer :: i, angmom
!  
      double_factorial = 1
!      
      do i = angmom, 0, -2
!      
         if(i.eq.0.or.i.eq.1) then
!      
           double_factorial = double_factorial
!      
         else
!      
           double_factorial = double_factorial * i
!      
        endif
!     
      enddo
!
    end function double_factorial
!
!
    function delta(i,j) result(delta_ij)
!!
!!    Delta    
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Returns the Kronecker delta as a real number:
!!
!!      delta_ij = 1 if i .eq. j
!!      delta_ij = 0 if i .ne. j
!!
      implicit none 
!
      real(dp) :: delta_ij 
!
      integer, intent(in) :: i, j
!
      if (i == j) then 
!
        delta_ij = one
!
      else 
!
        delta_ij = zero
!
      endif
!
    end function delta
!
!
end module math_utilities
