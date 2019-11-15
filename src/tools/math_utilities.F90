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
    use global_out, only : output
!
   implicit none
!
contains
!
   recursive function factorial_recursive(arg) result(f)
!!
!!    Recursive implementation of faculty operation
!!
!!    Written by Sarai D. Folkestad, Dec 2018
!!
      implicit none
!
      integer, intent(in) :: arg
      integer :: f
!
      f = 1
!
      if (arg == 0 ) then
!
         f = 1
         return
!
      elseif (arg .lt. 0) then
!
         call output%error_msg('factorial operation not defined for negative numbers!')
!
      else
!
         f = arg*factorial_recursive(arg - 1)
!
      endif
!
   end function factorial_recursive
!
!
   function factorial(arg) result(f)
!!
!!    Factorial operation
!!
!!    Written by Sarai D. Folkestad, Dec 2018
!!
      implicit none
!
      integer, intent(in) :: arg
      integer :: f
!
      if (arg .lt. 0) then
!
         call output%error_msg('factorial operation not defined for negative numbers!')
!
      endif
!
      f = factorial_recursive(arg)
!
   end function factorial
!
!
   recursive function double_factorial_recursive(arg) result(f)
!!
!!    Recursive implementation of double factorial operation
!!
!!    Written by Andreas Skeidsvoll, Jul 2019
!
      implicit none
!
      integer, intent(in) :: arg
      integer :: f
!
      f = 1
!
      if ((arg == 0) .or. (arg == -1)) then
!
         f = 1
         return
!
      elseif (arg .lt. -1) then
!
         call output%error_msg('double factorial operation not defined for numbers less than -1!')
!
      else
!
         f = arg*double_factorial_recursive(arg - 2)
!
      endif
!
   end function double_factorial_recursive
!
!
   function double_factorial(arg) result(f)
!!
!!    Double factorial operation
!!
!!    Written by Andreas Skeidsvoll, Jul 2019
!
      implicit none
!
      integer, intent(in) :: arg
      integer :: f
!
      if (arg .lt. -1) then
!
         call output%error_msg('double factorial operation not defined for numbers less than -1!')
!
      endif
!
      f = double_factorial_recursive(arg)
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
   function binomial(n, k) result(b)
!!
!!    Binomial coefficient
!!
!!    Written by Sarai D. Folkestad, Dec 2018
!!
!
      implicit none 
!
      integer, intent(in) :: n, k 
      integer :: b
!
      if (k .gt. n) then
!
         call output%error_msg('k ((i0)) is greater than n ((i0)) in binomial function.', &
                               ints=[k,n])
!
      endif
!
      if (k == 0) then
!
         b = 1
         return
!
      endif
!
      b = factorial(n)/(factorial(k)*factorial(n-k))
!
   end function binomial
!
!
end module math_utilities
