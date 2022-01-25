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
   function cross_product_R3(r1, r2) result(r3)
!!
!!    Cross product R^3
!!    Written by Sarai D. Folkestad Oct 2020
!!
      implicit none
!
      real(dp), dimension(3), intent(in) :: r1, r2
      real(dp), dimension(3) :: r3
!
      r3(1) =   r1(2)*r2(3) - r1(3)*r2(2)   !   y1 * z2 - y2 * z1
      r3(2) = - r1(1)*r2(3) + r1(3)*r2(1)   ! - x1 * z2 + x2 * z1
      r3(3) =   r1(1)*r2(2) - r1(2)*r2(1)   !   y1 * z2 - y2 * z1
!
   end function cross_product_R3
!
!
   function dot_R3(r1, r2) result(r1_dot_r2)
!!
!!    Dot product R^3
!!    Written by Sarai D. Folkestad Oct 2020
!!
      implicit none
!
      real(dp), dimension(3), intent(in) :: r1, r2
!
      real(dp) :: r1_dot_r2
!
      r1_dot_r2 =  r1(1) * r2(1) &
                 + r1(2) * r2(2) &
                 + r1(3) * r2(3)
!
   end function dot_R3
!
!
   function norm_R3(r) result(r_norm)
!!
!!    Norm R^3
!!    Written by Sarai D. Folkestad Oct 2020
!!
      implicit none
!
      real(dp), dimension(3), intent(in) :: r
!
      real(dp) :: r_norm
!
      r_norm =  sqrt(r(1)**2 + r(2)**2 + r(3)**2)
!
   end function norm_R3
!
!
   elemental function gcd(a,b) result(d)
!!
!!    Greatest common divisor
!!
!!    Written by Rolf H. Myhre, April 2021
!!
!!    Uses the Eucledean Algorithm
!!    to compute the greatest common divisor of a and b
!!    See https://en.wikipedia.org/wiki/Euclidean_algorithm for details
!!
      implicit none
!
      integer, intent(in) :: a,b
!
      integer :: tmp1, tmp2
      integer :: d
!
      d    = max(abs(a),abs(b))
      tmp1 = min(abs(a),abs(b))
!
      do while (tmp1 .ne. 0)
         tmp2 = tmp1
         tmp1 = mod(d, tmp2)
         d    = tmp2
      enddo
!
   end function gcd
!
!
   elemental function lcm(a,b) result(m)
!!
!!    Least common multiple
!!
!!    Written by Rolf H. Myhre, April 2021
!!
!!    Computes the least common multiple of a and b
!!    See https://en.wikipedia.org/wiki/Least_common_multiple for details
!!
      implicit none
!
      integer, intent(in) :: a,b
!
      integer :: m
!
      m = (abs(a)/gcd(a,b))*abs(b)
!
   end function lcm
!
!
   pure recursive function solve_n1s1_n2s2(s1, s2, x) result(n1)
!!
!!    Solve n1*s1 - n2*s2
!!
!!    Written by Rolf H. Myhre, May 2021
!!
!!    Returns an n1 solving the equation
!!    n1*s1 - n2*s2 = x
!!    if they exists.
!!    If they do not exist, the function returns garbage.
!!    Solutions exist iff x/gcd(s1,s2) is an integer.
!!
!!    The algorithm works by by exploiting that s1 = q*s2 + r
!!    Inserting this in the original expression, we get
!!    n1'*s2 - n2'*r = x
!!    with n1' = (n1*q - n2) and n2' = -n1
!!    This is repeated until x is dividable by r and n1 is 0.
!!
      implicit none
!
      integer, intent(in) :: s1, s2, x
      integer :: n1, n1p
!
      if (mod(x,s2) .eq. 0) then
         n1 = 0
      else
         n1p = solve_n1s1_n2s2(s2, mod(s1,s2), x)
         n1 = -(n1p*s2 - x)/mod(s1,s2)
      end if
!
   end function solve_n1s1_n2s2
!
!
   pure recursive function least_positive_n1s1_n2s2(s1, s2, x) result(n1)
!!
!!    Least positive n1*s1 - n2*s2
!!
!!    Written by Rolf H. Myhre, May 2021
!!
!!    Solve n1*s1 - n2*s2 = x
!!
!!    and return the least positive n1 corresponding to a positive n2
!!
      implicit none
!
      integer, intent(in) :: s1, s2, x
      integer :: n1, n2, l, c1, c2
!
      l = lcm(s1,s2)
      c1 = l/s1
      c2 = l/s2
!
      n1 = solve_n1s1_n2s2(s1, s2, x)
!
      n1 = n1 - (n1/c1)*c1
      if (n1 .lt. 0) n1 = n1 + c1
!
      n2 = (n1*s1 - x)/s2
!
      if (n2 .lt. 0) then
         c2 = l/s2
         n1 = n1 - (n2/c2)*c1
         if (mod(n2,c2) .ne. 0) n1 = n1 + c1
      endif
!
   end function least_positive_n1s1_n2s2
!
!
   pure function get_n_values_greater_sum(dim_, values, sum_) result(n)
!!
!!    Get n values greater sum
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Get the n first values that add up to sum or more
!!
      implicit none 
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_), intent(in) :: values
!
      real(dp), intent(in) :: sum_
!
      integer  :: n, p
      real(dp) :: local_sum
!
      local_sum = zero
!
      do p = 1, dim_
!
         if (local_sum .ge. sum_) exit
         local_sum = local_sum + values(p)
!
      end do
!
      n = p - 1
!
   end function get_n_values_greater_sum
!
!
end module math_utilities
