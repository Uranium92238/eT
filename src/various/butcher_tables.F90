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
module butcher_tables
!
!!
!! Butcher tables module
!! Written by Andreas Skeidsvoll, Oct 2019
!!
!! Contains parameters that give the values of Butcher tables (used in Runge-Kutta methods).
!! Butcher tables specify the given Runge-Kutta method, and are given on the form
!!
!!    c1 | a11 a12 ... a1n
!!    c2 | a21 a22 ... a2n
!!     : |  :   :       :
!!    cn | an1 an2 ... ann
!!       -----------------
!!         b1  b2  ... bn
!!
!
   use kinds
!
   implicit none
!
!  2nd order Gauss-Legendre table
!
   real(dp), parameter :: a11_gl2 = 0.5_dp,   &
                          b1_gl2 = 1.0_dp,    &
                          c1_gl2 = 0.5_dp
!
!  4th order Gauss-Legendre table
!
   real(dp), parameter :: a11_gl4 = 0.25_dp,                          &
                          a12_gl4 = 0.25_dp - sqrt(3.0_dp)/6.0_dp,    &
                          a21_gl4 = 0.25_dp + sqrt(3.0_dp)/6.0_dp,    &
                          a22_gl4 = 0.25_dp,                          &
                          b1_gl4 = 0.5_dp,                            &
                          b2_gl4 = 0.5_dp,                            &
                          c1_gl4 = 0.5_dp - sqrt(3.0_dp)/6.0_dp,      &
                          c2_gl4 = 0.5_dp + sqrt(3.0_dp)/6.0_dp
!
!  6th order Gauss-Legendre table
!
   real(dp), parameter :: a11_gl6 = 5.0_dp/36.0_dp,                          &
                          a12_gl6 = 2.0_dp/9.0_dp - 1.0_dp/sqrt(15.0_dp),    &
                          a13_gl6 = 5.0_dp/36.0_dp - sqrt(15.0_dp)/30.0_dp,  &
                          a21_gl6 = 5.0_dp/36.0_dp + sqrt(15.0_dp)/24.0_dp,  &
                          a22_gl6 = 2.0_dp/9.0_dp,                           &
                          a23_gl6 = 5.0_dp/36.0_dp - sqrt(15.0_dp)/24.0_dp,  &
                          a31_gl6 = 5.0_dp/36.0_dp + sqrt(15.0_dp)/30.0_dp,  &
                          a32_gl6 = 2.0_dp/9.0_dp + 1.0_dp/sqrt(15.0_dp),    &
                          a33_gl6 = 5.0_dp/36.0_dp,                          &
                          b1_gl6 = 5.0_dp/18.0_dp,                           &
                          b2_gl6 = 4.0_dp/9.0_dp,                            &
                          b3_gl6 = 5.0_dp/18.0_dp,                           &
                          c1_gl6 = 0.5_dp - 0.1_dp*sqrt(15.0_dp),            &
                          c2_gl6 = 0.5_dp,                                   &
                          c3_gl6 = 0.5_dp + 0.1_dp*sqrt(15.0_dp)
!
end module butcher_tables