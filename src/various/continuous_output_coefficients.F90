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
module continuous_output_coefficients
!
!!
!! Continuous output coefficients module
!! Written by Andreas Skeidsvoll, Oct 2019
!!
!! Contains parameters that give the values of coefficients used to accelerate the convergence
!! of equations in implicit integration methods, see Section VIII.6.1 (A) in Hairer et. al.,
!! Geometric Numerical Integration.
!!
!
   use kinds
!
   implicit none
!
!  2nd order Gauss-Legendre coefficient
!
   real(dp), parameter :: beta11_gl2 = 1.5_dp
!
!  4th order Gauss-Legendre coefficients
!
   real(dp), parameter :: beta11_gl4 = 1.25_dp - 0.5_dp*sqrt(3.0_dp), &
                          beta12_gl4 = 0.25_dp + 1.0_dp/sqrt(3.0_dp), &
                          beta21_gl4 = 0.25_dp - 1.0_dp/sqrt(3.0_dp), &
                          beta22_gl4 = 1.25_dp + 0.5_dp*sqrt(3.0_dp)
!
!  6th order Gauss-Legendre coefficients
!
   real(dp), parameter :: beta11_gl6 = 9.0_dp/4.0_dp - sqrt(15.0_dp)/2.0_dp, &
                          beta12_gl6 = -2.0_dp + 3.0_dp*sqrt(0.6_dp),        &
                          beta13_gl6 = 1.25_dp - sqrt(0.6_dp),               &
                          beta21_gl6 = 1.25_dp - sqrt(15.0_dp)/8.0_dp,       &
                          beta22_gl6 = -1.0_dp,                              &
                          beta23_gl6 = 1.25_dp + sqrt(15.0_dp)/8.0_dp,       &
                          beta31_gl6 = 1.25_dp + sqrt(0.6_dp),               &
                          beta32_gl6 = -2.0_dp - 3.0_dp*sqrt(0.6_dp),        &
                          beta33_gl6 = 9.0_dp/4.0_dp + sqrt(15.0_dp)/2.0_dp
!
end module continuous_output_coefficients