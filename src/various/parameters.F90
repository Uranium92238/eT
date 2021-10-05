!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module parameters
!
!!
!!    Parameters module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Defines a set of useful parameters, such as integers defined
!!    as double precision numbers (to avoid issues associated with
!!    integer*double precision, as in, e.g., 2*g_aibj-g_ajbi).
!!
   use kinds
!
   implicit none
!
!  Version
!
   integer, parameter :: major_version = 1
   integer, parameter :: minor_version = 5
   integer, parameter :: patch_version = 0
!
   character(len=27), parameter :: version_name = "Furia (development)"
!
!  Integers
!
   real(dp), parameter :: zero   =  0.0D0
   real(dp), parameter :: one    =  1.0D0
   real(dp), parameter :: two    =  2.0D0
   real(dp), parameter :: three  =  3.0D0
   real(dp), parameter :: four   =  4.0D0
   real(dp), parameter :: five   =  5.0D0
   real(dp), parameter :: six    =  6.0D0
   real(dp), parameter :: seven  =  7.0D0
   real(dp), parameter :: eight  =  8.0D0
   real(dp), parameter :: nine   =  9.0D0
   real(dp), parameter :: ten    = 10.0D0
   real(dp), parameter :: eleven = 11.0D0
!
!  Fractions
!
   real(dp), parameter :: half    = one/two
   real(dp), parameter :: third   = one/three
   real(dp), parameter :: quarter = one/four
   real(dp), parameter :: eighth  = one/eight
   real(dp), parameter :: tenth   = one/ten
!
!  Pi
!
   real(dp), parameter :: pi = four * atan(one)
   real(dp), parameter :: radians_to_degrees = 180.0d0/pi
   real(dp), parameter :: degrees_to_radians = pi/180.0d0
!
!  Square roots
!
   real(dp), parameter :: sqrt_3 = sqrt(three)
   real(dp), parameter :: sqrt_5 = sqrt(five)
!
   real(dp), parameter :: inv_sqrt_3  = one/sqrt_3
   real(dp), parameter :: inv_sqrt_5  = one/sqrt_5
   real(dp), parameter :: inv_sqrt_7  = one/sqrt(seven)
   real(dp), parameter :: inv_sqrt_9  = one/sqrt(nine)
   real(dp), parameter :: inv_sqrt_11 = one/sqrt(eleven)
!
!  Complex integers
!
   complex(dp), parameter :: zero_complex  = cmplx(zero, zero, dp)
   complex(dp), parameter :: one_complex   = cmplx(one, zero, dp)
   complex(dp), parameter :: two_complex   = cmplx(two, zero, dp)
   complex(dp), parameter :: three_complex = cmplx(three, zero, dp)
   complex(dp), parameter :: four_complex  = cmplx(four, zero, dp)
   complex(dp), parameter :: five_complex  = cmplx(five, zero, dp)
   complex(dp), parameter :: six_complex   = cmplx(six, zero, dp)
   complex(dp), parameter :: seven_complex = cmplx(seven, zero, dp)
   complex(dp), parameter :: eight_complex = cmplx(eight, zero, dp)
   complex(dp), parameter :: ten_complex   = cmplx(ten, zero, dp)
!
!  Complex fractions
!
   complex(dp), parameter :: half_complex    = cmplx(half, zero, dp)
   complex(dp), parameter :: third_complex   = cmplx(third, zero, dp)
   complex(dp), parameter :: quarter_complex = cmplx(quarter, zero, dp)
   complex(dp), parameter :: eighth_complex  = cmplx(eighth, zero, dp)
   complex(dp), parameter :: tenth_complex   = cmplx(tenth, zero, dp)
!
!  Conversion factors
!
   real(dp), parameter :: bohr_to_angstrom = 0.52917721092D0      ! 2010 CODATA
   real(dp), parameter :: angstrom_to_bohr = one/bohr_to_angstrom
   real(dp), parameter :: Hartree_to_eV    = 27.21138602D0        ! 2014 CODATA
                                                                  ! (https://physics.nist.gov/cgi-bin/cuu/Value?threv)
   real(dp), parameter :: Hartree_to_kcalmol = 627.509474D0       ! 2019 wikipedia
   real(dp), parameter :: au_to_debye        = 2.541746473D0      ! 2018 CODATA
!
   integer, parameter :: int32_mo_limit = 215
!
!
contains
!
!
   pure function get_conversion_factor(from, to) result(conversion_factor)
!!
!!    Get conversion factor
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns the conversion factor needed to convert from units
!!    'from' to units 'to' (angstrom/bohr).
!!
!!    Returns zero if from or two are not recognized units.
!!
      implicit none
!
      character(len=*), intent(in) :: from, to
!
      real(dp) :: conversion_factor
!
      if (trim(from) .eq. 'bohr' .and. trim(to) .eq. 'angstrom') then
!
         conversion_factor = bohr_to_angstrom
!
      elseif (trim(from) .eq. 'angstrom' .and. trim(to) .eq. 'bohr') then
!
         conversion_factor = angstrom_to_bohr
!
      else
!
         conversion_factor = one
!
      endif
!
      if (trim(from) .ne. 'angstrom' .and. trim(from) .ne. 'bohr' .or. &
          trim(to)   .ne. 'angstrom' .and. trim(to)   .ne. 'bohr') then
!
         conversion_factor = zero
!
      endif
!
   end function get_conversion_factor
!
!
   pure function get_units_label(units) result(label)
!!
!!    Get units label
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns 'a.u.' if units = 'bohr'
!!    Returns 'angstrom' if units = 'angstrom'
!!    Returns 'unknown' otherwise
!!
      implicit none
!
      character(len=*), intent(in) :: units
!
      character(len=100) :: label
!
      if (trim(units) == 'angstrom') then
!
         label = 'angstrom'
!
      elseif (trim(units) == 'bohr') then
!
         label = 'a.u.'
!
      else
!
         label = 'unknown'
!
      endif
!
   end function get_units_label
!
!
end module parameters
