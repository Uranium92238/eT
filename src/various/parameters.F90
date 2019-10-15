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
!  Integers 
!
   real(dp), parameter :: zero  = 0.0D0
   real(dp), parameter :: one   = 1.0D0
   real(dp), parameter :: two   = 2.0D0
   real(dp), parameter :: three = 3.0D0
   real(dp), parameter :: four  = 4.0D0
   real(dp), parameter :: five  = 5.0D0
   real(dp), parameter :: six   = 6.0D0
   real(dp), parameter :: seven = 7.0D0
   real(dp), parameter :: eight = 8.0D0
!
!  Fractions
!
   real(dp), parameter :: half           = 0.5D0
   real(dp), parameter :: quarter        = 0.25D0
   real(dp), parameter :: one_over_eight = one/eight
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
!
!  Complex fractions
!
   complex(dp), parameter :: half_complex           = cmplx(half, zero, dp)
   complex(dp), parameter :: quarter_complex        = cmplx(quarter, zero, dp)
   complex(dp), parameter :: one_over_eight_complex = cmplx(one_over_eight, zero, dp)
!
!  Conversion factors
!
   real(dp), parameter :: bohr_to_angstrom = 0.52917721092D0      ! 2010 CODATA
   real(dp), parameter :: angstrom_to_bohr = one/bohr_to_angstrom
   real(dp), parameter :: Hartree_to_eV    = 27.21138602D0        ! 2014 CODATA
                                                                  ! (https://physics.nist.gov/cgi-bin/cuu/Value?threv)
   real(dp), parameter :: Hartree_to_kcalmol = 627.509474D0       ! 2019 wikipedia
!
end module parameters
