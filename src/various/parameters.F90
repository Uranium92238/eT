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
!  Integers and half-integers
!
   real(dp), parameter :: zero  = 0.0D0
   real(dp), parameter :: half  = 0.5D0
   real(dp), parameter :: one   = 1.0D0
   real(dp), parameter :: two   = 2.0D0
   real(dp), parameter :: three = 3.0D0
   real(dp), parameter :: four  = 4.0D0
   real(dp), parameter :: five  = 5.0D0
   real(dp), parameter :: six   = 6.0D0
!
!  Conversion factors
!
   real(dp), parameter :: bohr_to_angstrom = 0.52917721067D0 ! 2014 CODATA
   real(dp), parameter :: angstrom_to_bohr = one/bohr_to_angstrom
!
end module parameters
