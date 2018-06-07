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
   real(dp), parameter :: zero  = 0.0D0
   real(dp), parameter :: half  = 0.5D0
   real(dp), parameter :: two   = 2.0D0
   real(dp), parameter :: three = 3.0D0
   real(dp), parameter :: four  = 4.0D0
   real(dp), parameter :: five  = 5.0D0
   real(dp), parameter :: six   = 6.0D0
!
end module parameters
