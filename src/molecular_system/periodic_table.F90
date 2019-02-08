module periodic_table
!
!!
!!    Periodic table module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!
!
   use kinds
!
   implicit none
!
   integer, parameter :: size_periodic_table = 20
   character(len=2), dimension(size_periodic_table) :: atomic_symbol = &
   ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', &
    'S ',  'Cl', 'Ar', 'K ', 'Ca']
!
end module periodic_table
