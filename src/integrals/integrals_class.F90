module integrals_class
!
!!
!!    Integrals class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   implicit none
!
   type :: integrals
!
   contains
!
      procedure :: get_ao_xy => get_ao_xy_integrals
!
   end type integrals
!
contains
!
   subroutine get_ao_xy_integrals(int)
!!
!!    Get AO XY integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(integrals) :: int
!
   end subroutine get_ao_xy_integrals
!
end module integrals_class
