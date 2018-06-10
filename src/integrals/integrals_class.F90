module integrals_class
!
!!
!!    Integrals class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinetic
   use file_class
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
      real(kind=8), dimension(2,1) :: g
!
      g(1,1) = 1.2D0
      g(2,1) = 1.3D0
!
      write(output%unit,*) 'Hello from fortran, g is now: ', g
!
      call get_ao_xy_kinetics(g)
!
      write(output%unit,*) 'Back in fortran, g is now: ', g
!
   end subroutine get_ao_xy_integrals
!
end module integrals_class
