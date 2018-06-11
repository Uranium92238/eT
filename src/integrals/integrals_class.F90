module integrals_class
!
!!
!!    Integrals class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use h_xy
   use file_class
   use memory_manager_class
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
   subroutine get_ao_xy_integrals(int, n_ao)
!!
!!    Get AO XY integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(integrals) :: int
!
      integer(i15) :: n_ao
!
      real(kind=8), dimension(:,:), allocatable :: h
!
      integer(i15) :: i = 0, j = 0
!
      call mem%alloc(h, n_ao, n_ao)
      h = zero
!
    !  write(output%unit,*) 'Hello from fortran, g is now: ', g
!
      call get_ao_xy_kinetics(h)
!
      do i = 1, n_ao
    !     do j = 1, n_ao
            write(output%unit,*) 'i i h_ii', i, i, h(i, i)
      !   enddo
      enddo
!
    !  write(output%unit,*) 'Back in fortran, g is now: ', g
!
   end subroutine get_ao_xy_integrals
!
end module integrals_class
