module integrals_class
!
!!
!!    Integrals class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use h_xy ! For one electron integrals, h_xy
   use s_xy ! For one electron overlaps, s_xy
!
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
      real(kind=8), dimension(:,:), allocatable :: h, s
!
      real(kind=8), dimension(:,:), allocatable :: temp
!
      integer(i15) :: i = 0, j = 0
!
!     Get h_xy
!
      call mem%alloc(h, n_ao, n_ao)
      h = zero
!
      call get_ao_xy(h)
!
!     Get s_xy
!
      call mem%alloc(s, n_ao, n_ao)
      s = zero
!
      call get_ao_s_xy(s)
!
!     Print norms and diagonal elements
!
      do i = 1, n_ao
         write(output%unit, *) 'i s_ii h_ii', i, s(i,i), h(i,i)
      enddo
!
   end subroutine get_ao_xy_integrals
!
end module integrals_class
