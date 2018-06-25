module integrals_class
!
!!
!!    Integrals class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
!  Fortran interfaces to C++ routines
!
   use h_xy
   use s_xy
   use g_xyzw
   use L_xyzw
   use libint_initialization
!
!  Disk & memory class modules
!
   use file_class
   use memory_manager_class
!
   implicit none
!
!  Class definition
!
   type :: integrals
!
!     No attributes yet
!
   contains
!
      procedure :: get_ao_h_xy   => get_ao_h_xy_integrals   ! h_αβ
      procedure :: get_ao_s_xy   => get_ao_s_xy_integrals   ! s_αβ
      procedure :: get_ao_g_xyzw => get_ao_g_xyzw_integrals ! g_αβγδ
!
   end type integrals
!
!
contains
!
!
   subroutine get_ao_h_xy_integrals(int, h)
!!
!!    Get h_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integrals in the array h.
!!
      implicit none
!
      class(integrals) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: h
!
      call get_ao_h_xy(h)
!
   end subroutine get_ao_h_xy_integrals
!
!
   subroutine get_ao_s_xy_integrals(int, s)
!!
!!    Get s_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the s_αβ integrals in the array s.
!!
      implicit none
!
      class(integrals) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: s
!
      call get_ao_s_xy(s)
!
   end subroutine get_ao_s_xy_integrals
!
!
   subroutine get_ao_g_xyzw_integrals(int, g)
!!
!!    Get g_αβγδ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integrals in the array g.
!!
      implicit none
!
      class(integrals) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: g
!
      call get_ao_g_xyzw(g)
!
   end subroutine get_ao_g_xyzw_integrals
!
!
end module integrals_class
