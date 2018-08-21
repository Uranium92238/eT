module ao_integral_manager_class
!
!!
!!    Integral_manager class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
!  Fortran interfaces to C++ routines
!
   use h_xy
   use s_xy
   use g_wxyz
!
   !use libint_initialization
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
   type :: ao_integral_manager
!
!     No attributes yet
!
   contains
!
      procedure :: get_ao_h_xy   => get_ao_h_xy_ao_integral_manager   ! h_αβ
      procedure :: get_ao_s_xy   => get_ao_s_xy_ao_integral_manager   ! s_αβ
      procedure :: get_ao_g_wxyz => get_ao_g_wxyz_ao_integral_manager ! g_αβγδ
      procedure :: get_ao_g_wxyz_epsilon => get_ao_g_wxyz_epsilon_ao_integral_manager ! g_αβγδ
!
   end type ao_integral_manager
!
!
contains
!
!
   subroutine get_ao_h_xy_ao_integral_manager(int, h)
!!
!!    Get h_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral_manager in the array h.
!!
      implicit none
!
      class(ao_integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: h
!
      call get_ao_h_xy(h)
!
   end subroutine get_ao_h_xy_ao_integral_manager
!
!
   subroutine get_ao_s_xy_ao_integral_manager(int, s)
!!
!!    Get s_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the s_αβ integral_manager in the array s.
!!
      implicit none
!
      class(ao_integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: s
!
      call get_ao_s_xy(s)
!
   end subroutine get_ao_s_xy_ao_integral_manager
!
!
   subroutine get_ao_g_wxyz_ao_integral_manager(int, g, s1, s2, s3, s4)
!!
!!    Get g_αβγδ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral_manager in the array g.
!!
!!    s1 is first shell index and so on.
!!
      implicit none
!
      class(ao_integral_manager), intent(in) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: g
!
      integer(kind=8), intent(in) :: s1, s2, s3, s4
!
      call get_ao_g_wxyz(g, s1, s2, s3, s4)
!
   end subroutine get_ao_g_wxyz_ao_integral_manager
!
!
   subroutine get_ao_g_wxyz_epsilon_ao_integral_manager(int, g, s1, s2, s3, s4, epsilon, thread, skip, &
                                                            n1, n2, n3, n4)
!!
!!    Get g_αβγδ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral_manager in the array g.
!!
!!    s1 is first shell index and so on.
!!
      implicit none
!
      class(ao_integral_manager), intent(in) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: g
!
      real(dp), intent(in) :: epsilon 
!
      integer(kind=8), intent(in) :: s1, s2, s3, s4, thread, n1, n2, n3, n4 
      integer(kind=8) :: skip 
!
      call get_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, epsilon, thread, skip, n1, n2, n3, n4)
!
   end subroutine get_ao_g_wxyz_epsilon_ao_integral_manager
!
!
end module ao_integral_manager_class

