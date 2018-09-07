module ao_integral_tool_class
!
!!
!!    Integral_tool class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
!  Fortran interfaces to C++ routines
!
   use h_xy
   use s_xy
   use g_wxyz
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
   type :: ao_integral_tool
!
!     No attributes yet
!
   contains
!
!     These are to be replaced with get routines in the WF - we shall not do 
!     do loops on the C++ side, only ask for arrays to handle Fortran side 
!
      procedure :: get_ao_h_xy           => get_ao_h_xy_ao_integral_tool            ! h_αβ
!
      procedure :: construct_ao_h_wx           => construct_ao_h_wx_ao_integral_tool            ! h_αβ
      procedure :: construct_ao_s_wx           => construct_ao_s_wx_ao_integral_tool            ! h_αβ
      procedure :: construct_ao_g_wxyz         => construct_ao_g_wxyz_ao_integral_tool          ! g_αβγδ
      procedure :: construct_ao_g_wxyz_epsilon => construct_ao_g_wxyz_epsilon_ao_integral_tool  ! g_αβγδ
!
   end type ao_integral_tool
!
!
contains
!
!
   subroutine get_ao_h_xy_ao_integral_tool(int, h)
!!
!!    Get h_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral_tool in the array h.
!!
      implicit none
!
      class(ao_integral_tool) :: int
!
      real(dp), dimension(:,:), intent(inout) :: h
!
      call get_ao_h_xy(h)
!
   end subroutine get_ao_h_xy_ao_integral_tool
!
!
   subroutine construct_ao_h_wx_ao_integral_tool(int, h, s1, s2)
!!
!!    Construct h_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral_tool in the array h.
!!
      implicit none
!
      class(ao_integral_tool) :: int
!
      real(dp), dimension(:,:), intent(inout) :: h
      integer(i15), intent(in) :: s1, s2
!
      call get_ao_h_xy_sp(h, s1, s2)
!
   end subroutine construct_ao_h_wx_ao_integral_tool
!
!
   subroutine construct_ao_s_wx_ao_integral_tool(int, s, s1, s2)
!!
!!    Construct s_αβ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral_tool in the array h.
!!
      implicit none
!
      class(ao_integral_tool) :: int
!
      real(dp), dimension(:,:), intent(inout) :: s
      integer(i15), intent(in) :: s1, s2
!
      call construct_ao_s_wx(s, s1, s2) 
!
   end subroutine construct_ao_s_wx_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_ao_integral_tool(int, g, s1, s2, s3, s4)
!!
!!    Construct g_αβγδ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral_tool in the array g. s1-s4 are 
!!    the shells that alpha, beta, gamma and delta belong to. 
!!
      implicit none
!
      class(ao_integral_tool), intent(in) :: int
!
      real(dp), dimension(:,:), intent(inout) :: g
!
      integer(dp), intent(in) :: s1, s2, s3, s4
!
      call get_ao_g_wxyz(g, s1, s2, s3, s4)
!
   end subroutine construct_ao_g_wxyz_ao_integral_tool
!
!
   subroutine construct_ao_g_wxyz_epsilon_ao_integral_tool(int, g, s1, s2, s3, s4, epsilon, thread, skip, &
                                                            n1, n2, n3, n4)
!!
!!    Construct g_αβγδ epsilon 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral_tool in the array g. s1-s4 are 
!!    the shells that alpha, beta, gamma and delta belong to. 
!!
!!    This is the most efficient routine to calculate these 
!!    integrals, due mostly to the precision control parameter 
!!    (epsilon) but also because we avoid computing information
!!    that might already be available (such as thread ID and 
!!    the size of each of the shells, n1-n4). 
!!
!!       - To get thread, use omp_get_thread_num()
!!       - To get n1-n4, see the shells array of the wavefunction's system object
!!       - Skip is an integer that is either 0 or 1 on exit, where a 0 means 
!!         Libint decided to skip calculating the integrals. In other words,
!!         be sure to zero g if skip is 1 and this is necessary for your 
!!         application! The reason for this integer is that it is sometimes
!!         not necessary to actually zero g (which can be a relevant penalty).
!!       - In order to determine epsilon, one useful approach is consider the 
!!         desired accuracy epsilon' = 1.0D-14 of Y, where e.g.
!!
!!             Y_αβ = Y_αβ + g_αβγδ X_γδ.
!!
!!         If X_γδ is on the order 1.0D-5, then to get 1.0D-14 accuracy in g_αβγδ X_γδ
!!         requires only that g_αβγδ is accurate to 1.0D-9 (the error is multiplied 
!!         by X_γδ). The Libint integral can be much faster for large epsilons,
!!         but care should be taken when dynamically changing epsilon.
!!
      implicit none
!
      class(ao_integral_tool), intent(in) :: int
!
      real(dp), dimension(:,:), intent(inout) :: g
!
      real(dp), intent(in) :: epsilon 
!
      integer(i15), intent(in) :: s1, s2, s3, s4, thread, n1, n2, n3, n4 
      integer(i15) :: skip 
!
      call get_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, epsilon, thread, skip, n1, n2, n3, n4)
!
   end subroutine construct_ao_g_wxyz_epsilon_ao_integral_tool
!
!
end module ao_integral_tool_class

