!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module bfgs_tool_class
!
!!
!!    BFGS solver tool class module
!!    Written by Eirik F. Kjønstad, 2019
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
   use array_utilities
!
!
   type :: bfgs_tool
!
      private 
      integer  :: iteration 
      integer  :: n_parameters
!
      real(dp), dimension(:,:), allocatable :: Hessian
!
      real(dp), dimension(:), allocatable :: prev_g 
      real(dp), dimension(:), allocatable :: prev_x 
!
   contains
!
      procedure, public :: get_direction => get_direction_bfgs_tool
      procedure, public :: set_previous_geometry_and_gradient => set_previous_geometry_and_gradient_bfgs_tool
      procedure, public :: update_hessian => update_hessian_bfgs_tool
!
   end type bfgs_tool
!
!
   interface bfgs_tool 
!
      procedure :: new_bfgs_tool
!
   end interface bfgs_tool 
!
!
contains
!
!
   function new_bfgs_tool(n_parameters) result(bfgs)
!!
!!    New BFGS tool
!!    Written by Eirik F. Kjønstad, 2019
!!
      type(bfgs_tool) :: bfgs
!
      integer, intent(in) :: n_parameters
!  
      integer :: k 
!
      bfgs%n_parameters = n_parameters
      bfgs%iteration    = 0
!
      call mem%alloc(bfgs%prev_g, bfgs%n_parameters)
      call mem%alloc(bfgs%prev_x, bfgs%n_parameters)
!
      bfgs%prev_g = zero 
      bfgs%prev_x = zero 
!
      call mem%alloc(bfgs%Hessian, bfgs%n_parameters, bfgs%n_parameters)
!
      bfgs%Hessian = zero 
!
      do k = 1, bfgs%n_parameters
!
         bfgs%Hessian(k,k) = one
!
      enddo
!
   end function new_bfgs_tool
!
!
   subroutine get_direction_bfgs_tool(bfgs, g, d)
!!
!!    Get direction 
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Solves H d = - g and returns d 
!!
      implicit none 
!
      class(bfgs_tool), intent(in) :: bfgs 
!
      real(dp), dimension(bfgs%n_parameters), intent(in) :: g 
!
      real(dp), dimension(bfgs%n_parameters) :: d 
!
      integer, dimension(bfgs%n_parameters) :: ipiv 
      real(dp), dimension(bfgs%n_parameters, bfgs%n_parameters) :: H 
!
      integer :: info 
!
      H = bfgs%Hessian
      d = -g
!
      call dgesv(bfgs%n_parameters,    &
                  1,                   & ! Number of RHS 
                  H,                   &
                  bfgs%n_parameters,   &
                  ipiv,                &
                  d,                   &
                  bfgs%n_parameters,   &
                  info)
!
      if (info /= 0) call output%error_msg('Could not solve linear equation in BFGS.')
!
   end subroutine get_direction_bfgs_tool
!
!
   subroutine set_previous_geometry_and_gradient_bfgs_tool(bfgs, x, g)
!!
!!    Set previous geometry and gradient 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(bfgs_tool) :: bfgs 
!
      real(dp), dimension(bfgs%n_parameters), intent(in) :: x, g 
!
      bfgs%prev_x = x 
      bfgs%prev_g = g 
!
   end subroutine set_previous_geometry_and_gradient_bfgs_tool
!
!
   subroutine update_hessian_bfgs_tool(bfgs, x, g)
!!
!!    Update Hessian 
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    g:  current gradient 
!!    x:  current geometry on entry, next geometry on exit 
!!
!!    Direction for line-search in iteration k: 
!!
!!         d_k = - H_k^-1 g_k     
!!
!!    The Hessian is updated according to 
!!    
!!       H_k+1 = H_k - (H_k s_k)(H_k s_k)^T / (s_k^T H_k s_k) + (y_k y_k^T) / y_k^T s_k 
!!             = H_k - (z_k)(z_k)^T / (s_k^T z_k) + (y_k y_k^T) / y_k^T s_k
!!
!!    where we have let 
!!
!!       s_k = x_k - x_k-1
!!       y_k = g_k - g_k-1 
!!       z_k = H_k s_k 
!!
      implicit none 
!
      class(bfgs_tool), intent(inout) :: bfgs 
!
      real(dp), dimension(bfgs%n_parameters) :: x 
      real(dp), dimension(bfgs%n_parameters) :: g  
!
      integer, dimension(bfgs%n_parameters) :: ipiv 
!
      real(dp), dimension(bfgs%n_parameters) :: s, y, z
      real(dp), dimension(bfgs%n_parameters, bfgs%n_parameters) :: H, yyT, zzT 
!
      real(dp) :: yTs, zTs, sTy, ddot
!
      integer :: info 
!
      bfgs%iteration = bfgs%iteration + 1
!
      if (bfgs%iteration == 1) return 
!
!     Compute s_k and y_k from the k+1- and k-th geometries and gradients,
!     and then compute z_k = H_k s_k 
!
      s = x - bfgs%prev_x 
      y = g - bfgs%prev_g 
!
      sTy = ddot(bfgs%n_parameters, s, 1, y, 1)
      if (sTy < 0) y = -y
!
      call dgemm('N', 'N',             &
                  bfgs%n_parameters,   &
                  1,                   &
                  bfgs%n_parameters,   &
                  one,                 &
                  bfgs%Hessian,        &
                  bfgs%n_parameters,   &
                  s,                   &
                  bfgs%n_parameters,   &
                  zero,                &
                  z,                   &
                  bfgs%n_parameters)
!     
!     Compute the outer products yyT and zzT
!
      yyT = zero 
      call dger(bfgs%n_parameters,     &
                  bfgs%n_parameters,   &
                  one,                 &
                  y,                   &
                  1,                   &
                  y,                   &
                  1,                   &
                  yyT,                 &
                  bfgs%n_parameters)
!
      zzT = zero 
      call dger(bfgs%n_parameters,     &
                  bfgs%n_parameters,   &
                  one,                 &
                  z,                   &
                  1,                   &
                  z,                   &
                  1,                   &
                  zzT,                 &
                  bfgs%n_parameters)
!
!     Compute the inner products yTs and zTs 
!
      yTs = ddot(bfgs%n_parameters, y, 1, s, 1)
      zTs = ddot(bfgs%n_parameters, z, 1, s, 1)
!
      call daxpy(bfgs%n_parameters**2, -one/zTs, zzT, 1, bfgs%Hessian, 1)
      call daxpy(bfgs%n_parameters**2, one/yTs, yyT, 1, bfgs%Hessian, 1)
!
      bfgs%prev_x = x
      bfgs%prev_g = g
!
   end subroutine update_hessian_bfgs_tool
!
!
end module bfgs_tool_class
