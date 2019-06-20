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
!!    Usage:
!!
!!    bfgs = bfgs_tool(dim_g, max_step_length)
!!
!!    if (.not. converged) then 
!! 
!!       call bfgs%udpate_hessian(x, g)      
!!       call bfgs%get_step(g, d)
!!       x = x + d 
!!
!!    endif      
!!
!!    Notation:
!!
!!       x : parameters (geometry)
!!       g : gradient 
!!       d : BFGS step (i.e., solution d of H d = -g, where H has been level-shifted)
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
      integer  :: iteration    ! Starts at 0, increments when calling "update"
      integer  :: n_parameters ! Length of gradient 
      real(dp) :: max_step     ! Maximum acceptable step length (in 2-norm)
!
      real(dp), dimension(:), allocatable :: prev_g ! Prev. gradient 
      real(dp), dimension(:), allocatable :: prev_x ! Prev. geometry ('parameters', more generally)
!
      real(dp), dimension(:,:), allocatable :: Hessian ! BFGS Hessian estimate
!
   contains
!
      procedure, public :: get_step                => get_step_bfgs_tool
      procedure, public :: update_hessian          => update_hessian_bfgs_tool
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
   function new_bfgs_tool(n_parameters, max_step) result(bfgs)
!!
!!    New BFGS tool
!!    Written by Eirik F. Kjønstad, 2019
!!
      type(bfgs_tool) :: bfgs
!
      integer, intent(in)  :: n_parameters
      real(dp), intent(in) :: max_step
!  
      integer :: k 
!
      bfgs%n_parameters = n_parameters
      bfgs%max_step     = max_step
!
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
   subroutine get_step_bfgs_tool(bfgs, g, d)
!!
!!    Get step 
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Solves H d = - g and returns d. 
!!
!!    H is level shifted using the augmented RF approach, where  
!!    the level shift is given by the lowest eigenvalue of the 
!!    augmented Hessian (H, g; g^T 0).
!!
      implicit none 
!
      class(bfgs_tool), intent(in) :: bfgs 
!
      real(dp), dimension(bfgs%n_parameters), intent(in) :: g 
!
      real(dp), dimension(bfgs%n_parameters) :: d 
!
      real(dp), dimension(bfgs%n_parameters + 1, bfgs%n_parameters + 1) :: aug_H 
      real(dp), dimension(bfgs%n_parameters+1) :: eigvals 
!
      integer :: info 
      real(dp), dimension(4*(bfgs%n_parameters+1)) :: work 
!
      real(dp) :: norm_d
!
!     Set up rational function (RF) augmented Hessian 
!
      aug_H = zero 
      aug_H(1:bfgs%n_parameters, 1:bfgs%n_parameters) = bfgs%Hessian(:,:)
      aug_H(1:bfgs%n_parameters, bfgs%n_parameters + 1) = g(:)
      aug_H(bfgs%n_parameters + 1, 1:bfgs%n_parameters) = g(:)
!
!     Get lowest eigenvalue (level shift) and eigenvector (d, step)
!     of the augmented Hessian 
!
      call dsyev('V', 'L',                   &
                  bfgs%n_parameters+1,       &
                  aug_H,                     & ! aug_H on entry, eigenvectors on exit
                  bfgs%n_parameters+1,       &
                  eigvals,                   &
                  work,                      &
                  4*(bfgs%n_parameters+1),   &
                  info)
!
      if (info .ne. 0) call output%error_msg('Could not solve eigenvalue equation in BFGS tool.')
!
      call output%printf('Level shift: (f19.12)', reals=[eigvals(1)], fs='(/t3,a)')
!
      d = aug_H(1:bfgs%n_parameters,1)/aug_H(bfgs%n_parameters+1,1)
!
!     Scale the vector to the boundary of the trust region (max step)
!     if the d vector is too long 
!
      norm_d = get_l2_norm(d, bfgs%n_parameters)
      if (norm_d > bfgs%max_step) then
!
         call output%printf('BFGS step exceeds max. Scaling down to max_step.')
         d = d*(bfgs%max_step/norm_d)
!
      endif
!
   end subroutine get_step_bfgs_tool
!
!
   subroutine update_hessian_bfgs_tool(bfgs, x, g)
!!
!!    Update Hessian 
!!    Written by Eirik F. Kjønstad, 2019
!!
!!       x:  current geometry
!!       g:  current gradient 
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
!!    For k = 1, the routine does not update the Hessian but keeps a copy 
!!    of the geometry and gradient for the next iteration. For k > 1, the 
!!    routine also updates the Hessian.
!!
      implicit none 
!
      class(bfgs_tool), intent(inout) :: bfgs 
!
      real(dp), dimension(bfgs%n_parameters) :: x 
      real(dp), dimension(bfgs%n_parameters) :: g  
!
      real(dp), dimension(bfgs%n_parameters) :: s, y, z
      real(dp), dimension(bfgs%n_parameters, bfgs%n_parameters) :: yyT, zzT
!
      real(dp) :: yTs, zTs, ddot
!
      bfgs%iteration = bfgs%iteration + 1
!
!     Compute s_k and y_k from the k+1- and k-th geometries and gradients
!
      s = x - bfgs%prev_x 
      y = g - bfgs%prev_g 
!
!     Keep a copy, for next time, of the current geometry and gradient 
!
      bfgs%prev_x = x
      bfgs%prev_g = g
!
      if (bfgs%iteration == 1) then 
!
         call output%printf('First iteration: no update of the Hessian', fs='(/t3,a)')
         return
!
      endif 
!
!     Compute z_k = H_k s_k 
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
!     Update the Hessian
!
      call daxpy(bfgs%n_parameters**2, -one/zTs, zzT, 1, bfgs%Hessian, 1)
      call daxpy(bfgs%n_parameters**2, one/yTs, yyT, 1, bfgs%Hessian, 1)
!
   end subroutine update_hessian_bfgs_tool
!
!
end module bfgs_tool_class
