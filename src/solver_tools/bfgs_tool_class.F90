!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   use parameters
   use global_out, only: output
   use memory_manager_class, only: mem
   use array_utilities, only: zero_array, copy_and_scale, get_l2_norm
!
!
   type :: bfgs_tool
!
      private
      integer  :: iteration    ! Starts at 0, increments when calling "update"
      integer  :: n_parameters ! Length of gradient
      real(dp) :: max_step     ! Maximum acceptable step length (in 2-norm)
!
      real(dp), dimension(:), allocatable, private :: previous_gradient
      real(dp), dimension(:,:), allocatable, private :: hessian
!
   contains
!
      procedure, public :: get_step                => get_step_bfgs_tool
      procedure, public :: update_hessian          => update_hessian_bfgs_tool
      procedure, public :: initialize_arrays       => initialize_arrays_bfgs_tool
!
      procedure, public :: get_hessian &
                        => get_hessian_bfgs_tool
!
      procedure, public :: set_hessian &
                        => set_hessian_bfgs_tool
!
      procedure, public :: set_initial_hessian_diagonal &
                        => set_initial_hessian_diagonal_bfgs_tool
!
      final :: destructor_bfgs_tool
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
      implicit none
!
      type(bfgs_tool) :: bfgs
!
      integer, intent(in)  :: n_parameters
      real(dp), intent(in) :: max_step
!
      bfgs%n_parameters = n_parameters
      bfgs%max_step     = max_step
!
      bfgs%iteration    = 0
!
   end function new_bfgs_tool
!
!
   subroutine initialize_arrays_bfgs_tool(bfgs)
!!
!!    Initialize arrays
!!    Written by Eirik F. Kjønstad, Jan 2020
!!
!!    Allocates and initializes the previous g and x (gradient, parameters),
!!    as well as the hessian.
!!
      implicit none
!
      class(bfgs_tool) :: bfgs
!
      integer :: k
!
      call mem%alloc(bfgs%previous_gradient, bfgs%n_parameters)
!
      call zero_array(bfgs%previous_gradient, bfgs%n_parameters)
!
      call mem%alloc(bfgs%hessian, bfgs%n_parameters, bfgs%n_parameters)
!
      call zero_array(bfgs%hessian, bfgs%n_parameters**2)
!
!$omp parallel do private(k)
      do k = 1, bfgs%n_parameters
!
         bfgs%hessian(k,k) = one
!
      enddo
!$omp end parallel do
!
   end subroutine initialize_arrays_bfgs_tool
!
!
   subroutine get_hessian_bfgs_tool(bfgs, H)
!!
!!    Get hessian
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(bfgs_tool), intent(in) :: bfgs
!
      real(dp), dimension(bfgs%n_parameters, bfgs%n_parameters), intent(out) :: H
!
      call dcopy(bfgs%n_parameters**2, bfgs%hessian, 1, H, 1)
!
   end subroutine get_hessian_bfgs_tool
!
!
   subroutine set_hessian_bfgs_tool(bfgs, H)
!!
!!    Set hessian
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(bfgs_tool), intent(inout) :: bfgs
!
      real(dp), dimension(bfgs%n_parameters, bfgs%n_parameters), intent(in) :: H
!
      call dcopy(bfgs%n_parameters**2, H, 1, bfgs%hessian, 1)
!
   end subroutine set_hessian_bfgs_tool
!
!
   subroutine destructor_bfgs_tool(bfgs)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(bfgs_tool) :: bfgs
!
      if (allocated(bfgs%hessian)) call mem%dealloc(bfgs%hessian, bfgs%n_parameters, bfgs%n_parameters)
      if (allocated(bfgs%previous_gradient)) call mem%dealloc(bfgs%previous_gradient, bfgs%n_parameters)
!
   end subroutine destructor_bfgs_tool
!
!
   subroutine get_step_bfgs_tool(bfgs, g, d)
!!
!!    Get step
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Solves H d = - g and returns d.
!!
!!    H is level shifted using the augmented rational function
!!    (RF) approach, where the level shift is given by the lowest
!!    eigenvalue of the augmented hessian (H, g; g^T 0).
!!
      implicit none
!
      class(bfgs_tool), intent(in) :: bfgs
!
      real(dp), dimension(bfgs%n_parameters), intent(in) :: g
!
      real(dp), dimension(bfgs%n_parameters) :: d
!
      real(dp), dimension(:,:), allocatable    :: aug_H
      real(dp), dimension(bfgs%n_parameters+1) :: eigvals
!
      integer :: info
      real(dp), dimension(4*(bfgs%n_parameters+1)) :: work
!
      real(dp) :: norm_d
!
      integer :: i, j
!
!     Set up rational function (RF) augmented hessian
!
      call mem%alloc(aug_H, bfgs%n_parameters + 1, bfgs%n_parameters + 1)
!
!$omp parallel do private(i,j)
      do i = 1, bfgs%n_parameters
         do j = 1, bfgs%n_parameters
!
            aug_H(i,j) = bfgs%hessian(i,j)
!
         end do
      end do
!$omp end parallel do
!
!$omp parallel do private(j)
      do j = 1, bfgs%n_parameters
!
         aug_H(j, bfgs%n_parameters + 1) = g(j)
         aug_H(bfgs%n_parameters + 1, j) = g(j)
!
      end do
!$omp end parallel do
!
      aug_H(bfgs%n_parameters + 1, bfgs%n_parameters + 1) = zero
!
!     Get lowest eigenvalue (level shift) and eigenvector (d, step)
!     of the augmented hessian
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
      if (info .ne. 0) then
         call output%error_msg('Could not solve eigenvalue equation in BFGS tool.' // &
                              ' "Dsyev" finished with info: (i0)', ints=[info])
      end if
!
      call output%printf('n', 'Rational function level shift: (f19.12)', &
                           reals=[eigvals(1)], fs='(/t6,a)')
!
!     Set d equal to the (scaled) first eigenvector
!
      call copy_and_scale(one/aug_H(bfgs%n_parameters+1,1), &
                          aug_H(1:bfgs%n_parameters,1),     &
                          d,                                &
                          bfgs%n_parameters)
!
      call mem%dealloc(aug_H, bfgs%n_parameters + 1, bfgs%n_parameters + 1)
!
!     Scale the vector to the boundary of the trust region (max step)
!     if the d vector is too long
!
      norm_d = get_l2_norm(d, bfgs%n_parameters)
      if (norm_d > bfgs%max_step) then
!
         call output%printf('n', 'BFGS step exceeds max. Scaling down to max_step.', fs='(/t6,a)')
         call dscal(bfgs%n_parameters, bfgs%max_step/norm_d, d, 1)
!
      endif
!
   end subroutine get_step_bfgs_tool
!
!
   subroutine update_hessian_bfgs_tool(bfgs, dx, g)
!!
!!    Update hessian
!!    Written by Eirik F. Kjønstad, 2019
!!
!!       x:  current geometry
!!       g:  current gradient
!!
!!    The hessian is updated according to
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
!!    For k = 1, the routine does not update the hessian but keeps a copy
!!    of the geometry and gradient for the next iteration. For k > 1, the
!!    routine also updates the hessian.
!!
      implicit none
!
      class(bfgs_tool), intent(inout) :: bfgs
!
      real(dp), dimension(bfgs%n_parameters) :: dx
      real(dp), dimension(bfgs%n_parameters) :: g
!
      real(dp), dimension(bfgs%n_parameters) :: s, y, z
      real(dp), dimension(:,:), allocatable  :: yyT, zzT
!
      real(dp) :: yTs, zTs, ddot
!
      bfgs%iteration = bfgs%iteration + 1
!
!     Compute s_k and y_k from the k+1- and k-th geometries and gradients
!
      call dcopy(bfgs%n_parameters, dx, 1, s, 1)
!
      call dcopy(bfgs%n_parameters, g, 1, y, 1)
      call daxpy(bfgs%n_parameters, -one, bfgs%previous_gradient, 1, y, 1)
!
!     Keep a copy, for next time, of the current geometry and gradient
!
      call dcopy(bfgs%n_parameters, g, 1, bfgs%previous_gradient, 1)
!
      if (bfgs%iteration == 1) then
!
         call output%printf('n', 'First iteration: no update of the hessian', fs='(/t6,a)')
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
                  bfgs%hessian,        &
                  bfgs%n_parameters,   &
                  s,                   &
                  bfgs%n_parameters,   &
                  zero,                &
                  z,                   &
                  bfgs%n_parameters)
!
!     Compute the outer products yyT and zzT
!
      call mem%alloc(yyT, bfgs%n_parameters, bfgs%n_parameters)
      call mem%alloc(zzT, bfgs%n_parameters, bfgs%n_parameters)
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
!     Update the hessian
!
      call daxpy(bfgs%n_parameters**2, -one/zTs, zzT, 1, bfgs%hessian, 1)
      call daxpy(bfgs%n_parameters**2, one/yTs, yyT, 1, bfgs%hessian, 1)
!
      call mem%dealloc(yyT, bfgs%n_parameters, bfgs%n_parameters)
      call mem%dealloc(zzT, bfgs%n_parameters, bfgs%n_parameters)
!
   end subroutine update_hessian_bfgs_tool
!
!
   subroutine set_initial_hessian_diagonal_bfgs_tool(bfgs, d)
!!
!!    Set initial hessian diagonal
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(bfgs_tool) :: bfgs
!
      real(dp), dimension(bfgs%n_parameters), intent(in) :: d
!
      integer :: k
!
      call zero_array(bfgs%hessian, bfgs%n_parameters**2)
!
!$omp parallel do private(k)
      do k = 1, bfgs%n_parameters
!
         bfgs%hessian(k,k) = d(k)
!
      enddo
!$omp end parallel do
!
   end subroutine set_initial_hessian_diagonal_bfgs_tool
!
!
end module bfgs_tool_class
