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
module rk4_cc_propagation_class
!
!!
!! rk4 CC propagation solver class
!! Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!! Makes a Runge-Kutta 4 solver by letting the propagation solver do Runge-Kutta 4 steps.
!
   use parameters
   use memory_manager_class, only: mem
   use ccs_class, only: ccs
   use electric_field_class, only: electric_field
   use cc_propagation_class, only: cc_propagation
!
   implicit none
!
   type, extends(cc_propagation) :: rk4_cc_propagation
!
   contains
!
      procedure :: step => rk4_step
!
   end type rk4_cc_propagation
!
!
   interface rk4_cc_propagation
!
      procedure :: new_rk4_cc_propagation
!
   end interface rk4_cc_propagation
!
!
contains
!
!
   function new_rk4_cc_propagation(wf) result(solver)
!!
!!    New rk4 CC propagation
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf
      type(rk4_cc_propagation) :: solver
!
      solver%tag = 'rk4'
      solver%name_ = 'Runge-Kutta 4 coupled cluster propagation solver'
!
      call solver%new_cc_propagation(wf)
!
   end function new_rk4_cc_propagation
!
!
   subroutine rk4_step(solver, wf, field, ti, dt, ui, uf, n)
!!
!!    rk4 step
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Takes one forward rk4 step
!!
      implicit none
!
      class(rk4_cc_propagation), intent(inout) :: solver
      class(ccs), intent(inout) :: wf
      class(electric_field), intent(inout) :: field
!
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: ti
!
      integer, intent(in) :: n
!
      complex(dp), dimension(n), intent(in) :: ui
      complex(dp), dimension(n), intent(out) :: uf
!
      complex(dp), dimension(:), allocatable :: u1, u2, u3, ddt_ui, ddt_u1, ddt_u2, ddt_u3
!
      real(dp) :: t1, t2, t3
!
!     Construct derivative of the first stage
!
      call solver%update_field_and_wavefunction(wf, field, ti, ui)
!
      call mem%alloc(ddt_ui, n)
      call wf%construct_complex_time_derivative(ddt_ui)
!
!     Construct derivative of the second stage
!
      t1 = ti + dt*half
!
!     u1 = ui + dt*half*ddt_ui
!
      call mem%alloc(u1, n)
      call zcopy(n, ui, 1, u1, 1)
      call zaxpy(n, cmplx(dt*half, zero, dp), ddt_ui, 1, u1, 1)
!
      call solver%update_field_and_wavefunction(wf, field, t1, u1)
      call mem%dealloc(u1, n)
!
      call mem%alloc(ddt_u1, n)
      call wf%construct_complex_time_derivative(ddt_u1)
!
!     Construct derivative of the third stage
!
      t2 = ti + dt*half
!
!     u2 = ui + dt*half*ddt_u1
!
      call mem%alloc(u2, n)
      call zcopy(n, ui, 1, u2, 1)
      call zaxpy(n, cmplx(dt*half, zero, dp), ddt_u1, 1, u2, 1)
!
      call solver%update_field_and_wavefunction(wf, field, t2, u2)
      call mem%dealloc(u2, n)
!
      call mem%alloc(ddt_u2, n)
      call wf%construct_complex_time_derivative(ddt_u2)
!
!
!     Construct derivative of the fourth stage
!
      t3 = ti + dt
!
!     u3 = ui + dt*ddt_u2
!
      call mem%alloc(u3, n)
      call zcopy(n, ui, 1, u3, 1)
      call zaxpy(n, cmplx(dt, zero, dp), ddt_u2, 1, u3, 1)
!
      call solver%update_field_and_wavefunction(wf, field, t3, u3)
      call mem%dealloc(u3, n)
!
      call mem%alloc(ddt_u3, n)
      call wf%construct_complex_time_derivative(ddt_u3)
!
!     Use the derivative of the four stages to estimate the solution uf at time t = ti + dt.
!
!     uf = ui + dt/six*(ddt_ui + two*ddt_u1 + two*ddt_u2 + ddt_u3)
!
      call zcopy(n, ddt_ui, 1, uf, 1)
      call zaxpy(n, two_complex, ddt_u1, 1, uf, 1)
      call zaxpy(n, two_complex, ddt_u2, 1, uf, 1)
      call zaxpy(n, one_complex, ddt_u3, 1, uf, 1)
!
      call mem%dealloc(ddt_ui, n)
      call mem%dealloc(ddt_u1, n)
      call mem%dealloc(ddt_u2, n)
      call mem%dealloc(ddt_u3, n)
!
      call zdscal(n, dt/six, uf, 1)
      call zaxpy(n, one_complex, ui, 1, uf, 1)
!
   end subroutine rk4_step
!
!
end module rk4_cc_propagation_class