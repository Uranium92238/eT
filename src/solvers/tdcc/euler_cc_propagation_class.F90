!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module euler_cc_propagation_class
!
!!
!! euler CC propagation solver class
!! Written by Andreas Skeidsvoll, Nov 2018
!!
!! Makes an Euler solver by letting the propagation solver do Euler steps.
!
   use parameters
   use memory_manager_class, only: mem
   use ccs_class, only: ccs
   use electric_field_class, only: electric_field
   use cc_propagation_class, only: cc_propagation
!
   implicit none
!
   type, extends(cc_propagation) :: euler_cc_propagation
!
   contains
!
      procedure :: step => euler_step
!
   end type euler_cc_propagation
!
!
   interface euler_cc_propagation
!
      procedure :: new_euler_cc_propagation
!
   end interface euler_cc_propagation
!
!
contains
!
!
   function new_euler_cc_propagation(wf) result(solver)
!!
!!    New Euler CC propagation
!!    Written by Andreas Skeidsvoll, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
      type(euler_cc_propagation) :: solver
!
      solver%tag = 'euler'
      solver%name_ = 'Euler coupled cluster propagation solver'
!
      call solver%new_cc_propagation(wf)
!
   end function new_euler_cc_propagation
!
!
   subroutine euler_step(solver, wf, field, ti, dt, ui, uf, n)
!!
!!    Euler step
!!    Written by Andreas Skeidsvoll, Nov 2018
!!
!!    Takes one forward Euler step
!!
      implicit none
!
      class(euler_cc_propagation), intent(inout) :: solver
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
      complex(dp), dimension(:), allocatable :: ddt_ui
!
      call solver%update_field_and_wavefunction(wf, field, ti, ui)
!
      call mem%alloc(ddt_ui, n)
      call wf%construct_complex_time_derivative(ddt_ui)
!
!     Use the derivative of the single stage to estimate the solution uf at t = ti + dt.
!
!     uf = ui + dt*ddt_ui
!
      call zcopy(n, ddt_ui, 1, uf, 1)
      call mem%dealloc(ddt_ui, n)
!
      call zdscal(n, dt, uf, 1)
      call zaxpy(n, one_complex, ui, 1, uf, 1)
!
   end subroutine euler_step
!
!
end module euler_cc_propagation_class