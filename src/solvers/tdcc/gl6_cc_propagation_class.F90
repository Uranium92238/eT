!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module gl6_cc_propagation_class
!
!!
!! gl6 CC propagation solver class
!! Written by Andreas Skeidsvoll, Feb 2019
!!
!! Makes an Gauss-Legendre 6 solver by letting the propagation solver do Gauss-Legendre 6 steps.
!!
!! Initial guesses are used to accelerate the solution of equations taking part in implicit
!! integration methods. The guesses at given time step are generated by using the solution from the
!! previous time step. The method used is described in Section VIII.6.1 (A) in Hairer et. al.,
!! Geometric Numerical Integration.
!!
!
   use parameters
   use butcher_tables
   use continuous_output_coefficients
   use memory_manager_class, only: mem
   use ccs_class, only: ccs
   use electric_field_class, only: electric_field
   use cc_propagation_class, only: cc_propagation
!
   implicit none
!
   type, extends(cc_propagation) :: gl6_cc_propagation
!
!     Initial guesses used to accelerate the convergence of equations in implicit methods, see
!     Section VIII.6.1 (A) in Hairer et. al., Geometric Numerical Integration.
!
      complex(dp), dimension(:), allocatable :: z1_guess, z2_guess, z3_guess
!
   contains
!
      procedure :: step => gl6_step
!
      procedure :: Initializations => initializations_gl6_cc_propagation
      procedure :: cleanup => cleanup_gl6_cc_propagation
!
   end type gl6_cc_propagation
!
!
   interface gl6_cc_propagation
!
      procedure :: new_gl6_cc_propagation
!
   end interface gl6_cc_propagation
!
!
contains
!
!
   function new_gl6_cc_propagation(wf) result(solver)
!!
!!    New gl6 CC propagation
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf
      type(gl6_cc_propagation) :: solver
!
      solver%tag = 'gl6'
      solver%name_ = 'Gauss-Legendre 6 coupled cluster propagation solver'
!
      call solver%new_cc_propagation(wf)
!
   end function new_gl6_cc_propagation
!
!
   subroutine cleanup_gl6_cc_propagation(solver, wf)
!!
!!    Destructor
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Deallocates arrays.
!!
      implicit none
!
      class(gl6_cc_propagation), intent(inout) :: solver
      class(ccs), intent(inout) :: wf
!
      if (solver%energy_output             &
          .or. solver%dipole_moment_output &
          .or. solver%density_matrix_output) call wf%destruct_gs_density_complex()
!
      call mem%dealloc(solver%z1_guess, solver%vector_length)
      call mem%dealloc(solver%z2_guess, solver%vector_length)
      call mem%dealloc(solver%z3_guess, solver%vector_length)
!
   end subroutine cleanup_gl6_cc_propagation
!
!
   subroutine initializations_gl6_cc_propagation(solver)
!!
!!    Initializations
!!    Written by Andreas Skedsvoll, Sep 2019
!!
!!    Allocates z1, z2, z3 guess and sets them to zero.
!!
!!    Moved from constructor, Eirik F. Kjønstad, Jan 2020.
!!
      implicit none
!
      class(gl6_cc_propagation), intent(inout) :: solver
!
      call mem%alloc(solver%z1_guess, solver%vector_length)
      call mem%alloc(solver%z2_guess, solver%vector_length)
      call mem%alloc(solver%z3_guess, solver%vector_length)
!
      solver%z1_guess = zero_complex
      solver%z2_guess = zero_complex
      solver%z3_guess = zero_complex
!
   end subroutine initializations_gl6_cc_propagation
!
!
   subroutine gl6_step(solver, wf, field, ti, dt, ui, uf, n)
!!
!!    gl6 step
!!    Written by Andreas Skeidsvoll, Feb 2019
!!
!!    Takes one forward gl6 step
!!
      use global_out, only: output
!
      implicit none
!
      class(gl6_cc_propagation), intent(inout) :: solver
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
      complex(dp), dimension(:), allocatable :: ddt_u1, ddt_u2, ddt_u3, z1, z2, z3, z1_next, &
                                                z2_next, z3_next
!
      integer :: i
!
      logical :: converged
!
      call mem%alloc(ddt_u1, n)
      call mem%alloc(ddt_u2, n)
      call mem%alloc(ddt_u3, n)
      call mem%alloc(z1, n)
      call mem%alloc(z2, n)
      call mem%alloc(z3, n)
      call mem%alloc(z1_next, n)
      call mem%alloc(z2_next, n)
      call mem%alloc(z3_next, n)
!
!     Get guess of z vectors
!
      call zcopy(n, solver%z1_guess, 1, z1, 1)
      call zcopy(n, solver%z2_guess, 1, z2, 1)
      call zcopy(n, solver%z3_guess, 1, z3, 1)
!
      converged = .false.
!
!     Iterate until convergence
!
      i = 0
!
      do
!
!        Construct derivative of the first stage
!
         call solver%update_field_and_wavefunction(wf, field, ti+c1_gl6*dt, ui+z1)
         call wf%construct_complex_time_derivative(ddt_u1)
!
!        Construct derivative of the second stage
!
         call solver%update_field_and_wavefunction(wf, field, ti+c2_gl6*dt, ui+z2)
         call wf%construct_complex_time_derivative(ddt_u2)
!
!        Construct derivative of the third stage
!
         call solver%update_field_and_wavefunction(wf, field, ti+c3_gl6*dt, ui+z3)
         call wf%construct_complex_time_derivative(ddt_u3)
!
         if (converged) exit
         i = i + 1
!
!        Construct first z vector
!
!        z1_next = dt*a11_gl6*ddt_u1 + dt*a12_gl6*ddt_u2 + dt*a13_gl6*ddt_u3
!
         call zcopy(n, ddt_u1, 1, z1_next, 1)
         call zdscal(n, dt*a11_gl6, z1_next, 1)
         call zaxpy(n, cmplx(dt*a12_gl6, zero, dp), ddt_u2, 1, z1_next, 1)
         call zaxpy(n, cmplx(dt*a13_gl6, zero, dp), ddt_u3, 1, z1_next, 1)
!
!        Construct second z vector
!
!        z2_next = dt*a21_gl6*ddt_u1 + dt*a22_gl6*ddt_u2 + dt*a23_gl6*ddt_u3
!
         call zcopy(n, ddt_u1, 1, z2_next, 1)
         call zdscal(n, dt*a21_gl6, z2_next, 1)
         call zaxpy(n, cmplx(dt*a22_gl6, zero, dp), ddt_u2, 1, z2_next, 1)
         call zaxpy(n, cmplx(dt*a23_gl6, zero, dp), ddt_u3, 1, z2_next, 1)
!
!        Construct third z vector
!
!        z3_next = dt*a31_gl6*ddt_u1 + dt*a32_gl6*ddt_u2 + dt*a33_gl6*ddt_u3
!
         call zcopy(n, ddt_u1, 1, z3_next, 1)
         call zdscal(n, dt*a31_gl6, z3_next, 1)
         call zaxpy(n, cmplx(dt*a32_gl6, zero, dp), ddt_u2, 1, z3_next, 1)
         call zaxpy(n, cmplx(dt*a33_gl6, zero, dp), ddt_u3, 1, z3_next, 1)
!
!        Determine whether convergence is reached or not
!
         if ((norm2(real(z1_next-z1))  .le. solver%implicit_threshold) &
             .and. (norm2(aimag(z1_next-z1)) .le. solver%implicit_threshold) &
             .and. (norm2(real(z2_next-z2))  .le. solver%implicit_threshold) &
             .and. (norm2(aimag(z2_next-z2)) .le. solver%implicit_threshold) &
             .and. (norm2(real(z3_next-z3))  .le. solver%implicit_threshold) &
             .and. (norm2(aimag(z3_next-z3)) .le. solver%implicit_threshold)) then
!
            converged = .true.
!
         endif
!
!        Update z vectors
!
         call zcopy(n, z1_next, 1, z1, 1)
         call zcopy(n, z2_next, 1, z2, 1)
         call zcopy(n, z3_next, 1, z3, 1)
!
      enddo
!
      call output%printf('n', 'GL6 iterations: (i6)', ints=[i], fs='(t3,a)')
!
      call mem%dealloc(z1, n)
      call mem%dealloc(z2, n)
      call mem%dealloc(z3, n)
      call mem%dealloc(z1_next, n)
      call mem%dealloc(z2_next, n)
      call mem%dealloc(z3_next, n)
!
!     Use the derivative of the three stages to estimate the solution uf at time t = ti + dt.
!
!     uf = ui + dt*b1_gl6*ddt_u3 + dt*b2_gl6*dvec4
!
      call zcopy(n, ui, 1, uf,1)
      call zaxpy(n, cmplx(dt*b1_gl6, zero, dp), ddt_u1, 1, uf, 1)
      call zaxpy(n, cmplx(dt*b2_gl6, zero, dp), ddt_u2, 1, uf, 1)
      call zaxpy(n, cmplx(dt*b3_gl6, zero, dp), ddt_u3, 1, uf, 1)
!
!     Calculate guess of next z vectors
!
!     solver%z1_guess = ui - uf + dt*beta11_gl6*ddt_u1 + dt*beta12_gl6*ddt_u2
!     + dt*beta13_gl6*ddt_u3
!
      call zcopy(n, ui, 1, solver%z1_guess, 1)
      call zaxpy(n, -one_complex, uf, 1, solver%z1_guess, 1)
      call zaxpy(n, cmplx(dt*beta11_gl6, zero, dp), ddt_u1, 1, solver%z1_guess, 1)
      call zaxpy(n, cmplx(dt*beta12_gl6, zero, dp), ddt_u2, 1, solver%z1_guess, 1)
      call zaxpy(n, cmplx(dt*beta13_gl6, zero, dp), ddt_u3, 1, solver%z1_guess, 1)
!
!     solver%z2_guess = ui - uf + dt*beta21_gl6*ddt_u1 + dt*beta22_gl6*ddt_u2
!     + dt*beta23_gl6*ddt_u3
!
      call zcopy(n, ui, 1, solver%z2_guess, 1)
      call zaxpy(n, -one_complex, uf, 1, solver%z2_guess, 1)
      call zaxpy(n, cmplx(dt*beta21_gl6, zero, dp), ddt_u1, 1, solver%z2_guess, 1)
      call zaxpy(n, cmplx(dt*beta22_gl6, zero, dp), ddt_u2, 1, solver%z2_guess, 1)
      call zaxpy(n, cmplx(dt*beta23_gl6, zero, dp), ddt_u3, 1, solver%z2_guess, 1)
!
!     solver%z3_guess = ui - uf + dt*beta31_gl6*ddt_u1 + dt*beta32_gl6*ddt_u2
!     + dt*beta33_gl6*ddt_u3
!
      call zcopy(n, ui, 1, solver%z3_guess, 1)
      call zaxpy(n, -one_complex, uf, 1, solver%z3_guess, 1)
      call zaxpy(n, cmplx(dt*beta31_gl6, zero, dp), ddt_u1, 1, solver%z3_guess, 1)
      call zaxpy(n, cmplx(dt*beta32_gl6, zero, dp), ddt_u2, 1, solver%z3_guess, 1)
      call zaxpy(n, cmplx(dt*beta33_gl6, zero, dp), ddt_u3, 1, solver%z3_guess, 1)
!
      call mem%dealloc(ddt_u1, n)
      call mem%dealloc(ddt_u2, n)
      call mem%dealloc(ddt_u3, n)
!
   end subroutine gl6_step
!
!
end module gl6_cc_propagation_class
