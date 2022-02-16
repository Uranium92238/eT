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
module ri_task_class
!
!!
!! RI task class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use parameters
   use cc_task_class, only: cc_task
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(cc_task) :: ri_task
!
   contains
!
      procedure, public :: execute &
                        => execute_ri_task
!
   end type ri_task
!
!
   interface ri_task
!
      procedure :: new_ri_task
!
   end interface ri_task
!
!
contains
!
!
   function new_ri_task() result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(ri_task) :: this
!
      this%name_ = 'Constructing Cholesky vectors for the RI-ERI approximation'
!
   end function new_ri_task
!
!
   subroutine execute_ri_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad, Feb 2021
!!
      use eri_ri_class, only: eri_ri
      use global_in, only: input
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(ri_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      type(eri_ri), allocatable :: ri_solver
!
      character(len=200) :: ri_basis_set
!
      call this%print_header()
      call this%start_timer()
!
      call input%get_keyword('ri', 'integrals', ri_basis_set)
!
      ri_solver = eri_ri(ri_basis_set, wf%ao%get_libint_epsilon())
      call ri_solver%initialize()
!
      call ri_solver%run(wf%ao)
!
      call wf%integral_preparations(ri_solver%get_n_J())
!
      call ri_solver%construct_cholesky_mo_vectors(wf%ao, wf%ao%n, wf%n_mo, &
                                                   wf%orbital_coefficients, wf%L_mo)
!
      call ri_solver%cleanup()
!
      call wf%L_t1%set_equal_to(wf%L_mo)
      call wf%mo_preparations()
!
      call this%end_timer()
!
   end subroutine execute_ri_task
!
!
end module ri_task_class
