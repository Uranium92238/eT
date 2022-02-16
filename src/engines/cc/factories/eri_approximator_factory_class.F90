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
module eri_approximator_factory_class
!
!!
!! ERI approximator factory class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use cc_task_class, only: cc_task
   use ri_task_class, only: ri_task
   use cholesky_decomposition_task_class, only: cholesky_decomposition_task
!
   implicit none
!
!
   type :: eri_approximator_factory
!
      logical :: ri
!
   contains
!
      procedure, public :: create
!
   end type eri_approximator_factory
!
   interface  eri_approximator_factory
!
      procedure :: new_eri_approximator_factory
!
   end interface  eri_approximator_factory
!
contains
!
!
   function new_eri_approximator_factory() result(this)
!!
!!    New eri approximator
!!    Written by Alexander C. Paul, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      type(eri_approximator_factory) :: this
!
      this%ri = input%is_keyword_present('ri', 'integrals')
!
   end function new_eri_approximator_factory
!
!
   function create(this) result(eri_task)
!!
!!    Create
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      class(eri_approximator_factory), intent(in) :: this
!
      class(cc_task), allocatable :: eri_task
!
      if (this%ri) then
!
         eri_task = ri_task()
!
      else
!
         eri_task = cholesky_decomposition_task()
!
      end if
!
   end function create
!
!
end module eri_approximator_factory_class
