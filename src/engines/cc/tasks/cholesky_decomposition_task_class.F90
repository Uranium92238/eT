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
module cholesky_decomposition_task_class
!
!!
!! Cholesky decomposition task class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use parameters
   use cc_task_class, only: cc_task
   use ccs_class, only: ccs
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(cc_task) :: cholesky_decomposition_task
!
   contains
!
      procedure, public :: execute &
                        => execute_cholesky_decomposition_task
!
   end type cholesky_decomposition_task
!
!
   interface cholesky_decomposition_task
!
      procedure :: new_cholesky_decomposition_task
!
   end interface cholesky_decomposition_task
!
!
contains
!
!
   function new_cholesky_decomposition_task() result(this)
!!
!!    New
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      type(cholesky_decomposition_task) :: this
!
      this%name_ = 'Cholesky-decomposing electron repulsion integrals'
!
   end function new_cholesky_decomposition_task
!
!
   subroutine execute_cholesky_decomposition_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Feb 2020
!!
!!    Performs Cholesky decomposition of the electron repulsion integral matrix.
!!
!!    For reduced space coupled cluster calculations (frozen HF or MLHF)
!!    the MO screening can be used to reduce the number of Cholesky vectors,
!!    as the accuracy of the integrals in the active MO basis,
!!    rather than the AO basis, is targeted. For details, see
!!    Folkestad, S. D., Kjønstad, E. F., and Koch, H., JCP, 150(19), 194112 (2019).
!!
      use eri_cd_class, only: eri_cd
      use global_in, only: input
      use global_out, only: output
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(cholesky_decomposition_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      type(eri_cd), allocatable :: eri_cholesky_solver
!
      logical :: do_MO_screening
!
      real(dp), dimension(:,:), allocatable, target :: screening_vector
!
      call this%print_header()
      call this%start_timer()
!
      do_MO_screening = input%is_keyword_present('mo screening', 'solver cholesky')
!
      eri_cholesky_solver = eri_cd(wf%ao)
!
      if (do_MO_screening) then
!
         call output%printf('m', 'Using the MO screening for the Cholesky decomposition', &
                        fs='(/t3,a)')
!
         call mem%alloc(screening_vector, wf%ao%n, wf%ao%n)
         call wf%construct_MO_screening_for_cd(screening_vector)
!
         call eri_cholesky_solver%run(wf%ao, screening_vector)   ! Do the Cholesky decomposition
                                                                 ! using the MO screening
!
         call mem%dealloc(screening_vector, wf%ao%n, wf%ao%n)
!
      else
!
         call eri_cholesky_solver%run(wf%ao) ! Do the Cholesky decomposition
!
      endif
!
      call eri_cholesky_solver%diagonal_test(wf%ao)  ! Determine the largest
                                                     ! deviation in the ERI matrix
!
      call wf%integral_preparations(eri_cholesky_solver%get_n_cholesky())
!
      call eri_cholesky_solver%construct_cholesky_mo_vectors(wf%ao, wf%ao%n, wf%n_mo, &
                                                   wf%orbital_coefficients, wf%L_mo)
!
      call eri_cholesky_solver%cleanup()
      call wf%L_t1%set_equal_to(wf%L_mo)
      call wf%mo_preparations()
!
      call this%end_timer()
!
   end subroutine execute_cholesky_decomposition_task
!
!
end module cholesky_decomposition_task_class
