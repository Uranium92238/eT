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
module scf_solver_factory_class
!
!!
!!    SCF solver factory class
!!    Written by Sarai D. Folkestad, May 2021
!!
!

   use kinds
   use global_out,                   only: output
   use hf_class,                     only: hf
   use convergence_tool_class,       only: convergence_tool
   use scf_solver_class,             only: scf_solver
!
   implicit none
!
   type:: scf_solver_factory
!
      integer                 :: max_iterations
      type(convergence_tool)  :: convergence_checker
      character(len=200)      :: acceleration_type
!
   contains
!
      procedure, public :: create => create_scf_solver_factory
!
      procedure, private :: read_settings
!
   end type scf_solver_factory
!
   interface scf_solver_factory
!
      procedure :: new_scf_solver_factory
!
   end interface scf_solver_factory
!
contains
!
!
   function new_scf_solver_factory(max_iterations,       &
                                   energy_threshold,     &
                                   acceleration_type) result(this)
!!
!!    New scf solver factory
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      type(scf_solver_factory) :: this
!
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: energy_threshold
      character(len=*), intent(in), optional :: acceleration_type
!
      call this%read_settings()
!
      if (present(max_iterations)) this%max_iterations = max_iterations
!
      if (present(acceleration_type)) this%acceleration_type = acceleration_type
!
      if (present(energy_threshold)) &
         call this%convergence_checker%set_energy_threshold(energy_threshold)
!
   end function new_scf_solver_factory
!
!
   subroutine create_scf_solver_factory(this, wf, solver, restart, skip)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(scf_solver_factory),                  intent(inout)   :: this
      class(hf),                                  intent(inout)   :: wf
      class(scf_solver), allocatable,             intent(out)     :: solver
!
      logical, intent(in) :: restart, skip
!
      call this%convergence_checker%set_residual_threshold(wf%gradient_threshold)
!
      solver = scf_solver(restart             = restart,                      &
                          acceleration_type   = this%acceleration_type,       &
                          max_iterations      = this%max_iterations,          &
                          skip                = skip,                         &
                          dim_                = wf%n_mo,                      &
                          n_equations         = wf%n_densities,               &
                          gradient_dimension  = wf%packed_gradient_dimension, &
                          convergence_checker = this%convergence_checker)
!
!
   end subroutine create_scf_solver_factory
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad, May 2021
!!
!
      use global_in, only: input
!
      implicit none
!
      class(scf_solver_factory),   intent(inout)  :: this
!
      real(dp) :: energy_threshold
!
      character(len=200) :: algorithm
!
      this%convergence_checker = convergence_tool(1.0d-7, 1.0d-7, energy_convergence=.false.)
!
      if (input%is_keyword_present('energy threshold', 'solver scf')) then
         call input%get_keyword('energy threshold',  &
                                'solver scf',        &
                                energy_threshold)
!
         call this%convergence_checker%set_energy_threshold(energy_threshold)
!
      endif
!
      this%max_iterations = 100
      call input%get_keyword('max iterations',      &
                             'solver scf',          &
                             this%max_iterations)
!
   algorithm = 'scf-diis'
   call input%get_keyword('algorithm', 'solver scf', algorithm)
!
   if (trim(algorithm) .ne. 'scf-diis'    .and. &
       trim(algorithm) .ne. 'scf'         .and. &
       trim(algorithm) .ne. 'mo-scf-diis') then
!
      call output%error_msg('did not recognize SCF algorithm')
!
   endif
!
   this%acceleration_type = 'none'
   if (trim(algorithm) == 'scf-diis' .or. &
      trim(algorithm) == 'mo-scf-diis') this%acceleration_type = 'diis'
!
   end subroutine read_settings
!
end module scf_solver_factory_class
