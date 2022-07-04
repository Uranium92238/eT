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
module hf_geoopt_task_class
!
!!
!! HF geometry optimization task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use hf_task_class,                  only: hf_task
   use hf_class,                       only: hf
   use bfgs_solver_class,              only: bfgs_solver
   use hf_geoopt_solver_factory_class, only: hf_geoopt_solver_factory
!
   implicit none
!
   type, extends(hf_task) :: hf_geoopt_task
!
      class(bfgs_solver), allocatable, private :: solver
      type(hf_geoopt_solver_factory), private  :: solver_factory
!
   contains
!
      procedure, public :: execute &
                        => execute_hf_geoopt_task
!
      procedure, private, nopass :: sanity_checks
!
   end type hf_geoopt_task
!
!
   interface hf_geoopt_task
!
      procedure :: new_hf_geoopt_task
!
   end interface hf_geoopt_task
!
!
contains
!
!
   function new_hf_geoopt_task() result(this)
!!
!!    New HF amplitudes task
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(hf_geoopt_task) :: this
!
      this%name_ = 'Determining optimal HF geometry'
!
   end function new_hf_geoopt_task
!
!
   subroutine execute_hf_geoopt_task(this, wf)
!!
!!    Execute
!!    Written by Alexander C. Paul, May 2022
!!
      use redundant_internal_coords_class, only: redundant_internal_coords
!
      implicit none
!
      class(hf_geoopt_task), intent(inout) :: this
      class(hf), target, intent(inout) :: wf
!
      type(redundant_internal_coords), allocatable :: internals
!
      call this%print_header()
      call this%start_timer()
!
      call this%sanity_checks(wf)
!
      internals = redundant_internal_coords(3*wf%n_atomic_centers)
      call wf%initialize_redundant_internal_coordinates(internals)
!
      this%solver_factory = hf_geoopt_solver_factory()
      call this%solver_factory%create(wf, internals, this%solver)
!
      call this%solver%initialize()
      call this%solver%run()
      call this%solver%finalize()
!
      call this%end_timer()
!
   end subroutine execute_hf_geoopt_task
!
!
   subroutine sanity_checks(wf)
!!
!!    Sanity checks
!!    Written by Alexander C. Paul, May 2022
!!
      use global_out, only: output
!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      if (wf%ao%has_ghost_atoms()) &
         call output%warning_msg("Ghosts are experimental in geometry optimization.")
!
      if (wf%embedded) &
         call output%error_msg('geometry optimization with embedding is not supported')
!
   end subroutine sanity_checks
!
!
end module hf_geoopt_task_class
