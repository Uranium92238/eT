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
module hf_geoopt_solver_factory_class
!
!!
!! HF geoopt solver factory class
!! Written by Alexander C. Paul, May 2022
!!
!
   use parameters
   use hf_class,          only: hf
   use bfgs_solver_class, only: bfgs_solver
!
   implicit none
!
   type:: hf_geoopt_solver_factory
!
      integer :: max_iterations
      real(dp) :: max_step, energy_threshold, gradient_threshold
!
   contains
!
      procedure, public :: create => create_hf_geoopt_solver_factory
!
      procedure, private :: read_settings
!
   end type hf_geoopt_solver_factory
!
   interface hf_geoopt_solver_factory
!
      procedure :: new_hf_geoopt_solver_factory
!
   end interface hf_geoopt_solver_factory
!
contains
!
!
   function new_hf_geoopt_solver_factory() result(this)
!!
!!    New geoopt solver factory
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      type(hf_geoopt_solver_factory) :: this
!
      call this%read_settings()
!
   end function new_hf_geoopt_solver_factory
!
!
   subroutine create_hf_geoopt_solver_factory(this, wf, internals, solver)
!!
!!    Create
!!    Written by Alexander C. Paul, May 2022
!!
      use hf_energy_function_class, only: hf_energy_function
      use redundant_internal_coords_class, only: redundant_internal_coords
!
      implicit none
!
      class(hf_geoopt_solver_factory), intent(inout) :: this
      class(hf),                       intent(inout) :: wf
      type(redundant_internal_coords), intent(inout) :: internals
      class(bfgs_solver), allocatable, intent(out)   :: solver
!
      type(hf_energy_function), allocatable :: energy_function
!
      energy_function = hf_energy_function(wf)
!
      solver = bfgs_solver(energy_function,          &
                           internals,                &
                           this%max_iterations,      &
                           this%max_step,            &
                           this%energy_threshold,    &
                           this%gradient_threshold)
!
   end subroutine create_hf_geoopt_solver_factory
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad, May 2021
!!
      use global_in, only: input
      use global_out, only: output
!
      implicit none
!
      class(hf_geoopt_solver_factory), intent(inout)  :: this
      character(len=200) :: algorithm
!
      algorithm = "bfgs"
!
      call input%get_keyword('algorithm', 'solver scf geoopt', algorithm)
!
      if (trim(algorithm) == 'bfgs') then
!
         this%max_iterations = 100
         call input%get_keyword('max iterations', 'solver scf geoopt', this%max_iterations)
!
         this%max_step = half
         call input%get_keyword('max step', 'solver scf geoopt', this%max_step)
!
         this%energy_threshold = 1.0d-6
         call input%get_keyword('energy threshold', 'solver scf geoopt', this%energy_threshold)
!
         this%gradient_threshold = 3.0d-4
         call input%get_keyword('gradient threshold', 'solver scf geoopt', this%gradient_threshold)
!
      else
!
         call output%error_msg('did not recognize hf geoopt algorithm: (a0)', &
                                chars=[algorithm])
!
      endif
!
   end subroutine read_settings
!
end module hf_geoopt_solver_factory_class
