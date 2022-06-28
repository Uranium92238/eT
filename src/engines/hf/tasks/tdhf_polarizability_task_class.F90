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
module tdhf_polarizability_task_class
!
!!
!! TDHF polarizability task class
!! Written by Sarai D. Folkestad, May 2022
!!
!
   use parameters
   use memory_manager_class, only: mem
!
   use hf_class,      only: hf
   use hf_task_class, only: hf_task
!
   use global_in, only: input
   use global_out, only: output
!
   implicit none
!
   type, extends(hf_task) :: tdhf_polarizability_task
!
      real(dp), dimension(:,:,:), allocatable :: polarizability
      real(dp), dimension(:),     allocatable :: frequencies
!
      integer :: n_frequencies, max_dim_red, max_iterations
!
      real(dp) :: threshold
!
   contains
!
      procedure, public :: execute &
                        => execute_tdhf_polarizability_task
!
      procedure, private :: create_polarizability_solver
      procedure, private :: read_polarizability_settings
      procedure, private :: compute_polarizability
!
   end type tdhf_polarizability_task
!
!
   interface tdhf_polarizability_task
!
      procedure :: new_tdhf_polarizability_task
!
   end interface tdhf_polarizability_task
!
!
contains
!
!
   function new_tdhf_polarizability_task() result(this)
!!
!!    New
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      type(tdhf_polarizability_task) :: this
!
      this%name_ = 'Determining TDHF polarizabilities'
!
   end function new_tdhf_polarizability_task
!
!
   subroutine execute_tdhf_polarizability_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(tdhf_polarizability_task), intent(inout) :: this
      class(hf), target, intent(inout) :: wf
!
      integer :: f
!
      call this%print_header()
      call this%start_timer()
!
      call this%read_polarizability_settings()
!
      call mem%alloc(this%polarizability, 3, 3, this%n_frequencies)
!
      call this%compute_polarizability(wf)
!
      do f = 1, this%n_frequencies
!
         call wf%print_polarizability(this%polarizability(:,:,f), this%frequencies(f))
!
      enddo
!
      call mem%dealloc(this%polarizability, 3, 3, this%n_frequencies)
      call mem%dealloc(this%frequencies, this%n_frequencies)
!
      call this%end_timer()
!
   end subroutine execute_tdhf_polarizability_task
!
!
   subroutine compute_polarizability(this, wf)
!!
!!    Compute polarizability
!!    Written by Sarai D. Folkestad, May 2021
!!
      use linear_davidson_solver_class, only: linear_davidson_solver
!
      implicit none

      class(tdhf_polarizability_task), intent(inout) :: this
      class(hf),                       intent(inout) :: wf
!
      class(linear_davidson_solver), allocatable :: solver
!
      real(dp), dimension(:,:), allocatable :: X, V
!
      real(dp) :: ddot
!
      integer :: k, l, f
!
      call mem%alloc(X, 2*wf%n_v*wf%n_o, 3)
      call mem%alloc(V, 2*wf%n_v*wf%n_o, 3)
!
      do f = 1, this%n_frequencies
!
         if (.not. input%is_keyword_present('print iterations', 'solver tdhf response')) &
            call output%mute()
!
         do k = 1, 3
!
            call this%create_polarizability_solver(wf, solver, k, f)
!
            call solver%run()
            call solver%get_solution(X(:,k), 1)
            call wf%get_dipole_gradient(V(:,k), k)
!
            call solver%cleanup()
!
         enddo
!
         if (.not. input%is_keyword_present('print iterations', 'solver tdhf response')) &
            call output%unmute()
!
         do k = 1, 3
            do l = 1, 3
!
               this%polarizability(k, l, f) = ddot(2*wf%n_o*wf%n_v, X(1,l), 1, V(1,k), 1)
!
            enddo
         enddo
!
      enddo
!
      call mem%dealloc(X, 2*wf%n_v*wf%n_o, 3)
      call mem%dealloc(V, 2*wf%n_v*wf%n_o, 3)
!
   end subroutine compute_polarizability
!
!
   subroutine read_polarizability_settings(this)
!!
!!    Read polarizability settings
!!    Written by Sarai D. Folkestad, May 2022
!!
      implicit none
!
      class(tdhf_polarizability_task), intent(inout) :: this
!
      if (input%is_keyword_present('frequencies', 'solver tdhf response')) then
!
         this%n_frequencies = input%get_n_elements_for_keyword('frequencies', 'solver tdhf response')
         call mem%alloc(this%frequencies, this%n_frequencies)
         call input%get_array_for_keyword('frequencies', 'solver tdhf response', this%n_frequencies, this%frequencies)
!
      else
!
         this%n_frequencies = 1
         call mem%alloc(this%frequencies, this%n_frequencies)
         this%frequencies(1) = zero
!
      endif
!
      this%threshold = 1.0d-3
      call input%get_keyword('residual threshold', 'solver tdhf response', this%threshold)
!
      this%max_dim_red = 100
      call input%get_keyword('max reduced dimension', 'solver tdhf response', this%max_dim_red)
!
      this%max_iterations = 100
      call input%get_keyword('max iterations', 'solver tdhf response', this%max_iterations)
!
   end subroutine read_polarizability_settings
!
!
   subroutine create_polarizability_solver(this, wf, solver, k, f)
!!
!!    Create polarizability solver
!!    Written by Sarai D. Folkestad, May 2022
!!
      use linear_davidson_solver_class,                     only: linear_davidson_solver
      use tdhf_response_start_vector_tool_class,            only: tdhf_response_start_vector_tool
      use linear_davidson_tool_class,                       only: linear_davidson_tool
      use tdhf_rhs_tool_class,                              only: tdhf_rhs_tool
      use linear_davidson_single_solution_print_tool_class, only: linear_davidson_single_solution_print_tool
      use null_linear_storage_tool_class,                   only: null_linear_storage_tool
      use rpa_response_preconditioner_getter_class,         only: rpa_response_preconditioner_getter
      use rpa_response_transformation_class ,               only: rpa_response_transformation
!
      implicit none
!
      class(hf),                                   intent(inout)  :: wf
      integer,                                     intent(in)     :: k, f
      class(linear_davidson_solver),  allocatable, intent(inout)  :: solver
      class(tdhf_polarizability_task), intent(in)                 :: this
!
      class(linear_davidson_tool),                          allocatable :: davidson
      class(tdhf_response_start_vector_tool),               allocatable :: start_vector
      class(tdhf_rhs_tool),                                 allocatable :: rhs
      class(linear_davidson_single_solution_print_tool),    allocatable :: printer
      class(null_linear_storage_tool),                      allocatable :: storer
      class(rpa_response_preconditioner_getter),            allocatable :: preconditioner
      class(rpa_response_transformation),                   allocatable :: transformer
!
      start_vector   = tdhf_response_start_vector_tool(wf, k)
      rhs            = tdhf_rhs_tool(wf, k)
      storer         = null_linear_storage_tool()
      preconditioner = rpa_response_preconditioner_getter(wf, this%frequencies(f))
!
      davidson = linear_davidson_tool(name_          = 'tdhf_response', &
                                   n_parameters      = 2*wf%n_v*wf%n_o, &
                                   lindep_threshold  = 1.0D-11,         &
                                   max_dim_red       = this%max_dim_red,&
                                   n_equations       = 1,               &
                                   records_in_memory = .false.)
!
      printer = linear_davidson_single_solution_print_tool()
      call printer%print_banner()
!
      transformer = rpa_response_transformation(wf, this%frequencies(f))
!
      solver = linear_davidson_solver(transformer             = transformer,           &
                                        davidson              = davidson,              &
                                        storer                = storer,                &
                                        start_vector          = start_vector,          &
                                        preconditioner        = preconditioner,        &
                                        rhs_getter            = rhs,                   &
                                        printer               = printer,               &
                                        n_rhs                 = 1,                     &
                                        n_solutions           = 1,                     &
                                        max_iterations        = this%max_iterations,   &
                                        residual_threshold    = this%threshold,        &
                                        frequencies           = [zero])
!
   end subroutine create_polarizability_solver
!
!
end module tdhf_polarizability_task_class
