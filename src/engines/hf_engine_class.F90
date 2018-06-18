module hf_engine_class
!
!!
!!		HF engine class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use diis_solver_class
   use file_class
   use hf_class
   use disk_manager_class
!
   implicit none
!
   type :: hf_engine
!
      real(dp) :: energy_threshold   = 1.0D-6
      real(dp) :: residual_threshold = 1.0D-6
!
      integer(i15) :: max_iterations = 100
!
      integer(i15) :: n_parameters = 0
      integer(i15) :: n_equations  = 0
!
      logical :: restart
!
   contains
!
      procedure :: initialize => initialize_hf_engine
      procedure :: solve      => solve_hf_engine
      procedure :: finalize   => finalize_hf_engine
!
   end type hf_engine
!
!
contains
!
!
   subroutine initialize_hf_engine(engine, wf)
!!
!!    Initialize HF engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_engine) :: engine
!
      class(hf) :: wf
!
!     Get number of parameters and equations to solve for
!
      engine%n_parameters = wf%get_n_hf_parameters()
      engine%n_equations  = wf%get_n_hf_equations()
!
   end subroutine initialize_hf_engine
!
!
   subroutine solve_hf_engine(engine, wf)
!!
!!    Solve
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Routine that solves the reference state equations associated with a
!!    wavefunction
!!
      implicit none
!
      class(hf_engine) :: engine
!
      class(hf) :: wf
!
      type(diis) :: solver
!
      logical :: converged = .false.
!
      integer(i15) :: iteration = 1
!
      real(dp), dimension(:,:), allocatable :: X ! Parameters
      real(dp), dimension(:,:), allocatable :: O ! Equations
!
!     Initialize engine (read thresholds, restart, etc., from file,
!     but also ask the wavefunction for the number of parameters to solve
!     for and related information)
!
      call engine%initialize(wf)
!
!     Solve the equations
!
      call solver%init('hf_diis', engine%n_parameters, engine%n_equations)
!
      call mem%alloc(X, engine%n_parameters, 1)
      call mem%alloc(O, engine%n_equations, 1)
!
      X = zero
      O = zero
!
      call wf%set_initial_hf_parameters(X)
!
      iteration = 1
!
      do while (.not. converged .and. iteration .lt. 100)
!
         call wf%calculate_hf_equations(O, X)
      !   call solver%update(O, X)
!
         if (iteration .eq. 1) then
!
       !     call wf%get_ao_density_from_mo_coefficients(X)
!
         endif
!
         call wf%get_ao_density_from_mo_coefficients(X)
!
         iteration = iteration + 1
!
      enddo
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call engine%finalize()
!
   end subroutine solve_hf_engine
!
!
   subroutine finalize_hf_engine(engine)
!!
!! 	Finalize SCF engine
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_engine) :: engine
!
!     A dummy routine for now
!
   end subroutine finalize_hf_engine
!
!
end module hf_engine_class
