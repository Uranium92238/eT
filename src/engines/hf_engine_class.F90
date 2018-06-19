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
      real(dp) :: error, ddot
!
      integer(i15) :: iteration = 1
!
      real(dp), dimension(:,:), allocatable :: D ! Parameters, D_αβ
      real(dp), dimension(:,:), allocatable :: F ! Equations, F_ia
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
      call mem%alloc(D, engine%n_parameters, 1)
      call mem%alloc(F, engine%n_equations, 1)
!
      D = zero
      F = zero
!
!     Get an initial AO density and initial MO coefficients for loop
!
      call wf%initialize_ao_density()
      call wf%set_ao_density_to_soad_guess() ! D^AO = D^SOAD
!
      call wf%initialize_ao_fock()
      call wf%construct_ao_fock() ! From current D^AO
!
      call wf%initialize_mo_coefficients()
      call wf%solve_roothan_hall() ! F^AO C = S C e to get new MOs C
!
      call wf%construct_ao_density() ! Construct AO density from C
      call wf%get_ao_density(D)
      call wf%construct_ao_fock()    ! Update the AO Fock
!
!     Iterative solution loop
!
      iteration = 1
      converged = .false.
!
      do while (.not. converged .and. iteration .le. engine%max_iterations)
!
!        Calculate the occ-vir block from the current AO density
!        and MO coefficients
!
         F = zero
         call wf%get_hf_equations(F)
!
!        Check the error
!
         error = ddot(engine%n_equations, F, 1, F, 1)
         error = sqrt(error)
!
         if (error .lt. 1.0D-6) then
!
            converged = .true.
!
         else
!
            call wf%calculate_hf_energy()
!
            write(output%unit, '(/t3,a11,i3)')     'Iteration: ', iteration
            write(output%unit, '(t3,a11,f17.12)')  'Error:     ', error
            write(output%unit, '(t3,a11,f17.12/)') 'Energy:    ', wf%hf_energy
            flush(output%unit)
!
!           Get an averaged density matrix by DIIS
!
            call solver%update(F, D)
!
!           Construct the AO fock matrix
!
            call wf%set_ao_density(D)
!
            call wf%construct_ao_fock()
!
            call wf%solve_roothan_hall() ! Updates the MO coefficients
!
!           Get the new AO density and new AO Fock matrix
!
            call wf%construct_ao_density()
            call wf%construct_ao_fock()
!
            call wf%get_ao_density(D)
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call wf%destruct_mo_coefficients()
      call wf%destruct_ao_density()
      call wf%destruct_ao_fock()
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
   end subroutine finalize_hf_engine
!
!
end module hf_engine_class
