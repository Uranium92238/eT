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
   use libint_initialization
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
      integer(i15) :: diis_dimension = 8
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
      engine%n_parameters = (wf%n_ao)*(wf%n_ao + 1)/2
      engine%n_equations  = (wf%n_o)*(wf%n_v)
!
!     Set AO density to superposition of atomic densities (SAD)
!
      call wf%initialize_ao_density()
      call wf%set_ao_density_to_sad()
!
   end subroutine initialize_hf_engine
!
!
   subroutine solve_hf_engine(engine, wf)
!!
!!    Solve 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    This routine solves the Roothan-Hall equations in each macro-iteration,
!!    using the resulting density and Fock matrix together with previous densities
!!    and errors to make a new effective density matrix.
!!
      implicit none
!
      class(hf_engine) :: engine
!
      class(hf) :: wf
!
      type(diis) :: diis_manager
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_residual
!
      real(dp) :: max_grad, energy, prev_energy
!
      real(dp) :: ddot, omp_get_wtime
!
      real(dp) :: t0, t1
!
      integer(i15) :: iteration = 1, ao = 0
!
      real(dp), dimension(:,:), allocatable :: D ! Parameters, D_αβ
      real(dp), dimension(:,:), allocatable :: F ! Equations, F_ia
!
      integer(i15) :: n_s 
!
      real(dp), dimension(:,:), allocatable :: eri_deg
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz
!
!     Print engine banner
!
      write(output%unit, '(/t3,a)') ':: Direct-integral Hartree-Fock engine'
      write(output%unit, '(t3,a/)') ':: E. F. Kjønstad, S. D. Folkestad, 2018'
      flush(output%unit)
!
!     Initialize engine (read thresholds, restart, etc., from file,
!     but also ask the wavefunction for the number of parameters to solve
!     for and related information)
!
      call engine%initialize(wf)
!
!     Construct screening vectors, as well as a degeneracy vector, 
!     needed to construct AO Fock efficiently
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s, n_s)
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, n_s)
!  
      call mem%alloc(eri_deg, n_s**2, n_s**2)
      call wf%determine_degeneracy(eri_deg, n_s)
!
!     Solve the equations
!
      call diis_manager%init('hf_diis', engine%n_parameters, engine%n_equations, engine%diis_dimension)
!
      call mem%alloc(D, engine%n_parameters, 1)
      call mem%alloc(F, engine%n_equations, 1)
!
      D = zero
      F = zero
!
!     Get an initial AO density and initial MO coefficients for loop
!
      call wf%initialize_orbital_energies()
      call wf%initialize_ao_fock()
!
      call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
      call wf%initialize_mo_coefficients()
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
!
      call wf%solve_roothan_hall() ! F^AO C = S C e to get new MOs C
!
      call wf%construct_ao_density() ! Construct AO density from C
      call wf%get_ao_density(D)
      call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
!     Iterative solution loop
!
      iteration   = 1
      converged   = .false.
!
      converged_energy   = .false.
      converged_residual = .false.
!
      prev_energy = zero
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)           Max(gradient) '
      write(output%unit, '(t3,a)') '---------------------------------------------------'
!
      do while (.not. converged .and. iteration .le. engine%max_iterations)
!
!        Calculate the occ-vir block, or gradient, from the current density
!
         F = zero
         call wf%get_fock_ov(F)
!
!        Determine the maximum element of the gradient
!
         max_grad = get_abs_max(F, (wf%n_o)*(wf%n_v))
!
!        Set current energy
!
         energy = wf%hf_energy
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,f17.12)') iteration, wf%hf_energy, max_grad
         flush(output%unit)
!
!        Test for convergence:
!
         converged_energy   = abs(energy-prev_energy) .lt. engine%energy_threshold
         converged_residual = max_grad                .lt. engine%residual_threshold
!
         converged = converged_residual .and. converged_energy
!
         if (converged) then
!
!
            write(output%unit, '(t3,a)') '---------------------------------------------------'
            write(output%unit, '(/t3,a13,i3,a12/)') 'Converged in ', iteration, ' iterations!'
!
         else
!
!           Get an optimal averaged density matrix from DIIS
!
            call diis_manager%update(F, D)
!
!           Construct the AO fock matrix from it
!
            call wf%set_ao_density(D)
            call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
!           Solve the Roothan-Hall equation and update the AO density 
!           and Fock matrix from the solution
!
            call wf%solve_roothan_hall() 
!
            call wf%construct_ao_density()
            call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
            call wf%get_ao_density(D) ! For next iteration
!
            prev_energy = wf%hf_energy
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(sp_eri_schwarz, n_s, n_s)
      call mem%dealloc(eri_deg, n_s**2, n_s**2)
!
      call wf%destruct_mo_coefficients()
      call wf%destruct_ao_density()
      call wf%destruct_ao_fock()
      call wf%destruct_ao_overlap()
      call wf%destruct_orbital_energies()
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call engine%finalize()
      call diis_manager%finalize()
!
      if (.not. converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
!
      endif 
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
