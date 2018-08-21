module scf_diis_solver_class
!
!!
!!		Self-consistent field DIIS solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    This solver employs a direct inversion of the iterative
!!    subspace (DIIS) algorithm to get an effective density, based
!!    on previous AO densities and gradients, and from this density
!!    it diagonalizes the Roothan-Hall equation, giving a new density
!!    and gradient, Fock matrix and energy. The process is repeated
!!    until convergence.
!!
!
   use kinds
   use hf_solver_class
   use diis_tool_class
   use file_class
   use hf_class
   use disk_manager_class
   use libint_initialization
!
   implicit none
!
   type, extends(hf_solver) :: scf_diis_solver
!
      integer(i15) :: diis_dimension = 8
!
      integer(i15) :: n_parameters = 0
      integer(i15) :: n_equations  = 0
!
   contains
!
      procedure :: initialize           => initialize_scf_diis_solver
      procedure :: run                  => run_scf_diis_solver
      procedure :: finalize             => finalize_scf_diis_solver
!
      procedure :: print_banner         => print_banner_scf_diis_solver
!
   end type scf_diis_solver
!
!
contains
!
!
   subroutine initialize_scf_diis_solver(solver, wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
!     Get number of parameters and equations to solve for
!
      solver%n_parameters = (wf%n_ao)*(wf%n_ao + 1)/2
      solver%n_equations  = (wf%n_o)*(wf%n_v)
!
!     Set AO density to superposition of atomic densities (SAD)
!
      call wf%initialize_ao_density()
      call wf%set_ao_density_to_sad()
!
   end subroutine initialize_scf_diis_solver
!
!
   subroutine run_scf_diis_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    This routine solves the Roothan-Hall equations in each macro-iteration,
!!    using the resulting density and Fock matrix together with previous densities
!!    and errors to make a new effective density matrix.
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
      type(diis_tool) :: diis_manager
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_residual
!
      real(dp) :: max_grad, energy, prev_energy
!
      real(dp) :: ddot
!
      integer(i15) :: iteration
!
      real(dp), dimension(:,:), allocatable :: D ! Parameters, D_αβ
      real(dp), dimension(:,:), allocatable :: F ! Equations, F_ia
!
      integer(i15) :: n_s 
!
      real(dp), dimension(:,:), allocatable     :: sp_eri_schwarz
      integer(i15), dimension(:,:), allocatable :: sp_eri_schwarz_list
!
!     Print solver banner
!
      call solver%print_banner()
!
!     :: Initialize solver 
!
      call solver%initialize(wf)
!
!     :: Cholesky decompose the AO overlap matrix 
!
!     This determines P and L such that P^T S P = L L^T, where P is the permutation 
!     matrix, and L the AO Cholesky matrix, both of which are members of the
!     solver object.
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call solver%decompose_ao_overlap(wf) 
!
!     :: Construct screening vectors, as well as a degeneracy vector, 
!     used to construct AO Fock efficiently
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 1)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     :: Initialize the DIIS object,
!     which is used to get new effective densities from previous densities and gradients
!
      call diis_manager%init('hf_diis', solver%n_parameters, solver%n_equations, solver%diis_dimension)
!
!     :: Allocate and the density and gradient, D and F
!
      call mem%alloc(D, solver%n_parameters, 1)
      call mem%alloc(F, solver%n_equations, 1)
!
      D = zero
      F = zero
!
!     :: Get an initial AO density and initial MO coefficients
!     by solving the Roothan-Hall equations for the SAD density
!
      call wf%initialize_ao_fock()
      call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
      call wf%initialize_mo_coefficients()
!
      call solver%do_roothan_hall(wf)
!
      call wf%construct_ao_density() ! Construct AO density from C
      call wf%get_ao_density(D)
      call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     Iterative solution loop
!
      iteration = 1
      converged = .false.
!
      converged_energy   = .false.
      converged_residual = .false.
!
      prev_energy = zero
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)           Max(gradient) '
      write(output%unit, '(t3,a)') '---------------------------------------------------'
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)
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
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_residual = max_grad                .lt. solver%residual_threshold
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
            call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!           Solve the Roothan-Hall equation and update the AO density 
!           and Fock matrix from the solution
!
            call solver%do_roothan_hall(wf) 
!
            call wf%construct_ao_density()
            call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
            call wf%get_ao_density(D)
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
!
      call wf%destruct_mo_coefficients()
      call wf%destruct_ao_density()
      call wf%destruct_ao_fock()
      call wf%destruct_ao_overlap()
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call solver%finalize(wf)
      call diis_manager%finalize()
!
      call mem%dealloc(D, solver%n_parameters, 1)
      call mem%dealloc(F, solver%n_equations, 1)
!
      if (.not. converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
!
      endif 
!
   end subroutine run_scf_diis_solver
!
!
   subroutine finalize_scf_diis_solver(solver, wf)
!!
!! 	Finalize
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
   end subroutine finalize_scf_diis_solver
!
!
   subroutine print_banner_scf_diis_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: Direct-integral self-consistent field DIIS Hartree-Fock solver'
      write(output%unit, '(t3,a/)') ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(t3,a)')  'This solver employs a direct inversion of the iterative'
      write(output%unit, '(t3,a)')  'subspace (DIIS) algorithm to get an effective density, based'
      write(output%unit, '(t3,a)')  'on previous AO densities and gradients, and from this density'
      write(output%unit, '(t3,a)')  'it diagonalizes the Roothan-Hall equation, giving a new density'
      write(output%unit, '(t3,a)')  'and gradient, Fock matrix and energy. The process is repeated'
      write(output%unit, '(t3,a/)') 'until convergence.'
      flush(output%unit)
!
   end subroutine print_banner_scf_diis_solver
!
!
end module scf_diis_solver_class
