module hf_solver_class
!
!!
!!		HF solver class module
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
   use diis_solver_class
   use file_class
   use hf_class
   use disk_manager_class
   use libint_initialization
!
   implicit none
!
   type :: hf_solver
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
      real(dp), dimension(:,:), allocatable :: cholesky_ao_overlap ! The non-zero leading rank-by-rank 
                                                                   ! block L' in S = L L^T = (L', L'^T 0; 0, 0) 
!
      real(dp), dimension(:,:), allocatable :: permutation_matrix  ! Permutation matrix 
!
   contains
!
      procedure :: initialize           => initialize_hf_solver
      procedure :: run                  => run_hf_solver
      procedure :: finalize             => finalize_hf_solver
!
      procedure :: print_banner         => print_banner_hf_solver
!
      procedure :: decompose_ao_overlap => decompose_ao_overlap_hf_solver
      procedure :: do_roothan_hall      => do_roothan_hall_hf_solver
!
   end type hf_solver
!
!
contains
!
!
   subroutine initialize_hf_solver(solver, wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_solver) :: solver
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
   end subroutine initialize_hf_solver
!
!
   subroutine run_hf_solver(solver, wf)
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
      class(hf_solver) :: solver
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
      real(dp) :: ddot
!
      integer(i15) :: iteration
!
      real(dp), dimension(:,:), allocatable :: D ! Parameters, D_αβ
      real(dp), dimension(:,:), allocatable :: F ! Equations, F_ia
!
      integer(i15) :: n_s 
!
      real(dp), dimension(:,:), allocatable :: eri_deg
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz
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
      call mem%alloc(sp_eri_schwarz, n_s, n_s)
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, n_s)
!  
      call mem%alloc(eri_deg, n_s**2, n_s**2)
      call wf%determine_degeneracy(eri_deg, n_s)
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
      call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
      call wf%initialize_mo_coefficients()
!
      call solver%do_roothan_hall(wf)
!
      call wf%construct_ao_density() ! Construct AO density from C
      call wf%get_ao_density(D)
      call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
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
            call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
!           Solve the Roothan-Hall equation and update the AO density 
!           and Fock matrix from the solution
!
            call solver%do_roothan_hall(wf) 
!
            call wf%construct_ao_density()
            call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
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
      call mem%dealloc(eri_deg, n_s**2, n_s**2)
!
      call wf%destruct_mo_coefficients()
      call wf%destruct_ao_density()
      call wf%destruct_ao_fock()
      call wf%destruct_ao_overlap()
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call solver%finalize()
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
   end subroutine run_hf_solver
!
!
   subroutine finalize_hf_solver(solver)
!!
!! 	Finalize
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_solver) :: solver
!
   end subroutine finalize_hf_solver
!
!
   subroutine decompose_ao_overlap_hf_solver(solver, wf)
!!
!!    Decompose AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_solver) :: solver
!
      class(hf) :: wf
!
      integer(kind=4), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(:, :), allocatable :: L
!
      real(dp) :: lin_dep_threshold = 1.0D-3
!
      integer(i15) :: i, j
!
      allocate(used_diag(wf%n_ao, 1))
      used_diag = 0
!
      call mem%alloc(L, wf%n_ao, wf%n_ao) ! Full Cholesky vector
      L = zero
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, wf%n_so, &
                                                 lin_dep_threshold, used_diag)
!
      call mem%alloc(solver%cholesky_ao_overlap, wf%n_so, wf%n_so) 
      solver%cholesky_ao_overlap(:,:) = L(1:wf%n_so, 1:wf%n_so)
!
      call mem%dealloc(L, wf%n_ao, wf%n_ao)
!
!     Make permutation matrix P
!
      call mem%alloc(solver%permutation_matrix, wf%n_ao, wf%n_so)
!
      solver%permutation_matrix = zero
!
      do j = 1, wf%n_so
!
         solver%permutation_matrix(used_diag(j, 1), j) = one
!
      enddo
!
      deallocate(used_diag)
!
   end subroutine decompose_ao_overlap_hf_solver
!
!
   subroutine do_roothan_hall_hf_solver(solver, wf)
!!
!!    Do Roothan Hall
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves P^T F P (P^T C) = P^T S P P^T C e = L L^T (P^T C) e
!!
      implicit none
!
      class(hf_solver) :: solver 
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: metric 
      real(dp), dimension(:,:), allocatable :: ao_fock 
      real(dp), dimension(:,:), allocatable :: orbital_energies 
      real(dp), dimension(:,:), allocatable :: FP 
!
      real(dp) :: ddot, norm
!
      integer(i15) :: info = 0
!
      call mem%alloc(metric, wf%n_so, wf%n_so)
!
      call dgemm('N','T',                     &
                  wf%n_so,                    &
                  wf%n_so,                    &
                  wf%n_so,                    &
                  one,                        &
                  solver%cholesky_ao_overlap, & 
                  wf%n_so,                    &
                  solver%cholesky_ao_overlap, &
                  wf%n_so,                    &
                  zero,                       &
                  metric,                     & ! metric = L L^T
                  wf%n_so)
!
!     Allocate reduced space matrices 
!
      call mem%alloc(ao_fock, wf%n_so, wf%n_so)
      call mem%alloc(orbital_energies, wf%n_so, 1)
!
!     Construct reduced space Fock matrix, F' = P^T F P,
!     which is to be diagonalized over the metric L L^T 
!
      call mem%alloc(FP, wf%n_ao, wf%n_so)
!
      call dgemm('N','N',                    &
                  wf%n_ao,                   &
                  wf%n_so,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%ao_fock,                &
                  wf%n_ao,                   &
                  solver%permutation_matrix, &
                  wf%n_ao,                   &
                  zero,                      &  
                  FP,                        &
                  wf%n_ao)
!
      call dgemm('T','N',                    &
                  wf%n_so,                   &
                  wf%n_so,                   &
                  wf%n_ao,                   &
                  one,                       &
                  solver%permutation_matrix, &
                  wf%n_ao,                   &
                  FP,                        &
                  wf%n_ao,                   &
                  zero,                      &
                  ao_fock,                   & ! F' = P^T F P
                  wf%n_so)   
!
      call mem%dealloc(FP, wf%n_so, wf%n_ao)   
!
!     Solve F'C' = L L^T C' e
!
      info = 0
!
      call mem%alloc(work, 4*wf%n_so, 1)
      work = zero
!
      call dsygv(1, 'V', 'L',       &
                  wf%n_so,          &
                  ao_fock,          & ! ao_fock on entry, orbital coefficients on exit
                  wf%n_so,          &
                  metric,           &
                  wf%n_so,          &
                  orbital_energies, &
                  work,             &
                  4*(wf%n_so),      &
                  info)
!
      call mem%dealloc(metric, wf%n_so, wf%n_so)
      call mem%dealloc(work, 4*wf%n_so, 1)
      call mem%dealloc(orbital_energies, wf%n_so, 1)
!
      if (info .ne. 0) then 
!
         write(output%unit, '(/t3,a/)') 'Error: could not solve Roothan-Hall equations.'
         stop
!
      endif
!
!     Transform back the solutions to original basis, C = P (P^T C) = P C'
!
      wf%orbital_coefficients = zero
!
      call dgemm('N','N',                    &
                  wf%n_ao,                   &
                  wf%n_so,                   &
                  wf%n_so,                   &
                  one,                       &
                  solver%permutation_matrix, &
                  wf%n_ao,                   &
                  ao_fock,                   & ! orbital coefficients 
                  wf%n_so,                   &
                  zero,                      &
                  wf%orbital_coefficients,   &
                  wf%n_ao)
!
      call mem%dealloc(ao_fock, wf%n_so, wf%n_so)
!
   end subroutine do_roothan_hall_hf_solver
!
!
   subroutine print_banner_hf_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(hf_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: Direct-integral Hartree-Fock solver'
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
   end subroutine print_banner_hf_solver
!
!
end module hf_solver_class
