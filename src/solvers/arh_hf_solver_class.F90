module arh_hf_solver_class
!
!!
!!		Augmented Roothan-Hall HF solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use diis_solver_class
   use file_class
   use hf_class
   use hf_solver_class
   use array_utilities
   use disk_manager_class
   use libint_initialization
!
   implicit none
!
   type, extends(hf_solver) :: arh_hf_solver
!
      integer(i15) :: max_micro_iterations = 250
!
      real(dp) :: purification_threshold           = 1.0D-8
      real(dp) :: relative_micro_threshold         = 1.0D-3 
      real(dp) :: relative_shifted_micro_threshold = 1.0D-3 
!
      real(dp) :: trust_radius                     = 0.50D0 
      real(dp) :: relative_trust_radius_threshold  = 0.10D0
!
      real(dp) :: rotation_norm_threshold          = 0.2D0
!
      integer(i15) :: diis_dimension = 10
!
      real(dp), dimension(:,:), allocatable :: trace_matrix        ! Trace matrix, T_ij = Tr (D_in D_jn), i, j < n
!
      type(file) :: RH_gradients_file 
      type(file) :: AO_densities_file 
!
      integer(i15) :: history = 50
      integer(i15) :: current_index
!
   contains
!
      procedure :: initialize  => initialize_arh_hf_solver
      procedure :: run         => run_arh_hf_solver
      procedure :: finalize    => finalize_arh_hf_solver
!
      procedure, private :: print_banner                                => print_banner_arh_hf_solver
!
      procedure, private :: solve_aug_Newton_equation                   => solve_aug_Newton_equation_arh_hf_solver
      procedure, private :: solve_level_shifted_aug_Newton_equation     => solve_level_shifted_aug_Newton_equation_arh_hf_solver
!
      procedure, private :: construct_stationary_roothan_hall_condition => construct_stationary_roothan_hall_condition_arh_hf_solver
      procedure, private :: add_augmented_Roothan_Hall_contribution     => add_augmented_Roothan_Hall_contribution_arh_hf_solver
!
      procedure, private :: construct_trace_matrix                      => construct_trace_matrix_arh_hf_solver
!
      procedure, private :: rotate_and_purify                           => rotate_and_purify_arh_hf_solver
      procedure, private :: construct_and_pack_gradient                 => construct_and_pack_gradient_arh_hf_solver

!
   end type arh_hf_solver
!
!
contains
!
!
   subroutine initialize_arh_hf_solver(solver, wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(arh_hf_solver) :: solver
!
      class(hf) :: wf
!
      call wf%initialize_ao_density()
      call wf%set_ao_density_to_sad()
!
   end subroutine initialize_arh_hf_solver
!
!
   subroutine run_arh_hf_solver(solver, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(arh_hf_solver) :: solver
!
      class(hf) :: wf
!
      type(diis) :: diis_solver
!
      real(dp), dimension(:,:), allocatable :: Xr       ! Full rotation matrix, antisymmetric, preconditioned 
      real(dp), dimension(:,:), allocatable :: Xr_pck   ! Packed variant of X 
!
      real(dp), dimension(:,:), allocatable :: RHC     ! Roothan-Hall stationary condition, full, preconditioned 
      real(dp), dimension(:,:), allocatable :: RHC_pck ! Packed variant of RHC 
!
      real(dp), dimension(:,:), allocatable :: X       ! Full (lin.dep.) rotation matrix, antisymmetric
      real(dp), dimension(:,:), allocatable :: X_pck   ! Full (lin.dep.) rotation matrix, antisymmetric
!
      real(dp), dimension(:,:), allocatable :: Po, Pv  ! Projection matrices
!
      real(dp), dimension(:,:), allocatable :: VT      ! Preconditioner, now = LT in S = L L^T, but may be changed
      real(dp), dimension(:,:), allocatable :: inv_VT
!
      real(dp), dimension(:,:), allocatable :: S       ! Preconditioned S 
!
      real(dp), dimension(:,:), allocatable :: H       ! Full (lin.dep.) Roothan-Hall Hessian
      real(dp), dimension(:,:), allocatable :: G       ! Full (lin.dep.) Roothan-Hall gradient
!
      real(dp), dimension(:,:), allocatable :: Hr      ! Roothan-Hall Hessian, preconditioned
      real(dp), dimension(:,:), allocatable :: Gr      ! Roothan-Hall gradient, preconditioned 
!
      real(dp) :: norm_X 
!
      logical :: transpose_left
!
      real(dp) :: beta, coulomb_thr, exchange_thr
!
      integer(i15) :: iteration = 1
!
      real(dp) :: start_timer, end_timer, omp_get_wtime
!
      logical :: converged          = .false.
      logical :: converged_residual = .false.
      logical :: converged_energy   = .false.
!
      real(dp) :: energy 
      real(dp) :: prev_energy 
!
      real(dp) :: max_grad, alpha, level_shift
!
      real(dp) :: ddot
!
      integer(i15) :: n_s, i 
!
      real(dp) :: trace_of_ao_density
!
      real(dp), dimension(:,:), allocatable :: eri_deg
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Initialize solver (thresholds, etc., and wf initialization)
!
      call solver%initialize(wf)
!
!     Initialize gradient and AO densities file 
!
      call solver%RH_gradients_file%init('RH_gradients', 'sequential', 'unformatted')
      call solver%AO_densities_file%init('AO_densities', 'sequential', 'unformatted')
!
!     Construct screening vectors, as well as a degeneracy vector, 
!     needed to construct AO Fock efficiently
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s, n_s)
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, n_s)
!
!     Construct initial AO Fock from the SOAD density
!
      call wf%initialize_ao_fock()
!
      trace_of_ao_density = zero 
      do i = 1, wf%n_ao 
         trace_of_ao_density = trace_of_ao_density + wf%ao_density(i,i)
      enddo
      !write(output%unit, *) 'Trace of ao density: ', trace_of_ao_density
!
      coulomb_thr  = 0.5D0
      exchange_thr = 0.5D0
!
      start_timer = omp_get_wtime()
      call wf%construct_ao_fock(sp_eri_schwarz, n_s, coulomb_thr, exchange_thr) 
      end_timer = omp_get_wtime()
     ! write(output%unit, *) 'Time to construct AO Fock from SAD: ', end_timer-start_timer
     ! flush(output%unit)
      prev_energy = wf%hf_energy
!
!     Construct AO overlap matrix, Cholesky decompose it,
!     then precondition (making it the identity matrix for this preconditioner, V)
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call solver%decompose_ao_overlap(wf) 
!
!     :: Solve Roothan Hall once - using the SOAD guess - to get a decent N-representable 
!     AO density on which to start doing conjugate gradient
!
      call wf%initialize_mo_coefficients()
      call solver%do_roothan_hall(wf)
!
!     From the obtained MO coefficients, update the AO density and from it the AO Fock matrix
!
      call wf%construct_ao_density()
      call wf%destruct_mo_coefficients()
!
      trace_of_ao_density = zero 
      do i = 1, wf%n_ao 
         trace_of_ao_density = trace_of_ao_density + wf%ao_density(i,i)
      enddo
      !write(output%unit, *) 'Trace of ao density: ', trace_of_ao_density
!
      start_timer = omp_get_wtime()
      call wf%construct_ao_fock(sp_eri_schwarz, n_s)  
      end_timer = omp_get_wtime()
     ! write(output%unit, *) 'Time to construct AO Fock: ', end_timer-start_timer 
!
!     :: Construct precondition matrices, used to transform H and G prior to solving the Newton equation 
!
      call mem%alloc(VT, wf%n_so, wf%n_so)
      call mem%alloc(inv_VT, wf%n_so, wf%n_so)
!
      VT = transpose(solver%cholesky_ao_overlap)
!
      call inv_lower_tri(inv_VT, solver%cholesky_ao_overlap, wf%n_so)
      inv_VT = transpose(inv_VT)
!
      call mem%dealloc(solver%cholesky_ao_overlap, wf%n_so, wf%n_so)    
!
!     :: Construct the preconditioned S matrix, which is the identity matrix for the 
!     particular case of preconditioner chosen here 
!
      call mem%alloc(S, wf%n_so, wf%n_so) ! P^T S P = V V^T 
!
      call dgemm('T','N',  &
                  wf%n_so, &
                  wf%n_so, &
                  wf%n_so, &
                  one,     &
                  VT,      &
                  wf%n_so, &
                  VT,      &
                  wf%n_so, &
                  zero,    &
                  S,       &
                  wf%n_so)
!
      call sandwich(S, inv_VT, inv_VT, wf%n_so) ! S <- (V-T)^T S V-T = I
!
!     :: Allocate and zero projection matrices on occupied and virtual spaces
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)  
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)  
!
      Po = zero
      Pv = zero
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)        Max(grad.)    ΔE (a.u.)'
      write(output%unit, '(t3,a)') '----------------------------------------------------------'
!
      iteration = 1
!
      converged          = .false.
      converged_energy   = .false.
      converged_residual = .false.
!
      call mem%alloc(G, wf%n_ao, wf%n_ao)
      call mem%alloc(H, wf%n_ao, wf%n_ao)
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)
!
!        Construct the projection matrices Po and Pv
!        from the current AO density
!
         call wf%construct_projection_matrices(Po, Pv)
!
!        Use the projections to construct the gradient G, then pack it in 
!
         call wf%construct_roothan_hall_gradient(G, Po, Pv)
!
!        Make the reduced space G (where linearly dependent directions have been removed),
!        and determine its maximum (absolute) element 
!
         call mem%alloc(Gr, wf%n_so, wf%n_so) 
         call symmetric_sandwich(Gr, G, solver%permutation_matrix, wf%n_ao, wf%n_so)

         max_grad = get_abs_max(Gr, wf%n_so*wf%n_so)
!
!        Set current energy
!
         energy = wf%hf_energy
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e10.4,4x,e10.4)') iteration, wf%hf_energy, &
                           max_grad, abs(wf%hf_energy-prev_energy)
         flush(output%unit)
!
!        Test for convergence:
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_residual = max_grad                .lt. solver%residual_threshold
!
         converged = converged_energy .and. converged_residual 
!
         if (converged) then ! Done, hooray
!
            write(output%unit, '(t3,a)') '----------------------------------------------------------'
            write(output%unit, '(/t3,a13,i3,a12/)') 'Converged in ', iteration, ' iterations!'
!
         else ! Rotate density matrix
!
!           :: Construct the Hessian H (one-electron terms), transform it to the linearly independent
!           basis, then precondition both it and G; i.e., replace Y by V-1 Y V-T for Y = G and H
!
            call wf%construct_roothan_hall_hessian(H, Po, Pv)
!
            call mem%alloc(Hr, wf%n_so, wf%n_so)
            call symmetric_sandwich(Hr, H, solver%permutation_matrix, wf%n_ao, wf%n_so)
!
            call sandwich(Gr, inv_VT, inv_VT, wf%n_so)
            call sandwich(Hr, inv_VT, inv_VT, wf%n_so)
!
!           Write current RH gradient and AO density to file 
!
            solver%current_index = iteration - &
                  (solver%history-1)*((iteration-1)/(solver%history-1))
!
            if (solver%current_index .eq. 1) then 
!
               call disk%open_file(solver%RH_gradients_file, 'readwrite', 'rewind')
               call disk%open_file(solver%AO_densities_file, 'readwrite', 'rewind')
!
            else
!
               call disk%open_file(solver%RH_gradients_file, 'readwrite', 'append')
               call disk%open_file(solver%AO_densities_file, 'readwrite', 'append')  
!
            endif 
!
            write(solver%RH_gradients_file%unit) Gr 
            write(solver%AO_densities_file%unit) wf%ao_density 
!
            call disk%close_file(solver%RH_gradients_file)
            call disk%close_file(solver%AO_densities_file)
!
!           Construct trace matrix, T_ij = Tr (D_in D_jn), whose inverse 
!           is needed to solve the augmented Newton equations 
!
            call solver%construct_trace_matrix(wf, iteration)
!
!           :: Perform micro-iterations to get a direction X in which to rotate
!           (solves the equation of RH gradient equal to zero to first order in X)
!
            call mem%alloc(Xr, wf%n_so, wf%n_so)              ! Full antisymmetric rotation matrix
            call mem%alloc(Xr_pck, packed_size(wf%n_so-1), 1) ! Packed rotation matrix (strictly lower triangular part)
!
            Xr     = zero
            Xr_pck = zero
!
            call solver%solve_aug_Newton_equation(wf, Xr_pck, Hr, Gr, S, max_grad)
            call squareup_anti(Xr_pck, Xr, wf%n_so)
!
!           :: Solve the augmented Hessian (i.e., level shifted Newton equation) if
!           the converged step X is longer than the trust radius within which the Roothan-Hall
!           energy expansion can (presumably) be trusted to second order in X
!
            norm_X = sqrt(ddot(packed_size(wf%n_so-1), Xr_pck, 1, Xr_pck, 1))
            level_shift = zero 
!
            if (norm_X .gt. solver%trust_radius .and. abs(norm_X-solver%trust_radius) .gt. & 
                        solver%relative_trust_radius_threshold*solver%trust_radius) then
!
               write(output%unit, '(/t6,a/,t6,a)') 'Step outside trust region. Will attempt to', &
                                                   'solve the augmented Hessian equation.'
               flush(output%unit)
!
               call solver%solve_level_shifted_aug_Newton_equation(wf, Xr_pck, Hr, Gr, &
                                                         S, max_grad, norm_X, level_shift)
               call squareup_anti(Xr_pck, Xr, wf%n_so)
! 
            endif
!          
            call mem%dealloc(Xr_pck, packed_size(wf%n_so-1), 1) 
!
            call mem%dealloc(Gr, wf%n_so, wf%n_so) 
            call mem%dealloc(Hr, wf%n_so, wf%n_so)    
!
!           :: Convert the converged direction X from the basis of the preconditioned system (X' = V^T X V)
!           back to the original system (X -> V-T X V-1), then transform the vector back to full (lin.dep.) 
!           space, packing it into the current direction vector s
!
            transpose_left = .false.
            call sandwich(Xr, inv_VT, inv_VT, wf%n_so, transpose_left) 
! 
            call mem%alloc(X, wf%n_ao, wf%n_ao)
            call symmetric_sandwich_right(X, Xr, solver%permutation_matrix, wf%n_ao, wf%n_so)         
            call mem%dealloc(Xr, wf%n_so, wf%n_so) 
!
            call mem%alloc(X_pck, packed_size(wf%n_ao-1), 1)
            call packin_anti(X_pck, X, wf%n_ao) 
            call mem%dealloc(X, wf%n_ao, wf%n_ao)
!
            norm_X = sqrt(ddot(packed_size(wf%n_ao-1), X_pck, 1, X_pck, 1))
!
!           :: Rotate the AO density and purify it           
!
            call solver%rotate_and_purify(wf, X_pck, one, norm_X)
            call mem%dealloc(X_pck, packed_size(wf%n_ao-1), 1)
!
!           :: Construct AO Fock (and energy) with the new rotated density,
!           and set previous gradient and step direction to current, in preparation 
!           for next conjugate gradient iteration
!
            trace_of_ao_density = zero 
            do i = 1, wf%n_ao 
               trace_of_ao_density = trace_of_ao_density + wf%ao_density(i,i)
            enddo
           ! write(output%unit, *) 'Trace of ao density: ', trace_of_ao_density
!
            coulomb_thr  = max(1.0D-10, max_grad)
            exchange_thr = max(1.0D-8, max_grad)
!
            prev_energy = wf%hf_energy
            start_timer = omp_get_wtime()
            call wf%construct_ao_fock(sp_eri_schwarz, n_s, coulomb_thr, exchange_thr)
            end_timer = omp_get_wtime()
            !write(output%unit, *) 'Time to construct AO Fock: ', end_timer-start_timer 
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(G, wf%n_ao, wf%n_ao)
      call mem%dealloc(H, wf%n_ao, wf%n_ao)
!
      call wf%destruct_ao_density()
      call wf%destruct_ao_fock()
      call wf%destruct_ao_overlap()
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(S, wf%n_so, wf%n_so)
!
      call mem%dealloc(sp_eri_schwarz, n_s, n_s)
!
      call mem%dealloc(VT, wf%n_so, wf%n_so)
      call mem%dealloc(inv_VT, wf%n_so, wf%n_so)
!
!     Initialize solver (make final deallocations, and other stuff)
!
      call solver%finalize(wf)
!
      if (.not. converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
!
      endif 
!
   end subroutine run_arh_hf_solver
!
!
   subroutine finalize_arh_hf_solver(solver, wf)
!!
!! 	Finalize
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(arh_hf_solver) :: solver
!
      class(hf) :: wf
!
   end subroutine finalize_arh_hf_solver
!
!
   subroutine print_banner_arh_hf_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(arh_hf_solver) :: solver
!
      write(output%unit, '(/t3,a)') ':: Direct-integral augmented Roothan-Hall Hartree-Fock solver'
      write(output%unit, '(t3,a)')  ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(/t3,a)') 'This solver uses a trust-region augmented Roothan-Hall method (TR-ARH)'
      write(output%unit, '(t3,a)')  'to converge the Hartree-Fock energy. See Høst et al. for a detailed'
      write(output%unit, '(t3,a/)') 'description of the algorithm [J. Chem. Phys. 129, 124106 (2008)].'
      flush(output%unit)
!
   end subroutine print_banner_arh_hf_solver
!
!
   subroutine rotate_and_purify_arh_hf_solver(solver, wf, X_pck, kappa, norm_X)
!!
!!    Rotate and purify
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Rotates the AO density by the antisymmetric rotation matrix X,
!!    which is sent to the routine as a packed antisymmetric matrix,
!!    and then purifies the resulting density (making it idempotent)
!!
!!    The parameter kappa is the prefactor on the rotation X.
!!    We rotate 2^k times by 2^-k to make the errors in the 
!!    rotation smaller, where k is a parameter determined by the 
!!    rotation norm threshold.
!!
      implicit none 
!
      class(arh_hf_solver), intent(in) :: solver
!
      class(hf) :: wf 
!
      real(dp), dimension((wf%n_ao-1)*(wf%n_ao)/2, 1), intent(in) :: X_pck
!
      real(dp), intent(in) :: kappa, norm_X 
!
      integer(i15) :: k
!
      integer(i15) :: i 
!
      real(dp), dimension(:,:), allocatable :: X 
!
!     Square up and scale rotation vector by kappa
!
      call mem%alloc(X, wf%n_ao, wf%n_ao)
      call squareup_anti(X_pck, X, wf%n_ao)
!
      X = kappa*X 
!
!     Determine k factor given the norm of X and the 
!     threshold on acceptable rotation lengths
!
      k = 0
      do while ((two**(-k))*norm_X .gt. solver%rotation_norm_threshold)
!
         k = k + 1
!
      enddo
!
!     Do the 2^-k X rotation, then purify, 2^k times 
!
      X = (two**(-k))*X 
!
      do i = 1, 2**k 
!
         call wf%rotate_ao_density(X)
         call wf%purify_ao_density(solver%purification_threshold)
!
      enddo
!
      call mem%dealloc(X, wf%n_ao, wf%n_ao)
!
   end subroutine rotate_and_purify_arh_hf_solver
!
!
   subroutine construct_and_pack_gradient_arh_hf_solver(solver, wf, G_pck, Po, Pv)
!!
!!    Construct and pack Roothan-Hall gradient
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall gradient G, then packs it into G_pck.
!!
      implicit none 
!
      class(arh_hf_solver), intent(in) :: solver 
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension((wf%n_ao-1)*(wf%n_ao)/2, 1) :: G_pck 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po  
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv  
!
      real(dp), dimension(:,:), allocatable :: G 
!
      call mem%alloc(G, wf%n_ao, wf%n_ao)
!
      call wf%construct_roothan_hall_gradient(G, Po, Pv)
      call packin_anti(G_pck, G, wf%n_ao)
!
      call mem%dealloc(G, wf%n_ao, wf%n_ao)
!
   end subroutine construct_and_pack_gradient_arh_hf_solver
!
!
   subroutine solve_aug_Newton_equation_arh_hf_solver(solver, wf, X_pck, H, G, S, max_grad)
!!
!!    Solve Newton equation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the Newton equation using a DIIS algorithm and no 
!!    level shift:
!!
!!       H X S - S X H + correction = -G.
!!
!!    On exit, X_pck is the solution X in packed antisymmetric form. If the 
!!    final X exceeds the trust  radius, a routine is called to solve the level 
!!    shifted equations instead.
!!
!!    The convergence threshold is max_grad*relative_micro_threshold (it becomes
!!    tighter as the maximum element of the gradient decreases) on the norm of 
!!    the equation error.
!!
      implicit none 
!
      class(arh_hf_solver), intent(in) :: solver 
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_so-1)*(wf%n_so)/2, 1) :: X_pck
!
!
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: H
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: G
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: S
!
      real(dp), intent(in) :: max_grad ! Maximum element of the gradient G
!
      real(dp), dimension(:,:), allocatable :: X 
!
      real(dp), dimension(:,:), allocatable :: RHC     ! H X S - S X H + G (= 0 on convergence)
      real(dp), dimension(:,:), allocatable :: RHC_pck ! Packed variant
!
      logical      :: micro_iterate
      integer(i15) :: micro_iteration
      real(dp)     :: micro_error 
!
      real(dp) :: ddot
!
      integer(i15) :: i, j
!
      type(diis) :: diis_solver 
!
!     Perform micro-iterations to get a direction X in which to rotate
!     (solves the equation of RH gradient equal to zero to first order in X)
!
      call mem%alloc(X, wf%n_so, wf%n_so)
!
      X = -G
      call packin_anti(X_pck, X, wf%n_so)
!
      micro_iteration = 1
      micro_iterate = .true.
      call diis_solver%init('hf_newton',                &
                              packed_size(wf%n_so - 1), &
                              packed_size(wf%n_so - 1), &
                              solver%diis_dimension)
!
      call mem%alloc(RHC, wf%n_so, wf%n_so)
      call mem%alloc(RHC_pck, packed_size(wf%n_so-1), 1)
!
      do while (micro_iterate .and. micro_iteration .le. solver%max_micro_iterations)
!
!        Construct stationary Roothan-Hall condition
!
         RHC = zero
         call solver%construct_stationary_Roothan_Hall_condition(wf, RHC, H, X, G, S)
!
!        Add second order correction 
!
         call solver%add_augmented_Roothan_Hall_contribution(wf, RHC, X, G)
!
!        Precondition the Roothan-Hall condition by the inverse of the linearized H matrix
!
         do i = 1, wf%n_so
            do j = 1, wf%n_so
!
               RHC(i,j) = RHC(i,j)/(H(i,i) + H(j,j))
!
            enddo
         enddo
!
!        Packin the strictly lower triangular part of the preconditioned Roothan-Hall condition
!
         RHC_pck = zero
         call packin_anti(RHC_pck, RHC, wf%n_so)
!
!        Compute the error (the L2 norm of the condition)
!
         micro_error = sqrt(ddot(packed_size(wf%n_so-1), RHC_pck, 1, RHC_pck, 1))
!
!        Ask DIIS for updated (packed) rotation parameters X
!
         X_pck = X_pck + RHC_pck
         call diis_solver%update(RHC_pck, X_pck) ! Solution placed in X_pck on exit 
!
!        Squareup anti-symmetric rotation matrix gotten from DIIS
!
         X = zero
         call squareup_anti(X_pck, X, wf%n_so)
!
         micro_iteration = micro_iteration + 1
!
         if (micro_error .lt. (solver%relative_micro_threshold)*max_grad) then
!
             micro_iterate = .false.
!
         endif
!
      enddo  
!
      call diis_solver%finalize()
!
      call mem%dealloc(X, wf%n_so, wf%n_so)
!
      call mem%dealloc(RHC, wf%n_so, wf%n_so)
      call mem%dealloc(RHC_pck, packed_size(wf%n_so-1), 1)
!
      if (micro_iterate) then 
!
         write(output%unit, '(/t3,a/)') 'Error: the Newton equations did not converge in the maximum number of micro-iterations.'
         stop
!
      endif
!
   end subroutine solve_aug_Newton_equation_arh_hf_solver
!
!
   subroutine construct_stationary_roothan_hall_condition_arh_hf_solver(solver, wf, RHC, H, X, G, S, level_shift)
!!
!!    Construct stationary Roothan-Hall condition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Sets
!!
!!       RHC = (H- mu S) X S + S X (H- mu S) + G, 
!!
!!    where mu is the level shift. Note that if similarity transformed H, S, and G 
!!    are used (Y <- V-1 Y V-T), then the iterated solution X' = V^T X V, from which 
!!    the actual X is easilyextractable.
!!
      implicit none
!
      class(arh_hf_solver), intent(in) :: solver 
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_so, wf%n_so) :: RHC
!
      real(dp), dimension(wf%n_so, wf%n_so) :: H
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: X
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: G
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: S
!
      real(dp), intent(in), optional :: level_shift
!
      real(dp), dimension(:,:), allocatable :: tmp
!
!     Construct tmp = X S => tmp^T = S^T X^T = S X^T = - S X
!
      call mem%alloc(tmp, wf%n_so, wf%n_so)
!
      call dgemm('N', 'N',       &
                  wf%n_so,       &
                  wf%n_so,       &
                  wf%n_so,       &
                  one,           &
                  X,             &
                  wf%n_so,       &
                  S,             &
                  wf%n_so,       &
                  zero,          &
                  tmp,           &
                  wf%n_so)
!
!     RHC = H X S = H tmp
!
      if (present(level_shift)) then 
!
         H = H - level_shift*S
!
      endif 
!
      call dgemm('N', 'N',       &
                  wf%n_so,       &
                  wf%n_so,       &
                  wf%n_so,       &
                  one,           &
                  H,             &
                  wf%n_so,       &
                  tmp,           &
                  wf%n_so,       &
                  zero,          &
                  RHC,           &
                  wf%n_so)
!
!     RHC = RHC - (- S X) H = RHC - tmp^T H = H X S + S X H
!
      call dgemm('T', 'N',       &
                  wf%n_so,       &
                  wf%n_so,       &
                  wf%n_so,       &
                  -one,          &
                  tmp,           &
                  wf%n_so,       &
                  H,             &
                  wf%n_so,       &
                  one,           &
                  RHC,           &
                  wf%n_so)
!
      if (present(level_shift)) then 
!
         H = H + level_shift*S
!
      endif 
!
      call mem%dealloc(tmp, wf%n_so, wf%n_so)
!
!     RHC = RHC + G = H X S + S X H + G
!
      RHC = RHC + G
!
   end subroutine construct_stationary_roothan_hall_condition_arh_hf_solver
!
!
   subroutine add_augmented_Roothan_Hall_contribution_arh_hf_solver(solver, wf, RHC, X, G)
!!
!!    Add augmented Roothan-Hall contribution 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
!!    Performs
!!
!!       RHC = RHC + sum_ij (G_i - G_n) T_ij^-1 E_j,
!!
!!    where T_ij is the trace matrix (of dimension n-1 x n-1) and E_j = Tr D_jn [D_n, X]
!!
      implicit none
!
      class(arh_hf_solver) :: solver 
!
      class(hf) :: wf 
!
      real(dp), dimension(wf%n_so, wf%n_so) :: RHC 
      real(dp), dimension(wf%n_so, wf%n_so) :: X 
      real(dp), dimension(wf%n_so, wf%n_so) :: G 
!
!
      real(dp), dimension(:, :), allocatable :: Xf
!
      real(dp), dimension(:, :), allocatable :: E 
      real(dp), dimension(:, :), allocatable :: DnX   ! [D_n, X]
      real(dp), dimension(:, :), allocatable :: D_j_n ! D_j - D_n
      real(dp), dimension(:, :), allocatable :: G_i_n ! G_i - G_n
!
      real(dp), dimension(:,:), allocatable :: inv_trace_matrix  
!
      real(dp) :: ddot
!
      integer(i15) :: i, j 
!
      if (solver%current_index .eq. 1) return ! Nothing to do in this case, i.e. no extra terms 
!
!     Pack out X, since it must be dotted with the AO density, which 
!     we have in full space 
!
      call mem%alloc(Xf, wf%n_ao, wf%n_ao)
      call symmetric_sandwich_right(Xf, X, solver%permutation_matrix, wf%n_ao, wf%n_so)
!
      call mem%alloc(E, solver%current_index - 1, 1)
      E = zero 
!
      call mem%alloc(DnX, wf%n_ao, wf%n_ao)
      call commute(wf%ao_density, Xf, DnX, wf%n_ao) ! DnX = [Dn, X]
!
      call disk%open_file(solver%AO_densities_file, 'read', 'rewind')
!
      call mem%alloc(D_j_n, wf%n_ao, wf%n_ao)
!
      do j = 1, solver%current_index - 1
!
         read(solver%AO_densities_file%unit) D_j_n
         D_j_n = D_j_n - wf%ao_density ! D_j - D_n
!
         E(j, 1) = ddot(wf%n_ao**2, D_j_n, 1, DnX, 1)
!
      enddo 
!
      call disk%close_file(solver%AO_densities_file)
!
      call mem%dealloc(D_j_n, wf%n_ao, wf%n_ao)
      call mem%dealloc(DnX, wf%n_ao, wf%n_ao)
!
      call mem%alloc(inv_trace_matrix, solver%current_index - 1, solver%current_index - 1)
      call inv(inv_trace_matrix, solver%trace_matrix, solver%current_index - 1)
!
      call disk%open_file(solver%RH_gradients_file, 'read', 'rewind')
!
      call mem%alloc(G_i_n, wf%n_ao, wf%n_ao)
!
      do i = 1, solver%current_index - 1
!
         read(solver%RH_gradients_file%unit) G_i_n
         G_i_n = G_i_n - G
!
         do j = 1, solver%current_index - 1
!
            RHC = RHC + (one/four)*G_i_n*inv_trace_matrix(i, j)*E(j, 1)
!
         enddo
      enddo
!
      call disk%close_file(solver%RH_gradients_file)
!
      call mem%dealloc(G_i_n, wf%n_ao, wf%n_ao)
      call mem%dealloc(inv_trace_matrix, solver%current_index - 1, solver%current_index - 1)
      call mem%dealloc(E, solver%current_index - 1, 1)
      call mem%dealloc(inv_trace_matrix, solver%current_index - 1, solver%current_index - 1)
!
   end subroutine add_augmented_Roothan_Hall_contribution_arh_hf_solver
!
!
   subroutine solve_level_shifted_aug_Newton_equation_arh_hf_solver(solver, wf, X_pck, H, G, S, max_grad, norm_X, level_shift)
!!
!!    Solve level shifted Newton equation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the level shifted Newton equation, written as the eigenvalue problem
!!    of the augmented Hessian (eigenvalue mu). This equation has two components,
!!
!!       H X S - S X H - mu S + gamma G = 0, (*)
!!
!!    and
!!
!!       gamma G^T X - mu = 0 (§),  ... note that G and X are here considered vectors with AO-pair index
!!
!!    where mu is the level shift and gamma is the parameter of the Hessian, H = H(gamma). 
!!    On exit, X_pck is the solution X in packed antisymmetric form. Gamma is updated 
!!    by fixed-point iteration — gamma = [trust radius / norm(X)]*gamma — until the norm 
!!    of X is within relative_trust_radius_threshold * trust radius (standard is within 
!!    10% of the trust radius).
!!
!!    The convergence threshold is max_grad*relative_micro_threshold (it becomes
!!    tighter as the maximum element of the gradient decreases) on the norm of 
!!    the equation error (*) for each given gamma.
!!
      implicit none 
!
      class(arh_hf_solver), intent(in) :: solver 
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_so-1)*(wf%n_so)/2, 1) :: X_pck
!
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: H
      real(dp), dimension(wf%n_so, wf%n_so)             :: G ! G is not changed, but it is scaled by gamma and then restored 
      real(dp), dimension(wf%n_so, wf%n_so), intent(in) :: S
!
      real(dp), intent(in) :: max_grad                       ! Maximum element of the gradient G
      real(dp)             :: norm_X                         ! Norm of the X vector (on input, from guess X; on exit, the final X norm)
!
      real(dp), dimension(:,:), allocatable :: X 
!
      real(dp), dimension(:,:), allocatable :: RHC           ! (*) H X S - S X H - mu S + gamma G (= 0 on convergence)
      real(dp), dimension(:,:), allocatable :: RHC_pck       ! Packed variant
!
      real(dp), dimension(:,:), allocatable :: G_pck         ! Packed variant of G
!
      real(dp), dimension(:,:), allocatable :: dz            ! The DIIS error vector (Newton eq. & level shift equation)
      real(dp), dimension(:,:), allocatable :: z_dz          ! The DIIS error + parameter vector
!
      logical      :: micro_iterate
      integer(i15) :: micro_iteration
      real(dp)     :: micro_error 
!
      real(dp) :: level_shift, gamma
!
      real(dp) :: projection_equation                        ! (§) gamma G^T X - mu (=0 on convergence)
!
      real(dp) :: ddot
!
      integer(i15) :: i, j
!
      type(diis) :: diis_solver 
!
      call mem%alloc(dz, packed_size(wf%n_so - 1) + 1, 1)
      call mem%alloc(z_dz, packed_size(wf%n_so - 1) + 1, 1)
!
      call mem%alloc(RHC, wf%n_so, wf%n_so)
      call mem%alloc(RHC_pck, packed_size(wf%n_so-1), 1)
!
      call mem%alloc(G_pck, packed_size(wf%n_so-1), 1)
      call packin_anti(G_pck, G, wf%n_so)
!
      call mem%alloc(X, wf%n_so, wf%n_so)
      call squareup_anti(X_pck, X, wf%n_so) 
!
      gamma = one
      do while (abs(norm_X-solver%trust_radius) .gt. solver%relative_trust_radius_threshold*solver%trust_radius)
!
         level_shift = zero
         gamma = (solver%trust_radius/norm_X)*gamma ! Guess for gamma given current norm of X
!
         call packin_anti(X_pck, X, wf%n_so)   
!
         micro_iteration = 1
         micro_iterate = .true.
         call diis_solver%init('level_shifted_hf_newton',      &
                                 packed_size(wf%n_so - 1) + 1, & 
                                 packed_size(wf%n_so - 1) + 1, &
                                 solver%diis_dimension)
!
         G = gamma*G
!
         do while (micro_iterate .and. micro_iteration .le. solver%max_micro_iterations)
!
!           Construct stationary Roothan-Hall condition
!
            RHC = zero
            call solver%construct_stationary_Roothan_Hall_condition(wf, RHC, H, X, G, S, level_shift)
!
!           Add second order correction 
!
            call solver%add_augmented_Roothan_Hall_contribution(wf, RHC, X, G)
!
!           Precondition the Roothan-Hall condition by the inverse of the linearized H matrix
!
            do i = 1, wf%n_so
               do j = 1, wf%n_so
!
                  RHC(i,j) = RHC(i,j)/(H(i,i) + H(j,j) - level_shift)
!
               enddo
            enddo
!
!           Packin the strictly lower triangular part of the preconditioned Roothan-Hall condition
!
            RHC_pck = zero
            call packin_anti(RHC_pck, RHC, wf%n_so)
!
!           Compute the projection equation, gamma * G^T X - level_shift
!
            projection_equation = ddot(packed_size(wf%n_so - 1), G_pck, 1, X_pck, 1) ! G^T X
            projection_equation = projection_equation*gamma - level_shift
!
!           Compute the error (the L2 norm of the condition)
!
            micro_error = sqrt(ddot(packed_size(wf%n_so-1), RHC_pck, 1, RHC_pck, 1) & 
                                       + projection_equation**2)                
!
!           Ask DIIS for updated (packed) rotation parameters X
!
            dz(1, 1) = projection_equation
            dz(2:(packed_size(wf%n_so - 1) + 1), 1) = RHC_pck(:, 1)
!
!           z_dz = z
!
            z_dz(1, 1) = level_shift
            z_dz(2:(packed_size(wf%n_so - 1) + 1), 1) = X_pck(:, 1)
!
!           z_dz = z_dz + dz = z + dz
!
            z_dz = z_dz + dz 
            call diis_solver%update(dz, z_dz) ! Updated z placed in z_dz on exit
!
            level_shift = z_dz(1, 1)
            X_pck(:, 1) = z_dz(2:(packed_size(wf%n_so - 1) + 1), 1)
!
!           Squareup anti-symmetric rotation matrix gotten from DIIS
!
            X = zero
            call squareup_anti(X_pck, X, wf%n_so)
!
            micro_iteration = micro_iteration + 1
!
            if (micro_error .lt. (solver%relative_shifted_micro_threshold)*max_grad) then
!
               micro_iterate = .false.
!
            endif
!
         enddo ! End of DIIS loop
!
         G = G/gamma
!
         call diis_solver%finalize()
!
         X = X/gamma
         call squareup_anti(X_pck, X, wf%n_so)
!
         norm_X = sqrt(ddot(packed_size(wf%n_so-1), X_pck, 1, X_pck, 1))
!
         write(output%unit, '(/t6,a28,i3)')     'Number of micro-iterations: ', micro_iteration - 1
         write(output%unit, '(t6,a28,f15.12)')  'Alpha:                      ', gamma ! Alpha is the name used in literature
         write(output%unit, '(t6,a28,f15.12/)') 'Level shift:                ', level_shift
         flush(output%unit)
      !   write(output%unit, '(t3,a28,f15.12/)') 'Rotation norm/trust_radius: ', abs(norm_X/solver%trust_radius)
!
         X = X*gamma ! restore for next it
!
      enddo ! End of gamma loop
!
      call mem%dealloc(dz, packed_size(wf%n_so - 1) + 1, 1)
      call mem%dealloc(z_dz, packed_size(wf%n_so - 1) + 1, 1)
!
      call mem%dealloc(RHC, wf%n_so, wf%n_so)
      call mem%dealloc(RHC_pck, packed_size(wf%n_so-1), 1)
!
      call mem%dealloc(G_pck, packed_size(wf%n_so-1), 1)
      call mem%dealloc(X, wf%n_so, wf%n_so)
!
      if (micro_iterate) then 
!
         write(output%unit, '(/t3,a,a/)') 'Error: the level shifted Newton equations did not converge ', & 
                                          'in the maximum number of micro-iterations.'
         stop
!
      endif
!
   end subroutine solve_level_shifted_aug_Newton_equation_arh_hf_solver
!
!
   subroutine construct_trace_matrix_arh_hf_solver(solver, wf, n)
!!
!!    Construct trace matrix 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the parts of T_ij = Tr (D_in D_jn) for iteration n,
!!    i.e. (i,j) = (n, :) and (:, n). Here, D_in = D_i - D_n, where D_n 
!!    is the current AO density - stored in wf%ao_density.
!!
!!    Called first time in iteration n = 2.
!!
      implicit none 
!
      class(arh_hf_solver) :: solver 
!
      class(hf) :: wf 
!
      integer(i15) :: n, j
!
      real(dp) :: ddot
!
      real(dp), dimension(:,:), allocatable :: trace_matrix_copy
!
      real(dp), dimension(:,:), allocatable :: D_nm1_n ! D_n-1, n 
      real(dp), dimension(:,:), allocatable :: D_j_n    ! D_j, n 
!
!     Preparations: allocate or expand trace matrix  
!
      if (solver%current_index .eq. 1) then ! Nothing to construct in first iteration
!
         if (n .ne. 1) call mem%dealloc(solver%trace_matrix, solver%history - 1, solver%history - 1)
         return 
!
      elseif (solver%current_index .eq. 2) then ! Just allocate and zero matrix 
!
         call mem%alloc(solver%trace_matrix, 1, 1)
         solver%trace_matrix = zero
!
      elseif (solver%current_index .gt. 2) then ! Reallocate and pad matrix with zeros to be filled
!
         call mem%alloc(trace_matrix_copy, solver%current_index-2, solver%current_index-2)
         trace_matrix_copy(:,:) = solver%trace_matrix(:,:)
!
         call mem%dealloc(solver%trace_matrix, solver%current_index-2, solver%current_index-2)
         call mem%alloc(solver%trace_matrix, solver%current_index-1, solver%current_index-1)
!
         solver%trace_matrix = zero 
         solver%trace_matrix(1:solver%current_index-2, 1:solver%current_index-2) = trace_matrix_copy(:,:)
         call mem%dealloc(trace_matrix_copy, solver%current_index-2, solver%current_index-2)
!
      endif
!
!     Read to calculate in D_n-1, n 
!
      call mem%alloc(D_nm1_n, wf%n_ao, wf%n_ao)
      call mem%alloc(D_j_n, wf%n_ao, wf%n_ao)
!
      call disk%open_file(solver%AO_densities_file, 'read', 'rewind')
!
      call solver%AO_densities_file%prepare_to_read_line(solver%current_index-1)
      read(solver%AO_densities_file%unit) D_nm1_n 
!
      D_nm1_n = D_nm1_n - wf%ao_density
!
      rewind(solver%AO_densities_file%unit)
!
!     Do dot-products
!
      do j = 1, solver%current_index - 1
!
         read(solver%AO_densities_file%unit) D_j_n
!
         D_j_n = D_j_n - wf%ao_density 
!
         solver%trace_matrix(solver%current_index-1, j) = (one/four)*ddot(wf%n_ao**2, D_nm1_n, 1, D_j_n, 1)
         solver%trace_matrix(j, solver%current_index-1) = solver%trace_matrix(solver%current_index-1, j)
!
      enddo 
!
      call disk%close_file(solver%AO_densities_file)
!
      call mem%dealloc(D_nm1_n, wf%n_ao, wf%n_ao)
      call mem%dealloc(D_j_n, wf%n_ao, wf%n_ao)
!
   end subroutine construct_trace_matrix_arh_hf_solver
!
!
end module arh_hf_solver_class
