module density_based_hf_engine_class
!
!!
!!		Density based HF engine class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use diis_solver_class
   use file_class
   use hf_class
   use array_utilities
   use disk_manager_class
   use libint_initialization
!
   implicit none
!
   type :: density_based_hf_engine
!
      integer(i15) :: max_iterations       = 150
      integer(i15) :: max_micro_iterations = 150
!
      real(dp) :: purification_threshold   = 1.0D-12
      real(dp) :: energy_threshold         = 1.0D-6
      real(dp) :: residual_threshold       = 1.0D-6
      real(dp) :: relative_micro_threshold = 1.0D-2         ! Newton equations treshold
      real(dp) :: line_search_threshold    = 1.0D-4         ! Projected gradient must be 1.0D-4 times smaller
                                                            ! than initial projected gradient (s^T g(alpha))
!
      real(dp) :: trust_radius                     = 0.1D0 
      real(dp) :: relative_trust_radius_threshold  = 0.1D0
      real(dp) :: relative_shifted_micro_threshold = 1.0D-3 ! Level shifted Newton equations threshold
!
      logical :: restart
!
   contains
!
      procedure :: initialize => initialize_density_based_hf_engine
      procedure :: solve      => solve_density_based_hf_engine
      procedure :: finalize   => finalize_density_based_hf_engine
!
      procedure :: print_banner => print_banner_density_based_hf_engine
!
      procedure :: rotate_and_purify           => rotate_and_purify_density_based_hf_engine
      procedure :: construct_and_pack_gradient => construct_and_pack_gradient_density_based_hf_engine
!
      procedure :: solve_Newton_equation               => solve_Newton_equation_density_based_hf_engine
      procedure :: solve_level_shifted_Newton_equation => solve_level_shifted_Newton_equation_density_based_hf_engine
!
   end type density_based_hf_engine
!
!
contains
!
!
   subroutine initialize_density_based_hf_engine(engine, wf)
!!
!!    Initialize density based HF engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(density_based_hf_engine) :: engine
!
      class(hf) :: wf
!
      integer(i15) :: ao
!
      real(dp), dimension(:,:), allocatable :: density_diagonal
!
      call wf%initialize_ao_density()
!
!     Set initial density to superposition of atomic densities (SOAD) guess
!
      call mem%alloc(density_diagonal, wf%n_ao, 1)
      call wf%system%SAD(wf%n_ao, density_diagonal)
!
      do ao = 1, wf%n_ao
!
         wf%ao_density(ao, ao) = density_diagonal(ao, 1)
!
      enddo
!
      call mem%dealloc(density_diagonal, wf%n_ao, 1)
!
   end subroutine initialize_density_based_hf_engine
!
!
   subroutine solve_density_based_hf_engine(engine, wf)
!!
!!    Solve
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the Hartree-Fock equations by minimizing the Roothan-Hall energy with respect
!!    to rotations of the AO density matrix (DMM), where a preconditioned conjugate-gradient
!!    (PCG) algorithm is used to determine search directions and a line search to determine
!!    the length of steps in the search direction. The uncorrected search direction, before
!!    using the conjugacy prefactor beta, is determined by solving the gradient equal to
!!    zero to first order in the rotation matrix, preconditioned by the Cholesky AO basis.
!!    These micro-iterations are performed with a direct inversion of the iterative subspace
!!    algorithm (DIIS).
!!
!!    The procedure is based heavily on the
!!
      implicit none
!
      class(density_based_hf_engine) :: engine
!
      class(hf) :: wf
!
      type(diis) :: diis_solver
!
      real(dp), dimension(:,:), allocatable :: X       ! Full rotation matrix, antisymmetric
      real(dp), dimension(:,:), allocatable :: X_pck   ! Packed rotation matrix: represents strictly lower
                                                       ! triangular part of X (excluding the diagonal)
!
      real(dp), dimension(:,:), allocatable :: RHC     ! Roothan-Hall stationary condition, full
      real(dp), dimension(:,:), allocatable :: RHC_pck ! Packed variant: represents strictly lower
                                                       ! triangular part of RHC
!
      real(dp), dimension(:,:), allocatable :: Po, Pv  ! Projection matrices
!
      real(dp), dimension(:,:), allocatable :: VT      ! Preconditioner, now = LT in S = L L^T, but may be changed
      real(dp), dimension(:,:), allocatable :: inv_VT
!
      real(dp), dimension(:,:), allocatable :: S       ! Preconditioned S (probably just the identity matrix)
      real(dp), dimension(:,:), allocatable :: H       ! Roothan-Hall Hessian
      real(dp), dimension(:,:), allocatable :: G       ! Roothan-Hall gradient
!
      real(dp), dimension(:,:), allocatable :: cur_s
      real(dp), dimension(:,:), allocatable :: prev_s
!
      real(dp), dimension(:,:), allocatable :: cur_g
      real(dp), dimension(:,:), allocatable :: prev_g
!
      real(dp), dimension(:,:), allocatable :: tmp
      real(dp), dimension(:,:), allocatable :: tmp_pck
!
      real(dp) :: norm_X 
!
      real(dp) :: p0, p1, pinit
      real(dp) :: alpha, alpha0, alpha1, alphainit, alpha_diff_length
!
      logical :: transpose_left
!
      real(dp) :: beta
      real(dp) :: beta_numerator
      real(dp) :: beta_denominator
!
      integer(i15) :: iteration = 1
!
      logical :: converged          = .false.
      logical :: converged_residual = .false.
      logical :: converged_energy   = .false.
!
      real(dp) :: energy 
      real(dp) :: prev_energy 
!
      real(dp) :: max_grad
!
      real(dp) :: ddot
!
!     Print engine banner
!
      call engine%print_banner()
!
!     Initialize engine (thresholds, etc., and wf initialization)
!
      call engine%initialize(wf)
!
!     :: Preparations for conjugate gradient loop
!
!     Construct initial AO Fock from the SOAD density
!
      call wf%initialize_ao_fock()
      call wf%construct_ao_fock() ! From current D^AO
!
!     Construct AO overlap matrix, Cholesky decompose it,
!     followed by preconditioning (making it the identity matrix
!     for this particular preconditioner - V)
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
!
      call mem%alloc(VT, wf%n_ao, wf%n_ao) ! Actual transformator
      call mem%alloc(inv_VT, wf%n_ao, wf%n_ao)
!
      VT     = zero
      inv_VT = zero
!
      call wf%decompose_ao_overlap(VT)
      VT = transpose(VT)
!
      call inv(inv_VT, VT, wf%n_ao)
!
      call mem%alloc(S, wf%n_ao, wf%n_ao)
!
      S = wf%ao_overlap
      call sandwich(S, inv_VT, inv_VT, wf%n_ao) ! S <- (V-T)^T S V-T = I
!
!     Solve Roothan Hall once - using the SOAD guess - to get a decent AO density
!     on  which to start the preconditioned conjugate gradient (PCG) algorithm
!
      call wf%initialize_orbital_energies()
      call wf%initialize_mo_coefficients()
      call wf%solve_roothan_hall() ! F^AO C = S C e to get new MOs C
!
!     Update the AO density and Fock matrices
!
      call wf%construct_ao_density() ! Construct AO density from C
      call wf%construct_ao_fock()    ! Update the AO Fock matrix
!
      call wf%destruct_mo_coefficients()
      call wf%destruct_orbital_energies()
!
!     Prepare for PCG loop by allocating matrices and zeroing them
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)            ! Projection matrices on orbital (Po) 
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)            ! and virtual (Pv) spaces
!
      Po = zero
      Pv = zero
!
      call mem%alloc(X, wf%n_ao, wf%n_ao)              ! Full antisymmetric rotation matrix
      call mem%alloc(X_pck, packed_size(wf%n_ao-1), 1) ! Packed rotation matrix (strictly lower triangular part)
!
      X     = zero
      X_pck = zero
!
      call mem%alloc(H, wf%n_ao, wf%n_ao)              ! Roothan-Hall Hessian
      call mem%alloc(G, wf%n_ao, wf%n_ao)              ! Roothan-Hall gradient
!
      H = zero
      G = zero
!
      call mem%alloc(cur_g, packed_size(wf%n_ao-1), 1)
      call mem%alloc(prev_g, packed_size(wf%n_ao-1), 1)
      call mem%alloc(prev_s, packed_size(wf%n_ao-1), 1)
      call mem%alloc(cur_s, packed_size(wf%n_ao-1), 1)
!
      cur_g  = zero
      prev_g = zero
!
      cur_s  = zero
      prev_s = zero
!
      call mem%alloc(tmp_pck, packed_size(wf%n_ao-1), 1)
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp     = zero
      tmp_pck = zero
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)           Max(gradient) '
      write(output%unit, '(t3,a)') '---------------------------------------------------'
!
      iteration = 1
      converged = .false.
!
      prev_energy = zero 
!
      do while (.not. converged .and. iteration .le. engine%max_iterations)
!
!        Construct the projection matrices Po and Pv
!
         call wf%construct_projection_matrices(Po, Pv)
!
!        Construct the Roothan-Hall Hessian H and gradient G
!
         call wf%construct_roothan_hall_hessian(H, Po, Pv)
         call wf%construct_roothan_hall_gradient(G, Po, Pv)
!
!        Precondition with V, i.e. replace matrix Y by V-1 Y V-T
!        for H and G, where S is already preconditioned before do loop
!
         call sandwich(G, inv_VT, inv_VT, wf%n_ao)
         call sandwich(H, inv_VT, inv_VT, wf%n_ao)
!
!        Pack-in gradient and determine its maximum (absolute) element
!
         call packin_anti(cur_g, G, wf%n_ao)
!
         max_grad = get_abs_max(cur_g, packed_size(wf%n_ao - 1))
!
!        Set current energy
!
         energy = wf%hf_energy
!
         converged_energy   = abs(energy-prev_energy) .lt. engine%energy_threshold
         converged_residual = max_grad                .lt. engine%residual_threshold
!
         converged = converged_energy .and. converged_residual 
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,f17.12)') iteration, wf%hf_energy, max_grad
!
         flush(output%unit)
!
!        Then test for convergence of residual
!
         if (converged) then ! Done!
!
            write(output%unit, '(t3,a)') '---------------------------------------------------'
            write(output%unit, '(/t3,a13,i3,a12/)') 'Converged in ', iteration, ' iterations!'
            converged = .true.
!
         else ! Determine next rotation of density matrix
!
!           :: Perform micro-iterations to get a direction X in which to rotate
!           (solves the equation of RH gradient equal to zero to first order in X)
!
            call engine%solve_Newton_equation(wf, X_pck, H, G, S, max_grad)
            call squareup_anti(X_pck, X, wf%n_ao)
!
!           :: Solve the augmented Hessian (i.e., level shifted Newton equation) if
!           the converged step X is longer than the trust radius where the Roothan-Hall
!           energy can be trusted to second order
!
            norm_X = sqrt(ddot(packed_size(wf%n_ao-1), X_pck, 1, X_pck, 1))
!
            if (norm_X .gt. engine%trust_radius) then
!
               write(output%unit, '(/t3,a/,t32, a)') 'Step outside trust region => solving the level shifted', &
                                                                                   'Newton equations'
!
               call engine%solve_level_shifted_Newton_equation(wf, X_pck, H, G, S, max_grad, norm_X)
               call squareup_anti(X_pck, X, wf%n_ao)
! 
            endif
!
!           :: Convert the converged direction X from the basis of the preconditioned system (X' = V^T X V)
!           back to the original system (X -> V-T X V-1); same with the gradient (G -> V G V^T)
!
            transpose_left = .false.
            call sandwich(X, inv_VT, inv_VT, wf%n_ao, transpose_left) 
!
            call packin_anti(cur_s, X, wf%n_ao) ! Save the transformed X as the uncorrected direction s
!
            call sandwich(G, VT, VT, wf%n_ao)
            call packin_anti(cur_g, G, wf%n_ao) ! Save the transformed G as the current G
!
!           :: Determine conjugacy factor beta, and adjust step direction
!
            beta_denominator = ddot(packed_size(wf%n_ao-1), prev_g, 1, prev_g, 1)
!
            if (beta_denominator .lt. 1.0D-15) then
!
               beta = zero
!
            else
!
               beta_numerator = ddot(packed_size(wf%n_ao-1), cur_g, 1, prev_g, 1) - &
                           ddot(packed_size(wf%n_ao-1), cur_g, 1, cur_g, 1)
!
               beta = beta_numerator/beta_denominator
!
            endif
!
            cur_s = cur_s - beta*prev_s
!
!           :: Perform line search to determine how far to step in the direction of s
!
!           Do linear extrapolation to get the alpha prefactor to use in front
!           of the conjugated s when rotating orbitals, k = alpha * s.
!
!           The interpolation is toward p(alpha) = s^T g(alpha) = 0, i.e.
!           the gradient is to be orthogonal to the search direction.
!
!
            alphainit = zero
            pinit     = ddot(packed_size(wf%n_ao-1), cur_s, 1, cur_g, 1) ! Projected gradient for current density
!
            alpha0 = one 
            call engine%rotate_and_purify(wf, cur_s, alpha0) 
            call wf%construct_projection_matrices(Po, Pv)
            call engine%construct_and_pack_gradient(wf, tmp_pck, Po, Pv)

            p0 = ddot(packed_size(wf%n_ao-1), cur_s, 1, tmp_pck, 1)
!
            p1    = p0
            alpha = alpha0
!
            do while (abs(p0/pinit) .gt. engine%line_search_threshold)
!
!              Adjust differentiation length,
!              where smaller lengths are used for smaller fractions
!
               alpha_diff_length = abs(p1/pinit)*(engine%line_search_threshold)

               alpha1 = alpha0 + alpha_diff_length
!
!              Rotate to get derivative
!
               call engine%rotate_and_purify(wf, cur_s, alpha1-alpha0)
               call wf%construct_projection_matrices(Po, Pv)
               call engine%construct_and_pack_gradient(wf, tmp_pck, Po, Pv)
!
               p1 = ddot(packed_size(wf%n_ao-1), cur_s, 1, tmp_pck, 1)
!
!              Do Newton step toward minimum
!
               alpha = alpha0 - p0/((p1-p0)/(alpha1-alpha0)) ! Newton step
!
!              Rotate to Newton updated alpha value
!
               call engine%rotate_and_purify(wf, cur_s, alpha-alpha1)
               call wf%construct_projection_matrices(Po, Pv)
               call engine%construct_and_pack_gradient(wf, tmp_pck, Po, Pv)
!
               p0 = ddot(packed_size(wf%n_ao-1), cur_s, 1, tmp_pck, 1)
               alpha0 = alpha
!
            enddo
!
            prev_energy = wf%hf_energy
!
!           :: Construct AO Fock (and energy) with the new rotated density
!
            call wf%construct_ao_fock()
!
!           :: Finally, set previous gradient and step direction to 
!           current, in preparation for next conjugate gradient iteration
!
            prev_g = cur_g
            prev_s = cur_s
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call wf%destruct_ao_density()
      call wf%destruct_ao_fock()
      call wf%destruct_ao_overlap()
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(X, wf%n_ao, wf%n_ao)
      call mem%dealloc(X_pck, packed_size(wf%n_ao-1), 1)
!
      call mem%dealloc(H, wf%n_ao, wf%n_ao)
      call mem%dealloc(G, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(cur_g, packed_size(wf%n_ao-1), 1)
      call mem%dealloc(prev_g, packed_size(wf%n_ao-1), 1)
      call mem%dealloc(prev_s, packed_size(wf%n_ao-1), 1)
      call mem%dealloc(cur_s, packed_size(wf%n_ao-1), 1)
!
      call mem%dealloc(tmp_pck, packed_size(wf%n_ao-1), 1)
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call engine%finalize()
!
   end subroutine solve_density_based_hf_engine
!
!
   subroutine finalize_density_based_hf_engine(engine)
!!
!! 	Finalize SCF engine
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(density_based_hf_engine) :: engine
!
   end subroutine finalize_density_based_hf_engine
!
!
   subroutine print_banner_density_based_hf_engine(engine)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(density_based_hf_engine) :: engine
!
      write(output%unit, '(/t3,a)') ':: Direct-integral density-based Hartree-Fock engine'
      write(output%unit, '(t3,a)')  ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(/t3,a)') 'This engine uses a preconditioned conjugate-gradient algorithm'
      write(output%unit, '(t3,a)')  'to minimize the gradient of the Roothan-Hall energy. In each step,'
      write(output%unit, '(t3,a)')  'the Newton equation is solved approximatively, giving a rotation of the'
      write(output%unit, '(t3,a)')  'density matrix. If a step exceeds the trust radius, an optimal step is'
      write(output%unit, '(t3,a/)') 'found on the boundary of the trust region (TR-DMM).'
!
   end subroutine print_banner_density_based_hf_engine
!
!
   subroutine rotate_and_purify_density_based_hf_engine(engine, wf, X_pck, kappa)
!!
!!    Rotate and purify
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Rotates the AO density by the antisymmetric rotation matrix X,
!!    which is sent to the routine as a packed antisymmetric matrix,
!!    and then purifies the resulting density (making it idempotent)
!!
!!    If the optional parameter kappa is passed to the routine, the 
!!    rotation is by kappa X instead of X.
!!
      implicit none 
!
      class(density_based_hf_engine), intent(in) :: engine
!
      class(hf) :: wf 
!
      real(dp), dimension((wf%n_ao-1)*(wf%n_ao)/2, 1), intent(in) :: X_pck
!
      real(dp), optional, intent(in) :: kappa
!
      real(dp), dimension(:,:), allocatable :: X 
!
      call mem%alloc(X, wf%n_ao, wf%n_ao)
      call squareup_anti(X_pck, X, wf%n_ao)
!
      if (present(kappa)) then 
!
         X = kappa*X
!
      endif
!
      call wf%rotate_ao_density(X)
      call wf%purify_ao_density(engine%purification_threshold)
!
      call mem%dealloc(X, wf%n_ao, wf%n_ao)
!
   end subroutine rotate_and_purify_density_based_hf_engine
!
!
   subroutine construct_and_pack_gradient_density_based_hf_engine(engine, wf, G_pck, Po, Pv)
!!
!!    Construct and pack Roothan-Hall gradient
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall gradient G, then packs it into G_pck.
!!
      implicit none 
!
      class(density_based_hf_engine), intent(in) :: engine 
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
   end subroutine construct_and_pack_gradient_density_based_hf_engine
!
!
   subroutine solve_Newton_equation_density_based_hf_engine(engine, wf, X_pck, H, G, S, max_grad)
!!
!!    Solve Newton equation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the Newton equation using a DIIS algorithm and no 
!!    level shift:
!!
!!       H X S - S X H = -G.
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
      class(density_based_hf_engine), intent(in) :: engine 
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_ao-1)*(wf%n_ao)/2, 1) :: X_pck
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: H
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: G
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: S
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
      call mem%alloc(X, wf%n_ao, wf%n_ao)
!
      X = -G
      call packin_anti(X_pck, X, wf%n_ao)
!
      micro_iteration = 1
      micro_iterate = .true.
      call diis_solver%init('hf_newton', packed_size(wf%n_ao - 1), 10)
!
      call mem%alloc(RHC, wf%n_ao, wf%n_ao)
      call mem%alloc(RHC_pck, packed_size(wf%n_ao-1), 1)
!
      do while (micro_iterate .and. micro_iteration .le. engine%max_micro_iterations)
!
!        Construct stationary Roothan-Hall condition
!
         RHC = zero
         call wf%construct_stationary_Roothan_Hall_condition(RHC, H, X, G, S)
!
!        Precondition the Roothan-Hall condition by the inverse of the linearized H matrix
!
         do i = 1, wf%n_ao
            do j = 1, wf%n_ao
!
               RHC(i,j) = RHC(i,j)/(H(i,i) + H(j,j))
!
            enddo
         enddo
!
!        Packin the strictly lower triangular part of the preconditioned Roothan-Hall condition
!
         RHC_pck = zero
         call packin_anti(RHC_pck, RHC, wf%n_ao)
!
!        Compute the error (the L2 norm of the condition)
!
         micro_error = sqrt(ddot(packed_size(wf%n_ao-1), RHC_pck, 1, RHC_pck, 1))
!
!        Ask DIIS for updated (packed) rotation parameters X
!
         X_pck = X_pck + RHC_pck
         call diis_solver%update(RHC_pck, X_pck)
         X_pck = RHC_pck
!
!        Squareup anti-symmetric rotation matrix gotten from DIIS
!
         X = zero
         call squareup_anti(X_pck, X, wf%n_ao)
!
         micro_iteration = micro_iteration + 1
!
         if (micro_error .lt. (engine%relative_micro_threshold)*max_grad) then
!
             micro_iterate = .false.
!
         endif
!
      enddo  
!
      call diis_solver%finalize()
!
      call mem%dealloc(X, wf%n_ao, wf%n_ao)
!
      if (micro_iterate) then 
!
         write(output%unit, '(/t3,a/)') 'Error: the Newton equations did not converge in the maximum number of micro-iterations.'
         stop
!
      endif
!
   end subroutine solve_Newton_equation_density_based_hf_engine    
!
!
   subroutine solve_level_shifted_Newton_equation_density_based_hf_engine(engine, wf, X_pck, H, G, S, max_grad, norm_X)
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
      class(density_based_hf_engine), intent(in) :: engine 
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_ao-1)*(wf%n_ao)/2, 1) :: X_pck
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: H
      real(dp), dimension(wf%n_ao, wf%n_ao)             :: G ! G is not changed, but it is scaled by gamma and then restored 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: S
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
      call mem%alloc(dz, packed_size(wf%n_ao - 1) + 1, 1)
      call mem%alloc(z_dz, packed_size(wf%n_ao - 1) + 1, 1)
!
      call mem%alloc(RHC, wf%n_ao, wf%n_ao)
      call mem%alloc(RHC_pck, packed_size(wf%n_ao-1), 1)
!
      call mem%alloc(G_pck, packed_size(wf%n_ao-1), 1)
      call packin_anti(G_pck, G, wf%n_ao)
!
      call mem%alloc(X, wf%n_ao, wf%n_ao)
      call squareup_anti(X_pck, X, wf%n_ao) 
!
      gamma = one
      do while (abs(norm_X-engine%trust_radius) .gt. engine%relative_trust_radius_threshold*engine%trust_radius)
!
         level_shift = zero
         gamma = (engine%trust_radius/norm_X)*gamma ! Guess for gamma given current norm of X
!
         call packin_anti(X_pck, X, wf%n_ao)   
!
         micro_iteration = 1
         micro_iterate = .true.
         call diis_solver%init('level_shifted_hf_newton', packed_size(wf%n_ao - 1) + 1, 10)
!
         G = gamma*G
!
         do while (micro_iterate .and. micro_iteration .le. engine%max_micro_iterations)
!
!           Construct stationary Roothan-Hall condition
!
            RHC = zero
            call wf%construct_stationary_Roothan_Hall_condition(RHC, H, X, G, S, level_shift)
!
!           Precondition the Roothan-Hall condition by the inverse of the linearized H matrix
!
            do i = 1, wf%n_ao
               do j = 1, wf%n_ao
!
                  RHC(i,j) = RHC(i,j)/(H(i,i) + H(j,j) - level_shift)
!
               enddo
            enddo
!
!           Packin the strictly lower triangular part of the preconditioned Roothan-Hall condition
!
            RHC_pck = zero
            call packin_anti(RHC_pck, RHC, wf%n_ao)
!
!           Compute the projection equation, gamma * G^T X - level_shift
!
            projection_equation = ddot(packed_size(wf%n_ao - 1), G_pck, 1, X_pck, 1) ! G^T X
            projection_equation = projection_equation*gamma - level_shift
!
!           Compute the error (the L2 norm of the condition)
!
            micro_error = sqrt(ddot(packed_size(wf%n_ao-1), RHC_pck, 1, RHC_pck, 1) & 
                                       + projection_equation**2)                
!
!           Ask DIIS for updated (packed) rotation parameters X
!
            dz(1, 1) = projection_equation
            dz(2:(packed_size(wf%n_ao - 1) + 1), 1) = RHC_pck(:, 1)
!
!           z_dz = z
!
            z_dz(1, 1) = level_shift
            z_dz(2:(packed_size(wf%n_ao - 1) + 1), 1) = X_pck(:, 1)
!
!           z_dz = z_dz + dz = z + dz
!
            z_dz = z_dz + dz 
            call diis_solver%update(dz, z_dz) ! Updated z placed in dz on exit
!
            level_shift = dz(1, 1)
            X_pck(:, 1) = dz(2:(packed_size(wf%n_ao - 1) + 1), 1)
!
!           Squareup anti-symmetric rotation matrix gotten from DIIS
!
            X = zero
            call squareup_anti(X_pck, X, wf%n_ao)
!
            micro_iteration = micro_iteration + 1
!
            if (micro_error .lt. (engine%relative_shifted_micro_threshold)*max_grad) then
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
         call squareup_anti(X_pck, X, wf%n_ao)
!
         norm_X = sqrt(ddot(packed_size(wf%n_ao-1), X_pck, 1, X_pck, 1))
!
        ! write(output%unit, '(/t3,a28,i3)')    'Number of micro-iterations: ', micro_iteration
         write(output%unit, '(/t6,a13,f15.12)') 'Gamma:       ', gamma
         write(output%unit, '(t6,a13,f15.12/)') 'Level shift: ', level_shift
         !write(output%unit, '(t3,a28,f15.12/)') 'Rotation norm/trust_radius: ', abs(norm_X/engine%trust_radius)
!
         X = X*gamma ! restore for next it
!
      enddo ! End of gamma loop
!
      call mem%dealloc(dz, packed_size(wf%n_ao - 1) + 1, 1)
      call mem%dealloc(z_dz, packed_size(wf%n_ao - 1) + 1, 1)
!
      call mem%dealloc(RHC, wf%n_ao, wf%n_ao)
      call mem%dealloc(RHC_pck, packed_size(wf%n_ao-1), 1)
!
      call mem%dealloc(G_pck, packed_size(wf%n_ao-1), 1)
      call mem%dealloc(X, wf%n_ao, wf%n_ao)
!
      if (micro_iterate) then 
!
         write(output%unit, '(/t3,a,a/)') 'Error: the level shifted Newton equations did not converge ', & 
                                          'in the maximum number of micro-iterations.'
         stop
!
      endif
!
   end subroutine solve_level_shifted_Newton_equation_density_based_hf_engine   
!
!
end module density_based_hf_engine_class
