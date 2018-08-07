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
      real(dp) :: purification_threshold = 1.0D-12
!
      integer(i15) :: max_iterations = 150
!
      real(dp) :: energy_threshold   = 1.0D-6
      real(dp) :: residual_threshold = 1.0D-6
!
      integer(i15) :: max_micro_iterations = 150
!
      real(dp) :: relative_micro_threshold = 1.0D-2  ! Error below 1.0D-2*(maximum gradient element)
!
      real(dp) :: line_search_threshold = 1.0D-4 ! Projected gradient must be 1.0D-4 times smaller
                                                 ! than initial projected gradient (s^T g(alpha))
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
      call wf%system%SOAD(wf%n_ao, density_diagonal)
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
      real(dp), dimension(:,:), allocatable :: X     ! Full rotation matrix, antisymmetric
      real(dp), dimension(:,:), allocatable :: X_pck ! Packed rotation matrix: represents strictly lower
                                                     ! triangular part of X (excluding the diagonal)
!
      real(dp), dimension(:,:), allocatable :: RHC     ! Roothan-Hall stationary condition, full
      real(dp), dimension(:,:), allocatable :: RHC_pck ! Packed variant: represents strictly lower
                                                       ! triangular part of RHC
!
      real(dp), dimension(:,:), allocatable :: Po, Pv  ! Projection matrices
!
      real(dp), dimension(:,:), allocatable :: H       ! Roothan-Hall Hessian
!
      real(dp), dimension(:,:), allocatable :: G       ! Roothan-Hall gradient
      real(dp), dimension(:,:), allocatable :: G_pck   ! Packed variant: represents strictly lower
                                                       ! triangular part of G
!
      real(dp), dimension(:,:), allocatable :: LT      ! L^T,  where S = L L^T
      real(dp), dimension(:,:), allocatable :: inv_LT  ! L^-T, where S = L L^T
!
      real(dp), dimension(:,:), allocatable :: VT      ! Preconditioner, now = LT, but may be changed
      real(dp), dimension(:,:), allocatable :: inv_VT
!
      real(dp), dimension(:,:), allocatable :: S       ! Preconditioned S (probably just the identity matrix)
!
      real(dp), dimension(:,:), allocatable :: cur_s
      real(dp), dimension(:,:), allocatable :: prev_s
      real(dp), dimension(:,:), allocatable :: cur_g
      real(dp), dimension(:,:), allocatable :: prev_g
!
      real(dp), dimension(:,:), allocatable :: tmp
      real(dp), dimension(:,:), allocatable :: tmp_pck
!
      real(dp) :: p0, p1, pinit
      real(dp) :: alpha, alpha0, alpha1, alphainit, alpha_diff_length
!
      real(dp) :: beta
      real(dp) :: beta_numerator
      real(dp) :: beta_denominator
!
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: iteration = 1
!
      real(dp)     :: micro_error
      integer(i15) :: micro_iteration = 1
      logical      :: micro_iterate = .false.
!
      logical :: converged = .false.
!
      real(dp) :: gradient_norm
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
      iteration = 1
      converged = .false.
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao) ! Projection matrices on orbital and virtual spaces
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
!
      Po = zero
      Pv = zero
!
      call mem%alloc(X, wf%n_ao, wf%n_ao)              ! Full antisymmetric rotation matrix
      call mem%alloc(X_pck, packed_size(wf%n_ao-1), 1) ! Packed rotation matrix (strictly lower triangular part)
!
      X = zero
      X_pck = zero
!
      call mem%alloc(RHC, wf%n_ao, wf%n_ao)              ! Full Roothan-Hall stationary condition matrix
      call mem%alloc(RHC_pck, packed_size(wf%n_ao-1), 1) ! Packed variant (strictly lower triangular part)
!
      RHC = zero
      RHC_pck = zero
!
      call mem%alloc(H, wf%n_ao, wf%n_ao)              ! Roothan-Hall Hessian
      call mem%alloc(G, wf%n_ao, wf%n_ao)              ! Roothan-Hall gradient
      call mem%alloc(G_pck, packed_size(wf%n_ao-1), 1) ! Packed gradient
!
      H = zero
!
      G = zero
      G_pck = zero
!
      call mem%alloc(cur_g, packed_size(wf%n_ao-1), 1)
      call mem%alloc(prev_g, packed_size(wf%n_ao-1), 1)
      call mem%alloc(prev_s, packed_size(wf%n_ao-1), 1)
      call mem%alloc(cur_s, packed_size(wf%n_ao-1), 1)
!
      call mem%alloc(tmp_pck, packed_size(wf%n_ao-1), 1)
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      cur_g = zero
      prev_g = zero
      cur_s = zero
      prev_s = zero
!
      tmp = zero
      tmp_pck = zero
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)           Max(gradient) '
      write(output%unit, '(t3,a)') '---------------------------------------------------'
!
      converged = .false.
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
!        Pack-in gradient and determine its maximum element,
!        which is used to check for convergence
!
         call packin_anti(G_pck, G, wf%n_ao)
!
         max_grad = zero
         do I = 1, packed_size(wf%n_ao - 1)
!
            if (abs(G_pck(I, 1)) .gt. max_grad) max_grad = abs(G_pck(I, 1))
!
         enddo
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,f17.12)') iteration, wf%hf_energy, max_grad
!
         flush(output%unit)
!
!        Then test for convergence of residual
!
         if (max_grad .lt. engine%residual_threshold) then ! Done!
!
            write(output%unit, '(t3,a)') '---------------------------------------------------'
            write(output%unit, '(/t3,a13,i3,a12/)') 'Converged in ', iteration, ' iterations!'
            converged = .true.
!
         else ! Determine next rotation of density matrix
!
!           Perform micro-iterations to get a direction X in which to rotate
!           (solves the equation of RH gradient equal to zero to first order in X)
!
            X = -G
            call packin_anti(X_pck, X, wf%n_ao)
!
            micro_iteration = 1
            micro_iterate = .true.
            call diis_solver%init('hf_newton', packed_size(wf%n_ao - 1), 25)
!
            do while (micro_iterate .and. micro_iteration .le. engine%max_micro_iterations)
!
!              Construct stationary Roothan-Hall condition
!
               RHC = zero
               call wf%construct_stationary_Roothan_Hall_condition(RHC, H, X, G, S)
!
!              Precondition the Roothan-Hall condition by the inverse of the linearized H matrix
!
               do i = 1, wf%n_ao
                  do j = 1, wf%n_ao
!
                     RHC(i,j) = RHC(i,j)/(H(i,i) + H(j,j))
!
                  enddo
               enddo
!
!              Packin the strictly lower triangular part of the preconditioned Roothan-Hall condition
!
               RHC_pck = zero
               call packin_anti(RHC_pck, RHC, wf%n_ao)
!
!              Compute the error (the L2 norm of the condition)
!
               micro_error = sqrt(ddot(packed_size(wf%n_ao-1), RHC_pck, 1, RHC_pck, 1))
!
!              Ask DIIS for updated (packed) rotation parameters X
!
               X_pck = X_pck + RHC_pck
               call diis_solver%update(RHC_pck, X_pck)
               X_pck = RHC_pck
!
!              Squareup anti-symmetric rotation matrix gotten from DIIS
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
            inv_VT = transpose(inv_VT) ! inv_V
            call sandwich(X, inv_VT, inv_VT, wf%n_ao) ! X <- V-T X' V, where X' = V^T X V
            call packin_anti(cur_s, X, wf%n_ao)
            inv_VT = transpose(inv_VT) ! inv_VT
!
            call wf%construct_roothan_hall_gradient(G, Po, Pv)
            call packin_anti(cur_g, G, wf%n_ao)
!
!           Compute the conjugacy factor beta (to add previous gradient contribution)
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
            alphainit = zero
            pinit = ddot(packed_size(wf%n_ao-1), cur_s, 1, cur_g, 1) ! Projected gradient for current density
!
            alpha0 = one ! Initial guess for optimal step length
            tmp_pck = alpha0*cur_s
            call squareup_anti(tmp_pck, tmp, wf%n_ao)
            call wf%rotate_ao_density(tmp)
            call wf%purify_ao_density(engine%purification_threshold)
            call wf%construct_projection_matrices(Po, Pv)
            call wf%construct_roothan_hall_gradient(tmp, Po, Pv)
            call packin_anti(tmp_pck, tmp, wf%n_ao)
            p0 = ddot(packed_size(wf%n_ao-1), cur_s, 1, tmp_pck, 1)
!
            p1 = p0 ! Initial value of projected gradient
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
               tmp_pck = (alpha1-alpha0)*cur_s
               call squareup_anti(tmp_pck, tmp, wf%n_ao)
               call wf%rotate_ao_density(tmp)
               call wf%purify_ao_density(engine%purification_threshold)
               call wf%construct_projection_matrices(Po, Pv)
               call wf%construct_roothan_hall_gradient(tmp, Po, Pv)
               call packin_anti(tmp_pck, tmp, wf%n_ao)
!
               p1 = ddot(packed_size(wf%n_ao-1), cur_s, 1, tmp_pck, 1)
!
!              Do Newton step toward minimum
!
               alpha = alpha0 - p0/((p1-p0)/(alpha1-alpha0)) ! Newton step
!
!              Rotate to Newton updated alpha value
!
               tmp_pck = (alpha-alpha1)*cur_s
               call squareup_anti(tmp_pck, tmp, wf%n_ao)
               call wf%rotate_ao_density(tmp)
               call wf%purify_ao_density(engine%purification_threshold)
               call wf%construct_projection_matrices(Po, Pv)
               call wf%construct_roothan_hall_gradient(tmp, Po, Pv)
               call packin_anti(tmp_pck, tmp, wf%n_ao)
!
               p0 = ddot(packed_size(wf%n_ao-1), cur_s, 1, tmp_pck, 1)
               alpha0 = alpha
!
            enddo
!
            prev_g = cur_g
            prev_s = cur_s
!
!           Construct AO Fock with the new rotated density
!
            call wf%construct_ao_fock()
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
      call mem%dealloc(G_pck, packed_size(wf%n_ao-1), 1)
!
      call mem%dealloc(RHC, wf%n_ao, wf%n_ao)
      call mem%dealloc(RHC_pck, packed_size(wf%n_ao-1), 1)
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
      write(output%unit, '(t3,a/)')  'to minimize the gradient of the Roothan-Hall energy.'
!
      write(output%unit, '(t3,a/)')  'Will soon include trust-radius to avoid premature large rotations.'
!
   end subroutine print_banner_density_based_hf_engine
!
!
end module density_based_hf_engine_class
