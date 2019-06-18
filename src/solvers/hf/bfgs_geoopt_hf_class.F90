module bfgs_geoopt_hf_class
!
!!
!!    BFGS geometry optimization for Hartree-Fock
!!    Written by Eirik F. Kjønstad, June 2019
!!
!
   use hf_class 
   use hf_engine_class
   use bfgs_tool_class
!
   implicit none
!
   type :: bfgs_geoopt_hf
!
      character(len=400) :: tag
      character(len=400) :: author
      character(len=400) :: description
!
      integer :: iteration 
!
      integer :: max_iterations
!
      real(dp) :: energy_threshold
      real(dp) :: gradient_threshold
!
      type(hf_engine) :: hf_gs_engine
!
      real(dp), dimension(:), allocatable :: energies, gradient_norms 
!
   contains
!
      procedure :: run                 => run_bfgs_geoopt_hf
      procedure :: cleanup             => cleanup_bfgs_geoopt_hf
!
      procedure :: line_search         => line_search_bfgs_geoopt_hf
!
      procedure :: read_settings       => read_settings_bfgs_geoopt_hf
      procedure :: print_banner        => print_banner_bfgs_geoopt_hf
      procedure :: print_summary       => print_summary_bfgs_geoopt_hf
      procedure :: determine_gradient  => determine_gradient_bfgs_geoopt_hf
!
   end type bfgs_geoopt_hf
!
!
   interface bfgs_geoopt_hf 
!
      procedure :: new_bfgs_geoopt_hf
!
   end interface bfgs_geoopt_hf
!
!
contains
!
!
   function new_bfgs_geoopt_hf() result(solver)
!!
!!    New BFGS geoopt HF 
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(bfgs_geoopt_hf) :: solver
!
!     Print solver banner
!
      solver%tag     = 'BFGS geometry optimization solver'
      solver%author  = 'Eirik F. Kjønstad, 2019'
!
      solver%description = 'A BFGS solver.'
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%max_iterations      = 250
      solver%energy_threshold    = 1.0d-6
      solver%gradient_threshold  = 1.0d-6
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call mem%alloc(solver%energies, solver%max_iterations)
      call mem%alloc(solver%gradient_norms, solver%max_iterations)
!
      solver%hf_gs_engine = hf_engine()
!
   end function new_bfgs_geoopt_hf
!
!
   subroutine read_settings_bfgs_geoopt_hf(solver)
!!
!!    Read settings 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
      implicit none 
!
      class(bfgs_geoopt_hf), intent(inout) :: solver 
!
      call input%get_keyword_in_section('gradient threshold', &
                                        'solver geometry optimization', solver%gradient_threshold)
!
      call input%get_keyword_in_section('energy threshold', &
                                        'solver geometry optimization', solver%gradient_threshold)
!
      call input%get_keyword_in_section('max iterations', &
                                        'solver geometry optimization', solver%max_iterations)
!
   end subroutine read_settings_bfgs_geoopt_hf
!
!
   function line_search_bfgs_geoopt_hf(solver, wf, descent_direction) result(alpha)
!!
!!    Line search 
!!    Written by Eirik F. Kjønstad, June 2019
!!
      implicit none 
!
      class(bfgs_geoopt_hf) :: solver 
!
      class(hf) :: wf 
!
      real(dp), dimension(3, wf%system%n_atoms) :: descent_direction
!
      real(dp) :: alpha 
!
      real(dp) :: gradient_norm_a, gradient_norm_b, gradient_norm_c, gradient_norm_d, a, b, c, d 
!
      real(dp), dimension(3, wf%system%n_atoms) :: geometry_0, geometry, gradient  
!
      real(dp), parameter :: phi = (one + sqrt(five))/two
!
!     Initialize by finding the first bracket [a,b]
!
      geometry_0 = wf%system%get_geometry()
!
      a = half
      geometry = geometry_0 + a*descent_direction
      gradient = solver%determine_gradient(wf, geometry)
      gradient_norm_a = get_l2_norm(gradient, 3*wf%system%n_atoms)
!
      b = two
      geometry = geometry_0 + b*descent_direction
      gradient = solver%determine_gradient(wf, geometry)
      gradient_norm_b = get_l2_norm(gradient, 3*wf%system%n_atoms)
! 
      c = b - (b - a)/phi 
      d = a + (b - a)/phi 
!
      do while (abs(c-d) >= 0.05d0)
!
         call output%printf('a = (f5.2), b = (f5.2)', reals=[a,b])
!
         geometry = geometry_0 + c*descent_direction
         gradient = solver%determine_gradient(wf, geometry)
         gradient_norm_c = get_l2_norm(gradient, 3*wf%system%n_atoms)
!
         geometry = geometry_0 + d*descent_direction
         gradient = solver%determine_gradient(wf, geometry)
         gradient_norm_d = get_l2_norm(gradient, 3*wf%system%n_atoms)
!
         call output%printf('gradient norm 0 = (f19.12)', reals=[solver%gradient_norms(solver%iteration)])
         call output%printf('gradient norm c = (f19.12)', reals=[gradient_norm_c])
         call output%printf('gradient norm d = (f19.12)', reals=[gradient_norm_d])
!
         if (gradient_norm_c < gradient_norm_d) then 
!
            b = d
!
         else
!
            a = c
!
         endif
!
         c = b - (b - a)/phi 
         d = a + (b - a)/phi 
!
      enddo  
!
      alpha = (a + b)/two
!
      write(output%unit, *) 'alpha:::', alpha 
   !   if (abs(alpha) < half) alpha = one
!
   end function line_search_bfgs_geoopt_hf
!    function line_search_bfgs_geoopt_hf(solver, wf, geometry, direction) result(alpha)
! !!
! !!    Line search 
! !!    Written by Eirik F. Kjønstad, June 2019
! !!
!       implicit none 
! !
!       class(bfgs_geoopt_hf) :: solver 
! !
!       class(hf) :: wf 
! !
!       real(dp) :: a, b, c, fa, fb, fc 
! !
!       real(dp), dimension(3, wf%system%n_atoms) :: x, g 
! !
! !     Compute three points 
! !
!       a = -one
!       x = geometry + a*direction
!       g = solver%determine_gradient(wf, x)
!       fa = get_l2_norm(g, 3*wf%system%n_atoms)
! !
!       b = zero 
!       fb = solver%gradient_norms(solver%iteration)
! !
!       c = one 
!       x = geometry + c*direction 
!       g = solver%determine_gradient(wf, x)
!       fc = get_l2_norm(g, 3*wf%system%n_atoms) 
! !
! !     Do cubic interpolation & determine minimum 
! !

! !
!    end function line_search_bfgs_geoopt_hf
!
!
   function determine_gradient_bfgs_geoopt_hf(solver, wf, geometry) result(gradient)
!!
!!    Determine gradient 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(bfgs_geoopt_hf) :: solver 
!  
      class(hf) :: wf 
!
      real(dp), dimension(3, wf%system%n_atoms) :: geometry, gradient 
!
      call wf%system%set_geometry(geometry)
!
      call wf%set_n_mo()                  ! Decomposes AO overlap (linear dependence)
      call solver%hf_gs_engine%ignite(wf) ! Converges HF orbitals/density
!
      call wf%construct_molecular_gradient(gradient)      
!
   end function determine_gradient_bfgs_geoopt_hf
!
!
   subroutine run_bfgs_geoopt_hf(solver, wf)
!!
!!    Run 
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(bfgs_geoopt_hf) :: solver
!
      class(hf) :: wf
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_gradient
!
      real(dp) :: energy, prev_energy, norm_gradient, alpha
!
      real(dp), dimension(3,wf%system%n_atoms) :: gradient
      real(dp), dimension(3,wf%system%n_atoms) :: bfgs_direction
      real(dp), dimension(3,wf%system%n_atoms) :: geometry
!
      type(bfgs_tool) :: bfgs 
!
      bfgs = bfgs_tool(3*wf%system%n_atoms)
!
      solver%iteration = 0
!
      converged          = .false.
      converged_energy   = .false.
      converged_gradient = .false.
!
      energy = zero
      prev_energy = zero
!
      geometry = wf%system%get_geometry()
!
      do while (.not. converged .and. solver%iteration <= solver%max_iterations)        
!
         solver%iteration = solver%iteration + 1
!
         gradient = solver%determine_gradient(wf, geometry)
!
         energy = wf%energy 
         norm_gradient = get_l2_norm(gradient, 3*wf%system%n_atoms)
!
         solver%energies(solver%iteration) = energy 
         solver%gradient_norms(solver%iteration) = norm_gradient
!
         converged_gradient = norm_gradient              <= solver%gradient_threshold
         converged_energy   = abs(energy - prev_energy)  <= solver%energy_threshold
!
         converged = converged_gradient .and. converged_energy
!
         if (converged) then
!
            write(output%unit, '(/t3,a,i0,a)') 'Converged geometry in ', solver%iteration, ' iterations!'
!
         else
!
!           Update Hessian 
!
            call bfgs%update_hessian(geometry, gradient)
!
!           Keep a copy of the current geometry and gradient (needed for BFGS)
!
            call bfgs%set_previous_geometry_and_gradient(geometry, gradient)
!
!           Get next search direction (according to the current BFGS Hessian)
!
            call bfgs%get_direction(gradient, bfgs_direction)
!
!           Do a line search in that direction,
!           which returns alpha (= how long the step should be)         
!
            alpha = solver%line_search(wf, bfgs_direction)
!
!           Update geometry
!
            geometry = geometry + alpha*bfgs_direction 
!
         endif
!
         call wf%system%print_geometry()
         prev_energy = energy 
!
      enddo
!
      if (.not. converged) then 
!
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         stop
!
      else
!
         call solver%print_summary()
!
      endif  
!
   end subroutine run_bfgs_geoopt_hf
!
!
   subroutine print_banner_bfgs_geoopt_hf(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(bfgs_geoopt_hf) :: solver 
!
      call output%long_string_print(solver%tag,'(//t3,a)',.true.)
      call output%long_string_print(solver%author,'(t3,a/)',.true.)
      call output%long_string_print(solver%description)
!
   end subroutine print_banner_bfgs_geoopt_hf
!
!
   subroutine print_summary_bfgs_geoopt_hf(solver)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(bfgs_geoopt_hf), intent(in) :: solver 
!
      integer :: iteration 
!
      write(output%unit, '(/t3,a)') '- Summary of geometry optimization iterations: '
!
      write(output%unit, '(/t6,a)') 'Iteration       Energy                Gradient norm       '
      write(output%unit, '(t6,a)')  '----------------------------------------------------------'
!
      do iteration = 1, solver%iteration 
!
         call output%printf('(i4)         (f19.12)     (e11.4)', &
                                    ints=[iteration], &
                                    reals=[solver%energies(iteration), solver%gradient_norms(iteration)], &
                                    fs='(t6,a)')
!
      enddo 
!
      write(output%unit, '(t6,a)')  '----------------------------------------------------------'
!
   end subroutine print_summary_bfgs_geoopt_hf
!
!
   subroutine cleanup_bfgs_geoopt_hf(solver)
!!
!!    Cleanup 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(bfgs_geoopt_hf) :: solver 
!
      call mem%dealloc(solver%energies, solver%max_iterations)
      call mem%dealloc(solver%gradient_norms, solver%max_iterations)
!
   end subroutine cleanup_bfgs_geoopt_hf
!
!
end module bfgs_geoopt_hf_class