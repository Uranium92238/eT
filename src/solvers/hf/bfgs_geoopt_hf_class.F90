module bfgs_geoopt_hf_class
!
!!
!!    BFGS geometry optimization Hartree-Fock solver
!!    Written by Eirik F. Kjønstad, June 2019
!!
!
   use hf_class 
   use hf_engine_class
   use bfgs_tool_class
   use diis_tool_class
!
   implicit none
!
   type :: bfgs_geoopt_hf
!
      character(len=400) :: tag
      character(len=400) :: author
      character(len=400) :: description
!
      logical :: restart 
!
      integer :: iteration 
!
      integer :: max_iterations
!
      real(dp) :: energy_threshold
      real(dp) :: gradient_threshold
      real(dp) :: max_step 
!
      type(hf_engine) :: hf_gs_engine
!
      real(dp), dimension(:), allocatable :: energies, gradient_maxs 
!
   contains
!
      procedure :: run                 => run_bfgs_geoopt_hf
      procedure :: cleanup             => cleanup_bfgs_geoopt_hf
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
   function new_bfgs_geoopt_hf(restart) result(solver)
!!
!!    New BFGS geoopt HF 
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(bfgs_geoopt_hf) :: solver
!
      logical, intent(in) :: restart 
!
!     Print solver banner
!
      solver%tag     = 'BFGS geometry optimization solver'
      solver%author  = 'Eirik F. Kjønstad, 2019'
!
      solver%description = 'Constructs an approximate Hessian using the BFGS algorithm &
                           &using the previous geometries and gradients. From the BFGS Hessian, &
                           &a level shift given by the rational function (RF) augmented Hessian &
                           &is applied. See J. Comput. Chem. 18: 1473-1483, 1997.'
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%max_iterations      = 250
      solver%energy_threshold    = 1.0d-4
      solver%gradient_threshold  = 1.0d-4
      solver%max_step            = 0.5d0
      solver%restart             = restart
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call mem%alloc(solver%energies, solver%max_iterations)
      call mem%alloc(solver%gradient_maxs, solver%max_iterations)
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
                                        'solver hf geoopt', solver%gradient_threshold)
!
      call input%get_keyword_in_section('energy threshold', &
                                        'solver hf geoopt', solver%energy_threshold)
!
      call input%get_keyword_in_section('max iterations', &
                                        'solver hf geoopt', solver%max_iterations)
!
      call input%get_keyword_in_section('max step', &
                                        'solver hf geoopt', solver%max_step)
!
   end subroutine read_settings_bfgs_geoopt_hf
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
!     Update geometry to the one requested
!
      call wf%system%set_geometry(geometry)
!
!     Re-decompose the AO overlap 
!
      call wf%set_n_mo() 
!
!     Attempt to converge HF orbitals/density, using restart 
!
      if (solver%restart .or. solver%iteration > 1) solver%hf_gs_engine%restart = .true.
      call solver%hf_gs_engine%ignite(wf) 
!
!     Compute gradient 
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
      real(dp) :: energy, prev_energy, max_gradient
!
      real(dp), dimension(3,wf%system%n_atoms) :: gradient
      real(dp), dimension(3,wf%system%n_atoms) :: bfgs_direction
      real(dp), dimension(3,wf%system%n_atoms) :: geometry
!
      type(bfgs_tool) :: bfgs 
!
      bfgs = bfgs_tool(3*wf%system%n_atoms, solver%max_step)
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
         max_gradient = get_abs_max(gradient, 3*wf%system%n_atoms)
!
         solver%energies(solver%iteration) = energy 
         solver%gradient_maxs(solver%iteration) = max_gradient
!
         converged_gradient = max_gradient               <= solver%gradient_threshold
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
            call bfgs%update_hessian(geometry, gradient)
            call bfgs%get_direction(gradient, bfgs_direction)
!
            geometry = geometry + bfgs_direction 
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
                                    reals=[solver%energies(iteration), solver%gradient_maxs(iteration)], &
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
      call mem%dealloc(solver%gradient_maxs, solver%max_iterations)
!
   end subroutine cleanup_bfgs_geoopt_hf
!
!
end module bfgs_geoopt_hf_class