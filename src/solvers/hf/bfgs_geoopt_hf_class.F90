module bfgs_geoopt_hf_class
!
!!
!!    BFGS geometry optimization Hartree-Fock solver
!!    Written by Eirik F. Kjønstad, June 2019
!!
!
   use parameters
!
   use global_out,            only: output
   use global_in,             only: input
   use memory_manager_class,  only: mem
!
   use hf_class,              only: hf
   use scf_diis_hf_class,     only: scf_diis_hf
   use bfgs_tool_class,       only: bfgs_tool
   use timings_class,         only: timings
!
   implicit none
!
   type :: bfgs_geoopt_hf
!
      character(len=400) :: name_
      character(len=400) :: tag
      character(len=400) :: note
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
      real(dp), dimension(:), allocatable :: energies, gradient_maxs 
!
      type(timings), allocatable :: timer 
!
   contains
!
      procedure :: run                 => run_bfgs_geoopt_hf
      procedure :: cleanup             => cleanup_bfgs_geoopt_hf
!
      procedure :: read_settings       => read_settings_bfgs_geoopt_hf
      procedure :: print_banner        => print_banner_bfgs_geoopt_hf
      procedure :: print_summary       => print_summary_bfgs_geoopt_hf
      procedure :: print_settings      => print_settings_bfgs_geoopt_hf
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
      solver%timer = timings('HF geometry optimization BFGS solver', pl='minimal')
      call solver%timer%turn_on()
!
!     Print solver banner
!
      solver%name_       = 'HF geometry optimization BFGS solver'
      solver%tag         = 'geometry optimization solver'
!
      solver%description = 'Constructs an approximate Hessian using the BFGS algorithm &
                           &using the previous geometries and gradients. From the BFGS Hessian, &
                           &a level shift given by the rational function (RF) augmented Hessian &
                           &is applied. See J. Comput. Chem. 18: 1473-1483, 1997. The HF gradient &
                           &is evaluated as described in J. Chem. Phys. 115.22 (2001): 10344-10352.'
!
      solver%note = 'the geometry optimization solver will &
                     &make successive calls to the SCF DIIS solver.'
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
!     Read & print settings (thresholds, etc.)
!
      call solver%read_settings()
      call solver%print_settings()
!
      allocate(solver%energies(solver%max_iterations))
      allocate(solver%gradient_maxs(solver%max_iterations))
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
                                        'solver scf geoopt', solver%gradient_threshold)
!
      call input%get_keyword_in_section('energy threshold', &
                                        'solver scf geoopt', solver%energy_threshold)
!
      call input%get_keyword_in_section('max iterations', &
                                        'solver scf geoopt', solver%max_iterations)
!
      call input%get_keyword_in_section('max step', &
                                        'solver scf geoopt', solver%max_step)
!
   end subroutine read_settings_bfgs_geoopt_hf
!
!
   subroutine print_settings_bfgs_geoopt_hf(solver)
!!
!!    Print settings
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(bfgs_geoopt_hf) :: solver
!
      call output%printf('m', '- BFGS geometry optimization settings:', fs='(/t3,a)')
!
      call output%printf('m', 'Gradient threshold: (e11.4)', &
                         reals=[solver%gradient_threshold], fs='(/t6,a)')
      call output%printf('m', 'Energy threshold:   (e11.4)', &
                         reals=[solver%energy_threshold], fs='(t6,a)')
!
      call output%printf('m', 'Max iterations:     (i11)', &
                         ints=[solver%max_iterations], fs='(/t6,a)')
      call output%printf('m', 'Max step size:      (e11.4)', &
                         reals=[solver%max_step], fs='(t6,a)')
!
   end subroutine print_settings_bfgs_geoopt_hf
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
      type(scf_diis_hf), allocatable :: hf_gs_solver
!
      logical :: restart
!
      if (solver%iteration > 1) then
!
!        Update geometry to the one requested
!
         call wf%system%set_geometry(geometry)
!
!        Re-decompose the AO overlap 
!
         call wf%destruct_pivot_matrix_ao_overlap()
         call wf%destruct_cholesky_ao_overlap()
!
         call wf%set_n_mo() 
!
!        Re-compute the AO h matrix 
!
         call wf%get_ao_h_wx(wf%ao_h)
!
!        Update Cauchy Schwarz screening list 
!
         call wf%construct_sp_eri_schwarz()
!
      endif
!
!     Attempt to converge HF orbitals/density
!
      restart = .false.
!
      if (solver%restart .or. solver%iteration > 1) restart = .true.
!
      hf_gs_solver = scf_diis_hf(wf, restart)
      call hf_gs_solver%run(wf)
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
!
      use array_utilities, only: get_abs_max
!
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
      real(dp), dimension(3,wf%system%n_atoms) :: step
      real(dp), dimension(3,wf%system%n_atoms) :: geometry
!
      type(bfgs_tool) :: bfgs
!
      type(timings), allocatable :: iteration_timer  
!
      bfgs = bfgs_tool(3*wf%system%n_atoms, solver%max_step)
      call bfgs%initialize_arrays()
!
      iteration_timer = timings('BFGS geoopt iteration time', pl='normal')
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
         call iteration_timer%turn_on()
!
         gradient = solver%determine_gradient(wf, geometry)
!
         energy = wf%energy 
         max_gradient = get_abs_max(gradient, 3*wf%system%n_atoms)
!
         call output%printf('n', 'Geometry optimization iteration (i0)', &
                            ints=[solver%iteration], fs='(/t3,a)')
         call output%printf('n', 'Absolute maximum of molecular gradient: (f17.12)', &
                            reals=[max_gradient], fs='(t3,a)')
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
            call output%printf('m', 'Geometry converged in (i0) iterations!', &
                               ints=[solver%iteration], fs='(/t3,a)')
!
         else
!
            call output%printf('m', 'Geometry not yet converged. Finding next geometry.', &
                               fs='(/t3,a)')
!
            call bfgs%update_hessian(geometry, gradient)
!
            call bfgs%get_step(gradient, step)
            geometry = geometry + step 
!
         endif
!
         call wf%system%print_geometry('angstrom')
         call wf%system%print_geometry('bohr')
!
         prev_energy = energy 
         call iteration_timer%turn_off()
         call iteration_timer%reset()
!
      enddo
!
      if (.not. converged) then 
!
         call output%error_msg('Was not able to converge the equations     &
                               &in the given number of maximum iterations.')
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
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none 
!
      class(bfgs_geoopt_hf) :: solver 
!
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('m', solver%description, ffs='(/t3,a)')
!
      call output%printf('m', 'Note: ' // solver%note, ffs='(/t3,a)')
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
      call output%printf('n', '- Summary of geometry optimization iterations: ', fs='(/t3,a)')
!
      call output%printf('n', 'Iteration       Energy                Gradient norm', &
                         fs='(/t6,a)')
      call output%print_separator('n', 58,'-', fs='(t6,a)')
!
      do iteration = 1, solver%iteration 
!
         call output%printf('n', '(i4)         (f19.12)     (e11.4)', &
                            ints=[iteration], &
                            reals=[solver%energies(iteration), &
                            solver%gradient_maxs(iteration)], fs='(t6,a)')
!
      enddo 
!
      call output%print_separator('n', 58, '-', fs='(t6,a)')
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
      call output%printf('m', '- Finished optimizing the HF geometry', fs='(/t3, a)')
!
      call solver%timer%turn_off()
!
   end subroutine cleanup_bfgs_geoopt_hf
!
!
end module bfgs_geoopt_hf_class
