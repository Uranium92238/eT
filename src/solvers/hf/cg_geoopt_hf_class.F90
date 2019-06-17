module cg_geoopt_hf_class
!
!!
!!    Geometry optimization solver using the conjugate gradient method for HF 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
!!    Supported wavefunctions: HF    
!!
!
   use hf_class 
   use hf_engine_class
   use conjugate_gradient_tool_class
!
   implicit none
!
   type :: cg_geoopt_hf
!
      character(len=400) :: tag
      character(len=400) :: author
      character(len=400) :: description
!
      integer :: max_iterations
!
      real(dp) :: energy_threshold
      real(dp) :: gradient_threshold
!
      real(dp), dimension(:), allocatable :: energies, gradient_norms 
!
   contains
!
      procedure :: run           => run_cg_geoopt_hf
      procedure :: cleanup       => cleanup_cg_geoopt_hf
!
      procedure :: read_settings => read_settings_cg_geoopt_hf
      procedure :: print_banner  => print_banner_cg_geoopt_hf
      procedure :: print_summary => print_summary_cg_geoopt_hf
!
   end type cg_geoopt_hf
!
!
   interface cg_geoopt_hf 
!
      procedure :: new_cg_geoopt_hf
!
   end interface cg_geoopt_hf
!
!
contains
!
!
   function new_cg_geoopt_hf() result(solver)
!!
!!    New CG geoopt HF 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(cg_geoopt_hf) :: solver
!
!     Print solver banner
!
      solver%tag     = 'Conjugate gradient geometry optimization solver'
      solver%author  = 'Åsmund H. Tveten, Eirik F. Kjønstad, 2019'
!
      solver%description = 'A solver that uses the Polak–Ribière conjugate gradient (CG) constant to find &
                           & a stationary geometry. Currently, the CG algorithm does not support a line &
                           & search.'
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%max_iterations      = 250
      solver%energy_threshold    = 1.0d-4
      solver%gradient_threshold  = 1.0d-4
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call mem%alloc(solver%energies, solver%max_iterations)
      call mem%alloc(solver%gradient_norms, solver%max_iterations)
!
   end function new_cg_geoopt_hf
!
!
   subroutine read_settings_cg_geoopt_hf(solver)
!!
!!    Read settings 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
      implicit none 
!
      class(cg_geoopt_hf), intent(inout) :: solver 
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
   end subroutine read_settings_cg_geoopt_hf
!
!
   subroutine run_cg_geoopt_hf(solver, wf)
!!
!!    Run 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(cg_geoopt_hf) :: solver
!
      class(hf) :: wf
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_gradient
!
      real(dp) :: energy, prev_energy, norm_gradient
!
      integer :: iteration
!
      real(dp), dimension(:,:), allocatable :: molecular_gradient
      real(dp), dimension(:,:), allocatable :: descent_direction
      real(dp), dimension(:,:), allocatable :: geometry
!
      type(hf_engine) :: ground_state_engine
!
      type(conjugate_gradient_tool) :: conjugate_gradient 
!
      conjugate_gradient   = conjugate_gradient_tool(3*wf%system%n_atoms)
      ground_state_engine  = hf_engine()
!
      iteration = 0
!
      converged          = .false.
      converged_energy   = .false.
      converged_gradient = .false.
!
      energy = zero
      prev_energy = zero
!
      call mem%alloc(molecular_gradient, 3, wf%system%n_atoms)
      call mem%alloc(descent_direction, 3, wf%system%n_atoms)
      call mem%alloc(geometry, 3, wf%system%n_atoms)
!
      molecular_gradient = zero 
      descent_direction = zero
!
      do while (.not. converged .and. iteration <= solver%max_iterations)        
!
         iteration = iteration + 1
!
         geometry = wf%system%get_geometry()
!
         call wf%set_n_mo()
         call ground_state_engine%ignite(wf)
!
         energy = wf%energy 
         call wf%construct_molecular_gradient(molecular_gradient)
!
         norm_gradient = get_l2_norm(molecular_gradient, 3*wf%system%n_atoms)
!
         solver%energies(iteration) = energy 
         solver%gradient_norms(iteration) = norm_gradient
!
         converged_gradient = norm_gradient              <= solver%gradient_threshold
         converged_energy   = abs(energy - prev_energy)  <= solver%energy_threshold
!
         converged = converged_gradient .and. converged_energy
!
         if (converged) then
!
            write(output%unit, '(/t3,a,i0,a)') 'Converged geometry in ', iteration, ' iterations!'
!
         else
!
            call conjugate_gradient%get_next_direction(molecular_gradient, descent_direction)
!
            geometry = geometry + descent_direction
            call wf%system%set_geometry(geometry)
!
         endif
!
         call wf%system%print_geometry()
         prev_energy = energy 
!
      enddo
!
      call mem%dealloc(molecular_gradient, 3, wf%system%n_atoms)
      call mem%dealloc(descent_direction, 3, wf%system%n_atoms)
      call mem%dealloc(geometry, 3, wf%system%n_atoms)
!
      if (.not. converged) then 
!
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         stop
!
      else
!
         call solver%print_summary(iteration)
!
      endif  
!
   end subroutine run_cg_geoopt_hf
!
!
   subroutine print_banner_cg_geoopt_hf(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(cg_geoopt_hf) :: solver 
!
      call output%long_string_print(solver%tag,'(//t3,a)',.true.)
      call output%long_string_print(solver%author,'(t3,a/)',.true.)
      call output%long_string_print(solver%description)
!
   end subroutine print_banner_cg_geoopt_hf
!
!
   subroutine print_summary_cg_geoopt_hf(solver, n_iterations)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(cg_geoopt_hf), intent(in) :: solver 
!
      integer, intent(in) :: n_iterations
!
      integer :: iteration 
!
      write(output%unit, '(/t3,a)') '- Summary of geometry optimization iterations: '
!
      write(output%unit, '(/t6,a)') 'Iteration       Energy                Gradient norm       '
      write(output%unit, '(t6,a)')  '----------------------------------------------------------'
!
      do iteration = 1, n_iterations
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
   end subroutine print_summary_cg_geoopt_hf
!
!
   subroutine cleanup_cg_geoopt_hf(solver)
!!
!!    Cleanup 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(cg_geoopt_hf) :: solver 
!
      call mem%dealloc(solver%energies, solver%max_iterations)
      call mem%dealloc(solver%gradient_norms, solver%max_iterations)
!
   end subroutine cleanup_cg_geoopt_hf
!
!
end module cg_geoopt_hf_class