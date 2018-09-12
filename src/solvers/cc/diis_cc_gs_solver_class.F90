module diis_cc_gs_solver_class
!
!!
!!		DIIS coupled cluster ground state solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    
!
   use kinds
   use file_class
   use ccs_class
   use diis_tool_class
!
   implicit none
!
   type :: diis_cc_gs_solver
!
      integer(i15) :: diis_dimension = 8
!
      integer(i15) :: max_iterations = 100
!
      real(dp) :: energy_threshold = 1.0d-6
      real(dp) :: omega_threshold  = 1.0d-6
!
      logical :: restart = .false.
!
   contains
!     
      procedure :: prepare        => prepare_diis_cc_gs_solver
      procedure :: run            => run_diis_cc_gs_solver
      procedure :: cleanup        => cleanup_diis_cc_gs_solver
!
      procedure :: print_banner   => print_banner_diis_cc_gs_solver
      procedure :: print_summary  => print_summary_diis_cc_gs_solver
!
      procedure :: read_settings  => read_settings_diis_cc_gs_solver
!
      procedure :: print_settings => print_settings_diis_cc_gs_solver
!
      procedure :: do_diagonal_precondition => do_diagonal_precondition_diis_cc_gs_solver
!
   end type diis_cc_gs_solver
!
!
contains
!
!
   subroutine prepare_diis_cc_gs_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_gs_solver) :: solver
!
      class(ccs) :: wf
!
!     Print solver banner
!
      call solver%print_banner()
!
      write(output%unit, '(/t3,a)') ':: Preparing DIIS ground state CC object.'
!
!     Read settings (thresholds, etc.)
!
      write(output%unit, '(/t3,a/)') 'Reading solver settings.'
!
      call solver%read_settings()
      call solver%print_settings()
!
!     Set the amplitudes to the initial guess or read if restart
!
      call wf%initialize_amplitudes()
!
      if (solver%restart) then
!
       !  call wf%read_amplitudes()
! 
      else
!
         call wf%set_initial_amplitudes_guess()
!
      endif
!
   end subroutine prepare_diis_cc_gs_solver
!
!
   subroutine print_settings_diis_cc_gs_solver(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_gs_solver) :: solver 
!
      write(output%unit, '(t6,a16,i2)') 'DIIS dimension: ', solver%diis_dimension
!
   end subroutine print_settings_diis_cc_gs_solver
!
!
   subroutine do_diagonal_precondition_diis_cc_gs_solver(solver, alpha, preconditioner, vector, n)
!!
!!    Do diagonal precondition 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Performs the following operation:
!!
!!       v(n) = alpha*(v(n)/preconditioner(n)).
!!
      implicit none 
!     
      class(diis_cc_gs_solver), intent(in) :: solver 
!
      integer(i15), intent(in) :: n
!
      real(dp), intent(in) :: alpha 
!
      real(dp), dimension(n, 1), intent(in)    :: preconditioner
      real(dp), dimension(n, 1), intent(inout) :: vector  
!
      integer(i15) :: I 
!
!$omp parallel do private(I)
      do I = 1, n 
!
         vector(I, 1) = alpha*vector(I, 1)/preconditioner(I, 1)
!
      enddo 
!$omp end parallel do
!
   end subroutine do_diagonal_precondition_diis_cc_gs_solver
!
!
   subroutine run_diis_cc_gs_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_gs_solver) :: solver
!
      class(ccs) :: wf
!
      type(diis_tool) :: diis_manager
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_omega
!
      real(dp) :: energy, prev_energy
      real(dp) :: omega_norm
!
      real(dp), dimension(:,:), allocatable :: omega 
      real(dp), dimension(:,:), allocatable :: amplitudes  
      real(dp), dimension(:,:), allocatable :: epsilon  
!
      integer(i15) :: iteration
!
      write(output%unit, '(/t3,a)') ':: Running DIIS ground state CC object'
!
      call diis_manager%init('cc_gs_diis', wf%n_amplitudes, wf%n_amplitudes, solver%diis_dimension)
!
      call mem%alloc(omega, wf%n_amplitudes, 1)
      call mem%alloc(amplitudes, wf%n_amplitudes, 1)
      call mem%alloc(epsilon, wf%n_amplitudes, 1)
!
      converged          = .false.
      converged_energy   = .false.
      converged_omega    = .false.
!
      write(output%unit, '(/t3,a)') 'Iteration    Energy (a.u.)        || Omega ||    Delta E (a.u.)'
      write(output%unit, '(t3,a)')  '---------------------------------------------------------------'
!
      prev_energy = zero
      iteration   = 1
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)         
!
!        Calculate the energy and error vector omega 
!
         call wf%calculate_energy()
         energy = wf%energy
!
         call wf%construct_fock()
         call wf%construct_omega(omega)
         omega_norm = get_l2_norm(omega, wf%n_amplitudes)
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e10.4,4x,e10.4)') iteration, wf%energy, &
                                          omega_norm, abs(wf%energy-prev_energy)
         flush(output%unit)
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_omega    = omega_norm              .lt. solver%omega_threshold
!
         converged = converged_omega .and. converged_energy
!
         if (converged) then
!
            write(output%unit, '(t3,a)')           '--------------------------------------------------------------'
            write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
!
            call solver%print_summary(wf)
!
         else
!
!           Precondition omega, shift amplitudes by preconditioned omega, 
!           then ask for the DIIS update of the amplitudes 
!
            call wf%get_orbital_differences(epsilon)
            call solver%do_diagonal_precondition(-one, epsilon, omega, wf%n_amplitudes)
!
            call wf%get_amplitudes(amplitudes)
            amplitudes = amplitudes + omega 
!
            call diis_manager%update(omega, amplitudes)
            call wf%set_amplitudes(amplitudes)
!
            prev_energy = energy 
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(omega, wf%n_amplitudes, 1)
      call mem%dealloc(amplitudes, wf%n_amplitudes, 1)
      call mem%dealloc(epsilon, wf%n_amplitudes, 1)
!
      call diis_manager%finalize()
!
      if (.not. converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         stop
!
      endif 
!
   end subroutine run_diis_cc_gs_solver
!
!
   subroutine cleanup_diis_cc_gs_solver(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_gs_solver) :: solver
!
      class(ccs) :: wf
!
!     Nothing here yet 
!
   end subroutine cleanup_diis_cc_gs_solver
!
!
   subroutine print_banner_diis_cc_gs_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(diis_cc_gs_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: DIIS coupled cluster ground state solver'
      write(output%unit, '(t3,a)')  ':: E. F. Kjønstad, S. D. Folkestad, 2018'

      flush(output%unit)
!
   end subroutine print_banner_diis_cc_gs_solver
!
!
   subroutine print_summary_diis_cc_gs_solver(solver, wf)
!!
!!    Print summary 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(diis_cc_gs_solver) :: solver 
!
      class(ccs) :: wf 
!
    !  call wf%print_wavefunction_summary()
!
   end subroutine print_summary_diis_cc_gs_solver
!
!
   subroutine read_settings_diis_cc_gs_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_gs_solver) :: solver 
!
! read stuff 
!
   end subroutine read_settings_diis_cc_gs_solver
!
!
end module diis_cc_gs_solver_class
