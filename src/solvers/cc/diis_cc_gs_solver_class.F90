!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module diis_cc_gs_solver_class
!
!!
!!		DIIS coupled cluster ground state solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
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
      character(len=100) :: tag = 'DIIS CC ground solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      character(len=500) :: description1 = 'A DIIS CC ground state amplitude equations solver. It uses &
                                           &an extrapolation of previous quasi-Newton perturbation theory &
                                           &estimates of the next amplitudes. See Helgaker et al., Molecular & 
                                           &Electronic Structure Theory, Chapter 13.'
!
      integer :: diis_dimension
!
      integer :: max_iterations 
!
      real(dp) :: energy_threshold
      real(dp) :: omega_threshold 
!
      type(file) :: restart_file
      logical    :: do_restart
!
   contains
!     
      procedure, nopass :: do_diagonal_precondition => do_diagonal_precondition_diis_cc_gs_solver
!
      procedure :: prepare                  => prepare_diis_cc_gs_solver
      procedure :: run                      => run_diis_cc_gs_solver
      procedure :: cleanup                  => cleanup_diis_cc_gs_solver
!
      procedure :: print_banner             => print_banner_diis_cc_gs_solver
      procedure :: read_settings            => read_settings_diis_cc_gs_solver
!
      procedure :: print_settings           => print_settings_diis_cc_gs_solver
      procedure, nopass :: print_summary    => print_summary_diis_cc_gs_solver
!
      procedure :: restart                  => restart_diis_cc_gs_solver
      procedure :: write_restart_file       => write_restart_file_diis_cc_gs_solver
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
!     Set standard settings 
!
      solver%diis_dimension   = 8 
      solver%max_iterations   = 100
      solver%energy_threshold = 1.0d-6
      solver%omega_threshold  = 1.0d-6
      solver%do_restart       = .false.
!
!     Read & print settings (thresholds, etc.)
!
      if (requested_section('cc ground state')) call solver%read_settings()
!
      call solver%print_settings()
!
!     Set the amplitudes to the initial guess or read if restart
!
      call wf%initialize_amplitudes()
!
!     Prepare restart information file 
!
      call solver%restart_file%init('diis_cc_gs_restart_info', 'sequential', 'formatted')
!
      if (solver%do_restart) then
!
         call solver%restart(wf)
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
   subroutine write_restart_file_diis_cc_gs_solver(solver, wf)
!!
!!    Write restart file 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2018
!!
      implicit none 
!
      class(diis_cc_gs_solver), intent(inout) :: solver
!
      class(ccs), intent(in) :: wf 
!
      call disk%open_file(solver%restart_file, 'write', 'rewind')
!
      write(solver%restart_file%unit, *) 'n_amplitudes'
      write(solver%restart_file%unit, *) wf%n_gs_amplitudes
!
      call disk%close_file(solver%restart_file)       
!
   end subroutine write_restart_file_diis_cc_gs_solver
!
!
   subroutine restart_diis_cc_gs_solver(solver, wf)
!!
!!    Restart 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(diis_cc_gs_solver), intent(inout) :: solver 
!
      class(ccs), intent(inout) :: wf 
!
      integer :: n_gs_amplitudes
!
      write(output%unit, '(/t6,a)') 'Requested restart. Reading amplitudes from file.'
!
!     Sanity checks 
!
      call disk%open_file(solver%restart_file, 'read', 'rewind')
!
      read(solver%restart_file%unit, *) ! Empty read to skip banner 
!
      read(solver%restart_file%unit, *) n_gs_amplitudes 
!
      call disk%close_file(solver%restart_file)
!
      if (n_gs_amplitudes .ne. wf%n_gs_amplitudes) call output%error_msg('Inconsistent dimensions on restart in DIIS-CC-GS.')
!
      call wf%read_amplitudes()
!
   end subroutine restart_diis_cc_gs_solver
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
      write(output%unit, '(/t3,a)')      '- DIIS CC ground state solver settings:'
!
      write(output%unit, '(/t6,a26,e9.2)') 'Omega threshold:          ', solver%omega_threshold
      write(output%unit, '(t6,a26,e9.2)')  'Energy threshold:         ', solver%energy_threshold

      write(output%unit, '(/t6,a26,i9)')   'DIIS dimension:           ', solver%diis_dimension
      write(output%unit, '(t6,a26,i9)')    'Max number of iterations: ', solver%max_iterations
!
      flush(output%unit)
!
   end subroutine print_settings_diis_cc_gs_solver
!
!
   subroutine do_diagonal_precondition_diis_cc_gs_solver(alpha, preconditioner, vector, n)
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
      integer, intent(in) :: n
!
      real(dp), intent(in) :: alpha 
!
      real(dp), dimension(n, 1), intent(in)    :: preconditioner
      real(dp), dimension(n, 1), intent(inout) :: vector  
!
      integer :: I 
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
      integer :: iteration
!
      call diis_manager%init('cc_gs_diis', wf%n_gs_amplitudes, wf%n_gs_amplitudes, solver%diis_dimension)
!
      call mem%alloc(omega, wf%n_gs_amplitudes, 1)
      call mem%alloc(amplitudes, wf%n_gs_amplitudes, 1)
      call mem%alloc(epsilon, wf%n_gs_amplitudes, 1)
!
      converged          = .false.
      converged_energy   = .false.
      converged_omega    = .false.
!
      write(output%unit, '(/t3,a)') 'Iteration    Energy (a.u.)        |omega|       Delta E (a.u.) '
      write(output%unit, '(t3,a)')  '---------------------------------------------------------------'
      flush(output%unit)
!
      prev_energy = zero
      iteration   = 1
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)         
!
!        Calculate the energy and error vector omega 
!
         call wf%integrals%write_t1_cholesky(wf%t1)
         call wf%integrals%can_we_keep_g_pqrs()
!
         call wf%construct_fock()
!
         call wf%calculate_energy()
         energy = wf%energy
!
         call wf%construct_omega(omega)
!
         omega_norm = get_l2_norm(omega, wf%n_gs_amplitudes)
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e11.4,4x,e11.4)') iteration, wf%energy, &
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
         if (iteration .eq. 1 .and. converged_omega) converged = .true. ! Exception to the rule
!
         if (converged) then
!
            write(output%unit, '(t3,a)')           '--------------------------------------------------------------'
            write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
!
            if (.not. converged_energy) then 
!
               write(output%unit, '(/t3,a,/t9,a)') 'Note: the omega vector converged in the first iteration,', &
                                                         'so the energy convergence has not been tested!'
!
            endif
!
         else
!
!           Precondition omega, shift amplitudes by preconditioned omega, 
!           then ask for the DIIS update of the amplitudes 
!
            call wf%get_gs_orbital_differences(epsilon, wf%n_gs_amplitudes)
            call solver%do_diagonal_precondition(-one, epsilon, omega, wf%n_gs_amplitudes)
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
      call mem%dealloc(omega, wf%n_gs_amplitudes, 1)
      call mem%dealloc(amplitudes, wf%n_gs_amplitudes, 1)
      call mem%dealloc(epsilon, wf%n_gs_amplitudes, 1)
!
      call diis_manager%finalize()
!
      if (.not. converged) then 
!   
         write(output%unit, '(t3,a)')   '---------------------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Warning: was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
!
      else
!
         call solver%print_summary(wf)
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
!     Write restart file & save amplitudes 
!
      call solver%write_restart_file(wf)
!
      call wf%save_amplitudes()
!
   end subroutine cleanup_diis_cc_gs_solver
!
!
   subroutine print_banner_diis_cc_gs_solver(solver)
!!
!!    Print banner
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none 
!
      class(diis_cc_gs_solver) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a)')
!
   end subroutine print_banner_diis_cc_gs_solver
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
      integer :: n_specs, i
!
      character(len=100) :: line
!
      call move_to_section('cc ground state', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:16) == 'omega threshold:' ) then
!
            read(line(17:100), *) solver%omega_threshold
!
         elseif (line(1:17) == 'energy threshold:' ) then
!
            read(line(18:100), *) solver%energy_threshold
!
         elseif (line(1:15) == 'diis dimension:' ) then
!
            read(line(16:100), *) solver%diis_dimension
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         elseif (trim(line) == 'restart') then
!
            solver%do_restart = .true.
!
         endif
!
      enddo
!
   end subroutine read_settings_diis_cc_gs_solver
!
!
   subroutine print_summary_diis_cc_gs_solver(wf)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp) :: t1_diagnostic 
!
      write(output%unit, '(/t3,a)') '- DIIS CC ground state solver summary:'
!
      write(output%unit, '(/t6,a33,f18.12)') 'Final ground state energy (a.u.):', wf%energy 
      call wf%print_dominant_amplitudes()
!
      t1_diagnostic = wf%get_t1_diagnostic() 
      write(output%unit, '(/t6,a32,f14.12)') 'T1 diagnostic (|T1|/sqrt(N_e)): ', t1_diagnostic
!
   end subroutine print_summary_diis_cc_gs_solver
!
!
end module diis_cc_gs_solver_class
