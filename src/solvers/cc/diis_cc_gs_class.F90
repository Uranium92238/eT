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
module diis_cc_gs_class
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
   type :: diis_cc_gs
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
      logical  :: restart
!
      type(timings) :: timer
!
   contains
!     
      procedure, nopass :: do_diagonal_precondition => do_diagonal_precondition_diis_cc_gs
!
      procedure :: run                      => run_diis_cc_gs
      procedure :: cleanup                  => cleanup_diis_cc_gs
!
      procedure :: print_banner             => print_banner_diis_cc_gs
      procedure :: read_settings            => read_settings_diis_cc_gs
!
      procedure :: print_settings           => print_settings_diis_cc_gs
      procedure, nopass :: print_summary    => print_summary_diis_cc_gs
!
   end type diis_cc_gs
!
!
   interface diis_cc_gs 
!
      procedure :: new_diis_cc_gs
!
   end interface diis_cc_gs 
!
!
contains
!
!
   function new_diis_cc_gs(wf) result(solver)
!!
!!    New DIIS CC GS 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(diis_cc_gs) :: solver
!
      class(ccs) :: wf
!
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' ground state')
      call solver%timer%turn_on()
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
      solver%restart          = .false.
!
!     Read & print settings (thresholds, etc.)
!
      call solver%read_settings()
      call solver%print_settings()
!
!     Set the amplitudes to the initial guess or read if restart
!
      call wf%initialize_amplitudes()
!
!     Prepare restart information file 
!
      if (solver%restart) then
!
         call wf%is_restart_safe('ground state')
!
         write(output%unit, '(/t3,a)') 'Requested restart. Reading in solution from file.'
!
         call wf%read_amplitudes()
         call wf%integrals%write_t1_cholesky(wf%t1) 
! 
      else
!
         call wf%integrals%write_t1_cholesky(wf%t1) 
         call wf%set_initial_amplitudes_guess()
!
      endif
!
   end function new_diis_cc_gs
!
!
   subroutine print_settings_diis_cc_gs(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_gs) :: solver 
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
   end subroutine print_settings_diis_cc_gs
!
!
   subroutine do_diagonal_precondition_diis_cc_gs(alpha, preconditioner, vector, n)
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
      real(dp), dimension(n), intent(in)    :: preconditioner
      real(dp), dimension(n), intent(inout) :: vector  
!
      integer :: I 
!
!$omp parallel do private(I)
      do I = 1, n 
!
         vector(I) = alpha*vector(I)/preconditioner(I)
!
      enddo 
!$omp end parallel do
!
   end subroutine do_diagonal_precondition_diis_cc_gs
!
!
   subroutine run_diis_cc_gs(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_gs) :: solver
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
      real(dp), dimension(:), allocatable :: omega 
      real(dp), dimension(:), allocatable :: amplitudes  
      real(dp), dimension(:), allocatable :: epsilon  
!
      integer :: iteration
!
      diis_manager = diis_tool('cc_gs_diis', wf%n_gs_amplitudes, wf%n_gs_amplitudes, solver%diis_dimension)
!
      call mem%alloc(omega, wf%n_gs_amplitudes)
      call mem%alloc(amplitudes, wf%n_gs_amplitudes)
      call mem%alloc(epsilon, wf%n_gs_amplitudes)
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
         flush(timing%unit)
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
!
            call wf%form_newton_raphson_t_estimate(amplitudes, omega)
!
            call diis_manager%update(omega, amplitudes)
            call wf%set_amplitudes(amplitudes)
!
            prev_energy = energy 
!
!           Compute the new T1 transformed Cholesky vectors,
!           and store in memory the entire ERI-T1 matrix if possible and necessary 
!
            call wf%integrals%write_t1_cholesky(wf%t1)
            if (wf%need_g_abcd()) call wf%integrals%can_we_keep_g_pqrs_t1()
!
         endif
!
         iteration = iteration + 1
!
!        Save amplitudes
!
         call wf%save_amplitudes()
!
      enddo
!
      call mem%dealloc(omega, wf%n_gs_amplitudes)
      call mem%dealloc(amplitudes, wf%n_gs_amplitudes)
      call mem%dealloc(epsilon, wf%n_gs_amplitudes)
!
      call diis_manager%cleanup()
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
   end subroutine run_diis_cc_gs
!
!
   subroutine cleanup_diis_cc_gs(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_gs) :: solver
      class(ccs) :: wf
!
!     Save amplitudes 
!
      call wf%save_amplitudes()
!
      call solver%timer%turn_off()
!
      write(output%unit, '(/t3, a)') '- Finished solving the ' // trim(convert_to_uppercase(wf%name_)) // &
                                       ' ground state equations'
!
      write(output%unit, '(/t6,a23,f20.5)')  'Total wall time (sec): ', solver%timer%get_elapsed_time('wall')
      write(output%unit, '(t6,a23,f20.5)')   'Total cpu time (sec):  ', solver%timer%get_elapsed_time('cpu')
!
   end subroutine cleanup_diis_cc_gs
!
!
   subroutine print_banner_diis_cc_gs(solver)
!!
!!    Print banner
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none 
!
      class(diis_cc_gs) :: solver 
!
      call output%long_string_print(solver%tag,'(//t3,a)',.true.)
      call output%long_string_print(solver%author,'(t3,a/)',.true.)
      call output%long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a)')
!
   end subroutine print_banner_diis_cc_gs
!
!
   subroutine read_settings_diis_cc_gs(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_gs) :: solver 
!
      call input%get_keyword_in_section('omega threshold', 'solver cc gs', solver%omega_threshold)
      call input%get_keyword_in_section('energy threshold', 'solver cc gs', solver%energy_threshold)
      call input%get_keyword_in_section('diis dimension', 'solver cc gs', solver%diis_dimension)
      call input%get_keyword_in_section('max iterations', 'solver cc gs', solver%max_iterations)
!
      if (input%requested_keyword_in_section('restart', 'solver cc gs')) solver%restart = .true.
!
   end subroutine read_settings_diis_cc_gs
!
!
   subroutine print_summary_diis_cc_gs(wf)
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
   end subroutine print_summary_diis_cc_gs
!
!
end module diis_cc_gs_class
