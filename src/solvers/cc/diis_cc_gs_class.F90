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
   use parameters
!
   use global_in, only : input
   use global_out, only : output
!
   use string_utilities, only : convert_to_uppercase
   use array_utilities, only : get_l2_norm
   use memory_manager_class, only : mem
   use ccs_class, only : ccs
   use diis_tool_class, only : diis_tool
   use timings_class, only : timings
   use precondition_tool_class, only : precondition_tool
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
      character(len=200) :: storage 
      logical :: restart, records_in_memory 
!
      type(timings) :: timer
!
      class(precondition_tool), allocatable :: preconditioner 
!
   contains
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
      real(dp), dimension(:), allocatable :: eps
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' ground state')
      call solver%timer%turn_on()
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%diis_dimension      = 8 
      solver%max_iterations      = 100
      solver%energy_threshold    = 1.0d-6
      solver%omega_threshold     = 1.0d-6
      solver%restart             = .false.
      solver%storage             = 'disk'
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
         call output%printf('Requested restart. Reading in solution from file.', fs='(/t3,a)', pl='minimal')
!
         call wf%read_amplitudes()
!
         call wf%integrals%write_t1_cholesky(wf%t1) 
         if(wf%need_g_abcd .and. wf%integrals%room_for_g_pqrs_t1()) &
            call wf%integrals%place_g_pqrs_t1_in_memory()
! 
      else
!
         call wf%integrals%write_t1_cholesky(wf%t1) 
         if(wf%need_g_abcd .and. wf%integrals%room_for_g_pqrs_t1()) &
            call wf%integrals%place_g_pqrs_t1_in_memory()
!
         call wf%set_initial_amplitudes_guess()
!
      endif
!
!     Determine whether to store records in memory or on file
!
      if (trim(solver%storage) == 'memory') then 
!
         solver%records_in_memory = .true.
!
      elseif (trim(solver%storage) == 'disk') then 
!
         solver%records_in_memory = .false.
!
      else 
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                 trim(solver%storage))
!
      endif 
!
!     Initialize preconditioner 
!
      call mem%alloc(eps, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(eps, wf%n_gs_amplitudes)
!
      solver%preconditioner = precondition_tool(eps, wf%n_gs_amplitudes)
!
      call mem%dealloc(eps, wf%n_gs_amplitudes)
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
      call output%printf('- DIIS CC ground state solver settings:', fs='(/t3,a)', pl='minimal')
!
      call output%printf('Omega threshold:          (e9.2)', reals=[solver%omega_threshold], fs='(/t6, a)', pl='minimal')
      call output%printf('Energy threshold:         (e9.2)', reals=[solver%energy_threshold], fs='(t6, a)', pl='minimal')
!
      call output%printf('DIIS dimension:           (i9)', ints=[solver%diis_dimension], fs='(/t6, a)', pl='minimal')
      call output%printf('Max number of iterations: (i9)', ints=[solver%max_iterations], fs='(t6, a)', pl='minimal')
!
   end subroutine print_settings_diis_cc_gs
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
      type(diis_tool) :: diis
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
!
      integer :: iteration
!
      diis = diis_tool('cc_gs_diis', wf%n_gs_amplitudes, wf%n_gs_amplitudes, &
               solver%records_in_memory, dimension_=solver%diis_dimension)
!
      call mem%alloc(omega, wf%n_gs_amplitudes)
      call mem%alloc(amplitudes, wf%n_gs_amplitudes)
!
      converged          = .false.
      converged_energy   = .false.
      converged_omega    = .false.
!
      call output%printf('Iteration    Energy (a.u.)        |omega|       Delta E (a.u.) ', fs='(/t3,a)', pl='normal')
      call output%print_separator('n', 63,'-')
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
         call output%printf('(i3)  (f25.12)    (e11.4)    (e11.4)', &
            ints=[iteration], &
            reals=[wf%energy, omega_norm, abs(wf%energy-prev_energy)], &
            fs='(t3, a)', pl='normal')
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
            call output%print_separator('n', 63,'-')
!
           call output%printf('Convergence criterion met in (i0) iterations!', ints=[iteration], fs='(/t3,a)', pl='normal') 
!
            if (.not. converged_energy) then 
!
!
               call output%printf('Note: the omega vector converged in the first iteration, &
                                 & so the energy convergence has not been tested!', &
                                 ffs='(/t3,a)', pl='normal')
!
            endif
!
         else
!
!           Precondition omega, shift amplitudes by preconditioned omega, 
!           then ask for the DIIS update of the amplitudes 
!
            call solver%preconditioner%do_(omega)
!
            call wf%get_amplitudes(amplitudes)
            call wf%form_newton_raphson_t_estimate(amplitudes, omega)
!
            call diis%update(omega, amplitudes)
            call wf%set_amplitudes(amplitudes)
!
            prev_energy = energy 
!
!           Compute the new T1 transformed Cholesky vectors,
!           and store in memory the entire ERI-T1 matrix if possible and necessary 
!
            call wf%integrals%write_t1_cholesky(wf%t1)
            if (wf%integrals%get_eri_t1_mem()) &
               call wf%integrals%update_g_pqrs_t1_in_memory()
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
!
      if (.not. converged) then 
!   
         call output%print_separator('n', 63,'-')
!
         call output%error_msg('Did not converge in the max number of iterations.')
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
      call output%printf('- Finished solving the (a0) ground state equations', & 
                        chars=[convert_to_uppercase(wf%name_)], &
                        fs='(/t3, a)', pl='minimal')
!
      call output%printf('Total wall time (sec): (f20.5)', &
            reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6, a)', pl='minimal')
      call output%printf('Total cpu time (sec):  (f20.5)', &
            reals=[solver%timer%get_elapsed_time('cpu')], fs='(t6, a)', pl='minimal')
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
      call input%get_keyword_in_section('storage', 'solver cc gs', solver%storage)
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
      call output%printf('- DIIS CC ground state solver summary:', fs='(/t3,a)', pl='minimal')
!
      call output%printf('Final ground state energy (a.u.): (f18.12)', reals=[wf%energy], fs='(/t6,a)', pl='minimal')
!
      call wf%print_dominant_amplitudes()
!
      t1_diagnostic = wf%get_t1_diagnostic() 
      call output%printf('T1 diagnostic (|T1|/sqrt(N_e)): (f14.12)', reals=[t1_diagnostic], fs='(/t6,a)', pl='minimal')
!
   end subroutine print_summary_diis_cc_gs
!
!
end module diis_cc_gs_class
