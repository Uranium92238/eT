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
module newton_raphson_cc_gs_class
!
!!
!!		Newton-Raphson coupled cluster ground state solver class module
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
   type :: newton_raphson_cc_gs
!
      character(len=100) :: tag = 'Newton-Raphson coupled cluster ground state solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2019'
!
      character(len=500) :: description1 = 'A Newton-Raphson CC ground state equations solver. Solves the &
                                             ground state equation using updates Δt based on the Newton equation, &
                                             A Δt = -Ω, where the A is the Jacobian and Ω the omega vector.'
!
      integer :: max_iterations
!
      real(dp) :: omega_threshold
      real(dp) :: energy_threshold
!
      logical :: restart
!
   contains
!     
      procedure, nopass :: cleanup          => cleanup_newton_raphson_cc_gs
      procedure :: prepare                  => prepare_newton_raphson_cc_gs
      procedure :: run                      => run_newton_raphson_cc_gs
!
      procedure :: read_settings            => read_settings_newton_raphson_cc_gs
      procedure :: print_banner             => print_banner_newton_raphson_cc_gs
      procedure :: print_settings           => print_settings_newton_raphson_cc_gs
!
   end type newton_raphson_cc_gs
!
!
contains
!
!
   subroutine prepare_newton_raphson_cc_gs(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(newton_raphson_cc_gs) :: solver
!
      class(ccs) :: wf
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%max_iterations   = 100
      solver%energy_threshold = 1.0d-6
      solver%omega_threshold  = 1.0d-6
      solver%restart       = .false.
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
   end subroutine prepare_newton_raphson_cc_gs
!
!
   subroutine print_settings_newton_raphson_cc_gs(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(newton_raphson_cc_gs) :: solver 
!
!
   end subroutine print_settings_newton_raphson_cc_gs
!
!
   subroutine run_newton_raphson_cc_gs(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(newton_raphson_cc_gs) :: solver
!
      class(ccs) :: wf
!
      logical :: converged_residual 
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:), allocatable :: epsilon  
      real(dp), dimension(:), allocatable :: omega   
!
      integer :: iteration
!
      integer :: micro_iteration
!
!     Get preconditioner for micro iterations 
!
      call mem%alloc(epsilon, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(epsilon, wf%n_gs_amplitudes)
!
      converged = .false.
!
      iteration = 0
      prev_energy = zero
!
      call mem%alloc(omega, wf%n_gs_amplitudes)
!
      do while (.not. converged .and. iteration .le. solver%max_iterations) 
!
         iteration = iteration + 1
!
!        Construct Fock, calculate energy, and construct omega 
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
!        Print energy, energy difference and residual, then test convergence 
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e11.4,4x,e11.4)') iteration, wf%energy, &
                                          omega_norm, abs(wf%energy-prev_energy)
         flush(output%unit)
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_omega    = omega_norm              .lt. solver%omega_threshold
!
         converged = converged_omega .and. converged_energy
!
         if (iteration .eq. 1 .and. converged_omega) converged = .true. ! Exception to the rule
!
!        If not converged, perform micro-iterations to get an estimate for the next amplitudes 
!
         call solver%do_micro_iterations(omega, )
!
      enddo
!
   end subroutine run_newton_raphson_cc_gs
!
!
   subroutine cleanup_newton_raphson_cc_gs(wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%save_amplitudes()  
!
   end subroutine cleanup_newton_raphson_cc_gs
!
!
   subroutine print_banner_newton_raphson_cc_gs(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(newton_raphson_cc_gs) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.)
!
   end subroutine print_banner_newton_raphson_cc_gs
!
!
   subroutine read_settings_newton_raphson_cc_gs(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(newton_raphson_cc_gs) :: solver    
!
   end subroutine read_settings_newton_raphson_cc_gs
!
!
end module newton_raphson_cc_gs_class
