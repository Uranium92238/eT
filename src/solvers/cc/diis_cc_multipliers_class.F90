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
module diis_cc_multipliers_class
!
!!
!!		DIIS coupled cluster multipliers solver class module
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
   type :: diis_cc_multipliers
!
      character(len=100) :: tag = 'Diis multipliers solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      character(len=500) :: description1 = 'A DIIS CC multiplier equations solver. It combines a quasi-Newton &
                                           &perturbation theory estimate of the next multipliers, using &
                                           &least square fitting to find an an optimal combination of &  
                                           &previous estimates such that the update is minimized.'
!
      character(len=500) :: description2 = 'See Helgaker et al., Molecular Electronic Structure Theory, &
                                           &Chapter 13, for the more details on this algorithm.'
!
      integer :: diis_dimension
!
      integer :: max_iterations
!
      real(dp) :: residual_threshold
!
      logical :: restart
!
   contains
!     
      procedure, nopass :: cleanup          => cleanup_diis_cc_multipliers
      procedure, nopass :: do_diagonal_precondition => do_diagonal_precondition_diis_cc_multipliers
!
      procedure :: prepare                  => prepare_diis_cc_multipliers
      procedure :: run                      => run_diis_cc_multipliers
!
      procedure :: read_settings            => read_settings_diis_cc_multipliers
!
      procedure :: print_banner             => print_banner_diis_cc_multipliers
      procedure, nopass :: print_summary    => print_summary_diis_cc_multipliers
!
      procedure :: print_settings           => print_settings_diis_cc_multipliers
!
   end type diis_cc_multipliers
!
!
contains
!
!
   subroutine prepare_diis_cc_multipliers(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set default settings
!
      solver%diis_dimension = 8
      solver%max_iterations = 100
!
      solver%residual_threshold  = 1.0d-6
!
      solver%restart = .false.
!
      call solver%read_settings()
!
      call solver%print_settings()
!
      call wf%construct_fock()
!
      call wf%initialize_multipliers()
!
   end subroutine prepare_diis_cc_multipliers
!
!
   subroutine print_settings_diis_cc_multipliers(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_multipliers) :: solver 
!
      write(output%unit, '(/t3,a)')      '- DIIS CC multipliers solver settings:'
!
      write(output%unit, '(/t6,a26,e9.2)') 'Residual threshold:       ', solver%residual_threshold

      write(output%unit, '(/t6,a26,i9)')   'DIIS dimension:           ', solver%diis_dimension
      write(output%unit, '(t6,a26,i9)')    'Max number of iterations: ', solver%max_iterations
!
   end subroutine print_settings_diis_cc_multipliers
!
!
   subroutine do_diagonal_precondition_diis_cc_multipliers(alpha, preconditioner, vector, n)
!!
!!    Do diagonal precondition 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Performs the following operation:
!!
!!       v(n) = alpha*(v(n)/preconditioner(n)).
!!
!!    This is now used in three solvers. I think the correct way to handle this 
!!    is not by inheritance but by generalizing this as a utility which can be 
!!    used by each solver irrespective of ancestry. - Eirik, Sep 2018
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
   end subroutine do_diagonal_precondition_diis_cc_multipliers
!
!
   subroutine run_diis_cc_multipliers(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
      type(diis_tool) :: diis_manager
!
      logical :: converged_residual 
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:), allocatable :: residual  
      real(dp), dimension(:), allocatable :: multipliers
      real(dp), dimension(:), allocatable :: epsilon  
!
      integer :: iteration
!
      call diis_manager%init('cc_multipliers_diis', wf%n_gs_amplitudes, wf%n_gs_amplitudes, solver%diis_dimension)
!
      call mem%alloc(residual, wf%n_gs_amplitudes)
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
      call mem%alloc(epsilon, wf%n_gs_amplitudes)
!
      call wf%get_gs_orbital_differences(epsilon, wf%n_gs_amplitudes)
!
      if (solver%restart) then 
!
         write(output%unit, '(/t3,a)') 'Requested restart. Reading multipliers from file.'
!
         call wf%read_multipliers()
         call wf%get_multipliers(multipliers) 
!
      else
!
         multipliers = zero
         call wf%set_multipliers(multipliers)
!
      endif
!
      converged_residual = .false.
!
      write(output%unit, '(/t3,a)') 'Iteration    Norm residual  '
      write(output%unit, '(t3,a)')  '----------------------------'
      flush(output%unit)
!
      iteration   = 1
!
      do while (.not. converged_residual .and. iteration .le. solver%max_iterations)         
!
!        Construct the multiplier equations 
!
         call wf%construct_multiplier_equation(residual)
         residual_norm = get_l2_norm(residual, wf%n_gs_amplitudes)
!
         write(output%unit, '(t3,i3,10x,e11.4)') iteration, residual_norm
         flush(output%unit)
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         converged_residual = residual_norm .lt. solver%residual_threshold
!
         if (converged_residual) then
!
            write(output%unit, '(t3,a)')           '----------------------------'
            write(output%unit, '(/t3,a29,i3,a12)') 'Convergence criterion met in ', iteration, ' iterations!'
!
         else
!
!           Precondition residual, shift multipliers by preconditioned residual, 
!           then ask for the DIIS update of the multipliers 
!
            call solver%do_diagonal_precondition(-one, epsilon, residual, wf%n_gs_amplitudes)
!
            call wf%get_multipliers(multipliers)
            multipliers = multipliers + residual 
!
            call diis_manager%update(residual, multipliers)
            call wf%set_multipliers(multipliers)
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call diis_manager%finalize()
!
      if (.not. converged_residual) then 
!   
         write(output%unit, '(t3,a)')   '---------------------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Warning: was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
!
      else
!
         call solver%print_summary(wf, multipliers)
!
      endif 
!
      flush(output%unit)
!
      call mem%dealloc(residual, wf%n_gs_amplitudes)
      call mem%dealloc(multipliers, wf%n_gs_amplitudes)
      call mem%dealloc(epsilon, wf%n_gs_amplitudes)
!
   end subroutine run_diis_cc_multipliers
!
!
   subroutine cleanup_diis_cc_multipliers(wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%save_multipliers()
!
   end subroutine cleanup_diis_cc_multipliers
!
!
   subroutine print_banner_diis_cc_multipliers(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(diis_cc_multipliers) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
      call long_string_print(solver%description2)
!
   end subroutine print_banner_diis_cc_multipliers
!
!
   subroutine print_summary_diis_cc_multipliers(wf, X)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: X
!
      write(output%unit, '(/t3,a)') '- Multipliers vector amplitudes:'
      flush(output%unit)      
!
      call wf%print_dominant_x_amplitudes(X, 'r')
!
   end subroutine print_summary_diis_cc_multipliers
!
!
   subroutine read_settings_diis_cc_multipliers(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_multipliers) :: solver 
!
      call input%get_keyword_in_section('threshold', 'solver cc multipliers', solver%residual_threshold)
      call input%get_keyword_in_section('max iterations', 'solver cc multipliers', solver%max_iterations)
!
      if (input%requested_keyword_in_section('restart', 'solver cc multipliers')) solver%restart = .true.    
!
   end subroutine read_settings_diis_cc_multipliers
!
!
end module diis_cc_multipliers_class
!
!  Testing multipliers against dalton:
!
!  Run Dalton, requesting transition moments (*CCOPA)
!
!  Copy the file CCL0_1__1 to the eT scratch folder
!
!  This file is sequential and unformatted and contains the following records:
!
!     1. Some integer and a string containing the CC-method used 
!
!     2. Multipliers singles-vector
!
!     3. Multipliers doubles-vector (for CCSD/CC2)
!  
