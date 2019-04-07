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
module cc_properties_class
!
!!
!!    Coupled cluster properties solver class module
!!    Written by Josefine H. Andersen, 2019
!!    
!
   use kinds
   use file_class
   use ccs_class
!
   implicit none
!
   type :: cc_properties
!
      character(len=100) :: tag = 'Coupled cluster properties solver'
      character(len=100) :: author = 'Josefine H. Andersen, 2019'
!
      character(len=500) :: description1 = 'A solver that calculates properties from coupled &
                                           &cluster excited state calculations'
!
      character(len=10)  :: es_type
      character(len=40)  :: operator_type
      character(len=40)  :: X
!
      logical    :: eom
      logical    :: linear_response
!
      character(len=1), dimension(3) :: component = ['X', 'Y', 'Z']
!
      real(dp) :: S
!
      real(dp), dimension(:), allocatable :: etaX
      real(dp), dimension(:), allocatable :: csiX
!
      real(dp) :: T_r, T_l
!
      integer :: n_singlet_states = 0
!
   contains
!
      procedure, non_overridable :: prepare          => prepare_cc_properties
      procedure, non_overridable :: run              => run_cc_properties
      procedure, non_overridable :: cleanup          => cleanup_cc_properties
!
      procedure :: intialize_variables               => intialize_variables_cc_properties
      procedure :: read_settings                     => read_settings_cc_properties
!
      procedure :: reset                             => reset_properties
!
      procedure :: print_banner                      => print_banner_cc_properties
      procedure :: print_settings                    => print_settings_cc_properties
      procedure :: print_summary                     => print_summary_cc_properties
!
   end type cc_properties
!
!
contains
!
!
   subroutine prepare_cc_properties(solver, wf)
!!
!!    Prepare 
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      class(ccs) :: wf
!
      call solver%print_banner()
!
!     Set default
!
      solver%eom             = .false.
      solver%linear_response = .false.
      solver%es_type         = 'valence'
      solver%operator_type   = 'dipole_length'
!
      call solver%read_settings()
      call solver%print_settings()
!
      call solver%intialize_variables(wf)
!
!     Read multipliers and prepare operator
!
      call wf%initialize_multipliers()
      call wf%read_multipliers()
!
      call wf%prepare_operator_pq(solver%operator_type)
!
   end subroutine prepare_cc_properties
!
!
   subroutine run_cc_properties(solver, wf) 
!!
!!    Run
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      class(ccs) :: wf
!
      integer :: i, n
!
      call solver%print_summary('header')
!      
!     Loop over components of X operator
!
      do i = 1, 3
!
         call solver%print_summary('top', i)
!         
         solver%X = solver%component(i) //'_'// trim(solver%operator_type) 
!
         call solver%reset()
!
         call wf%construct_etaX(solver%X, solver%etaX)      
!
         call wf%construct_csiX(solver%X, solver%csiX)      
!
!        Do EOM ?
!
         if (solver%eom) call wf%get_eom_contribution(solver%etaX, solver%csiX, solver%X)
!
!        Loop over excited states and calculate transition strength
!  
         do n = 1, solver%n_singlet_states
!
            call wf%calculate_transition_strength(solver%S, solver%etaX, &
                               solver%csiX, n, solver%T_l, solver%T_r)
!
            call solver%print_summary('results', n)
!            
         enddo
!         
         call solver%print_summary('bottom')
!
      enddo
!
   end subroutine run_cc_properties
!
!
   subroutine read_settings_cc_properties(solver)
!!
!!    Read cc excited state settings to find n_singlet_states
!!    and cc properties settings to determine operator
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      integer :: n_specs, n_ops, i
!
      character(len=100) :: line
!
      if (.not. requested_section('cc excited state')) then
!
         call output%error_msg('Number of excitations must be specified.')
!
      endif
!
      call move_to_section('cc excited state', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:15) == 'singlet states:' ) then
!
            read(line(16:100), *) solver%n_singlet_states
!
         elseif (line(1:15) == 'core excitation' ) then
!                 
            solver%es_type = 'core'
!            
         elseif (line(1:10) == 'valence excitation' ) then
!                 
            solver%es_type = 'valence'
!            
         endif
!
      enddo
!     
      if (.not. requested_section('cc properties')) then
!
         call output%error_msg('Operator type must be specified for properties')
!
      endif
!
      call move_to_section('cc properties', n_ops)
!
      do i = 1, n_ops
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:9) == 'operator:') then
!
            if (trim(line(11:100)) == 'dipole length') then
!
               solver%operator_type = 'dipole_length'
!
            else
!
               call output%error_msg('Operator must be specified.')
!
            endif
!
         elseif (line(1:3) == 'eom') then
!
            solver%eom = .true.
!
         elseif (line(1:15) == 'linear response') then
!
            solver%linear_response = .true.
!
            call output%error_msg('Linear response not implemented')
!
         endif
!
      enddo
!
   end subroutine read_settings_cc_properties
!
!
   subroutine reset_properties(solver)
!!
!!    Reset solver variables
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      solver%S    = zero
      solver%T_r  = zero
      solver%T_l  = zero
!
      solver%etaX = zero
      solver%csiX = zero
!
   end subroutine reset_properties
!
!
   subroutine cleanup_cc_properties(solver, wf) 
!!
!!    Cleanup
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      class(ccs) :: wf
!
      call wf%destruct_multipliers()
!
      if (allocated(solver%etaX)) &
         call mem%dealloc(solver%etaX, wf%n_es_amplitudes)
!
      if (allocated(solver%csiX)) &
         call mem%dealloc(solver%csiX, wf%n_es_amplitudes)
!
   end subroutine cleanup_cc_properties
!
!
   subroutine intialize_variables_cc_properties(solver, wf)
!!
!!    Allocate global variables
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      class(ccs) :: wf
!
      if (.not. allocated(solver%etaX)) &
         call mem%alloc(solver%etaX, wf%n_es_amplitudes)
!
      if (.not. allocated(solver%csiX)) &
         call mem%alloc(solver%csiX, wf%n_es_amplitudes)
!
   end subroutine intialize_variables_cc_properties
!
!
   subroutine print_banner_cc_properties(solver)
!!
!!    Print banner
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!              
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
!
   end subroutine print_banner_cc_properties
!
!
   subroutine print_settings_cc_properties(solver)
!!
!!    Print settings
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_properties) :: solver
!
      write(output%unit, '(t3,a)') '- Davidson CC properties solver settings:'
!
      write(output%unit,'(/t6,a26,a)') 'Operator:                 ', solver%operator_type
      write(output%unit,'(t6,a26,a)')  'Excitation:               ', solver%es_type 
      write(output%unit,'(t6,a,i1)')   'Number of singlet states: ', solver%n_singlet_states
!
   end subroutine print_settings_cc_properties
!
!
   subroutine print_summary_cc_properties(solver, output_type, i)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
      implicit  none
!
      class(cc_properties), intent(in) :: solver
!              
      character(len=*), intent(in) :: output_type
!
      integer, optional :: i
!
      if (output_type == 'header') then 
!              
         write(output%unit, '(/t3,a)') '- Property calculation results:'
!
      elseif (output_type == 'top') then
!              
         write(output%unit, '(/t6,a,a)') 'Operator component: ', solver%component(i)
!
         write(output%unit, '(/t6,a)')  '                Transition moments                             '
         write(output%unit, '(t6,a)')   '             -------------------------                         '
         write(output%unit, '(t6,a)')   'State         Right              Left                Strength  '
         write(output%unit, '(t6,a)')   '---------------------------------------------------------------'
!
      elseif (output_type == 'results') then
!
         write(output%unit, '(t6,i2,1x,f19.10,1x,f19.10,1x,f19.10)') i, solver%T_r, solver%T_l, solver%S
!
      elseif (output_type == 'bottom') then
!
         write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
      endif
!
   end subroutine print_summary_cc_properties
!
!
end module cc_properties_class
