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
      real(dp), dimension(:,:), allocatable :: etaX
      real(dp), dimension(:,:), allocatable :: csiX
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
      procedure :: do_eom_or_lr                      => do_eom_or_lr_properties
!
      procedure :: print_banner                      => print_banner_cc_properties
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
      call solver%print_summary(wf, 'header')
!      
!     Loop over components of X operator
!
      do i = 1, 3
!
         call solver%print_summary(wf, 'top', i)
!         
         solver%X = solver%component(i) //'_'// trim(solver%operator_type) 
!
         call solver%reset()
!
         call wf%construct_etaX(solver%X, solver%etaX)      
!
         call wf%construct_csiX(solver%X, solver%csiX)      
!
!        Do EOM or linear response
!
         call solver%do_eom_or_lr(wf)
!
!        Loop over excited states and calculate transition strength
!  
         do n = 1, solver%n_singlet_states
!
            call wf%calculate_transition_strength(solver%S, solver%etaX, &
                               solver%csiX, n, solver%T_l, solver%T_r)
!
            call solver%print_summary(wf, 'results', n)
!            
         enddo
!         
         call solver%print_summary(wf, 'bottom')
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
         call output%error_msg('number of excitations must be specified.')
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
         call output%error_msg('operator type must be specified for properties')
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
            read(line(10:100), *) solver%operator_type
!
         elseif (line(1:13) == 'dipole length') then
!
            solver%operator_type = 'dipole_length'
!        
         elseif (line(1:3) == 'eom') then
!
            solver%eom = .true.
!
         elseif (line(1:15) == 'linear response') then
!
            solver%linear_response = .true.
!
            call output%error_msg('Linear response not implemented for spectra calculations')
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
         call mem%dealloc(solver%etaX, wf%n_es_amplitudes, 1)
!
      if (allocated(solver%csiX)) &
         call mem%dealloc(solver%csiX, wf%n_es_amplitudes, 1)
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
         call mem%alloc(solver%etaX, wf%n_es_amplitudes, 1)
!
      if (.not. allocated(solver%csiX)) &
         call mem%alloc(solver%csiX, wf%n_es_amplitudes, 1)
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
   subroutine print_summary_cc_properties(solver, wf, output_type, i)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
      implicit  none
!
      class(cc_properties), intent(in) :: solver
!
      class(ccs), intent(in) :: wf
!              
      character(len=*), intent(in) :: output_type
!
      integer, optional :: i
!
      if (output_type == 'header') then 
!              
         write(output%unit, '(t3,a,a)') 'Type of excitation: ', solver%es_type
         write(output%unit, '(t3,a,a)') 'Type of operator: ', solver%operator_type
         write(output%unit, '(/t3,a)') '- Property calculation results:'
!
      elseif (output_type == 'top') then
!              
         write(output%unit, '(/t6,a,a)') 'Operator component: ', solver%component(i)
!
         write(output%unit, '(/t6,a)')  'State         etaX*R             L*csiX              Strength  '
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
   subroutine do_eom_or_lr_properties(solver, wf)
!!
!!    Calls construction of EOM or LR contribution.
!!    Based on wf, gives the correct input to wf%get_eom_contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(cc_properties), intent(inout) :: solver
!
      class(ccs), intent(in) :: wf
!
!     Do EOM or LR
!      
      if (solver%eom) then
!
            call wf%get_eom_contribution(solver%etaX, solver%csiX, solver%X)
!         
      elseif (solver%linear_response) then
!
         write(output%unit, '(t6,a)') 'Linear response has been selected but is not implemented. &
                                      & etaX will be calculated with no contribution '
!
      else
!
         write(output%unit, '(t6,a)') 'You have not specified EOM or LR. etaX &
                                      & will be calculated with no contributions'
      endif
!
   end subroutine do_eom_or_lr_properties
!
!
end module cc_properties_class
