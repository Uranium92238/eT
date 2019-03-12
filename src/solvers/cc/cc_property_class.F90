module cc_property_class
!
!!
!!    Coupled cluster property solver class module
!!    Written by Josefine H. Andersen, 2019
!!    
!
   use kinds
   use file_class
   use ccs_class
   use eigen_davidson_tool_class
!
   implicit none
!
   type :: cc_property
!
      character(len=100) :: tag = 'Coupled cluster property solver'
      character(len=100) :: author = 'Josefine H. Andersen, 2019'
!
      character(len=500) :: description1 = 'A solver that calculates spectral intensities from coupled &
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
      real(dp), dimension(:,:), allocatable :: l_vec_n
      real(dp), dimension(:,:), allocatable :: r_vec_n
!
      integer :: n_singlet_states = 0
!
   contains
!
      procedure, non_overridable :: prepare          => prepare_cc_property_solver
      procedure, non_overridable :: run              => run_cc_property_solver
      procedure, non_overridable :: cleanup          => cleanup_cc_property_solver
!
      procedure :: intialize_variables               => intialize_variables_cc_property_sover
      procedure :: read_settings                     => read_settings_cc_property_solver
!
      procedure :: reset                             => reset_property_solver
      procedure :: do_eom_or_lr                      => do_eom_or_lr_property_solver
!
      procedure :: print_banner                      => print_banner_cc_property_solver
      procedure :: print_summary                     => print_summary_cc_property_solver
!
   end type cc_property
!
!
contains
!
!
   subroutine prepare_cc_property_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
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
   end subroutine prepare_cc_property_solver
!
!
   subroutine run_cc_property_solver(solver, wf) 
!!
!!    Run
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
!
      class(ccs) :: wf
!
      integer :: i, n
!
      !real(dp), dimension(:,:), allocatable :: l_vec_n
      !real(dp), dimension(:,:), allocatable :: r_vec_n
!
!     Read multipliers and prepare operator
!
      call wf%initialize_multipliers()
      call wf%read_multipliers()
!
      call wf%prepare_operator_pq(solver%operator_type)
!
      call solver%print_summary(wf, 'header')
!      
!     Loop over components of X
!
      do i = 1, 3
!
      !   call solver%print_summary('top', i)
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
         call solver%print_summary(wf, 'top', i)
!
!        Loop over excited states
!  
         do n = 1, solver%n_singlet_states
!
            call wf%read_excited_state(solver%l_vec_n, n, 'left')
            call wf%read_excited_state(solver%r_vec_n, n, 'right')
!
            call wf%scale_left_excitation_vector(solver%l_vec_n, solver%r_vec_n)
!
            call wf%calculate_transition_strength(solver%S, solver%etaX, &
                               solver%csiX, solver%l_vec_n, solver%r_vec_n)
!
            call solver%print_summary(wf, 'results', n)
!            
         enddo
!         
         call solver%print_summary(wf, 'bottom')
!
      enddo
!
   end subroutine run_cc_property_solver
!
!
   subroutine read_settings_cc_property_solver(solver)
!!
!!    Read cc excited state settings to find n_singlet_states
!!    and cc properties settings to determine operator
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
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
         call output%error_msg('operator type must be specified for property calculations')
!
      endif
!
      call move_to_section('cc properties', n_ops)
!
!      if (n_ops .gt. 1) then
!
!         call output%error_msg('More than one operator specified')
!
!      endif
!      
      do i = 1, n_ops
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:13) == 'dipole length') then
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
   end subroutine read_settings_cc_property_solver
!
!
   subroutine reset_property_solver(solver)
!!
!!    Reset solver variables
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
!
      solver%S = zero
!
      solver%etaX = zero
      solver%csiX = zero
!
   end subroutine reset_property_solver
!
!
   subroutine cleanup_cc_property_solver(solver, wf) 
!!
!!    Cleanup
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
!
      class(ccs) :: wf
!
!     What to clean up?
      call wf%destruct_multipliers()
!
      if (allocated(solver%etaX)) &
         call mem%dealloc(solver%etaX, wf%n_es_amplitudes, 1)
!
      if (allocated(solver%csiX)) &
         call mem%dealloc(solver%csiX, wf%n_es_amplitudes, 1)
!
      if (allocated(solver%l_vec_n)) &
         call mem%dealloc(solver%l_vec_n, wf%n_es_amplitudes, 1)
!
      if (allocated(solver%r_vec_n)) &
         call mem%dealloc(solver%r_vec_n, wf%n_es_amplitudes, 1)
!
   end subroutine cleanup_cc_property_solver
!
!
   subroutine intialize_variables_cc_property_sover(solver, wf)
!!
!!    Allocate global variables
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
!
      class(ccs) :: wf
!
      if (.not. allocated(solver%etaX)) &
         call mem%alloc(solver%etaX, wf%n_es_amplitudes, 1)
!
      if (.not. allocated(solver%csiX)) &
         call mem%alloc(solver%csiX, wf%n_es_amplitudes, 1)
!
      if (.not. allocated(solver%l_vec_n)) &
         call mem%alloc(solver%l_vec_n, wf%n_es_amplitudes, 1)
!
      if (.not. allocated(solver%r_vec_n)) &
         call mem%alloc(solver%r_vec_n, wf%n_es_amplitudes, 1) 
!
   end subroutine intialize_variables_cc_property_sover
!
!
   subroutine print_banner_cc_property_solver(solver)
!!
!!    Print banner
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      !class(cc_property_solver) :: solver
      class(cc_property) :: solver
!              
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
!
   end subroutine print_banner_cc_property_solver
!
!
   subroutine print_summary_cc_property_solver(solver, wf, output_type, i)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
!!    The routine takes input, as the results are written out in the loops in run()
!!
      implicit  none
!
      !class(cc_property_solver), intent(in) :: solver
      class(cc_property), intent(in) :: solver
!
      class(ccs), intent(in) :: wf
!              
      character(len=*), intent(in) :: output_type
!
      integer, optional :: i
!
      real(dp) :: ddot
!
      if (output_type == 'header') then 
!              
         write(output%unit, '(t3,a,a)') 'Type of excitation: ', solver%es_type
         write(output%unit, '(t3,a,a)') 'Type of operator: ', solver%operator_type
         write(output%unit, '(/t3,a)') '- Property calculation results:'
!
      elseif (output_type == 'top') then
!              
        ! write(output%unit, '(t6,a)') '                                    '
         write(output%unit, '(/t6,a,a)') 'Operator component: ', solver%component(i)
        ! write(output%unit, '(/t6,a)') '                                    '
         write(output%unit, '(/t6,a)')  'State         etaX*R             L*csiX              Strength  '
         write(output%unit, '(t6,a)')   '---------------------------------------------------------------'
!
      elseif (output_type == 'results') then
!
         write(output%unit, '(t6,i2,1x,f19.10,1x,f19.10,1x,f19.10)') i, &
         ddot(wf%n_es_amplitudes, solver%etaX, 1, solver%r_vec_n, 1), &
         ddot(wf%n_es_amplitudes, solver%l_vec_n, 1, solver%csiX, 1), &
         solver%S
!
      elseif (output_type == 'bottom') then
!
         write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
      endif
!
   end subroutine print_summary_cc_property_solver
!
!
   subroutine do_eom_or_lr_property_solver(solver, wf)
!!
!!    Calls construction of EOM or LR contribution.
!!    Based on wf, gives the correct input to wf%get_eom_contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      !class(cc_property_solver), intent(inout) :: solver
      class(cc_property), intent(inout) :: solver
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
                                      & will be calculated with no contrubutions'
      endif
!
   end subroutine do_eom_or_lr_property_solver
!
!
end module cc_property_class
