module cc_property_solver_class
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
   type :: cc_property_solver 
!
      character(len=100) :: tag = 'Coupled cluster property solver'
      character(len=100) :: author = 'Josefine H. Andersen, 2019'
!
      character(len=500) :: description1 = 'A solver that calculates spectral intensities from coupled &
                                           &cluster excited state calculations'
!
      character(len=40)   :: X_operator
!
      real(dp), allocatable :: S
!
      real(dp), dimension(:,:), allocatable :: etaX
      real(dp), dimension(:,:), allocatable :: csiX
!
      integer :: n_singlet_states = 0
!
   contains
!
      procedure, non_overridable :: prepare          => prepare_cc_property_solver
      procedure, non_overridable :: run              => run_cc_property_solver
      procedure, non_overridable :: cleanup          => cleanup_cc_property_solver
!
      procedure :: read_settings                     => read_settings_cc_property_solver
!
      procedure :: print_banner                      => print_banner_cc_property_solver
      procedure :: print_summary                     => print_summary_cc_property_solver
!
   end type cc_property_solver
!
!
contains
!
!
   subroutine prepare_cc_property_solver(solver)
!!
!!    Prepare 
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_property_solver) :: solver
!
      call solver%print_banner()
!
      call solver%read_settings()
!
      solver%S = zero
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
      class(cc_property_solver) :: solver
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: etaX
      real(dp), dimension(:,:), allocatable :: csiX
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: Y
      real(dp), dimension(:,:), allocatable :: Z
!
      call mem%alloc(etaX, wf%n_amplitudes, 1)
      call mem%alloc(csiX, wf%n_amplitudes, 1)
      call mem%alloc(X, wf%n_amplitudes, 1)
      call mem%alloc(Y, wf%n_amplitudes, 1)
      call mem%alloc(Z, wf%n_amplitudes, 1)
!
!      call wf%get_transformed_dipole_operator(X, Y, Z)
      ! call wf%construct_etaX(wf, etaX, Xoperator)      
      ! call wf%construct_csiX(wf, csiX, Xoperator)      
!
      call solver%print_summary(wf)
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
      class(cc_property_solver) :: solver
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
         endif
!
      enddo
!     
      if (.not. requested_section('cc properties')) then
!
      call output%error_msg('Operator type must be specified for property calculations.')
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
      if (n_ops .gt. 1) then
!
         call output%error_msg('More than one operator specified')
!
      endif
!      
      read(input%unit, '(a100)') line
      line = remove_preceding_blanks(line)
!
      if (line(1:13) == 'dipole length') then
!
         solver%X_operator = 'dipole_length'
!         
      endif
!
   end subroutine read_settings_cc_property_solver
!
!
   subroutine cleanup_cc_property_solver(solver) 
!!
!!    Cleanup
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_property_solver) :: solver
!
!     What to clean up?
!
   end subroutine cleanup_cc_property_solver
!
!
   subroutine print_banner_cc_property_solver(solver)
!!
!!    Print banner
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_property_solver) :: solver
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
!
   end subroutine print_banner_cc_property_solver
!
!
   subroutine print_summary_cc_property_solver(solver, wf)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
      implicit  none
!
      class(cc_property_solver), intent(in) :: solver
!
      class(ccs), intent(in) :: wf
!
      integer :: state
!
!     Would be nice to have access to excitation energies as to not only
!     assign an intensity to electronic state number
!
      write(output%unit, '(/t3,a)') '- Summary of spectra calculation:'
!
      write(output%unit, '(/t6,a)') '                                    '
      write(output%unit, '(t6,a)')  'State        Excitation energy(Hartree)           Strength     '
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
      do state = 1, solver%n_singlet_states
!
         write(output%unit, '(t6,i2,14x,f19.12,4x,f19.12)') state, 'energy', solver%S
!
      enddo
!
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
   end subroutine print_summary_cc_property_solver
!
!
end module cc_property_solver_class
