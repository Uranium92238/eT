module cc_property_solver_class
!
!!
!!    Coupled cluster property solver class module
!!    Written by Josefine H. Andersen, 2019
!!
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
      character(len=500) :: description1 = 'A solver that calculates spectral  intesities from coupled &
                                           &cluster excited state calculations'
!
      real(dp), allocatable :: S
!
   contains
!
      procedure, non_overridable :: prepare          => prepare_cc_property_solver
      procedure, non_overridable :: run              => run_cc_property_solver
      procedure, non_overridable :: cleanup          => cleanup_cc_property_solver
!
      procedure :: print_banner                      => print_banner_cc_property_solver
      procedure :: print_summary                     => print_summary_cc_property_solver
!
      procedure :: construct_etaX                    => construct_etaX_cc_property_solver
      procedure :: construct_csiX                    => construct_csiX_cc_property_solver
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
      call read_settings()
!
      call print_settings()
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

!
   end subroutine run_cc_property_solver
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
      call solver%print_summary(solver,wf)
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
   subroutine print_summary_cc_property(solver, wf)
!!
!!    Print summary
!!    Written by Josefine H. Andersen
!!
      implicit  none
!
      class(cc_property_solver) :: solver
!
      class(ccs) :: wf
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
      do state = 1, wf%n_singlet_states
!
         write(output%unit, '(t6,i2,14x,f19.12,4x,f19.12)') state, 'energy', solver%get_S
!
      enddo
!
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!

!
   end subroutine print_summary_cc_property
!
!
   subroutine construct_etaX_cc_property_solver(solver, wf)
!!
!!    Construct left-hand-side vector etaX
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_property_solver) :: solver
!
      class(ccs) :: wf
!
!     Here call sequence of wf functions
!
   end subroutine construct_etaX_cc_property_solver
!
!
   subroutine construct_csiX_cc_property_solver(solver, wf) 
!!
!!    Construct right-hand-side vector csiX
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(cc_property_solver) :: solver
!
      class(ccs) :: wf
!
!     Here call sequence of wf functions
!
   end subroutine construct_csiX_cc_property_solver
!
!
end module cc_property_solver_class
