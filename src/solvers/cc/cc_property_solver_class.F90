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
      procedure :: print_banner                      => print_banner_cc_property_solver
      procedure :: print_summary                     => print_summary_cc_property_solver
!
      procedure :: transition_strength               => transition_strength_property_solver
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
!      call read_settings()
!
!      call print_settings()
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

      call solver%print_summary(wf)
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
   subroutine transition_strength_property_solver(solver, csiX, etaX, wf)
!!
!!    Calculate transition strength
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(cc_property_solver) :: solver
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(in) :: etaX
      real(dp), dimension(wf%n_amplitudes, 1), intent(in) :: csiX
!
      real(dp), dimension(:,:), allocatable :: T_left, T_right
!
      integer :: state
!
      call mem%alloc(T_left, wf%n_amplitudes, 1)
      call mem%alloc(T_right, wf%n_amplitudes, 1)
!
      do state = 1, solver%n_singlet_states
!
      ! get right and left excitation vector
      !
      ! calc dotproducts btwn exc. vectors and csiX/etaX
      ! T_left = ddot(wf%n_parameters, etaX, 1, right_j, 1)
      ! T_right = ddot(wf%n_parameters, left_j, 1, csiX, 1)
      !
      ! sum S over three components
      ! solver%S += T_left * T_right
!
      enddo
!
   end subroutine transition_strength_property_solver
!
!
end module cc_property_solver_class
