module davidson_cc_ip_solver_class
!
!!
!!    Davidson coupled cluster ionized state solver class module
!!    Written by Sarai D. Folkestad, 2018
!!
!!    
!
   use davidson_cc_es_solver_class
!
   implicit none
!
   type, extends(davidson_cc_es_solver) :: davidson_cc_ip_solver
!
   contains
!
      procedure :: print_banner           => print_banner_davidson_cc_ip_solver
    ! procedure :: print_summary  => print_summary_davidson_cc_ip_solver
!
      procedure :: read_settings          => read_settings_davidson_cc_ip_solver
!
    ! procedure :: print_settings => print_settings_davidson_cc_ip_solver
!
      procedure :: set_start_vectors      => set_start_vectors_davidson_cc_ip_solver
      procedure :: set_projection_vector  => set_projection_vector_davidson_cc_ip_solver
!
   end type davidson_cc_ip_solver
!
!
contains
!
!
   subroutine print_banner_davidson_cc_ip_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(davidson_cc_ip_solver) :: solver 
!
      write(output%unit, '(//t3,a)') ':: Davidson coupled cluster ionized state solver'
      write(output%unit, '(t3,a)')   ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(/t3,a)')  'A Davidson CVS solver that calculates core excitation energies and the'
      write(output%unit, '(t3,a)')   'corresponding right eigenvectors of the Jacobian matrix, A. The eigenvalue problem'
      write(output%unit, '(t3,a)')   'is solved in a reduced space, the dimension of which is expanded'
      write(output%unit, '(t3,a)')   'until the convergence criteria are met. Bath orbitals and projection is used'
      write(output%unit, '(t3,a)')   'to obtain ionized states.'
!
      write(output%unit, '(/t3,a)')  'A complete description of the Davidson algorithm can be found in'
      write(output%unit, '(t3,a)')   'E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      flush(output%unit)
!
   end subroutine print_banner_davidson_cc_ip_solver
!
!
   subroutine read_settings_davidson_cc_ip_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_ip_solver) :: solver 
!
      integer(i15) :: n_specs, i, j
!
      character(len=100) :: line
!
      call move_to_section('cc excited state', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:19) == 'residual threshold:' ) then
!
            read(line(20:100), *) solver%residual_threshold
!
         elseif (line(1:17) == 'energy threshold:' ) then
!
            read(line(18:100), *) solver%eigenvalue_threshold
!
         elseif (line(1:15) == 'singlet states:' ) then
!
            read(line(16:100), *) solver%n_singlet_states
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         elseif (trim(line) == 'restart') then
!
            solver%do_restart = .true.
!
         endif
!
      enddo
!
   end subroutine read_settings_davidson_cc_ip_solver
!
!
   subroutine set_start_vectors_davidson_cc_ip_solver(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cc_ip_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      integer(i15) :: trial, count_start_vecs, a, i, ai
!
      if (allocated(solver%start_vectors)) then
!
!        Initial trial vectors given on input
!
         call mem%alloc(c_i, wf%n_amplitudes, 1)
!
         c_i = zero
         c_i(solver%start_vectors(1, 1), 1) = one
!
         call davidson%write_trial(c_i, 'rewind')
!
         do trial = 2, solver%n_singlet_states
!
            c_i = zero
            c_i(solver%start_vectors(trial, 1), 1) = one
!
            call davidson%write_trial(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_amplitudes, 1)
!
      else
!
!        Initial trial vectors given by Koopman
!
         call mem%alloc(c_i, wf%n_amplitudes, 1)
!
         count_start_vecs = 0
         a = wf%n_v
!
         do i = wf%n_o, 1, -1
!
            ai = wf%n_v*(i - 1) + a
            count_start_vecs = count_start_vecs  + 1
!
            if (count_start_vecs .gt. solver%n_singlet_states) exit 
!
            c_i = zero
            c_i(ai, 1) = one
!
            if (count_start_vecs == 1) then
!
               call davidson%write_trial(c_i, 'rewind')
!
            else
!
               call davidson%write_trial(c_i, 'append')
!
            endif

!
         enddo
!
         call mem%dealloc(c_i, wf%n_amplitudes, 1)
!
      endif
!
   end subroutine set_start_vectors_davidson_cc_ip_solver
!
!
   subroutine set_projection_vector_davidson_cc_ip_solver(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cc_ip_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: projector
!
      call mem%alloc(projector, wf%n_amplitudes, 1)
!
      call wf%get_ip_projector(projector)
!
      call davidson%set_projector(projector)
      call mem%dealloc(projector, wf%n_amplitudes, 1)
!
   end subroutine set_projection_vector_davidson_cc_ip_solver
!
!
end module davidson_cc_ip_solver_class
