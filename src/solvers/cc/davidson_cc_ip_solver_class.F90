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

!
   contains
!
      procedure :: read_settings          => read_settings_davidson_cc_ip_solver
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
   subroutine read_settings_davidson_cc_ip_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_ip_solver) :: solver 
!
      integer :: n_specs, i
!
      character(len=100) :: line
!
      solver%tag = 'Davidson coupled cluster ionized state solver'
      solver%description1 = 'A Davidson CVS solver that calculates core ionization energies &
                            &and the corresponding right eigenvectors of the Jacobian matrix, & 
                            &A. The eigenvalue problem is solved in a reduced space, the &
                            &dimension of which is expanded until the convergence criteria & 
                            &are met. Bath orbitals and projection is used to obtain ionized &
                            &states.'
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
      integer :: trial, count_start_vecs, a, i, ai
!
      if (allocated(solver%start_vectors)) then
!
!        Initial trial vectors given on input
!
         call mem%alloc(c_i, wf%n_es_amplitudes, 1)
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
         call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
!
      else
!
!        Initial trial vectors given by Koopman
!
         call mem%alloc(c_i, wf%n_es_amplitudes, 1)
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
         call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
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
      call mem%alloc(projector, wf%n_es_amplitudes, 1)
!
      call wf%get_ip_projector(projector)
!
      call davidson%set_projector(projector)
      call mem%dealloc(projector, wf%n_es_amplitudes, 1)
!
      if (.false.) write(output%unit, *) solver%tag ! Hack to suppress unavoidable compiler warnings
!
   end subroutine set_projection_vector_davidson_cc_ip_solver
!
!
end module davidson_cc_ip_solver_class
