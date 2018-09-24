module davidson_cvs_cc_es_solver_class
!
!!
!!    Davidson CVS coupled cluster excited state solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    
!
   use davidson_cc_es_solver_class
!
   implicit none
!
   type, extends(davidson_cc_es_solver) :: davidson_cvs_cc_es_solver
!
      integer(i15) :: n_cores
!
      integer(i15), dimension(:,:), allocatable :: cores
!
   contains
!
      procedure :: print_banner => print_banner_davidson_cvs_cc_es_solver
      !procedure :: print_summary  => print_summary_davidson_cc_es_solver
!
      procedure :: read_settings  => read_settings_davidson_cvs_cc_es_solver
!
      !procedure :: print_settings => print_settings_davidson_cc_es_solver
!
      procedure :: set_start_vectors         => set_start_vectors_davidson_cvs_cc_es_solver
      !procedure :: set_precondition_vector   => set_precondition_vector_davidson_cc_es_solver
      !procedure :: set_projection_vector     => set_projection_vector_davidson_cc_es_solver
!
   end type davidson_cvs_cc_es_solver
!
!
contains
!
!
   subroutine print_banner_davidson_cvs_cc_es_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(davidson_cvs_cc_es_solver) :: solver 
!
      write(output%unit, '(//t3,a)') ':: Davidson core-valence separation coupled cluster excited state solver'
      write(output%unit, '(t3,a)')   ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(/t3,a)')  'A Davidson CVS solver that calculates core excitation energies and the'
      write(output%unit, '(t3,a)')   'corresponding right eigenvectors of the Jacobian matrix, A. The eigenvalue problem'
      write(output%unit, '(t3,a)')   'is solved in a reduced space, the dimension of which is expanded'
      write(output%unit, '(t3,a)')   'until the convergence criteria are met. In addition the CVS aproximation is used'
      write(output%unit, '(t3,a)')   'to obtain the core excitations'
!
      write(output%unit, '(/t3,a)')  'A complete description of the Davidson algorithm can be found in'
      write(output%unit, '(t3,a)')   'E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
      write(output%unit, '(/t3,a)')  'A description of the CVS approximation can be found in'
      write(output%unit, '(t3,a)')   ' S. Coriani & H. Koch, J. Chem. Phys. 143, 181103 (2015).'
!
      flush(output%unit)
!
   end subroutine print_banner_davidson_cvs_cc_es_solver
!
!
   subroutine read_settings_davidson_cvs_cc_es_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cvs_cc_es_solver) :: solver 
!
      integer(i15) :: n_specs, i, j
!
      character(len=100) :: line
!
      call move_to_section('excited state', n_specs)
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
         elseif (line(1:25) == 'number of singlet states:' ) then
!
            read(line(26:100), *) solver%n_singlet_states
!
         elseif (line(1:16) == 'core excitation:' ) then
!
            line = line(17:100)
            line = remove_preceding_blanks(line)
            solver%n_cores = 0
!
            do j = 1, 83
!
               if (line(j:j) .ne. ' ') solver%n_cores = solver%n_cores + 1
!
            enddo
!
            call mem%alloc_int(solver%cores, solver%n_cores, 1)
!
            read(line, *) solver%cores
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         elseif (trim(line) == 'restart') then
!
            solver%restart = .true.
!
         endif
!
      enddo
!
   end subroutine read_settings_davidson_cvs_cc_es_solver
!
!
!
!
   subroutine set_start_vectors_davidson_cvs_cc_es_solver(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cvs_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: orbital_differences
      real(dp), dimension(:,:), allocatable :: lowest_orbital_differences
!
      integer(i15), dimension(:,:), allocatable :: lowest_orbital_differences_index
!
      integer(i15) :: trial, i, j, k, first_ao_on_atom, last_ao_on_atom
!
      real(dp) :: mix_factor
!
   !  if (solver%n_cores .gt. solver%n_singlet_states) &
   !     call output%error_msg('number of roots requested should be equal or greater than the number of cores.')
!
   !  if (allocated(solver%start_vectors)) then
!
!  !     Initial trial vectors given on input
!
   !     call mem%alloc(c_i, wf%n_amplitudes, 1)
!
   !     c_i = zero
   !     c_i(solver%start_vectors(1, 1), 1) = one
!
   !     call davidson%write_trial(c_i, 'rewind')
!
   !     do trial = 2, solver%n_singlet_states
!
   !        c_i = zero
   !        c_i(solver%start_vectors(trial, 1), 1) = one
!
   !        call davidson%write_trial(c_i)
!
   !     enddo
!
   !     call mem%dealloc(c_i, wf%n_amplitudes, 1)
!
   !  else
!
!  !     Initial trial vectors given by Koopman
!
!  !     Calculate the mixing factor of equal mix 
!
   !     mix_factor = 1.0d0/(sqrt(real(solver%n_cores, kind=dp)))
!
!  !     Loop through the occupied MOs and determine if they are core mos
!
   !     n_MOs_found = 0
!
   !     call mem%alloc_int(core_MOs, solver%n_cores, 1)
!
   !     do i = 1, wf%n_o
!
   !        do j = 1, solver%n_cores
!
   !           first_ao_on_atom = wf%system%atoms(solver%cores(i))%shells(1)%first
   !           last_ao_on_atom = wf%system%atoms(solver%cores(i))%shells(wf%system%atoms(solver%cores(i))%n_shells)%last
!
   !           do k = first_ao_on_atom, last_ao_on_atom
!
   !              if (wf%orbital_coefficients(k, i) .ge. mix_factor) then
!
   !                 n_MOs_found = n_MOs_found + 1
!
   !                 if (n_MOs_found .gt. solver%n_cores) &
   !                          call output%error_msg('something went wrong in the selection of core MOs.')
!
   !                 core_MOs(n_MOs_found, 1) = i
!
   !              endif
!
   !           enddo
!
   !        enddo
!
   !     enddo
!
!  !     Calculate ai indices 
!
   !     call mem%alloc_int(ai_indices, solver%n_singlet_states, 1)
!
!  !     Set c(ai) = 1
!
   !     call mem%alloc(c_i, wf%n_amplitudes, 1)
!
   !     c_i = zero
!
   !     call davidson%write_trial(c_i, 'rewind')
!
   !     do trial = 1, solver%n_singlet_states
!
!
   !     enddo
!
   !     call mem%dealloc(c_i, wf%n_amplitudes, 1)
   !     call mem%dealloc_int(core_MOs, solver%n_cores, 1)
!
   !  endif
!
   end subroutine set_start_vectors_davidson_cvs_cc_es_solver
!
!
end module davidson_cvs_cc_es_solver_class
