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
      integer(i15), dimension(:,:), allocatable :: core_MOs
!
   contains
!
      procedure :: read_settings          => read_settings_davidson_cvs_cc_es_solver
!
      procedure :: set_start_vectors      => set_start_vectors_davidson_cvs_cc_es_solver
      procedure :: set_projection_vector  => set_projection_vector_davidson_cvs_cc_es_solver
!
      procedure :: initialize_core_MOs    => initialize_core_MOs_davidson_cvs_cc_es_solver
      procedure :: initialize_cores       => initialize_cores_davidson_cvs_cc_es_solver
!
      procedure :: destruct_core_MOs      => destruct_core_MOs_davidson_cvs_cc_es_solver
      procedure :: destruct_cores         => destruct_cores_davidson_cvs_cc_es_solver
!
   end type davidson_cvs_cc_es_solver
!
!
contains
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
      solver%tag = 'Davidson coupled cluster ionized state solver'
      solver%description1 = 'A Davidson CVS solver that calculates core excitation energies and the &
                            &corresponding right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                            &problem is solved in a reduced space, the dimension of which is expanded &
                            &until the convergence criteria are met. In addition the CVS aproximation is &
                            &used to obtain the core excitations'
      solver%description2 = 'A complete description of the Davidson algorithm can be found in &
                            &E. R. Davidson, J. Comput. Phys. 17, 87 (1975). &
                            &A description of the CVS approximation can be found in &
                            &S. Coriani & H. Koch, J. Chem. Phys. 143, 181103 (2015).'
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
            call solver%initialize_cores()
!
            read(line, *) solver%cores
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
   end subroutine read_settings_davidson_cvs_cc_es_solver
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
!
      integer(i15), dimension(:,:), allocatable :: ai_indices
!
      integer(i15) :: trial, core, i, j, k, l, a, first_ao_on_atom, last_ao_on_atom
      integer(i15) :: n_MOs_found, current_root
!
      real(dp) :: mix_factor
!
      logical :: all_selected, used
!
      if (solver%n_cores .gt. solver%n_singlet_states) &
         call output%error_msg('number of roots requested should be equal or greater than the number of cores.')
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
!        Calculate the mixing factor of equal mix 
!
         mix_factor = 1.0d0/(sqrt(real(solver%n_cores, kind=dp)))*0.9
!
!        Loop through the occupied MOs and determine if they are core mos
!
         n_MOs_found = 0
!
         call solver%initialize_core_MOs()
         solver%core_MOs = 0
!
         do j = 1, solver%n_cores
!
            first_ao_on_atom = wf%system%atoms(solver%cores(j, 1))%shells(1)%first
            last_ao_on_atom = wf%system%atoms(solver%cores(j, 1))%shells(wf%system%atoms(solver%cores(j, 1))%n_shells)%last
!
            do k = first_ao_on_atom, last_ao_on_atom

               do i = 1, wf%n_o
!
                  if (abs(wf%orbital_coefficients(k, i)) .ge. mix_factor) then
!
                     used =  .false.
!
                     do l = 1, n_MOs_found
!
                        if (solver%core_MOs(l, 1) == i) used = .true.
!
                     enddo
!
                     if (.not. used) then
!
                        n_MOs_found = n_MOs_found + 1
!
                        if (n_MOs_found .gt. solver%n_cores) &
                              call output%error_msg('something went wrong in the selection of core MOs.')
!

!
                        solver%core_MOs(n_MOs_found, 1) = i
                        exit
!
                     else
!
                        cycle
!
                     endif
!
                  endif
!
               enddo
!
            enddo
!
         enddo
!
!        Calculate ai indices 
!
         call mem%alloc(ai_indices, solver%n_singlet_states, 1)
!
         all_selected = .false.
         a =  0 
         current_root = 0
!
         do while (.not. all_selected)
!
            a = a + 1
!
            do core = 1, solver%n_cores
!
               i = solver%core_MOs(core, 1)
!
               current_root = current_root + 1
               ai_indices(current_root, 1) = wf%n_v*( i - 1) + a
!
               if (current_root .eq. solver%n_singlet_states) then
!
                  all_selected = .true.
                  exit
!
               endif
!
            enddo
!
         enddo
!
!        Set c(ai) = 1
!
         call mem%alloc(c_i, wf%n_es_amplitudes, 1)
!
         c_i = zero
!
         c_i(ai_indices(1, 1), 1) = one
!
         call davidson%write_trial(c_i, 'rewind')
!
         do trial = 2, solver%n_singlet_states
!
            c_i = zero
!
            c_i(ai_indices(trial, 1), 1) = one
!
            call davidson%write_trial(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes, 1)
!
      endif
!
   end subroutine set_start_vectors_davidson_cvs_cc_es_solver
!
!
   subroutine set_projection_vector_davidson_cvs_cc_es_solver(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cvs_cc_es_solver) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: projector
!
      call mem%alloc(projector, wf%n_es_amplitudes, 1)
!
      call wf%get_cvs_projector(projector, solver%n_cores, solver%core_MOs)
!
      call davidson%set_projector(projector)
      call mem%dealloc(projector, wf%n_es_amplitudes, 1)
!
   end subroutine set_projection_vector_davidson_cvs_cc_es_solver
!
!
   subroutine initialize_core_MOs_davidson_cvs_cc_es_solver(solver)
!!
!!    Initialize core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(davidson_cvs_cc_es_solver) :: solver
!
      if (.not. allocated(solver%core_MOs)) call mem%alloc(solver%core_MOs, solver%n_cores, 1)
!
   end subroutine initialize_core_MOs_davidson_cvs_cc_es_solver
!
!
   subroutine destruct_core_MOs_davidson_cvs_cc_es_solver(solver)
!!
!!    Destruct core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(davidson_cvs_cc_es_solver) :: solver
!
      if (allocated(solver%core_MOs)) call mem%dealloc(solver%core_MOs, solver%n_cores, 1)
!
   end subroutine destruct_core_MOs_davidson_cvs_cc_es_solver
!
!
   subroutine initialize_cores_davidson_cvs_cc_es_solver(solver)
!!
!!    Initialize cores
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(davidson_cvs_cc_es_solver) :: solver
!
      if (.not. allocated(solver%cores)) call mem%alloc(solver%cores, solver%n_cores, 1)
!
   end subroutine initialize_cores_davidson_cvs_cc_es_solver
!
!
   subroutine destruct_cores_davidson_cvs_cc_es_solver(solver)
!!
!!    Destruct cores
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(davidson_cvs_cc_es_solver) :: solver
!
      if (allocated(solver%cores)) call mem%dealloc(solver%cores, solver%n_cores, 1)
!
   end subroutine destruct_cores_davidson_cvs_cc_es_solver
!
!
end module davidson_cvs_cc_es_solver_class
