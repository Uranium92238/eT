submodule (cc2_class) excited_state
!
!!
!!    Excited state submodule (CC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Contains the following family of procedures of the CC2 class:
!!
!!    inititialize_excited_states: Initializes number of s2 amplitudes (n_s2am), and adds it n_parameters
!!    calculate_orbital_differences: Calculates the orbital differences, including the double excitation differences
!!                                   in the active CC2 spaces
!!    transform_trial_vectors: Transforms the new trial vectors. rho = Ac
!!
!!    Upper case indices are general indices, lower case indices are restricted
!!    to the CC2 orbital space.
!!
!
   implicit none
!
   logical :: debug   = .false.
   logical :: timings = .false.
!
!
contains
!
!
   module subroutine calculate_orbital_differences_cc2(wf, orbital_diff)
!!
!!    Calculate Orbital Differences (CC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!    Calculates orbital differences
!!
!!       1) ε_I^A = ε_A - ε_I
!!       2) ε_ij^ab = ε_a + ε_b - ε_i - ε_j (for active spaces only)
!!
!!    and puts them in orbital_diff, which is a vector of length n_parameters.
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
!
      integer(i15) :: A = 0, I = 0, b = 0, j = 0
      integer(i15) :: AI = 0, bj = 0
      integer(i15) :: aibj = 0
!
      do I = 1, wf%n_o
         do A = 1, wf%n_v
!
            AI = index_two(A, I, wf%n_v)
!
            orbital_diff(AI, 1) = wf%fock_diagonal(A + wf%n_o, 1) - wf%fock_diagonal(I, 1)
!
         enddo
      enddo
!
      do i = 1, wf%n_o
!
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do j = 1, wf%n_o
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  orbital_diff((wf%n_o)*(wf%n_v) + aibj, 1) &
                                                 = wf%fock_diagonal(wf%n_o + a , 1) &
                                                 - wf%fock_diagonal(i, 1) &
                                                 + wf%fock_diagonal(wf%n_o + b, 1) &
                                                 - wf%fock_diagonal(j, 1)
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine calculate_orbital_differences_cc2
!
!
   module subroutine transform_trial_vectors_cc2(wf, first_trial, last_trial)
!!
!!    Transformation of Trial Vectors (CC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Each trial vector in first_trial to last_trial is read from file and
!!    transformed before the transformed vector is written to file.
!!
!!    Singles and doubles part of the transformed vectors are written to
!!    the same record in file transformed_vec, record length is n_parameters long.
!!
      implicit none
!
      class(cc2) :: wf
!
      integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      real(dp), dimension(:,:), allocatable :: c_a_i
      real(dp), dimension(:,:), allocatable :: c_aibj
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
      integer(i15) :: trial = 0
!
!
!     Allocate c_a_i and c_aibj
!
      call wf%mem%alloc(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero
!
      call wf%mem%alloc(c_aibj, wf%n_s2am, 1)
      c_aibj = zero
!
!     Open trial vector- and transformed vector files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*wf%n_parameters, iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*wf%n_parameters, iostat=ioerror)
!
!     For each trial vector: read, transform and write
!
      do trial = first_trial, last_trial
!
         read(unit_trial_vecs, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
         if (wf%tasks%current == 'excited_state') then
!
            if (wf%excited_state_specifications%right) then
!
               call wf%jacobian_cc2_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_specifications%left) then
!
               write(unit_output,*)'Error: Jacobian transpose transformation not implemented for CC2'
               stop
!
            else
!
               write(unit_output,*) 'Error: Excited state task not recognized'
               stop
!
            endif
!
         elseif (wf%tasks%current == 'multipliers') then
!
               write(unit_output,*)'Error: multipliers not implemented for CC2'
               stop
!
         else
!
            write(unit_output,*) 'Error: Current task not recognized'
            stop
!
         endif
!
!        -::- Projections -::-
!
!        Test for core calculation
!
         if (wf%tasks%core_excited_state .or. wf%tasks%core_ionized_state) then
!
!           Project out contamination from valence contributions
!
            call wf%cvs_rho_a_i_projection(c_a_i)
            call wf%cvs_rho_aibj_projection(c_aibj)
!
         endif
!
!        Test for ionization calculation
!
         if (wf%tasks%ionized_state .or. wf%tasks%core_ionized_state) then
!
!           Project out contamination from regular excitations
!
            write(unit_output,*)'Error: Ionized state not implemented for CC2'
            stop
!
         endif
!
!
         write(unit_rho, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
      enddo

!
!     Close files
!
      close(unit_trial_vecs)
      close(unit_rho)
!
!     Deallocate c_a_i and c_aibj
!
      call wf%mem%dealloc(c_a_i, wf%n_v, wf%n_o)
      call wf%mem%dealloc(c_aibj, wf%n_s2am, 1)
!
   end subroutine transform_trial_vectors_cc2
!
!
   module subroutine cvs_residual_projection_cc2(wf, residual)
!!
!!    Residual projection (CC2),
!!    Written by Sarai D. Folkestad Aug. 2017
!!
      implicit none
!
      class(cc2) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital

      do i = 1, wf%n_o
!
         core_orbital = .false.
         do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
            if (i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) core_orbital = .true.
!
         enddo
!
         if (.not. core_orbital) then
            do a = 1, wf%n_v
               ai = index_two(a, i, wf%n_v)
               residual(ai, 1) = zero
            enddo
         endif
!
      enddo
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            core_orbital = .false.
            do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
               if ((i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) .or. &
                  (j .eq. wf%core_excited_state_specifications%index_core_mo(core, 1))) core_orbital = .true.
!
            enddo
!
            if (.not. core_orbital) then
               do a = 1, wf%n_v
                  do b = 1, wf%n_v
                     ai = index_two(a, i, wf%n_v)
                     bj = index_two(b, j, wf%n_v)
                     aibj = index_packed(ai, bj)

                     residual(wf%n_t1am + aibj, 1) = zero
                  enddo
               enddo
            endif
         enddo
      enddo
!
    end subroutine cvs_residual_projection_cc2
!
!
   module subroutine excited_state_preparations_cc2(wf)
!!
!!    Excited State Preparations (CC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.
!!
      class(cc2) :: wf
!
!     Change the number of parameters for the excited state problem
!
      wf%n_s2am = ((wf%n_v)*(wf%n_o))*((wf%n_v )*(wf%n_o)+1)/2
      wf%n_parameters = wf%n_parameters + wf%n_s2am
!
!     Store vvov-electronic repulsion integrals to file if there is space
!
      call wf%mem%alloc(wf%t1am, wf%n_v, wf%n_o)
      call wf%read_single_amplitudes
!
      call wf%store_t1_vv_ov_electronic_repulsion
!
!     Set current task to excited state calculation
!
      wf%tasks%current = 'excited_state'
!
!     Set filename for solution vectors
!
      if (wf%tasks%core_excited_state .or. wf%tasks%core_ionized_state) then   ! Core excitation
!
         if (wf%excited_state_specifications%right) then                         ! Right vectors
            wf%excited_state_specifications%solution_file = 'right_core'
         else                                                                    ! Left vectors
            write(unit_output,*)'Error: Jacobian transpose transformation not implemented for core excitations' ! S: should be able to get these with the same projections however so...
            stop
         endif
!
      else                                                                    ! Valence excitation
!
         if (wf%excited_state_specifications%left) then                          ! Right vectors
            wf%excited_state_specifications%solution_file = 'left_valence'
         else                                                                    ! Left vectors
            wf%excited_state_specifications%solution_file = 'right_valence'
         endif
!
      endif
!
      call deallocator(wf%t1am, wf%n_v, wf%n_o)
!
   end subroutine excited_state_preparations_cc2
!
!
   module subroutine summary_excited_state_info_cc2(wf, energies)
!!
!!    Summary of excited state info (CC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017 / May 2018
!!
!!    Prints the analysis of the excitation vectors
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%excited_state_specifications%n_singlet_states, 1) :: energies
!
      type(file) :: solution_file
!
      integer(i15) :: state = 0
!
      real(dp), dimension(:,:), allocatable :: solution_ai, solution_aibj
!
      real(dp) :: norm_sq_singles ! Note: the total vector is normalized in this case
!
!     Open solution vector file
!
      solution_file%name = wf%excited_state_specifications%solution_file
      call wf%disk%open_file(solution_file, 'unformatted', 'read', 'direct', dp*(wf%n_parameters))
!
!     Allocations
!
      call wf%mem%alloc(solution_ai, wf%n_t1am, 1)
      call wf%mem%alloc(solution_aibj, wf%n_s2am, 1)
!
      do state = 1, wf%excited_state_specifications%n_singlet_states
!
         write(unit_output,'(/t3,a30,i3,a1/)') 'Analysis of excitation vector ',state, ':'
         write(unit_output,'(t6, a, f14.8)')    'Excitation energy [a.u.]:   ', energies(state,1)
         write(unit_output,'(t6, a, f14.8/)')   'Excited state energy [a.u.]:', wf%energy + energies(state,1)
!
!        Read the solution
!
         solution_ai   = zero
         solution_aibj = zero
!
         read(solution_file%unit, rec=state) solution_ai, solution_aibj
!
!        Print dominant single & double excitations
!
         norm_sq_singles = dot_product(solution_ai, solution_ai, wf%n_t1am)
         call print_dominant_two_index(solution_ai, wf%n_v, wf%n_o, 'a', 'i')
         call print_dominant_four_index(solution_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o, &
                                                         'a', 'i', 'b', 'j')
         write(unit_output,'(/t6,a41,f14.12/)') &
               'Singles contribution to the full vector: ', norm_sq_singles
!
      enddo
!
!     Deallocations
!
      call wf%mem%dealloc(solution_ai, wf%n_t1am, 1)
      call wf%mem%dealloc(solution_aibj, wf%n_s2am, 1)
!
      call wf%disk%close_file(solution_file)
!
   end subroutine summary_excited_state_info_cc2
!
!
!
!
end submodule
