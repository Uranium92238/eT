!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
submodule (ccs_class) file_handling_ccs
!
!!
!!    File handling submodule
!!
!!    Gathers routines that save wavefunction parameters to file,
!!    and reads them from file, plus other routines related to the
!!    handling of the files that belong to the wavefunction.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_ground_state_files_ccs(wf)
!!
!!    Initialize ground state files
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      integer, dimension(:), allocatable :: amplitude_block_sizes
!
      call wf%get_gs_amplitude_block_sizes(amplitude_block_sizes)
!
      wf%t_storer = amplitude_file_storer('t', amplitude_block_sizes)
!
      wf%tbar_storer = amplitude_file_storer('tbar', amplitude_block_sizes)
!
   end subroutine initialize_ground_state_files_ccs
!
!
   module subroutine initialize_excited_state_files_ccs(wf)
!!
!!    Initialize excited state files
!!    Written by Alexander C. Paul, 2019-2020
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=5) :: file_name
      integer :: state
      integer, dimension(:), allocatable :: amplitude_block_sizes
!
      call wf%get_es_amplitude_block_sizes(amplitude_block_sizes)
!
      if (.not. allocated(wf%r_storers)) then
!
         allocate(wf%r_storers(wf%n_singlet_states))
!
         do state = 1, wf%n_singlet_states
!
            write(file_name,'(a,i3.3)') 'r_', state
            wf%r_storers(state) = amplitude_file_storer(trim(file_name), &
                                                   amplitude_block_sizes)
!
         end do
!
      end if
!
      if (.not. allocated(wf%l_storers)) then
!
         allocate(wf%l_storers(wf%n_singlet_states))
!
         do state = 1, wf%n_singlet_states
!
            write(file_name,'(a,i3.3)') 'l_', state
            wf%l_storers(state) = amplitude_file_storer(trim(file_name), &
                                                   amplitude_block_sizes)
!
         end do
!
      end if
!
   end subroutine initialize_excited_state_files_ccs
!
!
   module subroutine initialize_triplet_excited_state_files_ccs(wf)
!!
!!    Initialize triplet excited state files
!!    Written by Alexander C. Paul, 2019-2020
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=13) :: file_name
      integer :: state
      integer, dimension(:), allocatable :: amplitude_block_sizes
!
      call wf%get_triplet_es_amplitude_block_sizes(amplitude_block_sizes)
!
      if (.not. allocated(wf%triplet_r_storers)) then
!
         allocate(wf%triplet_r_storers(wf%n_triplet_states))
!
         do state = 1, wf%n_triplet_states
!
            write(file_name,'(a,i3.3)') 'r_triplet_', state
            wf%triplet_r_storers(state) = amplitude_file_storer(trim(file_name), &
                                                            amplitude_block_sizes)
!
         end do
!
      end if
!
      if (.not. allocated(wf%triplet_l_storers)) then
!
         allocate(wf%triplet_l_storers(wf%n_triplet_states))
!
         do state = 1, wf%n_triplet_states
!
            write(file_name,'(a,i3.3)') 'l_triplet_', state
            wf%triplet_l_storers(state) = amplitude_file_storer(trim(file_name), &
                                                   amplitude_block_sizes)
!
         end do
!
      end if
!
   end subroutine initialize_triplet_excited_state_files_ccs
!
!
   module subroutine save_amplitudes_ccs(wf)
!!
!!    Save amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Alexander C. Paul, 2017-2022
!!
!!    File format: energy n_t1, t1
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(:), allocatable :: t
!
      call mem%alloc(t, wf%n_gs_amplitudes)
!
      call wf%get_amplitudes(t)
      call wf%t_storer%save_(t, wf%energy)
!
      call mem%dealloc(t, wf%n_gs_amplitudes)
!
   end subroutine save_amplitudes_ccs
!
!
   module subroutine read_amplitudes_ccs(wf, read_n)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Alexander C. Paul, 2017-2022
!!
!!    read_n: optionally returns the number of amplitudes read.
!!            This is especially useful e.g. in CCSD to provide a start guess
!!            for the doubles if only singles were found on file.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(out), optional :: read_n
!
      real(dp), dimension(:), allocatable :: t
      real(dp) :: energy
!
      call mem%alloc(t, wf%n_gs_amplitudes)
!
      call wf%t_storer%read_(t, energy, read_n)
      call wf%set_amplitudes(t)
!
      call mem%dealloc(t, wf%n_gs_amplitudes)
!
   end subroutine read_amplitudes_ccs
!
!
   module subroutine save_multipliers_ccs(wf)
!!
!!    Save multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Alexander C. Paul, 2017-2022
!!
!!    File format: energy n_t1, t1bar
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(:), allocatable :: tbar
!
      call mem%alloc(tbar, wf%n_gs_amplitudes)
!
      call wf%get_multipliers(tbar)
      call wf%tbar_storer%save_(tbar, wf%energy)
!
      call mem%dealloc(tbar, wf%n_gs_amplitudes)
!
   end subroutine save_multipliers_ccs
!
!
   module subroutine read_multipliers_ccs(wf, read_n)
!!
!!    Read multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Alexander C. Paul, 2017-2022
!!
!!    read_n: optionally returns the number of multipliers read.
!!            This is especially useful e.g. in CCSD to provide a start guess
!!            for the doubles if only singles were found on file.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(out), optional :: read_n
!
      real(dp), dimension(:), allocatable :: tbar
      real(dp) :: energy
!
      call mem%alloc(tbar, wf%n_gs_amplitudes)
!
      call wf%tbar_storer%read_(tbar, energy, read_n)
      call wf%set_multipliers(tbar)
!
      call mem%dealloc(tbar, wf%n_gs_amplitudes)
!
   end subroutine read_multipliers_ccs
!
!
   module subroutine save_excited_state_ccs(wf, X, first, last, side, energies)
!!
!!    Save excited state
!!    Written by Eirik F. Kjønstad and Alexander C. Paul, 2019-2022
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: first, last ! first, last state number
!
      real(dp), dimension(wf%n_es_amplitudes, first:last), intent(in) :: X
!
      character(len=*), intent(in) :: side ! 'left' or 'right'
!
      real(dp), dimension(first:last), optional, intent(in) :: energies
!
      integer :: state
!
      if (trim(side) .eq. 'right') then
!
         do state = first, last
!
            call wf%r_storers(state)%save_(X(:,state), energies(state))
!
         end do
!
      else if (trim(side) .eq. 'left') then
!
         do state = first, last
!
            call wf%l_storers(state)%save_(X(:,state), energies(state))
!
         end do
!
      endif
!
   end subroutine save_excited_state_ccs
!
!
   module subroutine save_triplet_excited_state_ccs(wf, X, first, last, side, energies)
!!
!!    Save triplet excited state
!!    Written by Eirik F. Kjønstad and Alexander C. Paul, 2019-2022
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: first, last ! first, last state number
!
      real(dp), dimension(wf%n_es_amplitudes, first:last), intent(in) :: X
!
      character(len=*), intent(in) :: side ! 'left' or 'right'
!
      real(dp), dimension(first:last), optional, intent(in) :: energies
!
      integer :: state
!
      if (trim(side) .eq. 'right') then
!
         do state = first, last
!
            call wf%triplet_r_storers(state)%save_(X(:,state), energies(state))
!
         end do
!
      else if (trim(side) .eq. 'left') then
!
         do state = first, last
!
            call wf%triplet_l_storers(state)%save_(X(:,state), energies(state))
!
         end do
!
      endif
!
   end subroutine save_triplet_excited_state_ccs
!
!
   module subroutine read_excited_state_ccs(wf, X, first, last, side, energies)
!!
!!    Read excited state
!!    Written by Eirik F. Kjønstad and Alexander C. Paul, 2019-2022
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: first, last ! first, last state number
!
      real(dp), dimension(wf%n_es_amplitudes, first:last), intent(out) :: X
      real(dp), dimension(first:last), optional,           intent(out) :: energies
!
      character(len=*), intent(in) :: side ! 'left' or 'right'
!
      real(dp), dimension(first:last) :: temp_energies
!
      integer :: state
!
      if (trim(side) .eq. 'right') then
!
         do state = first, last
!
            call wf%r_storers(state)%read_(X(:,state), temp_energies(state))
!
         end do
!
      elseif (trim(side) .eq. 'left') then
!
         do state = first, last
!
            call wf%l_storers(state)%read_(X(:,state), temp_energies(state))
!
         end do
!
      else
!
         call output%error_msg('Tried to read an excited state, &
                               &but argument side not recognized: ' // side)
!
      endif
!
      if (present(energies)) then
         energies = temp_energies
      end if
!
   end subroutine read_excited_state_ccs
!
!
   module subroutine read_triplet_excited_state_ccs(wf, X, first, last, side, energies)
!!
!!    Read triplet excited state
!!    Written by Eirik F. Kjønstad and Alexander C. Paul, 2019-2022
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: first, last ! first, last state number
!
      real(dp), dimension(wf%n_triplet_amplitudes, first:last), intent(out) :: X
      real(dp), dimension(first:last), optional,           intent(out) :: energies
!
      character(len=*), intent(in) :: side ! 'left' or 'right'
!
      real(dp), dimension(first:last) :: temp_energies
!
      integer :: state
!
      if (trim(side) .eq. 'right') then
!
         do state = first, last
!
            call wf%triplet_r_storers(state)%read_(X(:,state), temp_energies(state))
!
         end do
!
      elseif (trim(side) .eq. 'left') then
!
         do state = first, last
!
            call wf%triplet_l_storers(state)%read_(X(:,state), temp_energies(state))
!
         end do
!
      else
!
         call output%error_msg('Tried to read an excited state, &
                               &but argument side not recognized: ' // side)
!
      endif
!
      if (present(energies)) then
         energies = temp_energies
      end if
!
   end subroutine read_triplet_excited_state_ccs
!
!

   module subroutine check_and_get_restart_vector_ccs(wf, vector, energy, n, side, found)
!!
!!    Check and get restart vector
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Inquires if the nth file of "side" exists and
!!    returns the corresponding vector if the file exists.
!!    Otherwise checks the "other side" and tries to read it.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
      logical,  intent(out) :: found
!
      integer,          intent(in)  :: n
      character(len=*), intent(in)  :: side
      character(len=40) :: file_name
!
      found = .false.
!
!     File_name cannot be uninitialized
      file_name = ""
!
      if (trim(side) == 'right') then
!
         if (wf%r_storers(n)%file_exists()) then
!
            call wf%get_restart_vector(wf%r_storers(n), vector, energy, 'right', 'right')
            file_name = wf%r_storers(n)%get_filename()
            found = .true.
!
         else
!
            if (wf%l_storers(n)%file_exists()) then
!
               call wf%get_restart_vector(wf%l_storers(n), vector, energy, 'left', 'right')
               file_name = wf%l_storers(n)%get_filename()
               found = .true.
!
            end if
!
         end if
!
      else if (trim(side) == 'left') then
!
         if (wf%l_storers(n)%file_exists()) then
!
            call wf%get_restart_vector(wf%l_storers(n), vector, energy, 'left', 'left')
            file_name = wf%l_storers(n)%get_filename()
            found = .true.
!
         else
!
            if (wf%r_storers(n)%file_exists()) then
!
               call wf%get_restart_vector(wf%r_storers(n), vector, energy, 'right', 'left')
               file_name = wf%r_storers(n)%get_filename()
               found = .true.
!
            end if
!
         end if
!
      else
!
         call output%error_msg('Trying to restart an excited state, but &
                               &argument side not recognized: ' // trim(side))
!
      end if
!
      if (found) then
!
         call output%printf('n', 'Restarting '// trim(side) //' vector (i0) from file (a0).', &
                           ints=[n], chars=[trim(file_name)], fs='(t6,a)')
!
      end if
!
   end subroutine check_and_get_restart_vector_ccs
!
!
   module subroutine get_restart_vector_ccs(wf, storer, vector, energy, restart_from, restart_to)
!!
!!    Get restart vector
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Gets start vector and energy from storer and
!!    handles the basis transformations according to:
!!
!!    restart from "right" to "left"
!!    L^a_i = 2R^a_i
!!
!!    restart from "left" to "right"
!!    R^a_i = 1/2 L^a_i
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      type(amplitude_file_storer), intent(inout) :: storer
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
!
      real(dp), intent(out) :: energy
!
      character(len=*), intent(in) :: restart_from, restart_to
!
      if (restart_to == restart_from) then
!
         call storer%read_(vector, energy)
!
      else if (restart_from == 'left' .and. restart_to == 'right') then
!
         call storer%read_(vector, energy)
         call dscal(wf%n_es_amplitudes, two, vector, 1)
!
      else if (restart_from == 'right' .and. restart_to == 'left') then
!
         call storer%read_(vector, energy)
         call dscal(wf%n_es_amplitudes, half, vector, 1)
!
      end if
!
   end subroutine get_restart_vector_ccs
!
!
   module subroutine save_tbar_intermediates_ccs(wf)
!!
!!    Save tbar intermediates multiplier equation
!!    Written by Alexander C. Paul, Aug 2019
!!
      use warning_suppressor
!
      implicit none
!
      class(ccs) :: wf
!
!     Suppress unused variable compiler warning for wf
      call do_nothing(wf)
!
   end subroutine save_tbar_intermediates_ccs
!
!
   module subroutine remove_parallel_states_ccs(wf, threshold, transformation)
!!
!!    Remove parallel states
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Find parallel left and right vectors and remove them
!!    from the corresponding files and wf variables.
!!
!!    NB: This routine may only be called once during the execution of eT,
!!        otherwise n_singlet_states is reduced twice.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), intent(in) :: threshold
!
      character(len=*), intent(in) :: transformation
!
      real(dp), dimension(:), allocatable :: left_energies, right_energies
!
      logical, dimension(:), allocatable :: parallel
!
      integer :: n_parallel_left, n_parallel_right, n_parallel
!
      character(len=20) :: warning_print
!
      n_parallel = 0
      call mem%alloc(parallel, wf%n_singlet_states)
!
      if (trim(transformation) == 'left' .or. trim(transformation) == 'both') then
!
         call mem%alloc(left_energies, wf%n_singlet_states)
         call dcopy(wf%n_singlet_states, wf%left_excitation_energies, 1, left_energies, 1)
!
         call wf%check_for_parallel_states('left', threshold, parallel)
!
         call wf%remove_parallel_states_from_file(parallel, n_parallel_left, &
                                                  left_energies, wf%l_storers, &
                                                  'left')
!
         deallocate(wf%l_storers)
         call wf%destruct_left_excitation_energies
!
      end if
!
      if (trim(transformation) == 'right' .or. trim(transformation) == 'both') then
!
         call mem%alloc(right_energies, wf%n_singlet_states)
         call dcopy(wf%n_singlet_states, wf%right_excitation_energies, 1, right_energies, 1)
!
         call wf%check_for_parallel_states('right', threshold, parallel)
!
         call wf%remove_parallel_states_from_file(parallel, n_parallel_right, &
                                                  right_energies, wf%r_storers, &
                                                  'right')
!
         deallocate(wf%r_storers)
         call wf%destruct_right_excitation_energies
!
      end if
!
      call mem%dealloc(parallel, wf%n_singlet_states)
!
!     Reinitialize excited state variables
!
      if (trim(transformation) == 'left') then
!
         n_parallel = n_parallel_left
         wf%n_singlet_states = wf%n_singlet_states - n_parallel
!
         call wf%initialize_left_excitation_energies()
         call dcopy(wf%n_singlet_states, left_energies, 1, wf%left_excitation_energies, 1)
         call mem%dealloc(left_energies, wf%n_singlet_states+n_parallel)
!
      else if (trim(transformation) == 'right') then
!
         n_parallel = n_parallel_right
         wf%n_singlet_states = wf%n_singlet_states - n_parallel
!
         call wf%initialize_right_excitation_energies()
         call dcopy(wf%n_singlet_states, right_energies, 1, wf%right_excitation_energies, 1)
         call mem%dealloc(right_energies, wf%n_singlet_states+n_parallel_right)
!
      else if (trim(transformation) == 'both') then
!
         if (n_parallel_left .ne. n_parallel_right) then
!
            call output%printf('n', 'Obtained different number of parallel left/right states.')
!
            call output%error_msg('Inconsistent number of left and right states.')
!
         end if
!
         n_parallel = n_parallel_right
         wf%n_singlet_states = wf%n_singlet_states - n_parallel
!
         call wf%initialize_left_excitation_energies()
         call wf%initialize_right_excitation_energies()
         call dcopy(wf%n_singlet_states, left_energies, 1, wf%left_excitation_energies, 1)
         call dcopy(wf%n_singlet_states, right_energies, 1, wf%right_excitation_energies, 1)
         call mem%dealloc(left_energies, wf%n_singlet_states+n_parallel)
         call mem%dealloc(right_energies, wf%n_singlet_states+n_parallel)
!
      end if
!
      call wf%initialize_excited_state_files()
!
      if (n_parallel > 0) then
!
         if (n_parallel == 1) then
            write(warning_print,'(i0,a)') n_parallel, ' parallel state'
         else
            write(warning_print,'(i0,a)') n_parallel, ' parallel states'
         end if
!
         call output%warning_msg('(a0) removed. Number of states reduced to: (i0)', &
                                 chars=[trim(warning_print)], &
                                 ints=[wf%n_singlet_states])
      end if
!
   end subroutine remove_parallel_states_ccs
!
!
   module subroutine remove_parallel_states_from_file_ccs(wf, parallel, n_parallel, &
                                                          energies, es_storers, side)
!!
!!    Remove parallel states from file
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Send in vector of logicals determining which states are parallel.
!!    If a parallel states is found among the non-parallel ones,
!!    the file with the parallel state is overwritten with the following
!!    non-parallel states. Afterwards the last n_parallel files are deleted.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      logical, dimension(wf%n_singlet_states), intent(in) :: parallel
!
      integer, intent(out) :: n_parallel
!
      real(dp), dimension(wf%n_singlet_states), intent(out) :: energies
!
      character(len=*), intent(in) :: side
!
      type(amplitude_file_storer), dimension(wf%n_singlet_states), intent(inout) :: es_storers
!
      real(dp), dimension(:), allocatable :: vector
!
      integer :: state
!
      n_parallel = 0
!
      call mem%alloc(vector, wf%n_es_amplitudes)
!
      do state = 1, wf%n_singlet_states
!
         if (parallel(state)) then
!
            n_parallel = n_parallel + 1
!
            call output%printf('v', 'Removed (a0) state (i0) as it was parallel &
                                &to another state.', ints=[state], chars=[trim(side)])
!
            cycle
!
         else if (n_parallel > 0) then
!
            call wf%read_excited_state(vector, state, state, side, &
                                       energies(state-n_parallel))
            call wf%save_excited_state(vector, state-n_parallel, state-n_parallel, &
                                       side, energies(state-n_parallel))
!
         end if
!
      end do
!
      call mem%dealloc(vector, wf%n_es_amplitudes)
!
      do state = wf%n_singlet_states - n_parallel + 1, wf%n_singlet_states
!
         call es_storers(state)%delete_file()
!
      end do
!
   end subroutine remove_parallel_states_from_file_ccs
!
!
   module subroutine get_density_for_plotting_ccs(wf, c_D_ct, density, state_p, &
                                                  state_q, file_read)
!!
!!    Get density for plotting
!!    Written by Alexander C. Paul, May 2021
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: c_D_ct
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
      integer, intent(in) :: state_p, state_q
!
      logical, intent(out) :: file_read
!
      character(len=10) :: file_name
      type(stream_file) :: density_file
!
      write(file_name, '(a, i3.3, a, i3.3)') 'dm_', state_p, '_', state_q
!
      density_file = stream_file(file_name)
!
      file_read = .false.
!
      if (density_file%exists()) then
!
         call density_file%open_()
         call density_file%read_(density, wf%n_mo**2)
!
         file_read = .true.
         call wf%add_t1_terms_and_transform(density, c_D_ct)
!
         call density_file%close_('keep')
!
      end if
!
   end subroutine get_density_for_plotting_ccs
!
!
end submodule file_handling_ccs
