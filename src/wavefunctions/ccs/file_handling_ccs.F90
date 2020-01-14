!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   module subroutine initialize_files_ccs(wf)
!!
!!    Initialize files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Initializes the wavefucntion files for wavefunction parameters.
!!
      class(ccs) :: wf 
!
      call wf%initialize_cc_files()
      call wf%initialize_ground_state_files()
!
   end subroutine initialize_files_ccs
!
!
   module subroutine initialize_ground_state_files_ccs(wf)
!!
!!    Initialize singles files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      class(ccs) :: wf 
!
      wf%t_file = sequential_file('t')
      wf%tbar_file = sequential_file('tbar')
!
   end subroutine initialize_ground_state_files_ccs
!
!
   module subroutine initialize_cc_files_ccs(wf)
!!
!!    Initialize singles files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      class(ccs) :: wf 
!
      wf%restart_file = sequential_file('cc_restart_file')
!
      wf%excitation_energies_file = sequential_file('excitation_energies')
!
   end subroutine initialize_cc_files_ccs
!
!
   module subroutine initialize_excited_state_files_ccs(wf, transformation)
!!
!!    Initialize files for excited state vectors and energies
!!    Written by Alexander C. Paul, Oct 2019 
!!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      if(trim(transformation) .eq. 'right') then
!
         wf%r_files = file_storer('r_', wf%n_es_amplitudes,       &
                                         wf%n_singlet_states,     &
                                         delete=.false.,          &
                                         direct_=.false.)
!
         wf%excitation_energies_file = sequential_file('excitation_energies')
!
         call wf%r_files%initialize_storer()
!
      else if(trim(transformation) .eq. 'left') then
!
         wf%l_files = file_storer('l_', wf%n_es_amplitudes,       &
                                         wf%n_singlet_states,     &
                                         delete=.false.,          &
                                         direct_=.false.)
!
         wf%excitation_energies_file = sequential_file('excitation_energies')
!
         call wf%l_files%initialize_storer()
!
      else if(trim(transformation) .eq. 'both') then 
!
         wf%r_files = file_storer('r_', wf%n_es_amplitudes,       &
                                         wf%n_singlet_states,     &
                                         delete=.false.,          &
                                         direct_=.false.)
!
         wf%l_files = file_storer('l_', wf%n_es_amplitudes,       &
                                         wf%n_singlet_states,     &
                                         delete=.false.,          &
                                         direct_=.false.)
!
         wf%excitation_energies_file = sequential_file('excitation_energies')
!
         call wf%l_files%initialize_storer()
         call wf%r_files%initialize_storer()
!
      else
!
         call output%error_msg('Tried to initialize files for excited states &
                              & but argument ' // trim(transformation) //    &
                              ' not recognized')
!
      end if
!
   end subroutine initialize_excited_state_files_ccs
!
!
   module subroutine read_singles_vector_ccs(wf, X, file_)
!!
!!    Read singles vector X from a file
!!    Written by Alexander C. Paul, Oct 2019 
!!
!!    Files are written with the singles part in the first record
!!    and the doubles part in the second. 
!!    Thus, read wf%n_o*wf%n_v elements.
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t1), intent(out) :: X 
!
      type(sequential_file), intent(inout) :: file_
!
      call file_%open_('read', 'rewind')
!
      call file_%read_(X, wf%n_t1)
!
      call file_%close_()
!
   end subroutine read_singles_vector_ccs
!
!
   module subroutine save_amplitudes_ccs(wf)
!!
!!    Save amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
      call wf%t_file%open_('write', 'rewind')
!
      call wf%t_file%write_(wf%t1, wf%n_t1)
!
      call wf%t_file%close_()
!
   end subroutine save_amplitudes_ccs
!
!
   module subroutine read_amplitudes_ccs(wf)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call wf%t_file%open_('read', 'rewind')
!
      call wf%t_file%read_(wf%t1, wf%n_t1)
!
      call wf%t_file%close_()
!
   end subroutine read_amplitudes_ccs
!
!
   module subroutine save_multipliers_ccs(wf)
!!
!!    Save multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
      call wf%tbar_file%open_('write', 'rewind')
!
      call wf%tbar_file%write_(wf%t1bar, wf%n_t1)
!
      call wf%tbar_file%close_()
!
   end subroutine save_multipliers_ccs
!
!
   module subroutine read_multipliers_ccs(wf)
!!
!!    Read multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
      call wf%tbar_file%open_('read', 'rewind')
!
      call wf%tbar_file%read_(wf%t1bar, wf%n_t1)
!
      call wf%tbar_file%close_()
!
   end subroutine read_multipliers_ccs
!
!
   module subroutine save_excited_state_ccs(wf, X, n, side)
!!
!!    Save excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019
!!    modified by Alexander C. Paul, Oct 2019
!!
!!    Saves an excited state to disk. 
!!    Since the solvers  keep these vectors in full length, 
!!    we save the vector in full length (n_es_amplitudes), 
!!
!!    Uses file_storer to distinguish different states
!!
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
      if (trim(side) .eq. 'right') then 
!
         call wf%r_files%set(x, n)
!
      elseif (trim(side) .eq. 'left') then 
!
         call wf%l_files%set(x, n)
!
      else
!
         call output%error_msg('Tried to save an excited state, &
                              & but argument side not recognized: ' // side)
!
      endif
!
   end subroutine save_excited_state_ccs
!
!
   module subroutine read_excited_state_ccs(wf, X, n, side)
!!
!!    Read excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019
!!    modified by Alexander C. Paul, Oct 2019
!!
!!    Reads an excited state from disk. Since this routine is used by 
!!    solvers, it returns the vector in the full space.
!!
!!    Uses file_storer to distinguish different states
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: X 
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
      if (trim(side) .eq. 'right') then 
!
         call wf%r_files%get(x, n)
!
      elseif (trim(side) .eq. 'left') then 
!
         call wf%l_files%get(x, n)
!
      else
!
         call output%error_msg('Tried to read an excited state, but argument side not recognized: ' // side)
!
      endif
!
   end subroutine read_excited_state_ccs
!
!
   module subroutine save_excitation_energies_ccs(wf, n_states, energies, r_or_l)
!!
!!    Save excitation energies 
!!    Written by Sarai D. Folkestad, Mar 2019 
!!
!!    Saves 'n_states' excitation energies to disk & in memory. 
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_states ! number of states
!
      real(dp), dimension(n_states), intent(in) :: energies
!
      character(len=*), intent(in) :: r_or_l 
!
      if (trim(r_or_l) .eq. 'right') then 
!
         call wf%initialize_right_excitation_energies()
         wf%right_excitation_energies = energies 
!
      elseif (trim(r_or_l) .eq. 'left') then 
!
         call wf%initialize_left_excitation_energies()
         wf%left_excitation_energies = energies 
!
      else
!
         call output%error_msg('Could not recognize transformation in save_excitation_energies_ccs')
!
      endif 
!
      call wf%excitation_energies_file%open_('write', 'rewind')
!
      call wf%excitation_energies_file%write_(n_states)
      call wf%excitation_energies_file%write_(energies, n_states)
!
      call wf%excitation_energies_file%close_()
!     
   end subroutine save_excitation_energies_ccs
!
!
   module subroutine read_excitation_energies_ccs(wf, n_states, energies)
!!
!!    Read excitation energies 
!!    Written by Sarai D. Folkestad, Mar 2019 
!!
!!    Reads excitation energies from file
!!
!!    Note: "n_states" gives the dimension of the array "energies".
!!          It should match the actual number of states on file,
!!          given by the first record, "n_states" before the routine is called
!!          by either reading the restart file or by calling the function 
!!          read_n_excitation_energies
!!
!!    n_states: number of states found on file
!!              Obtained by the number of existing records in the file storer
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_states
!
      real(dp), dimension(:), intent(out) :: energies
!
      integer :: dim_
!
      dim_ = size(energies)
!
      call wf%excitation_energies_file%open_('read', 'rewind')
!
!     First line of excitation_energies_file contains n_states
!
      call wf%excitation_energies_file%skip(1)
!
      if (dim_ .gt. n_states) then
!
         call wf%excitation_energies_file%read_(energies(1:n_states), n_states)
         call wf%excitation_energies_file%close_()
!
      else if (dim_ .le. n_states) then
!
         call wf%excitation_energies_file%read_(energies, dim_)
         call wf%excitation_energies_file%close_()
!
      endif
!     
   end subroutine read_excitation_energies_ccs
!
!
   integer module function get_n_excitation_energies_on_file_ccs(wf)
!!
!!    Read n excitation energies 
!!    Written by Sarai D. Folkestad, Mar 2019 
!!
!!    Reads and returns the number of excitation energies on file
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
      call wf%excitation_energies_file%open_('read', 'rewind')
!
      call wf%excitation_energies_file%read_(get_n_excitation_energies_on_file_ccs)
!
      call wf%excitation_energies_file%close_()
!     
   end function get_n_excitation_energies_on_file_ccs
!
!
   module function get_n_excited_states_on_file_ccs(wf, side) result(n_states)
!!
!!    Get number of excited states on file 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Figures out the number of excited states on file, 
!!    using the r1 and l1 files. This should be sufficient for 
!!    all coupled cluster models (i.e., it is most likely  
!!    unneccessary to overwrite this routine in descendants)
!!
      implicit none
!
      class(ccs) :: wf 
!
      character(len=*), intent(in) :: side 
!
      integer :: n_states
!
      n_states = 0
!
      if (trim(side) .eq. 'right') then
!
         n_states = wf%r_files%get_n_existing_records()
!
      else if (trim(side) .eq. 'left') then
!
         n_states = wf%l_files%get_n_existing_records()
!
      else 
!
         call output%error_msg('Tried to compute number of excited states. &
                              & Unrecognized _side_: ' // side)
!
      end if
!
   end function get_n_excited_states_on_file_ccs

!
!
   module subroutine save_tbar_intermediates_ccs(wf)
!!
!!    Save tbar intermediates multiplier equation
!!    Written by Alexander C. Paul, Aug 2019
!!
      implicit none
!
      class(ccs) :: wf
!
!     For now, do nothing.
!
      call output%printf('v', 'No intermediates to save in (a0)', chars=[trim(wf%name_)])
!
   end subroutine save_tbar_intermediates_ccs
!
!
end submodule file_handling_ccs
