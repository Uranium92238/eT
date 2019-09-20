!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
!!    File handling submodule (CCSD)
!!    Set up by Andreas Skeidsvoll, Aug 2019
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
      call wf%initialize_wavefunction_files()
      call wf%initialize_singles_files()
      call wf%initialize_cc_files()
!
   end subroutine initialize_files_ccs
!
!
   module subroutine initialize_singles_files_ccs(wf)
!!
!!    Initialize singles files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      class(ccs) :: wf 
!
      wf%t1_file = sequential_file('t1')
      wf%t1bar_file = sequential_file('t1bar')
!
      wf%r1_file = sequential_file('r1')
      wf%l1_file = sequential_file('l1')
!
   end subroutine initialize_singles_files_ccs
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
   module subroutine save_amplitudes_ccs(wf)
!!
!!    Save amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
      call wf%t1_file%open_('write', 'rewind')
!
      call wf%t1_file%write_(wf%t1, wf%n_t1)
!
      call wf%t1_file%close_()
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
      call wf%t1_file%open_('read', 'rewind')
!
      call wf%t1_file%read_(wf%t1, wf%n_t1)
!
      call wf%t1_file%close_()
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
      call wf%t1bar_file%open_('write', 'rewind')
!
      call wf%t1bar_file%write_(wf%t1bar, wf%n_t1)
!
      call wf%t1bar_file%close_()
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
      call wf%t1bar_file%open_('read', 'rewind')
!
      call wf%t1bar_file%read_(wf%t1bar, wf%n_t1)
!
      call wf%t1bar_file%close_()
!
   end subroutine read_multipliers_ccs
!
!
   module subroutine save_singles_vector_ccs(wf, X, n, file_)
!!
!!    Save singles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!    Writes singles vector "X" to the sequential
!!    and unformatted file "file_".
!!    
!!    NB! If n = 1, then the routine WILL REWIND the file before writing,
!!    thus DELETING every record in the file. For n >=2, we just append to
!!    the file. The purpose of this setup is that the files should be saved in 
!!    the correct order, from n = 1 to n = # states.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t1), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(sequential_file) :: file_
!
      call file_%open_('write', 'append')
!
      if (n .eq. 1) then
         call file_%rewind_()
      endif
!
      call file_%write_(X, wf%n_t1)
!
      call file_%close_()
!
   end subroutine save_singles_vector_ccs
!
!
   module subroutine read_singles_vector_ccs(wf, X, n, file_)
!!
!!    Read singles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!    Reads singles vector "X" from the "n"'th line
!!    of the sequential and unformatted file "file_".
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(out) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(sequential_file) :: file_
!
      call file_%open_('read', 'rewind')
!
      call file_%skip(n-1)
!
      call file_%read_(X, wf%n_t1)
!
      call file_%close_()
!
   end subroutine read_singles_vector_ccs
!
!
   module subroutine save_excited_state_ccs(wf, X, n, side)
!!
!!    Save excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Saves an excited state to disk. 
!!    Since the solvers  keep these vectors in full length, 
!!    we receive a vector in full length (n_es_amplitudes), 
!!    and then distribute the different parts of that vector 
!!    to singles, doubles, etc., files (if there are doubles, etc.).
!!
!!    NB! If n = 1, then the routine WILL REWIND the file before writing,
!!    thus DELETING every record in the file. For n >=2, we just append to
!!    the file. The purpose of this setup is that the files should be saved in 
!!    the correct order, from n = 1 to n = # states.
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
      if (trim(side) == 'right') then 
!
         call wf%save_singles_vector(X, n, wf%r1_file)
!
      elseif (trim(side) == 'left') then 
!
         call wf%save_singles_vector(X, n, wf%l1_file)
!
      else
!
         call output%error_msg('Tried to save an excited state, but argument side not recognized: ' // side)
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
!!
!!    Reads an excited state to disk. Since this routine is used by 
!!    solvers, it returns the vector in the full space. Thus, we open 
!!    files for singles, doubles, etc., paste them together, and return 
!!    the result in X.
!!
!!    NB! This will place the cursor of the file at position n + 1.
!!    Be cautious when using this in combination with writing to the files.
!!    We recommend to separate these tasks---write all states or read all
!!    states; don't mix if you can avoid it.
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
      if (trim(side) == 'right') then 
!
         call wf%read_singles_vector(X, n, wf%r1_file)
!
      elseif (trim(side) == 'left') then 
!
         call wf%read_singles_vector(X, n, wf%l1_file)
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
      if (trim(r_or_l) == 'right') then 
!
         call wf%initialize_right_excitation_energies()
         wf%right_excitation_energies = energies 
!
      elseif (trim(r_or_l) == 'left') then 
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
      implicit none 
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_states ! Obtained by reading the restart file 
!                                     ! or by calling read_n_excitation_energies 
!
      real(dp), dimension(n_states), intent(out) :: energies
!
      integer :: local_n_states
!
      call wf%excitation_energies_file%open_('read', 'rewind')
!
      call wf%excitation_energies_file%read_(local_n_states)
!
      if (local_n_states .ne. n_states) then
!
         call output%error_msg('Dimension of excited state array does not match what is on file.')
!
      endif
!
      call wf%excitation_energies_file%read_(energies, n_states)
!
      call wf%excitation_energies_file%close_()
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
      class(ccs) :: wf 
!
      character(len=*), intent(in) :: side 
!
      integer :: n_states 
!
      n_states = 0
!
      if (trim(side) == 'right') then 
!
         n_states = wf%r1_file%get_size()
         n_states = n_states/(dp*wf%n_t1)
!
      elseif (trim(side) == 'left') then 
!
         n_states = wf%l1_file%get_size()
         n_states = n_states/(dp*wf%n_t1)
!
      else
!
         call output%error_msg('Tried to compute number of excited states. Unrecognized _side_: ' // side)
!
      endif
!
   end function get_n_excited_states_on_file_ccs
!
!
   module subroutine save_tbar_intermediates_ccs(wf)
!!
!!    Save tbar intermediates multiplier equation
!!    Written by Alexander Paul, Aug 2019
!!
      implicit none
!
      class(ccs) :: wf
!
!     For now, do nothing.
!
      call output%printf('No intermediates to save in (a0)', chars=[trim(wf%name_)])
!
   end subroutine save_tbar_intermediates_ccs
!
end submodule file_handling_ccs
