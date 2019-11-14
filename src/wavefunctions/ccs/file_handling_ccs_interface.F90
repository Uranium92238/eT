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
   module subroutine initialize_files_ccs(wf)
!!
!!    Initialize files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Initializes the wavefucntion files for wavefunction parameters.
!!
      class(ccs) :: wf 
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
   end subroutine initialize_cc_files_ccs
!
!
   module subroutine initialize_excited_state_files_ccs(wf, transformation)
!!
!!    Initialize files for excited state vectors and energies
!!    Written by Alexander C. Paul, Oct 2019 
!!
      class(ccs), intent(inout) :: wf
      character(len=*), intent(in) :: transformation
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
      real(dp), dimension(wf%n_t1), intent(out) :: X 
      type(sequential_file), intent(inout) :: file_
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
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X 
      integer, intent(in) :: n ! state number 
      character(len=*), intent(in) :: side ! 'left' or 'right' 
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
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: X 
      integer, intent(in) :: n ! state number 
      character(len=*), intent(in) :: side ! 'left' or 'right' 
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
      integer, intent(in) :: n_states ! number of states
      real(dp), dimension(n_states), intent(in) :: energies
      character(len=*), intent(in) :: r_or_l 
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
      integer, intent(in) :: n_states ! Obtained by reading the restart file 
      real(dp), dimension(n_states), intent(out) :: energies
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
      character(len=*), intent(in) :: side 
      integer :: n_states
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
   end subroutine save_tbar_intermediates_ccs
