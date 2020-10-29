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
   module subroutine initialize_excited_state_files_ccs(wf)
!!
!!    Initialize files for excited state vectors and energies
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    Modified by Alexander C. Paul, May 2020: array of stream files
!!
      class(ccs), intent(inout) :: wf
!
   end subroutine initialize_excited_state_files_ccs
!
!
   module subroutine read_singles_vector_ccs(wf, file_, vector, read_n)
!!
!!    Read singles vector
!!    Written by Alexander C. Paul, Oct 2019 
!!
!!    read_n: optionally adds the number of amplitudes read to read_n
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_t1), intent(out) :: vector
      integer, intent(inout) :: read_n
!
   end subroutine read_singles_vector_ccs
!
!
   module subroutine save_singles_vector_ccs(wf, file_, vector)
!!
!!    Save singles vector
!!    Written by Alexander C. Paul, Oct 2019 
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_t1), intent(in) :: vector 
!
   end subroutine save_singles_vector_ccs
!
!
   module subroutine save_amplitudes_ccs(wf)
!!
!!    Save amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    File format: energy n_t1, t1
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
   end subroutine save_amplitudes_ccs
!
!
   module subroutine read_amplitudes_ccs(wf, read_n)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!    Adapted to return the number of read amplitdues if requested 
!!    by Alexander C. Paul, Oct 2020
!!
!!    read_n: optionally returns the number of amplitudes read. 
!!            This is especially useful e.g. in CCSD to provide a start guess 
!!            for the doubles if only singles were found on file.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      integer, intent(out), optional :: read_n
!
   end subroutine read_amplitudes_ccs
!
!
   module subroutine save_multipliers_ccs(wf)
!!
!!    Save multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    File format: energy n_t1, t1bar
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
!
   end subroutine save_multipliers_ccs
!
!
   module subroutine read_multipliers_ccs(wf, read_n)
!!
!!    Read multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!    Adapted to return the number of read multipliers if requested 
!!    by Alexander C. Paul, Oct 2020
!!
!!    read_n: optionally returns the number of multipliers read. 
!!            This is especially useful e.g. in CCSD to provide a start guess 
!!            for the doubles if only singles were found on file.
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
      integer, intent(out), optional :: read_n
!
   end subroutine read_multipliers_ccs
!
!
   module subroutine save_excited_state_ccs(wf, X, first, last, side, energies)
!!
!!    Save excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019
!!    Modified by Alexander C. Paul, Oct 2019
!!    Modified by Eirik F. Kjønstad, Mar 2020
!!
!!    Writes excited states in the columns of X to disk.
!!
!!    first: first state to write
!!    last:  last state to write
!!    side:  'right' or 'left' depending on right or left excited states 
!!
!!    Modified by Eirik F. Kjønstad, Mar 2020: made changes for direct stream,
!!                                             and added [first, last] range 
!!    Modified by Alexander C. Paul, May 2020: introduced array of stream files
!!    for the excited states
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
      integer, intent(in) :: first, last ! first, last state number 
      real(dp), dimension(wf%n_es_amplitudes, first:last), intent(in) :: X 
      character(len=*), intent(in) :: side ! 'left' or 'right'
      real(dp), dimension(first:last), optional, intent(in) :: energies
!
   end subroutine save_excited_state_ccs
!
!
   module subroutine save_excitation_vector_on_file_ccs(wf, file_, vector, energy)
!!
!!    Save excitation vector on file 
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Writes excitation vector o file structured as follows:
!!    excitation_energy, n_t1, X1
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: vector
      real(dp), intent(in) :: energy
!
   end subroutine save_excitation_vector_on_file_ccs
!
!
   module subroutine read_excited_state_ccs(wf, X, first, last, side, energies)
!!
!!    Read excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019
!!    modified by Alexander C. Paul, Oct 2019
!!    modified by Eirik F. Kjønstad, Mar 2020
!!
!!    Reads excited states from disk into the columns of X.
!!
!!    first: first state to read
!!    last:  last state to read
!!    side:  'right' or 'left' depending on right or left excited states 
!!
!!    Modified by Eirik F. Kjønstad, Mar 2020: made changes for direct stream,
!!                                             and added [first, last] range 
!!    Modified by Alexander C. Paul, May 2020: introduced array of stream files
!!    for the excited states.
!!    Modified by Alexander C. Paul, Sep 2020: Files now contain energy and n_es
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
      integer, intent(in) :: first, last ! first, last state number 
      real(dp), dimension(wf%n_es_amplitudes, first:last), intent(out) :: X
      real(dp), dimension(first:last), optional,           intent(out) :: energies
      character(len=*), intent(in) :: side ! 'left' or 'right'
!
   end subroutine read_excited_state_ccs
!
!
   module subroutine read_excitation_vector_file_ccs(wf, file_, vector, energy, read_n)
!!
!!    Read excitation vector file 
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Reads excitation vector from file structured as follows:
!!    excitation_energy, n_t1, X1
!!
!!    read_n: optionally returns the number of amplitudes read
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
      integer, intent(inout), optional :: read_n
!
   end subroutine read_excitation_vector_file_ccs
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
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
      logical,  intent(out) :: found
      integer,          intent(in) :: n
      character(len=*), intent(in) :: side
!
   end subroutine check_and_get_restart_vector_ccs
!
!
   module subroutine get_restart_vector_ccs(wf, file_, vector, energy)
!!
!!    Get restart vector
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Gets start vector and energy from file
!!
!!    Only a wrapper in CCS but overwritten for doubles
!!    to handle the restart from a pure singles vector.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
!
   end subroutine get_restart_vector_ccs
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
