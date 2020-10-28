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
   module subroutine save_amplitudes_mlccsd(wf)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
   end subroutine save_amplitudes_mlccsd
!
!
   module subroutine read_amplitudes_mlccsd(wf, read_n)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!    Adapted to return the number of read amplitdues if requested 
!!    by Alexander C. Paul, Oct 2020
!!
!!    read_n: returns the number of amplitudes read. 
!!            This is especially useful e.g. in CCSD to provide a start guess 
!!            for the doubles if only singles were found on file.
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
      integer, intent(out), optional :: read_n
!
   end subroutine read_amplitudes_mlccsd
!
!
   module subroutine read_excitation_vector_file_mlccsd(wf, file_, vector, energy, read_n)
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
      class(mlccsd), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
      integer, intent(inout), optional :: read_n
!
   end subroutine read_excitation_vector_file_mlccsd
!
!
   module subroutine save_excitation_vector_on_file_mlccsd(wf, file_, vector, energy)
!!
!!    Save excitation vector on file 
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Writes excitation vector o file structured as follows:
!!    excitation_energy, n_t1, X1
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: vector
      real(dp), intent(in) :: energy
!
   end subroutine save_excitation_vector_on_file_mlccsd
!
!
   module subroutine get_restart_vector_mlccsd(wf, file_, vector, energy)
!!
!!    Get restart vector
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Gets start vector and energy from file
!!    We set the doubles part to 0 if only singles amplitudes
!!    were found on file, 
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
!
   end subroutine get_restart_vector_mlccsd
