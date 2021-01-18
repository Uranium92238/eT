!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   module subroutine get_restart_vector_mlccsd(wf, file_, vector, energy, &
                                              restart_from, restart_to)
!!
!!    Get restart vector
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Gets start vector and energy from file and
!!    handles the basis transformations according to:
!!
!!    restart from "right" to "left"
!!    L^a_i = 2R^a_i
!!    L^ab_ij = 4R^ab_ij - 2R^ba_ij
!!
!!    restart from "left" to "right"
!!    R^a_i = 1/2 L^a_i
!!    R^ab_ij = 1/6 (2L^ab_ij + L^ba_ij)
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf 
      type(stream_file), intent(inout) :: file_
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
      character(len=*), intent(in) :: restart_from, restart_to
!
   end subroutine get_restart_vector_mlccsd
!
!
   module subroutine save_mlcc_orbitals_mlccsd(wf)
!!
!!    Save MLCC orbitals
!!    Written by Sarai D. Folkestad
!!
!!    File format:
!!
!!    1: n_ccs_o, n_ccs_v, n_cc2_o, n_cc2_v, n_ccsd_o, n_ccsd_v, orbitals 
!!    2: orbital energies
!!
      implicit none 
!
      class(mlccsd), intent(inout) :: wf 
!
   end subroutine save_mlcc_orbitals_mlccsd
!
!
   module subroutine read_mlcc_orbitals_mlccsd(wf)
!!
!!    Read MLCC orbitals
!!    Written by Sarai D. Folkestad
!!
!!    File format:
!!
!!    1: n_ccs_o, n_ccs_v, n_cc2_o, n_cc2_v, n_ccsd_o, n_ccsd_v, orbitals 
!!    2: orbital energies
!!
      implicit none 
!
      class(mlccsd), intent(inout) :: wf 
!
   end subroutine read_mlcc_orbitals_mlccsd
