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
submodule (ccsd_class) file_handling_ccsd
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
   module subroutine save_amplitudes_ccsd(wf)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      call wf%t_file%open_('write')
!
      call wf%t_file%write_(wf%energy)
      call wf%save_singles_vector(wf%t_file, wf%t1)
      call wf%save_doubles_vector(wf%t_file, wf%t2, wf%n_t2)
!
      call wf%t_file%close_
!
   end subroutine save_amplitudes_ccsd
!
!
   module subroutine read_amplitudes_ccsd(wf, read_n)
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
      class(ccsd), intent(inout) :: wf
!
      integer, intent(out), optional :: read_n
!
      integer :: n
!
      n = 0
!
      call wf%t_file%open_('read','rewind')
!
      call wf%read_singles_vector(wf%t_file, wf%t1, n)
      call wf%read_doubles_vector(wf%t_file, wf%t2, wf%n_t2, n)
!
      call wf%t_file%close_
!
      if (present(read_n)) read_n = n
!
   end subroutine read_amplitudes_ccsd
!
!
   module subroutine save_multipliers_ccsd(wf)
!!
!!    Save multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf 
!
      call wf%tbar_file%open_('write')
!
      call wf%tbar_file%write_(wf%energy)
      call wf%save_singles_vector(wf%tbar_file, wf%t1bar)
      call wf%save_doubles_vector(wf%tbar_file, wf%t2bar, wf%n_t2)
!
      call wf%tbar_file%close_
!
   end subroutine save_multipliers_ccsd
!
!
   module subroutine read_multipliers_ccsd(wf, read_n)
!!
!!    Read multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!    Adapted to return the number of read multipliers if requested 
!!    by Alexander C. Paul, Oct 2020
!!
!!    read_n: optionally returns the number of mutliplers read. 
!!            This is especially useful e.g. in CCSD to provide a start guess 
!!            for the doubles if only singles were found on file.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      integer, intent(out), optional :: read_n
!
      integer :: n
!
      n = 0
!
      call wf%tbar_file%open_('read','rewind')
!
      call wf%read_singles_vector(wf%tbar_file, wf%t1bar, n)
      call wf%read_doubles_vector(wf%tbar_file, wf%t2bar, wf%n_t2, n)
!
      call wf%tbar_file%close_
!
      if (present(read_n)) read_n = n
!
   end subroutine read_multipliers_ccsd
!
!
end submodule file_handling_ccsd 
