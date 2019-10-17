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
submodule (ccsd_class) file_handling_ccsd
!
!!
!!    File handling submodule (CCSD)
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
   module subroutine initialize_files_ccsd(wf)
!!
!!    Initialize files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Initializes the wavefucntion files for wavefunction parameters.
!!
      class(ccsd) :: wf 
!
      call wf%initialize_wavefunction_files()
      call wf%initialize_cc_files()
      call wf%initialize_singles_files()
      call wf%initialize_doubles_files()
!
   end subroutine initialize_files_ccsd
!
!
   module subroutine initialize_doubles_files_ccsd(wf)
!!
!!    Initialize doubles files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      class(ccsd) :: wf 
!
      wf%t2_file = sequential_file('t2')
      wf%t2bar_file = sequential_file('t2bar') 
!
      wf%r2_file = sequential_file('r2')
      wf%l2_file = sequential_file('l2') 
!
   end subroutine initialize_doubles_files_ccsd
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
      call wf%t1_file%open_('write', 'rewind')
      call wf%t2_file%open_('write', 'rewind')
!
      call wf%t1_file%write_(wf%t1, wf%n_t1)
      call wf%t2_file%write_(wf%t2, wf%n_t2)
!
      call wf%t1_file%close_()
      call wf%t2_file%close_()
!
   end subroutine save_amplitudes_ccsd
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
      call wf%t1bar_file%open_('write', 'rewind')
      call wf%t2bar_file%open_('write', 'rewind')
!
      call wf%t1bar_file%write_(wf%t1bar, wf%n_t1)
      call wf%t2bar_file%write_(wf%t2bar, wf%n_t2)
!
      call wf%t1bar_file%close_()
      call wf%t2bar_file%close_()
!
   end subroutine save_multipliers_ccsd
!
!
   module subroutine read_amplitudes_ccsd(wf)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      call wf%t1_file%open_('read', 'rewind')
      call wf%t2_file%open_('read', 'rewind')
!
      call wf%t1_file%read_(wf%t1, wf%n_t1)
      call wf%t2_file%read_(wf%t2, wf%n_t2)
!
      call wf%t1_file%close_()
      call wf%t2_file%close_()
!
   end subroutine read_amplitudes_ccsd
!
!
   module subroutine read_multipliers_ccsd(wf)
!!
!!    Read multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      call wf%t1bar_file%open_('read', 'rewind')
      call wf%t2bar_file%open_('read', 'rewind')
!
      call wf%t1bar_file%read_(wf%t1bar, wf%n_t1)
      call wf%t2bar_file%read_(wf%t2bar, wf%n_t2)
!
      call wf%t1bar_file%close_()
      call wf%t2bar_file%close_()
!
   end subroutine read_multipliers_ccsd
!
!
end submodule file_handling_ccsd 
