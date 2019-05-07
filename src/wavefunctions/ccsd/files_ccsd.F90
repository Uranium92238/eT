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
submodule (ccsd_class) files_ccsd
!
!!
!!    Files submodule (CCSD)
!!    Set up by Eirik F. Kjønstad, May 2019
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
      call wf%t2_file%init('t2', 'sequential', 'unformatted')
      call wf%t2bar_file%init('t2bar', 'sequential', 'unformatted')
      call wf%l2_file%init('l2', 'sequential', 'unformatted')
      call wf%r2_file%init('r2', 'sequential', 'unformatted')
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
      call disk%open_file(wf%t1_file, 'write')
      call disk%open_file(wf%t2_file, 'write')
!
      rewind(wf%t1_file%unit)
      rewind(wf%t2_file%unit)
!
      write(wf%t1_file%unit) wf%t1  
      write(wf%t2_file%unit) wf%t2
!
      call disk%close_file(wf%t1_file) 
      call disk%close_file(wf%t2_file) 
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
      call disk%open_file(wf%t1bar_file, 'write')
      call disk%open_file(wf%t2bar_file, 'write')
!
      rewind(wf%t1bar_file%unit)
      rewind(wf%t2bar_file%unit)
!
      write(wf%t1bar_file%unit) wf%t1bar  
      write(wf%t2bar_file%unit) wf%t2bar
!
      call disk%close_file(wf%t1bar_file) 
      call disk%close_file(wf%t2bar_file) 
!
   end subroutine save_multipliers_ccsd
!
!
   module subroutine save_excited_state_ccsd(wf, X, n, side)
!!
!!    Save excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Saves an excited state to disk. Since the solvers 
!!    keep these vectors in full length, we receive a vector 
!!    in full length (n_es_amplitudes), and then distribute 
!!    the different parts of that vector to singles, doubles, etc.,
!!    files (if there are doubles, etc.).
!!
!!    NB! If n = 1, then the routine WILL REWIND the files before writing,
!!    thus DELETING every record in the file. For n >=2, we just append to
!!    the file. The purpose of this setup is that the files should be saved in 
!!    the correct order, from n = 1 to n = # states. 
!!
      implicit none 
!
      class(ccsd), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
      if (trim(side) == 'right') then 
!
         call wf%save_singles_vector(X(1 : wf%n_t1), n, wf%r1_file)
         call wf%save_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%r2_file)
!
      elseif (trim(side) == 'left') then 
!
         call wf%save_singles_vector(X(1 : wf%n_t1), n, wf%l1_file)
         call wf%save_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%l2_file)
!
      else
!
         call output%error_msg('Tried to save an excited state, but argument side not recognized: ' // side)
!
      endif
!
   end subroutine save_excited_state_ccsd
!
!
   module subroutine save_doubles_vector_ccsd(wf, X, n, file_)
!!
!!    Save doubles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!    Writes doubles vector "X" to the sequential
!!    and unformatted file "file_".
!!    
!!    NB! If n = 1, then the routine WILL REWIND the file before writing,
!!    thus DELETING every record in the file. For n >=2, we just append to
!!    the file. The purpose of this setup is that the files should be saved in 
!!    the correct order, from n = 1 to n = # states.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(file) :: file_
!
      call disk%open_file(file_, 'write', 'append')
!
      if (n .eq. 1) rewind(file_%unit)
!
      write(file_%unit) X
!
      call disk%close_file(file_, 'keep')
!
   end subroutine save_doubles_vector_ccsd
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
      call disk%open_file(wf%t1_file, 'read', 'rewind')
      call disk%open_file(wf%t2_file, 'read', 'rewind')
!
      read(wf%t1_file%unit) wf%t1  
      read(wf%t2_file%unit) wf%t2
!
      call disk%close_file(wf%t1_file) 
      call disk%close_file(wf%t2_file) 
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
      call disk%open_file(wf%t1bar_file, 'read')
      call disk%open_file(wf%t2bar_file, 'read')
!
      rewind(wf%t1bar_file%unit)
      rewind(wf%t2bar_file%unit)
!
      read(wf%t1bar_file%unit) wf%t1bar  
      read(wf%t2bar_file%unit) wf%t2bar
!
      call disk%close_file(wf%t1bar_file) 
      call disk%close_file(wf%t2bar_file) 
!
   end subroutine read_multipliers_ccsd
!
!
   module subroutine read_excited_state_ccsd(wf, X, n, side)
!!
!!    Read excited state 
!!    Written by Sarai D. Fokestad, Mar 2019 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: X
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
      if (trim(side) == 'right') then
!
         call wf%read_singles_vector(X(1 : wf%n_t1), n, wf%r1_file)
         call wf%read_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%r2_file)
!
      elseif (trim(side) == 'left') then
!
         call wf%read_singles_vector(X(1 : wf%n_t1), n, wf%l1_file)
         call wf%read_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%l2_file)
!
      endif
!
   end subroutine read_excited_state_ccsd
!
!
   module subroutine read_doubles_vector_ccsd(wf, X, n, file_)
!!
!!    Read doubles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!    Reads doubles vector "X" from the "n"'th line
!!    of the sequential and unformatted file "file_".
!!
      implicit none 
!
      class(ccsd), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(out) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(file) :: file_
!
      call disk%open_file(file_, 'read')
!
      call file_%prepare_to_read_line(n)
!
      read(file_%unit) X
!
      call disk%close_file(file_, 'keep')
!
   end subroutine read_doubles_vector_ccsd
!
!
end submodule files_ccsd 
