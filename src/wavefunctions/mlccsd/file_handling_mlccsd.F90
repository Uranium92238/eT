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
submodule (mlccsd_class) file_handling_mlccsd
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
   module subroutine save_amplitudes_mlccsd(wf)
!!
!!    Save amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      call wf%t_file%open_('write', 'rewind')
!
      call wf%t_file%write_(wf%t1, wf%n_t1)
      call wf%t_file%write_(wf%t2, wf%n_t2)
!
      call wf%t_file%close_()
!
   end subroutine save_amplitudes_mlccsd
!
!
   module subroutine read_amplitudes_mlccsd(wf)
!!
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      integer :: records
!
      call wf%t_file%open_('read', 'rewind')
!
!     Check if doubles exist then read
!     Sinlges in first record, doubles in second
!
      records = wf%t_file%number_of_records()
!
      call wf%t_file%rewind_()
!
      call wf%t_file%read_(wf%t1, wf%n_t1)
!
      if (records .gt. 1) then
!
         call wf%t_file%read_(wf%t2, wf%n_t2)
!
      else
!
         call zero_array(wf%t2, wf%n_t2)
!
      end if
!
      call wf%t_file%close_()
!
   end subroutine read_amplitudes_mlccsd
!
!
end submodule file_handling_mlccsd
