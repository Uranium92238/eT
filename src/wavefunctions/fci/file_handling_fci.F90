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
submodule (fci_class) file_handling_fci
!
!!
!! File handling submodule
!!
!! Gathers routines that save wavefunction parameters to file,
!! and reads them from file, plus other routines related to the
!! handling of the files that belong to the wavefunction.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_fci_state_files_fci(wf)
!!
!!    Initialize FCI state files
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      character(len=21) :: file_name
      integer :: state
!
      allocate(wf%fci_files(wf%n_states))
!
      do state = 1, wf%n_states
!
         write(file_name,'(a,i3.3)') 'fci_eigenvectors_', state
         wf%fci_files(state) = stream_file(trim(file_name))
!
      end do
!
   end subroutine initialize_fci_state_files_fci
!
!
   module subroutine save_fci_state_fci(wf, vector, energy, state)
!!
!!    Save fci state
!!    Written by Enrico Ronca, 2022
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      real(dp), dimension(wf%n_determinants), intent(in) :: vector
!
      real(dp), intent(in) :: energy
!
      integer, intent(in) :: state
!
      integer(i64) :: vector_length
!
      vector_length = int(wf%n_determinants, kind=i64)
!
      call wf%fci_files(state)%open_('rewind')
!
      call wf%fci_files(state)%write_(energy)
      call wf%fci_files(state)%write_(vector_length)
      call wf%fci_files(state)%write_(vector, wf%n_determinants)
!
      call wf%fci_files(state)%close_
!
   end subroutine save_fci_state_fci
!
!
   module subroutine read_fci_state_fci(wf, vector, state)
!!
!!    Read FCI state
!!    Written By Enrico Ronca, 2022
!!
      implicit none
!
      class(fci), intent(inout) :: wf
      real(dp), dimension(wf%n_determinants), intent(out) :: vector
!
      integer, intent(in) :: state
      integer(i64) :: vector_length
!
      call wf%fci_files(state)%open_('rewind')
      call wf%fci_files(state)%read_(vector_length, dp + 1)
!
      if (int(vector_length) /= wf%n_determinants) &
            call output%error_msg('Wrong number of parameters in (a0): (i0)', &
                                  chars=[wf%fci_files(state)%get_name()], ints=[int(vector_length)])
!
      call wf%fci_files(state)%read_(vector, wf%n_determinants)
!
      call wf%fci_files(state)%close_
!
   end subroutine read_fci_state_fci
!
!
end submodule file_handling_fci
