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
!!    Read amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      call wf%save_singles_vector(wf%t_file, wf%t1)
      call wf%save_doubles_vector(wf%t_file, wf%t2)
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
!
      integer, intent(out), optional :: read_n
!
      integer :: n
!
      n = 0
!
      call wf%read_singles_vector(wf%t_file, wf%t1, n)
!
      call wf%read_doubles_vector(wf%t_file, wf%t2, n)
!
      if (present(read_n)) read_n = n
!
   end subroutine read_amplitudes_mlccsd
!
!
   module subroutine read_doubles_vector_mlccsd(wf, file_, vector, read_n)
!!
!!    Read doubles vector
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    File format:
!!    n_t1, t1, n_t2, t2
!!
!!    read_n: optionally returns the number of amplitudes read.
!!
      implicit none 
!
      class(mlccsd), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_t2), intent(out) :: vector
!
      integer, intent(inout), optional :: read_n
!
      integer :: n_t2, iostat
!
      call file_%open_('read')
!
      call file_%read_(n_t2, int_size + wf%n_t1*dp + 1, iostat)
!
      if (.not. is_iostat_end(iostat)) then
!
         if (n_t2 .ne. wf%n_t2) then
            call output%error_msg('Wrong number of doubles amplitudes in (a0)', &
                                 chars=[file_%get_name()])
         end if
!
         call file_%read_(vector, wf%n_t2)
!
         if (present(read_n)) read_n = read_n + n_t2
!
      end if
!
      call file_%close_()
!
   end subroutine read_doubles_vector_mlccsd
!
!
   module subroutine save_doubles_vector_mlccsd(wf, file_, vector)
!!
!!    Save doubles vector
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    File format:
!!    n_t1, t1, n_t2, t2
!!
      implicit none 
!
      class(mlccsd), intent(inout) :: wf 
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_t2), intent(in) :: vector 
!
      call file_%open_('write')
!
      call file_%write_(wf%n_t2, int_size + wf%n_t1*dp + 1)
      call file_%write_(vector, wf%n_t2)
!
      call file_%close_()
!
   end subroutine save_doubles_vector_mlccsd
!
!
end submodule file_handling_mlccsd
