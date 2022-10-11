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
module amplitude_file_storer_class
!
!!
!! Amplitude file storer class
!! Written by Alexander C. Paul, Mar 2022
!!
!! Tool that handles storage of amplitudes, e.g. coupled cluster
!! ground or excited state amplitudes on file.
!! The tool is set up with a block structure, such that the amplitudes
!! can be stored in blocks.
!!
!! For the t-amplitudes of CCSD the file structure is:
!! energy, n_t1, t1, n_t2, t2
!!
!
   use parameters
   use stream_file_class, only: stream_file
   use global_out, only: output
!
   implicit none
!
   type :: amplitude_file_storer
!
      integer, private :: n_blocks, n_amplitudes
      integer(i64), dimension(:), allocatable, private :: n_amplitudes_per_block
!
      type(stream_file), private :: file_
!
   contains
!
      procedure, public :: save_ => save_amplitude_file_storer
      procedure, public :: read_ => read_amplitude_file_storer
      procedure, public :: get_filename => get_filename_amplitude_file_storer
      procedure, public :: file_exists => file_exists_amplitude_file_storer
      procedure, public :: delete_file => delete_file_amplitude_file_storer
!
   end type amplitude_file_storer
!
!
   interface amplitude_file_storer
!
      procedure :: new_amplitude_file_storer
!
   end interface amplitude_file_storer
!
!
contains
!
!
   function new_amplitude_file_storer(file_name, n_amplitudes_per_block) result(this)
!!
!!    New
!!    Written by Alexander C. Paul, Mar 2022
!!
      implicit none
!
      type(amplitude_file_storer) :: this
!
      character(len=*), intent(in) :: file_name
!
      integer, dimension(:), intent(in) :: n_amplitudes_per_block
!
      integer :: i
!
      this%file_ = stream_file(file_name)
!
      this%n_blocks = size(n_amplitudes_per_block)
      allocate(this%n_amplitudes_per_block(this%n_blocks))
      this%n_amplitudes_per_block = int(n_amplitudes_per_block, kind=i64)
!
      this%n_amplitudes = 0
!
      do i = 1, this%n_blocks
         this%n_amplitudes = this%n_amplitudes + int(this%n_amplitudes_per_block(i))
      end do
!
   end function new_amplitude_file_storer
!
!
   subroutine save_amplitude_file_storer(this, vector, energy)
!!
!!    Save
!!    Written by Alexander C. Paul, Mar 2022
!!
      implicit none
!
      class(amplitude_file_storer), intent(inout) :: this
!
      real(dp), dimension(this%n_amplitudes), intent(in) :: vector
      real(dp), intent(in), optional :: energy
!
      real(dp) :: local_energy
!
      integer(i64) :: n, elements_written
      integer :: i
!
      local_energy = zero
      if (present(energy)) local_energy = energy
!
      call this%file_%open_('rewind')
!
      call this%file_%write_(local_energy)
!
      elements_written = 0
!
      do i = 1, this%n_blocks
!
         n = this%n_amplitudes_per_block(i)
!
         call this%file_%write_(n)
         call this%file_%write_(vector(elements_written+1:), int(n))
!
         elements_written = elements_written + n
!
      end do
!
      call this%file_%close_()
!
   end subroutine save_amplitude_file_storer
!
!
   subroutine read_amplitude_file_storer(this, vector, energy, read_n)
!!
!!    Read
!!    Written by Alexander C. Paul, Mar 2022
!!
      implicit none
!
      class(amplitude_file_storer), intent(inout) :: this
!
      real(dp), dimension(this%n_amplitudes), intent(out) :: vector
      real(dp), intent(out) :: energy
      integer, intent(out), optional :: read_n
!
      integer(i64) :: n, elements_read, n_from_file
      integer :: i, iostat
!
      call this%file_%open_('rewind')
!
      call this%file_%read_(energy)
!
      elements_read = 0
!
      do i = 1, this%n_blocks
!
         n = int(this%n_amplitudes_per_block(i),kind=i64)
!
         iostat = 0
         call this%file_%read_(n_from_file, status_=iostat)
!
         if (is_iostat_end(iostat)) exit
!
         if (n_from_file .ne. n) then
            call output%error_msg('Wrong number of singles amplitudes in (a0). &
                                  &Expected (i0), found (i0)', &
                                 chars=[this%file_%get_name()], &
                                 ints=[int(n), int(n_from_file)])
         end if
!
         call this%file_%read_(vector(elements_read+1:), int(n))
!
         elements_read = elements_read + n
!
      end do
!
      if (present(read_n)) read_n = int(elements_read)
!
      call this%file_%close_()
!
   end subroutine read_amplitude_file_storer
!
!
   function get_filename_amplitude_file_storer(this) result(filename)
!!
!!    Get filename
!!    Written by Alexander C. Paul, Mar 2022
!!
      implicit none
!
      class(amplitude_file_storer), intent(in) :: this
!
      character(len=:), allocatable :: filename
!
      filename = this%file_%get_name()
!
   end function get_filename_amplitude_file_storer
!
!
   function file_exists_amplitude_file_storer(this) result(exists)
!!
!!    File exists
!!    Written by Alexander C. Paul, Mar 2022
!!
      implicit none
!
      class(amplitude_file_storer), intent(in) :: this
!
      logical :: exists
!
      exists = this%file_%exists()
!
   end function file_exists_amplitude_file_storer
!
!
   subroutine delete_file_amplitude_file_storer(this)
!!
!!    Delete file
!!    Written by Alexander C. Paul, Mar 2022
!!
      implicit none
!
      class(amplitude_file_storer), intent(inout) :: this
!
      call this%file_%delete_()
!
   end subroutine delete_file_amplitude_file_storer
!
!
end module amplitude_file_storer_class
