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
module io_utilities
!
!!
!!    IO utilities module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!
   use kinds
   use direct_file_class, only : direct_file
   use output_file_class, only : output
   use batching_index_class, only : batching_index
!
   interface single_record_reader
!
      procedure   read_1_array_single_record,         &
                  read_1_array_single_record_batch,   &
                  read_2_arrays_single_record_batch,  &
                  read_3_arrays_single_record_batch
!
   end interface
!
!
   interface compound_record_reader
!
      procedure   read_1_array_compound_record_2batches,    &
                  read_2_arrays_compound_record_2batches,   &
                  read_3_arrays_compound_record_2batches,   &
                  read_4_arrays_compound_record_2batches,   &
                  read_1_array_compound_record_1batch
!
   end interface
!
!
   interface single_record_writer
!
      procedure   write_array_single_record_batch
      procedure   write_array_single_record
!
   end interface single_record_writer
!
!
   interface compound_record_writer
!
      procedure   write_array_compound_record_2batches, &
                  write_array_compound_record_1batch
!
   end interface compound_record_writer
!
!
contains
!
!
   subroutine long_string_print(string,format_string,colons,fformat_string,lformat_string,line_length)
!!
!!    Long string print
!!    Written by Rolf H. Myhre, Nov 2018
!!
!!    Prints a reasonbly formatted long string over several lines
!!    format_string: optional format string
!!    fformat_string: optional format string for first printed line, will be used for single line if present
!!    lformat_string: optional format string for last printed line
!!    line_length: optional argument for length printed lines, routine adds extra length to not split up words
!!
!!
      implicit none
!
      character(len=*), intent(in) :: string
      character(len=*), intent(in), optional :: format_string,fformat_string,lformat_string
      integer, intent(in), optional :: line_length
      logical, intent(in), optional :: colons
!
      character(90)  :: temp
      integer        :: l, l_left, lines, l_length
      integer        :: i,j, padd, printed 
      character(20)  :: fs,fstring,ffstring,lfstring
      logical        :: col 
!
!     Default line length
      l_length = 70
      if(present(line_length)) then
         l_length = line_length
      endif
!
!     Figure out the formatting
      fstring = '(t3,a)'
      if(present(format_string)) then
         fstring = format_string
      endif
!
      ffstring = fstring
      if(present(fformat_string)) then
         ffstring = fformat_string
      endif
!
      lfstring = fstring
      if(present(lformat_string)) then
         lfstring = lformat_string
      endif
!
!     First line format
      fs = ffstring
!
!     Fancy colons for banner headings
      col = .false.
      if(present(colons)) then
         col = colons
      endif
      if(col) then
         l_length = l_length - 3
      endif
!
      l = len_trim(string)      
      l_left = l
      lines = l/l_length + 1
      printed = 1
!
      do i = 1,lines
!
         if(i .ne. lines) then
!
            if(i .ne. 1) then
               fs = fstring
            endif
!
!           Add some extra padding to not split words
            do j = 1,18
               padd = j
               if(string(printed+l_length+j:printed+l_length+j) == ' ') then
                  exit
               endif
            enddo
!
!           Copy string to be printed. Add hyphen if word must be split
            if(padd == 18) then
               temp(1:l_length+padd+1) = string(printed:printed+l_length+padd)
               temp(l_length+padd+1:l_length+padd+1) = '-'
               printed = printed + l_length + padd
            else
               temp(1:l_length+padd+1) = string(printed:printed+l_length+padd+1)
               printed = printed + l_length + padd + 1
            endif
!
!           Print
            if(col) then
               write(output%unit, fs)  ':: '//temp(1:l_length+padd+1)
            else
               write(output%unit, fs)  temp(1:l_length+padd+1)
            endif
!
         else
!           Print the remaining string
            fs = lfstring
            if(col) then
               write(output%unit, fs)  ':: '//string(printed:l)
            else
               write(output%unit, fs)  string(printed:l)
            endif
            printed = l
         endif
!
      enddo
!
      flush(output%unit)
!
   end subroutine long_string_print
!
!
   subroutine read_1_array_single_record(dim_z, file_1, g_pqrz)
!!
!!    Read the direct access file "file_1" with record z into g_pqrz
!!    NB: It is assumed that z is sorted at the end
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      integer, intent(in) :: dim_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
!
      type(direct_file), intent(in) :: file_1
!
      integer :: ioerror
      integer :: z
!
      character(len=100) :: iom
!
      do z = 1, dim_z
!
         read(file_1%unit, rec=z, iostat=ioerror, iomsg=iom) g_pqrz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read file: ', trim(file_1%file_name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file:file')
         endif
!
      enddo
!
   end subroutine read_1_array_single_record
!
!
   subroutine read_1_array_single_record_batch(batch_z, file_1, g_pqrz)
!!
!!    Read parts of the direct access files "file" into g_pqrz for the current batch z
!!    NB: It is assumed that the batching index is sorted at the end
!!        and that the record is equal to the batching index
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
!
      type(direct_file), intent(in) :: file_1
!
      integer :: z, z_abs
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_1%reader(g_pqrz(:,:,:,z),z_abs)
!
      enddo
!
   end subroutine read_1_array_single_record_batch
!
!
   subroutine read_2_arrays_single_record_batch(batch_z, file_1, g_pqrz, file_2, g_stuz)
!!
!!    Read parts of the direct access files "file_1/2" 
!!    into g_pqrz/g_stuz for the current batch z
!!    NB: It is assumed that the batching index is sorted at the end
!!        and that the record is equal to the batching index
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_stuz
!
      type(direct_file), intent(in) :: file_1
      type(direct_file), intent(in) :: file_2
!
      integer :: z, z_abs
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_1%reader(g_pqrz(:,:,:,z),z_abs)
!
      enddo
!
!     Second file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_2%reader(g_stuz(:,:,:,z),z_abs)
!
      enddo
!
   end subroutine read_2_arrays_single_record_batch
!
!
   subroutine read_3_arrays_single_record_batch(batch_z, file_1, g_pqrz, file_2, g_stuz, file_3, g_vwxz)
!!
!!    Read parts of the direct access files "file_1/2" 
!!    into g_pqrz/g_stuz for the current batch z
!!    NB: It is assumed that the batching index is sorted at the end
!!        and that the record is equal to the batching index
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_stuz
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_vwxz
!
      type(direct_file), intent(in) :: file_1
      type(direct_file), intent(in) :: file_2
      type(direct_file), intent(in) :: file_3
!
      integer :: z, z_abs
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_1%reader(g_pqrz(:,:,:,z),z_abs)
!
      enddo
!
!     Second file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_2%reader(g_stuz(:,:,:,z),z_abs)
!
      enddo
!
!     Third file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_3%reader(g_vwxz(:,:,:,z),z_abs)
!
      enddo
!
   end subroutine read_3_arrays_single_record_batch
!
!
   subroutine read_1_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
!!
!!    Read parts of the direct access file "file_1" with records of zy into g_pqzy
!!    Read in batches of z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
!
      type(direct_file), intent(in) :: file_1
!
      integer :: record
      integer :: z, z_abs, y, y_abs
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_1%reader(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_1_array_compound_record_2batches
!
!
   subroutine read_1_array_compound_record_1batch(dim_z, batch_y, file_1, g_pqzy, switch)
!!
!!    Read parts of the direct access file "file_1" with records of zy into g_pqzy
!!    Reads z in full dimension y in batches
!!    NB: It is assumed that the indices zy are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Based on read_1_array_compound_record_2batches written by Alexander Paul and Rolf H. Myhre
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      integer, intent(in) :: dim_z
!
      type(batching_index), intent(inout) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
!
      type(direct_file), intent(inout) :: file_1
!
      logical, intent(in), optional :: switch
      logical :: switched
!
      type(batching_index) :: batch_z
!
!     Can't overload a function based on ordering alone, 
!     so optional keyword to revert order of y and z
!
      switched = .false.
      if(present(switch)) then
         switched = switch
      endif
!
!     Fake a batching_index with full dimensions and call 2batches_reader
!
      call batch_z%init(dim_z)
!
      batch_z%first = 1
      batch_z%last = dim_z
      batch_z%length = dim_z
      batch_z%max_length = dim_z
      batch_z%num_batches = 1
!
      if(.not. switched) then
         call read_1_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
      else
         call read_1_array_compound_record_2batches(batch_y, batch_z, file_1, g_pqzy)
      endif
!
   end subroutine read_1_array_compound_record_1batch
!
!
   subroutine read_2_arrays_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy, &
                                                   file_2, g_rszy)
!!
!!    Read parts of the direct access files "file_1/2" with records of zy into g_pqzy/g_rszy
!!    Reads in batches of z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_rszy
!
      type(direct_file), intent(in) :: file_1
      type(direct_file), intent(in) :: file_2
!
      integer :: record
      integer :: z, z_abs, y, y_abs
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_1%reader(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
!     Second file
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_2%reader(g_rszy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_2_arrays_compound_record_2batches
!
!
   subroutine read_3_arrays_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy, &
                                                   file_2, g_rszy, file_3, g_tuzy)
!!
!!    Read parts of the direct access files "file_1/2/3" with records of zy 
!!    into g_pqzy/g_rszy/g_tuzy
!!    Reads in batches of z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_rszy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_tuzy
!
      type(direct_file), intent(in) :: file_1
      type(direct_file), intent(in) :: file_2
      type(direct_file), intent(in) :: file_3
!
      integer :: record
      integer :: z, z_abs, y, y_abs
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_1%reader(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
!     Second file
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_2%reader(g_rszy(:,:,z,y), record)
!
         enddo
      enddo
!
!     Third file
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_3%reader(g_tuzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_3_arrays_compound_record_2batches
!
!
   subroutine read_4_arrays_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy, file_2,  &
                                                   g_rszy, file_3, g_tuzy, file_4, g_vwzy)
!!
!!    Read parts of the direct access files "file_1/2/3/4" with records of zy 
!!    into g_pqzy/g_rszy/g_tuzy/g_vwzy
!!    Reads in batches of z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_rszy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_tuzy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_vwzy
!
      type(direct_file), intent(in) :: file_1
      type(direct_file), intent(in) :: file_2
      type(direct_file), intent(in) :: file_3
      type(direct_file), intent(in) :: file_4
!
      integer :: record
      integer :: z, z_abs, y, y_abs
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_1%reader(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
!     Second file
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_2%reader(g_rszy(:,:,z,y), record)
!
         enddo
      enddo
!
!     Third file
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_3%reader(g_tuzy(:,:,z,y), record)
!
         enddo
      enddo
!
!     Fourth file
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = batch_z%index_dimension*(y_abs - 1) + z_abs
!
            call file_4%reader(g_vwzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_4_arrays_compound_record_2batches
!
!
   module subroutine write_array_single_record(dim_z, file_1, g_pqrz)
!!
!!    Writes an array to a direct access file "file_1" with records z
!!    The last index of the array (z) is the record number
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      integer, intent(in) :: dim_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
!
      type(direct_file), intent(in) :: file_1
!
      integer :: z
!
      do z = 1, dim_z
!
         call file_1%writer(g_pqrz(:,:,:,z), z)
!
      enddo
!
   end subroutine write_array_single_record
!
!
   module subroutine write_array_single_record_batch(batch_z, file_1, g_pqrz)
!!
!!    Writes an array to a direct access file "file_1" with records z
!!    The last index of the array (z) is batched over and is also the record number
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
!
      type(direct_file), intent(in) :: file_1
!
      integer :: z, z_abs
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_1%writer(g_pqrz(:,:,:,z), z_abs)
!
      enddo
!
   end subroutine write_array_single_record_batch
!
!
   module subroutine write_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
!!
!!    Writes an array to a direct access file "file_1" with records zy
!!
!!    NB: It is assumed that the array is split in batches of z and y 
!!        which are sorted in y,z order at the end of the array
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
!
      type(direct_file), intent(in) :: file_1
!
      integer :: record
      integer :: z, z_abs, y, y_abs
!
      do y = 1, batch_y%length
!
         y_abs = batch_y%first + y - 1
!
         do z = 1, batch_z%length
!
            z_abs = batch_z%first + z - 1
!
            record = (y_abs - 1) * batch_z%index_dimension + z_abs
!
            call file_1%writer(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine write_array_compound_record_2batches
!
!
   module subroutine write_array_compound_record_1batch(dim_z, batch_y, file_1, g_pqzy)
!!
!!    Writes an array to a direct access file "file_1" with records zy
!!
!!    NB: It is assumed that z is of full dimension and y is batched over
!!        z and y are sorted in z,y order at the end of the array
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      integer, intent(in) :: dim_z
!
      type(batching_index), intent(inout) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(inout) :: g_pqzy
!
      type(direct_file), intent(inout) :: file_1
!
      type(batching_index) :: batch_z
!
!     Fake a batching_index with full dimensions and call 2batches_writer
!
      call batch_z%init(dim_z)
!
      batch_z%first = 1
      batch_z%last = dim_z
      batch_z%length = dim_z
      batch_z%max_length = dim_z
      batch_z%num_batches = 1
!
      call write_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
!
   end subroutine write_array_compound_record_1batch
!
!
end module
