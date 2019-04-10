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
   use input_file_class
   use batching_index_class
!
   interface single_batch_reader
!
      procedure   read_1_4dim_array_with_1index_record, &
                  read_2_4dim_array_with_1index_record, &
                  read_3_4dim_array_with_1index_record
!
   end interface
!
   interface double_batch_reader
!
      procedure   read_1_4dim_array_with_2index_record, &
                  read_2_4dim_array_with_2index_record, &
                  read_3_4dim_array_with_2index_record, &
                  read_4_4dim_array_with_2index_record
!
   end interface
!
!
contains
!
   function remove_preceding_blanks(line)
!!  
!!     Remove preceding blanks
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017
!!  
!!     Removes white spaces before text from line
!!  
      implicit none
!
      character(len=100) :: line
!
      character(len=100) :: remove_preceding_blanks
!
      integer :: i = 0, length = 0
!
      remove_preceding_blanks = ' '
!
      do i = 1, 100
         if (line(i:i) == ' ') then
!
            continue
!
         else
!
            length = 100 - (i - 1)
            remove_preceding_blanks(1:length) = line(i:100)
            remove_preceding_blanks(length+1:100) = ' '
            return
!
         endif
      enddo
!
   end function remove_preceding_blanks
!
!
   logical function requested_section(string)
!!
!!
!!
      implicit none
!
      character(len=*) :: string
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      requested_section = .false.
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') then
            backspace(input%unit)
            return
         endif
!
         if (trim(line) == 'end ' // string) then
!
            requested_section = .true.
            return
! 
         endif
!
      enddo
!
   end function requested_section
!
!
   subroutine move_to_section(string, n_records)
!!
!!    Move to section,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Moves cursor to section given by string, and counts the number of records in the section
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      integer :: n_records
!
      character(len=100) :: line
!
      integer :: n_start, count_start, count_rec_end, count_rec_start, i
!
      rewind(input%unit)
!
      count_rec_end = 0
      n_start = 0
!
      count_start = 0
      count_rec_start = 0
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         count_rec_end = count_rec_end + 1
!
         if (trim(line) .eq. 'geometry') then
!
            backspace(input%unit)
            call output%error_msg('could not move to requested section: '// string)
!
         endif
!
         if (trim(line) == 'end ' // string) then        
!
            rewind(input%unit)
!
            do i = 1, count_rec_end - 1
!
               read(input%unit, '(a100)') line
               line = remove_preceding_blanks(line) 
!  
               if (trim(line) == string) then
!
                  n_start =  n_start + 1
!
               endif
!
            enddo 
!
            rewind(input%unit)
!
            do i = 1, count_rec_end - 1
!
              read(input%unit, '(a100)') line
              line = remove_preceding_blanks(line)
!  
              count_rec_start = count_rec_start + 1
!  
               if (trim(line) == string) then
!
                  count_start = count_start + 1
!
                  if (count_start == n_start) then
!
                     n_records = count_rec_end - count_rec_start - 1
!
                     return
!
                  endif
!
               endif
!
            enddo    
! 
         endif
!
      enddo
!
   end subroutine move_to_section
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
   module subroutine read_1_4dim_array_with_1index_record(batch_z, file_1, g_pqrz)
!!
!!    Read parts of the direct access files "file" into g_pqrz for the current batch z
!!    NB: It is assumed that the batching index is sorted at the end
!!        and that the record is equal to the batching index
!!
!!    Based on omega_cc3_vvv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
!
      type(file), intent(in) :: file_1
!
      integer :: ioerror
      integer :: z, z_abs
!
      character(len=100) :: iom
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         read(file_1%unit, rec=z_abs, iostat=ioerror, iomsg=iom) g_pqrz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
   end subroutine read_1_4dim_array_with_1index_record
!
!
   subroutine read_2_4dim_array_with_1index_record(batch_z, file_1, g_pqrz,  &
                                                      file_2, g_stuz)
!!
!!    Read parts of the direct access files "file_1/2" 
!!    into g_pqrz/g_stuz for the current batch z
!!    NB: It is assumed that the batching index is sorted at the end
!!        and that the record is equal to the batching index
!!
!!    Based on omega_cc3_vvv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_stuz
!
      type(file), intent(in) :: file_1
      type(file), intent(in) :: file_2
!
      integer :: ioerror
      integer :: z, z_abs
!
      character(len=100) :: iom
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         read(file_1%unit, rec=z_abs, iostat=ioerror, iomsg=iom) g_pqrz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!     Second file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         read(file_2%unit, rec=z_abs, iostat=ioerror, iomsg=iom) g_stuz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_2%name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
   end subroutine read_2_4dim_array_with_1index_record
!
!
   subroutine read_3_4dim_array_with_1index_record(batch_z, file_1, g_pqrz,  &
                                                      file_2, g_stuz, file_3, g_vwxz)
!!
!!    Read parts of the direct access files "file_1/2" 
!!    into g_pqrz/g_stuz for the current batch z
!!    NB: It is assumed that the batching index is sorted at the end
!!        and that the record is equal to the batching index
!!
!!    Based on omega_cc3_vvv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqrz
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_stuz
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_vwxz
!
      type(file), intent(in) :: file_1
      type(file), intent(in) :: file_2
      type(file), intent(in) :: file_3
!
      integer :: ioerror
      integer :: z, z_abs
!
      character(len=100) :: iom
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         read(file_1%unit, rec=z_abs, iostat=ioerror, iomsg=iom) g_pqrz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!     Second file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         read(file_2%unit, rec=z_abs, iostat=ioerror, iomsg=iom) g_stuz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_2%name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
!     Third file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         read(file_3%unit, rec=z_abs, iostat=ioerror, iomsg=iom) g_vwxz(:,:,:,z)
!
         if(ioerror .ne. 0) then
            write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_3%name)
            write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
            write(output%unit,'(t3,a)') trim(iom)
            call output%error_msg('Failed to read file')
         endif
!
      enddo
!
   end subroutine read_3_4dim_array_with_1index_record
!
!
   subroutine read_1_4dim_array_with_2index_record(batch_z, batch_y, file_1, g_pqzy)
!!
!!    Read parts of the direct access file "file_1" into g_pqzy for the current batches in z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Based on omega_cc3_ov_vv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
!
      type(file), intent(in) :: file_1
!
      integer :: ioerror, record
      integer :: z, z_abs, y, y_abs
!
      character(len=100) :: iom
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
            read(file_1%unit, rec=record, iostat=ioerror, iomsg=iom) g_pqzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
!
         enddo
      enddo
!
   end subroutine read_1_4dim_array_with_2index_record
!
!
   subroutine read_2_4dim_array_with_2index_record(batch_z, batch_y, file_1, g_pqzy,  &
                                                      file_2, g_rszy)
!!
!!    Read parts of the direct access files "file_1/2" into g_pqzy/g_rszy 
!!    for the current batches in z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Based on omega_cc3_ov_vv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_rszy
!
      type(file), intent(in) :: file_1
      type(file), intent(in) :: file_2
!
      integer :: ioerror, record
      integer :: z, z_abs, y, y_abs
!
      character(len=100) :: iom
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
            read(file_1%unit, rec=record, iostat=ioerror, iomsg=iom) g_pqzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
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
            read(file_2%unit, rec=record, iostat=ioerror, iomsg=iom) g_rszy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_2%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
!
         enddo
      enddo
!
   end subroutine read_2_4dim_array_with_2index_record
!
!
   subroutine read_3_4dim_array_with_2index_record(batch_z, batch_y, file_1, g_pqzy,  &
                                                      file_2, g_rszy, file_3, g_tuzy)
!!
!!    Read parts of the direct access files "file_1/2/3" into g_pqzy/g_rszy/g_tuzy
!!    for the current batches in z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Based on omega_cc3_ov_vv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
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
      type(file), intent(in) :: file_1
      type(file), intent(in) :: file_2
      type(file), intent(in) :: file_3
!
      integer :: ioerror, record
      integer :: z, z_abs, y, y_abs
!
      character(len=100) :: iom
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
            read(file_1%unit, rec=record, iostat=ioerror, iomsg=iom) g_pqzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
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
            read(file_2%unit, rec=record, iostat=ioerror, iomsg=iom) g_rszy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_2%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
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
            read(file_3%unit, rec=record, iostat=ioerror, iomsg=iom) g_tuzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_3%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
!
         enddo
      enddo
!
   end subroutine read_3_4dim_array_with_2index_record
!
!
   subroutine read_4_4dim_array_with_2index_record(batch_z, batch_y, file_1, g_pqzy,  &
                                                      file_2, g_rszy, file_3, g_tuzy,    &
                                                      file_4, g_vwzy)
!!
!!    Read parts of the direct access files "file_1/2/3/4" into g_pqzy/g_rszy/g_tuzy/g_vwzy 
!!    for the current batches in z and y
!!    NB: It is assumed that the batching indices are sorted at the end of the array
!!        in z,y order and that the record is equal to the compound index zy
!!
!!    Based on omega_cc3_ov_vv_reader_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, April 2019
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
      type(file), intent(in) :: file_1
      type(file), intent(in) :: file_2
      type(file), intent(in) :: file_3
      type(file), intent(in) :: file_4
!
      integer :: ioerror, record
      integer :: z, z_abs, y, y_abs
!
      character(len=100) :: iom
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
            read(file_1%unit, rec=record, iostat=ioerror, iomsg=iom) g_pqzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_1%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
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
            read(file_2%unit, rec=record, iostat=ioerror, iomsg=iom) g_rszy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_2%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
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
            read(file_3%unit, rec=record, iostat=ioerror, iomsg=iom) g_tuzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_3%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
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
            read(file_4%unit, rec=record, iostat=ioerror, iomsg=iom) g_vwzy(:,:,z,y)
!
            if(ioerror .ne. 0) then
               write(output%unit,'(t3,a,a)') 'Failed to read ', trim(file_4%name)
               write(output%unit,'(t3,a,i14)') 'Error code: ', ioerror
               write(output%unit,'(t3,a)') trim(iom)
               call output%error_msg('Failed to read file')
            endif
!
         enddo
      enddo
!
   end subroutine read_4_4dim_array_with_2index_record
!
!
end module
