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
module direct_stream_file_class
!
!!
!! Direct access stream file class module
!! Written by Rolf H. Myhre, Feb. 2020
!!
!! Stream file that acts as a Fortran direct file
!! Byte positions and read/write lengths are calculated based on record dimension
!! in new_direct_stream and records sent in to read and write.
!!
!! Because of various compiler bugs, specific wrappers like real_2 must be called
!! where the number refers to the rank of the array.
!!
!
   use kinds
   use abstract_stream_class, only : abstract_stream
   use global_out, only : output
!
   type, extends(abstract_stream) :: direct_stream_file
!
      integer, private  :: record_dim     ! Number of words per record
      integer, private  :: word_size      ! Size of a word, default is double precision
      integer, private  :: record_length  ! record_dim*word_size
!
   contains
!
!     Read routines
!
!     Real double precision read
!
      procedure, private :: read_1_real_dp_direct_stream_file
      procedure, private :: read_2_real_dp_direct_stream_file
      procedure, private :: read_3_real_dp_direct_stream_file
      procedure, private :: read_4_real_dp_direct_stream_file
!
!     Complex double precision read
!
      procedure, private :: read_1_complex_dp_direct_stream_file
      procedure, private :: read_2_complex_dp_direct_stream_file
      procedure, private :: read_3_complex_dp_direct_stream_file
      procedure, private :: read_4_complex_dp_direct_stream_file
!
!     Read generic
!
      generic, public :: read_ => read_1_real_dp_direct_stream_file,    &
                                  read_2_real_dp_direct_stream_file,    &
                                  read_3_real_dp_direct_stream_file,    &
                                  read_4_real_dp_direct_stream_file,    &
                                  read_1_complex_dp_direct_stream_file, &
                                  read_2_complex_dp_direct_stream_file, &
                                  read_3_complex_dp_direct_stream_file, &
                                  read_4_complex_dp_direct_stream_file

!
!     Specialized read routines
!
      procedure, public :: read_range &
                        => read_range_direct_stream_file
      procedure, public :: read_compound_full_batch &
                        => read_compound_full_batch_direct_stream_file
!
!     Write routines
!
!     Real double precision write
!
      procedure, private :: write_1_real_dp_direct_stream_file
      procedure, private :: write_2_real_dp_direct_stream_file
      procedure, private :: write_3_real_dp_direct_stream_file
      procedure, private :: write_4_real_dp_direct_stream_file
!
!     complex double precision write
!
      procedure, private :: write_1_complex_dp_direct_stream_file
      procedure, private :: write_2_complex_dp_direct_stream_file
      procedure, private :: write_3_complex_dp_direct_stream_file
      procedure, private :: write_4_complex_dp_direct_stream_file
!
!     Write generic
!
      generic, public :: write_ => write_1_real_dp_direct_stream_file,    &
                                   write_2_real_dp_direct_stream_file,    &
                                   write_3_real_dp_direct_stream_file,    &
                                   write_4_real_dp_direct_stream_file,    &
                                   write_1_complex_dp_direct_stream_file, &
                                   write_2_complex_dp_direct_stream_file, &
                                   write_3_complex_dp_direct_stream_file, &
                                   write_4_complex_dp_direct_stream_file
!
!     Specialized write routines
!
      procedure, public :: write_range &
                        => write_range_direct_stream_file
      procedure, public :: write_compound_batch_full &
                        => write_compound_batch_full_direct_stream_file
!
      procedure, public :: get_n_records &
                        => get_n_records_direct_stream_file
!
      procedure, private :: check_records_to_be_accessed
!
   end type direct_stream_file
!
   interface direct_stream_file
!
      procedure new_direct_stream_file
!
   end interface direct_stream_file
!
contains
!
!
   function new_direct_stream_file(name_, record_dim, word_size, status_) result(this)
!!
!!    new direct stream file
!!    Written by Rolf H. Myhre, May 2019
!!
!!    record_dim is number of words in each record
!!    word_size (optional) is the size of each word, default is double precision
!!    record length is record_dim*word_size
!!
      implicit none
!
      type(direct_stream_file) :: this
!
      character(len=*), intent(in)           :: name_
      integer, intent(in)                    :: record_dim
!
      integer, intent(in), optional          :: word_size
      character(len=*), intent(in), optional :: status_
!
      call this%initialize(name_, status_)
!
      this%word_size = dp
      if (present(word_size)) then
         if (word_size < 1) then
            call output%error_msg("Word size less than one for file (a0)", &
                                  chars=[this%get_name()])
         endif
         this%word_size = word_size
      endif
!
      if (record_dim < 1) then
         call output%error_msg("Record dimension less than one for file (a0)", &
                               chars=[this%get_name()])
      endif
!
      this%record_dim = record_dim
      this%record_length = record_dim * this%word_size
!
   end function new_direct_stream_file
!
!
   subroutine read_1_real_dp_direct_stream_file(this, array, first_rec, last_rec)
!!
!!    read 1 real dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    array :: real double precision array
!!
!!    first_rec :: first record to read
!!
!!    last_rec  :: last record to read
!!
!!    Calculates position to read in the underlying stream file based on
!!    first_rec and record_length
!!
!!    Calculates number of bytes to read based on first_rec, last_rec and this%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      real(dp), dimension((last_rec - first_rec + 1)*this%record_dim), intent(out) :: array
!
      integer :: position_, read_length, records
!
      records = last_rec - first_rec + 1
      read_length = this%record_dim * records
      position_ = (first_rec - 1) * this%record_length + 1
!
      call this%check_records_to_be_accessed(first_rec, last_rec, "read")
!
      call this%read_1_real_dp_abstract_stream(array, read_length, position_)
!
   end subroutine read_1_real_dp_direct_stream_file
!
!
!! Wrapper routines for read_1_real_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 2020
!
   subroutine read_2_real_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:), intent(out) :: array
!
      call this%read_1_real_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_2_real_dp_direct_stream_file
!
   subroutine read_3_real_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:), intent(out) :: array
!
      call this%read_1_real_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_3_real_dp_direct_stream_file
!
   subroutine read_4_real_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:,:), intent(out) :: array
!
      call this%read_1_real_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_4_real_dp_direct_stream_file
!
!
   subroutine read_1_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
!!
!!    read direct stream
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    array :: array of class(*) to dump file contents in.
!!             See below in select_type construct for accecpted types.
!!             GFortran does not currently accept c_loc of polymorphic types
!!
!!    first_rec :: The first record to read
!!
!!    last_rec  :: Optional, last record to read, default: first_rec
!!
!!    Calculates position to read in the underlying stream file based on
!!    first_rec and record_length
!!
!!    Calculates number of bytes to read based on first_rec, last_rec and this%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      complex(dp), dimension((last_rec - first_rec + 1)*this%record_dim), intent(out) :: array
!
      integer :: position_, read_length, records
!
      records = last_rec - first_rec + 1
      read_length = this%record_dim*records
      position_ = (first_rec - 1) * this%record_length + 1
!
      call this%check_records_to_be_accessed(first_rec, last_rec, "read")
!
      call this%read_1_complex_dp_abstract_stream(array, read_length, position_)
!
   end subroutine read_1_complex_dp_direct_stream_file
!
!
!! Wrapper routines for read_1_complex_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 2020
!
   subroutine read_2_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:), intent(out) :: array
!
      call this%read_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_2_complex_dp_direct_stream_file
!
   subroutine read_3_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:), intent(out) :: array
!
      call this%read_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_3_complex_dp_direct_stream_file
!
   subroutine read_4_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:,:), intent(out) :: array
!
      call this%read_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_4_complex_dp_direct_stream_file
!
!
   subroutine read_range_direct_stream_file(this, array, z_range)
!!
!!    read interval direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump file contents in.
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    z_range :: batching index
!!
!!    Wrapper for read_1 that calculates first and last record from
!!    a batching index
!!
      use range_class
!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      class(range_), intent(in) :: z_range
!
      real(dp), dimension(z_range%length) :: array
!
      call this%read_1_real_dp_direct_stream_file(array, z_range%first, z_range%get_last())
!
   end subroutine read_range_direct_stream_file
!
!
   subroutine read_compound_full_batch_direct_stream_file(this, array, dim_y, z_range)
!!
!!    read compound full batch direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump file contents in.
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    dim_y   :: integer, full dimension of first component of batching index
!!    z_range :: batching index
!!
!!    Reads a file whose records are represented by a compound index yz.
!!    Reads full y dimension and batches of z.
!!
      use range_class
!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      class(range_), intent(in) :: z_range
      integer, intent(in) :: dim_y
!
      real(dp), dimension(1) :: array
!
      call this%read_1_real_dp_direct_stream_file(array, &
                           (z_range%first-1)*dim_y + 1, &
                            z_range%get_last()*dim_y)
!
   end subroutine read_compound_full_batch_direct_stream_file
!
!
   subroutine write_1_real_dp_direct_stream_file(this, array, first_rec, last_rec)
!!
!!    write 1 real dp
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    array :: real double precision array
!!
!!    first_rec :: first record to read
!!
!!    last_rec  :: last record to read
!!
!!    Calculates position to write in the underlying stream file based on
!!    first_rec and record_length
!!
!!    Calculates number of bytes to write based on first_rec, last_rec and this%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      real(dp), dimension((last_rec - first_rec + 1)*this%record_dim), intent(in) :: array
!
      integer :: position_, write_length, records
!
      records = last_rec - first_rec + 1
      write_length = this%record_dim*records
      position_ = (first_rec - 1) * this%record_length + 1
!
      call this%check_records_to_be_accessed(first_rec, last_rec, "write")
!
      call this%write_1_real_dp_abstract_stream(array, write_length, position_)
!
   end subroutine write_1_real_dp_direct_stream_file
!
!
!! Wrapper routines for write_1_real_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 2020
!
   subroutine write_2_real_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:), intent(in) :: array
!
      call this%write_1_real_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_2_real_dp_direct_stream_file
!
   subroutine write_3_real_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:), intent(in) :: array
!
      call this%write_1_real_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_3_real_dp_direct_stream_file
!
   subroutine write_4_real_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:,:), intent(in) :: array
!
      call this%write_1_real_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_4_real_dp_direct_stream_file
!
!
   subroutine write_1_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
!!
!!    write 1 complex dp
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    array :: complex double precision array
!!
!!    first_rec :: first record to read
!!
!!    last_rec  :: last record to read
!!
!!    Calculates position to write in the underlying stream file based on
!!    first_rec and record_length
!!
!!    Calculates number of bytes to write based on first_rec, last_rec and this%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      complex(dp), dimension((last_rec - first_rec + 1)*this%record_dim), intent(in) :: array
!
      integer :: position_, write_length, records
!
      records = last_rec - first_rec + 1
      write_length = this%record_dim*records
      position_ = (first_rec - 1) * this%record_length + 1
!
      call this%check_records_to_be_accessed(first_rec, last_rec, "write")
!
      call this%write_1_complex_dp_abstract_stream(array, write_length, position_)
!
   end subroutine write_1_complex_dp_direct_stream_file
!
!
!! Wrapper routines for write_1_complex_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 2020
!
   subroutine write_2_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:), intent(in) :: array
!
      call this%write_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_2_complex_dp_direct_stream_file
!
   subroutine write_3_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:), intent(in) :: array
!
      call this%write_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_3_complex_dp_direct_stream_file
!
   subroutine write_4_complex_dp_direct_stream_file(this, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: this
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:,:), intent(in) :: array
!
      call this%write_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_4_complex_dp_direct_stream_file
!
!
   subroutine write_range_direct_stream_file(this, array, z_range)
!!
!!    write interval direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump to file
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    z_range :: batching index
!!
!!    Wrapper for write_1_real_dp that calculates first and last record from
!!    a batching index
!!
      use range_class
!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      class(range_), intent(in) :: z_range
!
      real(dp), dimension(1) :: array
!
      call this%write_1_real_dp_direct_stream_file(array, z_range%first, &
                                                              z_range%get_last())
!
   end subroutine write_range_direct_stream_file
!
!
   subroutine write_compound_batch_full_direct_stream_file(this, array, batch_y, dim_z)
!!
!!    write compound batch full direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump to file
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    batch_y :: batching index
!!    dim_z   :: integer, full dimension of second component of batching index
!!
!!    Writes a file whose records are represented by a compound index yz.
!!    Writes data in batches of y and full dimension of z.
!!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      class(batching_index), intent(in) :: batch_y
      integer, intent(in) :: dim_z
!
      real(dp),dimension(this%record_dim*batch_y%max_length, dim_z) :: array
!
      integer :: z
!
!     Check if we can do a single continuous write, else we have to loop
      if (batch_y%length .eq. batch_y%index_dimension) then
!
         call this%write_1_real_dp_direct_stream_file(array, 1, batch_y%index_dimension*dim_z)
!
      else
         do z = 1, dim_z
!
            call this%write_1_real_dp_direct_stream_file(array(:, z), &
                                  batch_y%index_dimension*(z-1) + batch_y%first, &
                                  batch_y%index_dimension*(z-1) + batch_y%get_last())
         enddo
      endif
!
   end subroutine write_compound_batch_full_direct_stream_file
!
!
   function get_n_records_direct_stream_file(this) result(n_records)
!!
!!    Get number of records
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Returns the number of records on the file.
!!
      implicit none
!
      class(direct_stream_file) :: this
!
      integer :: n_records
!
      integer :: file_size
!
      file_size = this%get_file_size()
!
      if (file_size .eq. -1) then
!
!        File size is -1 if get_file_size is not able to determine the size.
!        Assume that the number of records is zero in this case.
!
         n_records = 0
!
      else
!
         n_records = file_size/this%record_length
!
      endif
!
   end function get_n_records_direct_stream_file
!
!
   subroutine check_records_to_be_accessed(this, first_record, last_record, task)
!!
!!    Check records to be accessed
!!    Written by Alexander C. Paul, 2022
!!
!!    Sanity checks for first and last record when reading or writing
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: this
!
      integer, intent(in) :: first_record, last_record
      character(len=*), intent(in) :: task
!
      if (first_record < 1) then
         call output%error_msg('First record (i0) less than 1 for write in file ' // &
                              this%get_name(), chars=[task], ints=[first_record])
      endif
!
      if (last_record < 1) then
         call output%error_msg('Last record (i0) less than 1 for (a0)) in file' // &
                              this%get_name(), chars=[task], ints=[last_record])
      endif
!
      if (last_record < first_record) then
         call output%error_msg('Last record (i0) less than first record (i0) &
                               &for (a0) in file ' // this%get_name(), &
                               chars=[task], ints=[last_record, first_record])
      endif
!
   end subroutine check_records_to_be_accessed
!
!
end module direct_stream_file_class
