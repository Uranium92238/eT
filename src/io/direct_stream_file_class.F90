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
module direct_stream_file_class
!
!!
!!    Direct access stream file class module
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Stream file that acts as a Fortran direct file
!!    Byte positions and read/write lengths are calculated based on record dimension
!!    in new_direct_stream and records sent in to read and write.
!!
!!    Because of various compiler bugs, specific wrappers like real_2 must be called
!!    where the number refers to the rank of the array.
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
      procedure :: read_1_real_dp_direct_stream_file
      procedure :: read_2_real_dp_direct_stream_file
      procedure :: read_3_real_dp_direct_stream_file
      procedure :: read_4_real_dp_direct_stream_file
!
!     Complex double precision read
!
      procedure :: read_1_complex_dp_direct_stream_file
      procedure :: read_2_complex_dp_direct_stream_file
      procedure :: read_3_complex_dp_direct_stream_file
      procedure :: read_4_complex_dp_direct_stream_file
!
!     Read generic
!
      generic :: read_ => read_1_real_dp_direct_stream_file,&
                          read_2_real_dp_direct_stream_file,&
                          read_3_real_dp_direct_stream_file,&
                          read_4_real_dp_direct_stream_file,&
                          read_1_complex_dp_direct_stream_file,&
                          read_2_complex_dp_direct_stream_file,&
                          read_3_complex_dp_direct_stream_file,&
                          read_4_complex_dp_direct_stream_file

!
!     Specialized read routines
!
      procedure :: read_interval            => read_interval_direct_stream_file
      procedure :: read_compound_full_batch => read_compound_full_batch_direct_stream_file
      procedure :: read_compound_batch_full => read_compound_batch_full_direct_stream_file
      procedure :: read_compound            => read_compound_direct_stream_file
!
!     Write routines
!
!     Real double precision write
!
      procedure :: write_1_real_dp_direct_stream_file
      procedure :: write_2_real_dp_direct_stream_file
      procedure :: write_3_real_dp_direct_stream_file
      procedure :: write_4_real_dp_direct_stream_file
!
!     complex double precision write
!
      procedure :: write_1_complex_dp_direct_stream_file
      procedure :: write_2_complex_dp_direct_stream_file
      procedure :: write_3_complex_dp_direct_stream_file
      procedure :: write_4_complex_dp_direct_stream_file
!
!     Write generic
!
      generic :: write_ => write_1_real_dp_direct_stream_file,&
                           write_2_real_dp_direct_stream_file,&
                           write_3_real_dp_direct_stream_file,&
                           write_4_real_dp_direct_stream_file,&
                           write_1_complex_dp_direct_stream_file,&
                           write_2_complex_dp_direct_stream_file,&
                           write_3_complex_dp_direct_stream_file,&
                           write_4_complex_dp_direct_stream_file
!
      procedure :: write_interval            => write_interval_direct_stream_file
      procedure :: write_compound_full_batch => write_compound_full_batch_direct_stream_file
      procedure :: write_compound_batch_full => write_compound_batch_full_direct_stream_file
      procedure :: write_compound            => write_compound_direct_stream_file
!
      procedure :: get_n_records             => get_n_records_direct_stream_file
!
      final :: destructor
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
   function new_direct_stream_file(name_, rec_dim, w_size, status_) result(the_file)
!!
!!    Direct file constructer
!!    Writen by Rolf H. Myhre, May 2019
!!
!!    rec_dim is number of words in each record
!!    w_size (optional) is the size of each word, default is double precision
!!    record length is rec_dim*w_size
!!
      implicit none
!
      type(direct_stream_file), allocatable      :: the_file
!
      character(len=*), intent(in)           :: name_
      integer, intent(in)                    :: rec_dim
!
      integer, intent(in), optional          :: w_size
      character(len=*), intent(in), optional :: status_
!
      allocate(the_file)
!
      if (present(w_size)) then
         if (w_size .gt. 0) then
            the_file%word_size = w_size
         else
            call output%error_msg("Word size less than one for file (a0)", &
                                  chars=[trim(name_)])
         endif
      else
         the_file%word_size = dp
      endif
!
      call the_file%set_name(name_)
!
      if (rec_dim .lt. 1) then
         call output%error_msg("Record dimension less than one for file (a0)", &
                               chars=[trim(name_)])
      endif
!
      if (present(status_)) then
         call the_file%set_status(status_)
      else
         call the_file%set_status('unknown')
      endif
!
      the_file%record_dim = rec_dim
      the_file%record_length = rec_dim*the_file%word_size
!
   end function new_direct_stream_file
!
!
   subroutine destructor(the_file)
!!
!!    Destructor
!!    Written by Rolf H. Myhre, Feb. 2020
!!
      implicit none
!
      type(direct_stream_file) :: the_file
!
      if (the_file%get_open()) then
         call output%error_msg('Destructor for file (a0) called &
                               &while the file is still open', chars=[the_file%get_name()])
      endif
!
   end subroutine destructor
!
!
   subroutine read_1_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
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
!!    Calculates number of bytes to read based on first_rec, last_rec and the_file%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      real(dp), dimension((last_rec - first_rec + 1)*the_file%record_dim), intent(out) :: array
!
      integer :: position_, read_length, records
!
      records = last_rec - first_rec + 1
!
      if (first_rec .ge. 1) then
         position_ = (first_rec-1)*the_file%record_length + 1
      else
         call output%error_msg('Record (i0) less than 1 for read in file (a0)', &
                               chars=[the_file%get_name()], ints=[first_rec])
      endif
!
      if (records .ge. 1) then
         read_length = the_file%record_dim*records
      else
         call output%error_msg('Last record (i0) less than first record (i0) &
                               &for read in file (a0)', &
                               chars=[trim(the_file%get_name())], ints=[first_rec, last_rec])
      endif
!
      call the_file%read_real_dp(array, read_length, position_)
!
   end subroutine read_1_real_dp_direct_stream_file
!
!
!! Wrapper routines for read_1_real_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 20202
!
   subroutine read_2_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:), intent(out) :: array
!
      call the_file%read_1_real_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_2_real_dp_direct_stream_file
!
   subroutine read_3_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:), intent(out) :: array
!
      call the_file%read_1_real_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_3_real_dp_direct_stream_file
!
   subroutine read_4_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:,:), intent(out) :: array
!
      call the_file%read_1_real_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_4_real_dp_direct_stream_file
!
!
   subroutine read_1_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
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
!!    Calculates number of bytes to read based on first_rec, last_rec and the_file%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      complex(dp), dimension((last_rec - first_rec + 1)*the_file%record_dim), intent(out) :: array
!
      integer :: position_, read_length, records
!
      records = last_rec - first_rec + 1
!
      if (first_rec .ge. 1) then
         position_ = (first_rec-1)*the_file%record_length + 1
      else
         call output%error_msg('Record (i0) less than 1 for read in file (a0)', &
                               chars=[the_file%get_name()], ints=[first_rec])
      endif
!
      if (records .ge. 1) then
         read_length = the_file%record_dim*records
      else
         call output%error_msg('Last record (i0) less than first record (i0) &
                               &for read in file (a0)', &
                               chars=[trim(the_file%get_name())], ints=[first_rec, last_rec])
      endif
!
      call the_file%read_complex_dp(array, read_length, position_)
!
   end subroutine read_1_complex_dp_direct_stream_file
!
!
!! Wrapper routines for read_1_complex_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 20202
!
   subroutine read_2_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:), intent(out) :: array
!
      call the_file%read_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_2_complex_dp_direct_stream_file
!
   subroutine read_3_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:), intent(out) :: array
!
      call the_file%read_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_3_complex_dp_direct_stream_file
!
   subroutine read_4_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:,:), intent(out) :: array
!
      call the_file%read_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
!
   end subroutine read_4_complex_dp_direct_stream_file
!
!
   subroutine read_interval_direct_stream_file(the_file, array, batch_z)
!!
!!    read interval direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump file contents in.
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    batch_z :: batching index
!!
!!    Wrapper for read_1 that calculates first and last record from
!!    a batching index
!!
      use interval_class, only : interval
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(interval), intent(in) :: batch_z
!
      real(dp), dimension(the_file%record_dim*batch_z%length) :: array
!
      call the_file%read_1_real_dp_direct_stream_file(array, batch_z%first, batch_z%last)
!
   end subroutine read_interval_direct_stream_file
!
!
   subroutine read_compound_full_batch_direct_stream_file(the_file, array, dim_y, batch_z)
!!
!!    read compound full batch direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump file contents in.
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    dim_y   :: integer, full dimension of first component of batching index
!!    batch_z :: batching index
!!
!!    Wrapper for read_1 that loops over batch_z
!!    and reads all compound records in dim_y
!!
!!    Similar to read_compound, but takes the full y dimension
!!    for cases where y is not batched.
!!
      use interval_class, only : interval
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(interval), intent(in) :: batch_z
      integer, intent(in) :: dim_y
!
      real(dp), dimension(the_file%record_dim*dim_y*batch_z%length) :: array
!
      call the_file%read_1_real_dp_direct_stream_file(array, &
                           (batch_z%first-1)*dim_y + 1, &
                            batch_z%last*dim_y)
!
   end subroutine read_compound_full_batch_direct_stream_file
!
!
   subroutine read_compound_batch_full_direct_stream_file(the_file, array, batch_y, dim_z)
!!
!!    read compound batch full direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump file contents in
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    batch_y :: batching index
!!    dim_z   :: integer, full dimension of second component of batching index
!!
!!    Wrapper for read_1 that loops over dim_z
!!    and reads all compound records in batch_y
!!
!!    Note! There are very few safeguards in these routines and the compiler/runtime
!!          won't tell you if you mess up. Make very sure that you're batching
!!          dimensions are correct.
!!
!!    Note! This routine assumes that array is allocated something like (: ,batch_y%max_length, :)
!!          If array is reallocated in the batching loop and the second to last dimension
!!          changes, the routine will fail!
!!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(batching_index), intent(in) :: batch_y
      integer, intent(in) :: dim_z
!
      real(dp),dimension(the_file%record_dim*batch_y%max_length, dim_z) :: array
!
      integer :: z
!
!     Check if we can do a single continuous write, else we have to loop
      if (batch_y%length .eq. batch_y%index_dimension) then
!
         call the_file%read_1_real_dp_direct_stream_file(array, 1, batch_y%index_dimension*dim_z)
!
      else
         do z = 1, dim_z
!
            call the_file%read_1_real_dp_direct_stream_file(array(:, z), &
                                  batch_y%index_dimension*(z-1) + batch_y%first, &
                                  batch_y%index_dimension*(z-1) + batch_y%last)
         enddo
      endif
!
   end subroutine read_compound_batch_full_direct_stream_file
!
!
   subroutine read_compound_direct_stream_file(the_file, array, batch_y, batch_z)
!!
!!    read compound direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump file contents in.
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    batch_y :: batching index
!!    batch_z :: batching index
!!
!!    Wrapper for read_1_real_dp that loops over batch_z
!!    and reads all compound records in batch y
!!
!!    This routine is useful for coumpund records were the record number
!!    is calculated from two batched indexes
!!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(batching_index), intent(in) :: batch_y, batch_z
!
      real(dp), &
      dimension(the_file%record_dim*batch_y%max_length, batch_z%first:batch_z%last) :: array
!
      integer :: z
!
!     Check if we can do a single continuous read, else we have to loop
      if (batch_y%length .eq. batch_y%index_dimension) then
!
         call the_file%read_1_real_dp_direct_stream_file(array, &
                              batch_y%index_dimension*(batch_z%first-1) + 1, &
                              batch_y%index_dimension*batch_z%last)
      else
         do z = batch_z%first, batch_z%last
!
            call the_file%read_1_real_dp_direct_stream_file(array(:, z), &
                                 batch_y%index_dimension*(z-1) + batch_y%first, &
                                 batch_y%index_dimension*(z-1) + batch_y%last)
         enddo
      endif
!
   end subroutine read_compound_direct_stream_file
!
!
   subroutine write_1_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
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
!!    Calculates number of bytes to write based on first_rec, last_rec and the_file%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      real(dp), dimension((last_rec - first_rec + 1)*the_file%record_dim), intent(in) :: array
!
      integer :: position_, write_length, records
!
      records = last_rec - first_rec + 1
!
      if (first_rec .ge. 1) then
         position_ = (first_rec-1)*the_file%record_length + 1
      else
         call output%error_msg('Record (i0) less than 1 for write in file (a0)', &
                               chars=[the_file%get_name()], ints=[first_rec])
      endif
!
      if (records .ge. 1) then
         write_length = the_file%record_dim*records
      else
         call output%error_msg('Last record (i0) less than first record (i0) &
                               &for read in file (a0)', &
                               chars=[the_file%get_name()], ints=[first_rec, last_rec])
      endif
!
      call the_file%write_real_dp(array, write_length, position_)
!
   end subroutine write_1_real_dp_direct_stream_file
!
!
!! Wrapper routines for write_1_real_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 20202
!
   subroutine write_2_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:), intent(in) :: array
!
      call the_file%write_1_real_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_2_real_dp_direct_stream_file
!
   subroutine write_3_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:), intent(in) :: array
!
      call the_file%write_1_real_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_3_real_dp_direct_stream_file
!
   subroutine write_4_real_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      real(dp), dimension(:,:,:,:), intent(in) :: array
!
      call the_file%write_1_real_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_4_real_dp_direct_stream_file
!
!
   subroutine write_1_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
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
!!    Calculates number of bytes to write based on first_rec, last_rec and the_file%record_length
!!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      integer, intent(in) :: first_rec
      integer, intent(in) :: last_rec
!
      complex(dp), dimension((last_rec - first_rec + 1)*the_file%record_dim), intent(in) :: array
!
      integer :: position_, write_length, records
!
      records = last_rec - first_rec + 1
!
      if (first_rec .ge. 1) then
         position_ = (first_rec-1)*the_file%record_length + 1
      else
         call output%error_msg('Record (i0) less than 1 for write in file (a0)', &
                               chars=[the_file%get_name()], ints=[first_rec])
      endif
!
      if (records .ge. 1) then
         write_length = the_file%record_dim*records
      else
         call output%error_msg('Last record (i0) less than first record (i0) &
                               &for read in file (a0)', &
                               chars=[the_file%get_name()], ints=[first_rec, last_rec])
      endif
!
      call the_file%write_complex_dp(array, write_length, position_)
!
   end subroutine write_1_complex_dp_direct_stream_file
!
!
!! Wrapper routines for write_1_complex_dp_direct_stream_file
!! that accepts rank 2, 3, and 4 arrays
!!
!! Written by Rolf H. Myhre, Mar 20202
!
   subroutine write_2_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:), intent(in) :: array
!
      call the_file%write_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_2_complex_dp_direct_stream_file
!
   subroutine write_3_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:), intent(in) :: array
!
      call the_file%write_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_3_complex_dp_direct_stream_file
!
   subroutine write_4_complex_dp_direct_stream_file(the_file, array, first_rec, last_rec)
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
      integer, intent(in) :: first_rec, last_rec
      complex(dp), dimension(:,:,:,:), intent(in) :: array
!
      call the_file%write_1_complex_dp_direct_stream_file(array, first_rec, last_rec)
   end subroutine write_4_complex_dp_direct_stream_file
!
!
   subroutine write_interval_direct_stream_file(the_file, array, batch_z)
!!
!!    write interval direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump to file
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    batch_z :: batching index
!!
!!    Wrapper for write_1_real_dp that calculates first and last record from
!!    a batching index
!!
      use interval_class, only : interval
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(interval), intent(in) :: batch_z
!
      real(dp), dimension(the_file%record_dim*batch_z%length) :: array
!
      call the_file%write_1_real_dp_direct_stream_file(array, batch_z%first, batch_z%last)
!
   end subroutine write_interval_direct_stream_file
!
!
   subroutine write_compound_full_batch_direct_stream_file(the_file, array, dim_y, batch_z)
!!
!!    write compound full direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real(dp) to dump to file
!!             Only real(dp) allowed until GFortran figures out explicit shape
!!
!!    dim_y   :: integer, full dimension of first component of batching index
!!    batch_z :: batching index
!!
!!    Wrapper for write_1_real_dp that loops over batch_z
!!    and writes all compound records in dim_y
!!
!!    Similar to write_compound above, but takes the full y dimension
!!    for cases where y is not batched.
!!
      use interval_class, only : interval
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(interval), intent(in) :: batch_z
      integer, intent(in) :: dim_y
!
      real(dp), dimension(the_file%record_dim*dim_y*batch_z%length) :: array
!
      call the_file%write_1_real_dp_direct_stream_file(array, &
                                                       (batch_z%first-1)*dim_y + 1, &
                                                        batch_z%last*dim_y)
!
   end subroutine write_compound_full_batch_direct_stream_file
!
!
   subroutine write_compound_batch_full_direct_stream_file(the_file, array, batch_y, dim_z)
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
!!    Wrapper for write_1_real_dp that loops over dim_z
!!    and writes all compound records in batch_y
!!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(batching_index), intent(in) :: batch_y
      integer, intent(in) :: dim_z
!
      real(dp),dimension(the_file%record_dim*batch_y%max_length, dim_z) :: array
!
      integer :: z
!
!     Check if we can do a single continuous write, else we have to loop
      if (batch_y%length .eq. batch_y%index_dimension) then
!
         call the_file%write_1_real_dp_direct_stream_file(array, 1, batch_y%index_dimension*dim_z)
!
      else
         do z = 1, dim_z
!
            call the_file%write_1_real_dp_direct_stream_file(array(:, z), &
                                  batch_y%index_dimension*(z-1) + batch_y%first, &
                                  batch_y%index_dimension*(z-1) + batch_y%last)
         enddo
      endif
!
   end subroutine write_compound_batch_full_direct_stream_file
!
!
   subroutine write_compound_direct_stream_file(the_file, array, batch_y, batch_z)
!!
!!    write compound direct stream
!!    Written by Rolf H. Myhre and Alexander C. Paul, Feb. 2020
!!
!!    array :: real double precision array
!!
!!    batch_y :: batching index
!!    batch_z :: batching index
!!
!!    Wrapper for write_1_real_dp that loops over batch_z
!!    and writes all compound records in batch y
!!
!!    This routine is useful for coumpund records were the record number
!!    is calculated from two batched indexes
!!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(direct_stream_file), intent(in) :: the_file
!
      class(batching_index), intent(in) :: batch_y, batch_z
!
      real(dp), &
      dimension(the_file%record_dim*batch_y%max_length, batch_z%first:batch_z%last) :: array
!
      integer :: z
!
!     Check if we can do a single continuous write, else we have to loop
      if (batch_y%length .eq. batch_y%index_dimension) then
!
         call the_file%write_1_real_dp_direct_stream_file(array, &
                               batch_y%index_dimension*(batch_z%first-1) + 1, &
                               batch_y%index_dimension*batch_z%last)
      else
         do z = batch_z%first, batch_z%last
!
            call the_file%write_1_real_dp_direct_stream_file(array(:, z), &
                                  batch_y%index_dimension*(z-1) + batch_y%first, &
                                  batch_y%index_dimension*(z-1) + batch_y%last)
         enddo
      endif
!
!
   end subroutine write_compound_direct_stream_file
!
!
   function get_n_records_direct_stream_file(the_file) result(n_records)
!!
!!    Get number of existing records 
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Returns the number of records on the file.
!!
      implicit none 
!
      class(direct_stream_file) :: the_file
!
      integer :: n_records
!
      integer :: file_size 
!
      file_size = the_file%get_file_size()
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
         n_records = file_size/the_file%record_length
!
      endif 
!
   end function get_n_records_direct_stream_file
!
!
end module direct_stream_file_class
