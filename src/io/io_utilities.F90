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
   use global_out, only : output
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
                  read_1_array_compound_record_1batch,      &
                  read_1_array_compound_record_0batches
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
      procedure   write_array_compound_record_2batches,  &
                  write_array_compound_record_1batch,    &
                  write_array_compound_record_0batches
!
   end interface compound_record_writer
!
!
contains
!
!
   subroutine read_1_array_single_record(dim_z, file_1, g_pqrz)
!!
!!    Read one array with single index as record
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read the whole direct access file "file_1" with record z into g_pqrz
!!    NB: It is assumed that z is the last index
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
         call file_1%read_(g_pqrz(:,:,:,z),z)
!
      enddo
!
   end subroutine read_1_array_single_record
!
!
   subroutine read_1_array_single_record_batch(batch_z, file_1, g_pqrz)
!!
!!    Read one array with single index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access file "file_1" into g_pqrz for the current batch z
!!    NB: It is assumed that the batching index is the last index
!!        and that the record is equal to the batching index
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
         call file_1%read_(g_pqrz(:,:,:,z),z_abs)
!
      enddo
!
   end subroutine read_1_array_single_record_batch
!
!
   subroutine read_2_arrays_single_record_batch(batch_z, file_1, g_pqrz, file_2, g_stuz)
!!
!!    Read two arrays with single index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access files "file_1/2" for the current batch z
!!    NB: It is assumed that the batching index is the last index
!!        and that the record is equal to the batching index
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
         call file_1%read_(g_pqrz(:,:,:,z),z_abs)
!
      enddo
!
!     Second file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_2%read_(g_stuz(:,:,:,z),z_abs)
!
      enddo
!
   end subroutine read_2_arrays_single_record_batch
!
!
   subroutine read_3_arrays_single_record_batch(batch_z, file_1, g_pqrz, file_2, g_stuz, file_3, g_vwxz)
!!
!!    Read three arrays with single index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access files "file_1/2/3" for the current batch z
!!    NB: It is assumed that the batching index is the last index
!!        and that the record is equal to the batching index
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
         call file_1%read_(g_pqrz(:,:,:,z),z_abs)
!
      enddo
!
!     Second file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_2%read_(g_stuz(:,:,:,z),z_abs)
!
      enddo
!
!     Third file
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_3%read_(g_vwxz(:,:,:,z),z_abs)
!
      enddo
!
   end subroutine read_3_arrays_single_record_batch
!
!
   subroutine read_1_array_compound_record_0batches(dim_z, dim_y, file_1, g_pqzy)
!!
!!    Read one array with compound index as record
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read whole direct access file "file_1" into g_pqzy for the current batch z
!!    NB: It is assumed that the batching indices are the last indices of the
!!        array (z,y order) and that the record is equal to the compound index zy
!!
      implicit none
!
      integer, intent(in) :: dim_z
      integer, intent(in) :: dim_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_pqzy
!
      type(direct_file), intent(in) :: file_1
!
      integer :: record
      integer :: z, y
!
      do y = 1, dim_y
         do z = 1, dim_z
!
            record = dim_z*(y - 1) + z
!
            call file_1%read_(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_1_array_compound_record_0batches
!
!
   subroutine read_1_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
!!
!!    Read one array with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access file "file_1" with records of zy into g_pqzy
!!    NB: It is assumed that the batching indices are the last indices of the
!!        array (z,y order) and that the record is equal to the compound index zy
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
            call file_1%read_(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_1_array_compound_record_2batches
!
!
   subroutine read_1_array_compound_record_1batch(dim_z, batch_y, file_1, g_pqzy, switch)
!!
!!    Read one array with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access file "file_1" with records of zy
!!    Reads 1 index in full dimension depending on switch
!!    NB: It is assumed that the batching indices are the last indices of the
!!        array (z,y order) and that the record is equal to the compound index zy
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
      batch_z = batching_index(dim_z)
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
!!    Read two arrays with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access file "file_1/2" with records of zy
!!    NB: It is assumed that the batching indices are the last indices of the
!!        array (z,y order) and that the record is equal to the compound index zy
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
            call file_1%read_(g_pqzy(:,:,z,y), record)
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
            call file_2%read_(g_rszy(:,:,z,y), record)
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
!!    Read three arrays with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access file "file_1/2/3" with records of zy
!!    NB: It is assumed that the batching indices are the last indices of the
!!        arrays (z,y order) and that the record is equal to the compound index zy
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
            call file_1%read_(g_pqzy(:,:,z,y), record)
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
            call file_2%read_(g_rszy(:,:,z,y), record)
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
            call file_3%read_(g_tuzy(:,:,z,y), record)
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
!!    Read four arrays with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Read parts of the direct access file "file_1/2/3/4" with records of zy
!!    NB: It is assumed that the batching indices are the last indices of the
!!        arrays (z,y order) and that the record is equal to the compound index zy
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
            call file_1%read_(g_pqzy(:,:,z,y), record)
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
            call file_2%read_(g_rszy(:,:,z,y), record)
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
            call file_3%read_(g_tuzy(:,:,z,y), record)
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
            call file_4%read_(g_vwzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine read_4_arrays_compound_record_2batches
!
!
   subroutine write_array_single_record(dim_z, file_1, g_pqrz)
!!
!!    Write array with single index as record
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Writes parts of the direct access file "file_1" with record in z
!!    NB: It is assumed that the batching index is the last index
!!
      implicit none
!
      integer, intent(in) :: dim_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(in) :: g_pqrz
!
      type(direct_file), intent(in) :: file_1
!
      integer :: z
!
      do z = 1, dim_z
!
         call file_1%write_(g_pqrz(:,:,:,z), z)
!
      enddo
!
   end subroutine write_array_single_record
!
!
   subroutine write_array_single_record_batch(batch_z, file_1, g_pqrz)
!!
!!    Write array with single index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Writes parts of the direct access file "file_1" with record in z
!!    NB: It is assumed that the batching index is the last index
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
!
      real(dp), dimension(:,:,:,:), contiguous, intent(in) :: g_pqrz
!
      type(direct_file), intent(in) :: file_1
!
      integer :: z, z_abs
!
      do z = 1, batch_z%length
!
         z_abs = batch_z%first + z - 1
!
         call file_1%write_(g_pqrz(:,:,:,z), z_abs)
!
      enddo
!
   end subroutine write_array_single_record_batch
!
!
   subroutine write_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
!!
!!    Write array with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Writes parts of the direct access file "file_1" with records in z,y
!!    NB: It is assumed that the array is split in batches of z and y
!!        which are sorted (z,y order) at the end of the array
!!
      implicit none
!
      type(batching_index), intent(in) :: batch_z
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(in) :: g_pqzy
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
            call file_1%write_(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine write_array_compound_record_2batches
!
!
   subroutine write_array_compound_record_1batch(dim_z, batch_y, file_1, g_pqzy, reverse)
!!
!!    Write array with compound index as record (batched)
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Writes parts of the direct access file "file_1" with records in z,y
!!    NB: It is assumed that z is of full dimension and y is batched over
!!        z and y are sorted (z,y order) at the end of the array
!!
!!    Writes an array to a direct access file "file_1" with records zy
!!
      implicit none
!
      integer, intent(in) :: dim_z
!
      type(batching_index), intent(in) :: batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(in) :: g_pqzy
!
      type(direct_file), intent(in) :: file_1
!
      logical, optional, intent(in) :: reverse
!
      type(batching_index) :: batch_z
!
      logical :: rev
!
      if(present(reverse)) then
         rev = reverse
      else
         rev = .false.
      endif
!
!     Fake a batching_index with full dimensions and call 2batches_writer
!
      batch_z = batching_index(dim_z)
!
      batch_z%first = 1
      batch_z%last = dim_z
      batch_z%length = dim_z
      batch_z%max_length = dim_z
      batch_z%num_batches = 1
!
      if (.not. rev) then
         call write_array_compound_record_2batches(batch_z, batch_y, file_1, g_pqzy)
      else
         call write_array_compound_record_2batches(batch_y, batch_z, file_1, g_pqzy)
      endif
!
   end subroutine write_array_compound_record_1batch
!
!
   subroutine write_array_compound_record_0batches(dim_z, dim_y, file_1, g_pqzy)
!!
!!    Write array with compound index as record
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Writes whole direct access file "file_1" with records in z,y
!!    NB: It is assumed that z and y are sorted (z,y order) at the end of the array
!!
      implicit none
!
      integer, intent(in) :: dim_z
      integer, intent(in) :: dim_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(in) :: g_pqzy
!
      type(direct_file), intent(in) :: file_1
!
      integer :: record
      integer :: z, y
!
      do y = 1, dim_y
         do z = 1, dim_z
!
            record = (y - 1) * dim_z + z
!
            call file_1%write_(g_pqzy(:,:,z,y), record)
!
         enddo
      enddo
!
   end subroutine write_array_compound_record_0batches
!
!
end module
