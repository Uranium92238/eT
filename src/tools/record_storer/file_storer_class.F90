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
module file_storer_class
!
!!
!!    File storer class module
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Currently only for fixed-length real(dp) arrays.
!!
!!    Mimics the functionality of direct access files by handling
!!    the storage of a set of vectors (each vector corresponds to
!!    a record). Uses a direct stream file.
!!
!
   use kinds
   use record_storer_class, only: record_storer
   use direct_stream_file_class, only: direct_stream_file
   use memory_manager_class, only: mem
   use global_out, only: output
!
   type, extends(record_storer) :: file_storer
!
      logical :: delete 
!
      logical, dimension(:), allocatable :: written_to_record
!
      type(direct_stream_file), allocatable :: direct_file
!
   contains
!
      procedure :: get                    => get_file_storer
      procedure :: set                    => set_file_storer
!
      procedure :: open_                  => open_file_storer
      procedure :: close_                 => close_file_storer
!
      procedure :: delete_                => delete_file_storer
!
      final :: destructor
!
      procedure :: initialize_storer      => initialize_storer_file_storer
      procedure :: finalize_storer        => finalize_storer_file_storer
!
   end type file_storer
!
!
   interface file_storer
!
      procedure :: new_file_storer
!
   end interface file_storer
!
!
contains
!
!
   function new_file_storer(name_, record_dim, n_records, delete) result(storer)
!!
!!    File storer constructer
!!    Writen by Eirik F. Kjønstad, 2019
!!
!!    name_:      Name of sequential storer (used for sequential file names).
!!                The files are given the name "name_record_number".
!!
!!    record_dim: The length of each record (i.e. the dimensionality of the array
!!                you want to store in each record)
!!
!!    n_records:  The number of records (you can store arrays in
!!                records 1,2,3,..., n_records)
!!
!!    delete:     Delete the files where the records are stored
!!                when the storer is deallocated.
!!
      implicit none
!
      type(file_storer) :: storer
!
      character(len=*), intent(in) :: name_
      integer, intent(in) :: record_dim
      integer, intent(in) :: n_records
      logical, intent(in) :: delete
!
      storer%name_       = trim(name_)
      storer%record_dim  = record_dim
      storer%n_records   = n_records
      storer%delete      = delete
!
!     Initialize file 
!
      call output%printf('v', 'Using direct stream file to store records with name: ' &
                            // trim(storer%name_) // '.', fs='(/t3,a)')
!
      storer%direct_file = direct_stream_file(trim(storer%name_), storer%record_dim)
!
   end function new_file_storer
!
!
   subroutine get_file_storer(storer, x, n)
!!
!!    Get
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Reads record n and returns the result in the array x.
!!
      implicit none
!
      class(file_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(out) :: x
!
      integer, intent(in) :: n
!
      integer :: record
!
      record = storer%record_indices(n)
!
      call storer%direct_file%read_(x, record, record)
!
   end subroutine get_file_storer
!
!
   subroutine set_file_storer(storer, x, n)
!!
!!    Set
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Writes x to record n.
!!
      implicit none
!
      class(file_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(in) :: x
!
      integer, intent(in) :: n
!
      integer :: record
!
      record = storer%record_indices(n)
!
      call storer%direct_file%write_(x, record, record)
!
      storer%written_to_record(record) = .true.
!
   end subroutine set_file_storer
!
!
   subroutine open_file_storer(storer)
!!
!!    Open
!!    Written by Anders Hutcheson, 2019
!!
!!    Opens all files.
!!
      implicit none
!
      class(file_storer) :: storer
!
      call storer%direct_file%open_()
!
   end subroutine open_file_storer
!
!
   subroutine close_file_storer(storer)
!!
!!    Close
!!    Written by Anders Hutcheson, 2019
!!
!!    Closes all files.
!!
      implicit none
!
      class(file_storer) :: storer
!
      call storer%direct_file%close_()
!
   end subroutine close_file_storer
!
!
   subroutine delete_file_storer(storer)
!!
!!    Delete
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Deletes all the files in the array.
!!
      implicit none
!
      class(file_storer) :: storer
!
      if (any(storer%written_to_record)) &
         call storer%direct_file%delete_()
!
   end subroutine delete_file_storer
!
!
   subroutine initialize_storer_file_storer(storer)
!!
!!    Initialize storer 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Initializes indices for records.
!!    Opens file if direct.
!!
      implicit none 
!
      class(file_storer) :: storer 
!
      integer :: I 
!
      call output%printf('debug', 'Doing preparations for file storer (a0)', &
                         chars=[storer%name_], fs='(/t3,a)')
!
!     Set up logical array telling us record has been written to or not 
!     
      call mem%alloc(storer%written_to_record, storer%n_records)
!
      storer%written_to_record = .false.
!
!     Set up index array telling us which record is 
!     stored in which position
!
      call mem%alloc(storer%record_indices, storer%n_records)
!
      do I = 1, storer%n_records
!
         storer%record_indices(I) = I
!
      enddo
!
!     Open file 
!
      call storer%direct_file%open_('readwrite')
!
   end subroutine initialize_storer_file_storer
!
!
   subroutine finalize_storer_file_storer(storer)
!!
!!    Finalize storer 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Deallocates index array. 
!!
!!    Deletes file(s) if requested. 
!!
!!    If direct file, closes the file.     
!!
      implicit none 
!
      class(file_storer) :: storer 
!
      call output%printf('debug', 'Doing finalizations for file storer (a0)', &
                         chars=[storer%name_], fs='(/t3,a)')
!
      call mem%dealloc(storer%record_indices, storer%n_records)
!
      call storer%direct_file%close_('keep')
!
      if (storer%delete) call storer%delete_()
!
      call mem%dealloc(storer%written_to_record, storer%n_records)     
!
   end subroutine finalize_storer_file_storer
!
!
   subroutine destructor(storer)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(file_storer) :: storer
!
      deallocate(storer%direct_file)
!
   end subroutine destructor
!
!
end module file_storer_class
