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
module file_storer_class
!
!!
!!    File storer class module
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Currently only for unformatted fixed-length real(dp) arrays.
!!
!!    Mimics the functionality of direct access files by handling
!!    the storage of a set of vectors (each vector corresponds to
!!    a record). Stores records either in a direct file (when 
!!    recl < approx. 2 GB) or in an array of sequential files 
!!    (when recl > approx. 2GB)
!!
!
   use kinds
   use record_storer_class, only: record_storer
   use sequential_file_class, only: sequential_file
   use direct_file_class, only: direct_file
   use memory_manager_class, only: mem
   use global_out, only: output
!
   type, extends(record_storer) :: file_storer
!
      logical :: delete, direct_ 
!
      type(direct_file), allocatable :: direct_file
      type(sequential_file), dimension(:), allocatable :: sequential_files
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
      procedure :: get_n_existing_records => get_n_existing_records_file_storer
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
   function new_file_storer(name_, record_dim, n_records, delete, direct_) result(storer)
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
!!    direct:     (Optional) Whether to use a direct file. If not, the storer  
!!                will use an array of sequential files. Default is direct==.true.  
!!                for records that are shorter than 2 GB and .false. for records 
!!                that are longer.
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
      logical, intent(in), optional :: direct_ 
!
      integer :: I
      character(len=200) :: record_name
!
!     https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-record-length
!     recl limit: 2147483647 - minus the bytes for record overhead
!
      integer, parameter :: recl_limit = 2147000000 
!
      storer%name_       = trim(name_)
      storer%record_dim  = record_dim
      storer%n_records   = n_records
      storer%delete      = delete
!
!     Use direct file or sequential file array? 
!
      if (present(direct_)) then 
!
         storer%direct_ = direct_
!
      else 
!
         if (record_dim*dp > recl_limit) then 
!
            storer%direct_ = .false.
!
         else 
!
            storer%direct_ = .true.
!
         endif
!
      endif
!
!     Initialize file(s) 
!
      if (.not. storer%direct_) then 
!
!        Use array of sequential files 
!
         call output%printf('Using sequential file array to store records using the prefix: ' &
                              // trim(storer%name_) // '.', pl='verbose', fs='(/t3,a)')
!
         allocate(storer%sequential_files(storer%n_records))
!
         do I = 1, storer%n_records
!
            write(record_name, '(a, i3.3)') trim(storer%name_), I
            storer%sequential_files(I) = sequential_file(record_name)
!
         enddo
!
      else 
!
!        Use one direct file 
!
         call output%printf('Using direct file to store records with name: ' &
                              // trim(storer%name_) // '.', pl='verbose', fs='(/t3,a)')
!
         storer%direct_file = direct_file(trim(storer%name_), storer%record_dim)
!
      endif
!
   end function new_file_storer
!
!
   function get_n_existing_records_file_storer(storer) result(n_existing_records)
!!
!!    Get n existing records 
!!    Written by Rolf H. Myhre and Eirik F. Kjønstad, 2019
!!
      implicit none 
!
      class(file_storer) :: storer 
!
      integer :: I 
!
      integer :: n_existing_records
!
      if (storer%direct_) call output%error_msg('cannot ask for number of existing records '&
                                             // 'for storer using direct access file.')
!
      call storer%open_()
!  
      n_existing_records = 0
!
      do I = 1, storer%n_records
!
         if (storer%sequential_files(I)%number_of_records() .ge. 1) then 
!  
            n_existing_records = n_existing_records + 1
!
         endif
!
      enddo
!
      call storer%close_()
!
   end function get_n_existing_records_file_storer
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
      if (storer%direct_) then 
!
         call storer%direct_file%read_(x, record)
!
      else 
!
         call storer%sequential_files(record)%open_('read', 'rewind')
         call storer%sequential_files(record)%read_(x, storer%record_dim)
         call storer%sequential_files(record)%close_()
!
      endif
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
      if (storer%direct_) then 
!
         call storer%direct_file%write_(x, record)
!
      else
!
         call storer%sequential_files(record)%open_('write', 'rewind')
         call storer%sequential_files(record)%write_(x, storer%record_dim)
         call storer%sequential_files(record)%close_()
!
      endif
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
      integer :: I
!
      if (storer%direct_) then 
!
         call storer%direct_file%open_()
!
      else 
!
         do I = 1, storer%n_records
!
            call storer%sequential_files(I)%open_()
!
         enddo
!
      endif
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
      integer :: I
!
      if (storer%direct_) then 
!
         call storer%direct_file%close_()
!
      else 
!
         do I = 1, storer%n_records
!
            call storer%sequential_files(I)%close_()
!
         enddo
!
      endif
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
      integer :: I
!
      if (storer%direct_) then 
!
         call storer%direct_file%delete_()
!
      else 
!
         do I = 1, storer%n_records
!
            call storer%sequential_files(I)%delete_()
!
         enddo
!
      endif
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
      call output%printf('Doing preparations for file storer (a0)', pl='debug', &
                           chars=[storer%name_], fs='(/t3,a)')
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
!     Open file if direct 
!
      if (storer%direct_) call storer%direct_file%open_('readwrite')
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
      call output%printf('Doing finalizations for file storer (a0)', pl='debug', &
                           chars=[storer%name_], fs='(/t3,a)')
!
      call mem%dealloc(storer%record_indices, storer%n_records)
!
      if (storer%delete) then 
!
         call storer%delete_()
!
      else 
!
         if (storer%direct_) call storer%direct_file%close_('keep')
!
      endif      
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
      if (storer%direct_) then 
!
         deallocate(storer%direct_file)
!
      else 
!
         deallocate(storer%sequential_files)
!
      endif
!
   end subroutine destructor
!
!
end module file_storer_class
