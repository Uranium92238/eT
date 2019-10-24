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
module sequential_storer_class
!
!!
!!    Sequential storer class module
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Currently only for unformatted fixed-length real(dp) arrays.
!!
!!    Mimics the functionality of direct access files by handling
!!    the storage of a set of vectors (each vector corresponds to
!!    a record). Stores records in sequential files. Can be kept
!!    or deleted when deallocated (see constructor).
!!
!
   use kinds
   use record_storer_class, only: record_storer
   use sequential_file_class, only: sequential_file
   use memory_manager_class, only: mem
   use global_out, only: output
!
   type, extends(record_storer) :: sequential_storer
!
      logical :: delete
      integer :: n_preexisting_files
!
      type(sequential_file), dimension(:), allocatable :: files
!
   contains
!
      procedure :: get      => get_sequential_storer
      procedure :: set      => set_sequential_storer
!
      procedure :: open_     => open_sequential_storer
      procedure :: close_    => close_sequential_storer
!
      procedure :: delete_  => delete_sequential_storer
!
      final :: destructor
!
   end type sequential_storer
!
!
   interface sequential_storer
!
      procedure :: new_sequential_storer
!
   end interface sequential_storer
!
!
contains
!
!
   function new_sequential_storer(name_, record_dim, n_records, delete) result(storer)
!!
!!    Sequential storer constructer
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
      type(sequential_storer) :: storer
!
      character(len=*), intent(in) :: name_
      integer, intent(in) :: record_dim
      integer, intent(in) :: n_records
      logical, intent(in) :: delete
!
      integer :: I
      character(len=200) :: record_name
!
      storer%name_       = trim(name_)
      storer%record_dim  = record_dim
      storer%n_records   = n_records
      storer%delete      = delete
!
      storer%n_preexisting_files = 0
!
      call mem%alloc(storer%record_indices, storer%n_records)
      allocate(storer%files(storer%n_records))
!
      do I = 1, storer%n_records
!
         storer%record_indices(I) = I
!
         write(record_name, '(a, i3.3)') trim(storer%name_), I
         storer%files(I) = sequential_file(record_name)
!
         call storer%files(I)%open_('readwrite', 'rewind')
!
         if (storer%files(I)%number_of_records() .ge. 1) then
            storer%n_preexisting_files = storer%n_preexisting_files + 1
         endif
!
      enddo
!
   end function new_sequential_storer
!
!
   subroutine get_sequential_storer(storer, x, n)
!!
!!    Get
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Reads record n and returns the result in the array x.
!!
      implicit none
!
      class(sequential_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(out) :: x
!
      integer, intent(in) :: n
!
      integer :: record
!
      record = storer%record_indices(n)
!
      call storer%files(record)%rewind_()
      call storer%files(record)%read_(x, storer%record_dim)
!
   end subroutine get_sequential_storer
!
!
   subroutine set_sequential_storer(storer, x, n)
!!
!!    Set
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Writes x to record n.
!!
      implicit none
!
      class(sequential_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(in) :: x
!
      integer, intent(in) :: n
!
      integer :: record
!
      record = storer%record_indices(n)
!
      call storer%files(record)%rewind_()
      call storer%files(record)%write_(x, storer%record_dim)
!
   end subroutine set_sequential_storer
!
!
   subroutine open_sequential_storer(storer)
!!
!!    Open
!!    Written by Anders Hutcheson, 2019
!!
!!    Opens all files in the array.
!!
      implicit none
!
      class(sequential_storer) :: storer
!
      integer :: I
!
      do I = 1, storer%n_records
!
         call storer%files(I)%open_()
!
      enddo
!
   end subroutine open_sequential_storer
!
!
   subroutine close_sequential_storer(storer)
!!
!!    Close
!!    Written by Anders Hutcheson, 2019
!!
!!    Closes all files in the array
!!
      implicit none
!
      class(sequential_storer) :: storer
!
      integer :: I
!
      do I = 1, storer%n_records
!
         call storer%files(I)%close_()
!
      enddo
!
   end subroutine close_sequential_storer
!
!
   subroutine delete_sequential_storer(storer)
!!
!!    Delete
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Deletes all the files in the array.
!!
      implicit none
!
      class(sequential_storer) :: storer
!
      integer :: I
!
      do I = 1, storer%n_records
!
         call storer%files(I)%delete_()
!
      enddo
!
   end subroutine delete_sequential_storer
!
!
   subroutine destructor(storer)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(sequential_storer) :: storer
!
      integer :: I
!
      do I = 1, storer%n_records
!
         call storer%files(I)%close_()
!
      enddo
!
      if (storer%delete) call storer%delete_()
!
      call mem%dealloc(storer%record_indices, storer%n_records)
      deallocate(storer%files)
!
   end subroutine destructor
!
!
end module sequential_storer_class
