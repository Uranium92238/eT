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
module memory_storer_class
!
!!
!!    Memory storer class module
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Currently only for unformatted fixed-length real(dp) arrays.
!!
!!    Mimics the functionality of direct access files by handling 
!!    the storage of a set of vectors (each vector corresponds to 
!!    a record). Stores records in memory, and is lost when deallocated.
!!
!
   use kinds 
   use record_storer_class, only: record_storer 
   use memory_manager_class, only: mem   
   use global_out, only: output  
!
   type, extends(record_storer) :: memory_storer
!
      logical :: delete
!
      real(dp), dimension(:,:), allocatable :: array ! records are stored in columns of this array 
!
   contains
!
      procedure :: get => get_memory_storer
      procedure :: set => set_memory_storer
!
      procedure :: initialize_storer   => initialize_storer_memory_storer 
      procedure :: finalize_storer     => finalize_storer_memory_storer 
!
   end type memory_storer
!
!
   interface memory_storer
!
      procedure :: new_memory_storer
!
   end interface memory_storer
!
!
contains
!
!
   function new_memory_storer(name_, record_dim, n_records) result(storer)
!!
!!    Memory storer constructer
!!    Writen by Eirik F. Kjønstad, 2019
!!
!!    name_:      Name of sequential storer.
!!
!!    record_dim: The length of each record (i.e. the dimensionality of the array 
!!                you want to store in each record)
!!
!!    n_records:  The number of records (you can store arrays in 
!!                records 1,2,3,..., n_records)
!!
      implicit none
!
      type(memory_storer) :: storer
!
      character(len=*), intent(in) :: name_
      integer, intent(in) :: record_dim
      integer, intent(in) :: n_records 
!
      storer%name_       = trim(name_)
      storer%record_dim  = record_dim 
      storer%n_records   = n_records 
!
   end function new_memory_storer
!
!
   subroutine get_memory_storer(storer, x, n)
!!
!!    Get 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Reads record n and returns the result in the array x. 
!!
      implicit none 
!
      class(memory_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(out) :: x 
!
      integer, intent(in) :: n 
!
      integer :: record
!
      record = storer%record_indices(n)
!
      call dcopy(storer%record_dim, storer%array(:,record), 1, x, 1)
!
   end subroutine get_memory_storer
!
!
   subroutine set_memory_storer(storer, x, n)
!!
!!    Set 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Writes x to record n.
!!
      implicit none 
!
      class(memory_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(in) :: x 
!
      integer, intent(in) :: n 
!
      integer :: record 
!
      record = storer%record_indices(n)
!
      call dcopy(storer%record_dim, x, 1, storer%array(:,record), 1)
!
   end subroutine set_memory_storer
!
!
   subroutine initialize_storer_memory_storer(storer)
!!
!!    Initialize storer  
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Allocates arrays for storing and indices. 
!!
      implicit none 
!
      class(memory_storer) :: storer 
!
      integer :: I 
!
      call output%printf('Doing preparations for memory storer (a0)', pl='debug', &
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
!     Allocate array to hold records 
!
      call mem%alloc(storer%array, storer%record_dim, storer%n_records)
!
   end subroutine initialize_storer_memory_storer
!
!
   subroutine finalize_storer_memory_storer(storer)
!!
!!    Finalize storer  
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Deallocates arrays for storing and indices.
!!
      implicit none 
!
      class(memory_storer) :: storer 
!
      call output%printf('Doing finalizations for memory storer (a0)', pl='debug', &
                           chars=[storer%name_], fs='(/t3,a)')
!
      call mem%dealloc(storer%record_indices, storer%n_records)
      call mem%dealloc(storer%array, storer%record_dim, storer%n_records)
!
   end subroutine finalize_storer_memory_storer
!
!
end module memory_storer_class
