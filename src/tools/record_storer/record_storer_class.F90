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
module record_storer_class
!
!!
!!    Record storer class module
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Mimics the functionality of direct access files by handling 
!!    the storage of a set of vectors (each vector corresponds to 
!!    a record). In descendants, this is either stored in memory 
!!    or on disk. 
!!
!
   use kinds
!
   type, abstract :: record_storer
!
      character(len=255) :: name_ = 'no_name'
!
      integer :: n_records 
      integer :: record_dim
!
      integer, dimension(:), allocatable :: record_indices ! map from record to storage position
!
   contains
!
      procedure(get_record_storer), deferred :: get
      procedure(set_record_storer), deferred :: set
!
      procedure :: cycle_left => cycle_left_record_storer
!
   end type record_storer
!
!
   abstract interface 
!
      subroutine get_record_storer(storer, x, n)
!!
!!       Get 
!!       Written by Eirik F. Kjønstad, 2019 
!!
!!       Retrieves record n and places the result in x 
!!
         import :: record_storer, dp
!
         implicit none 
!
         class(record_storer) :: storer 
!
         integer, intent(in) :: n 
!
         real(dp), dimension(storer%record_dim), intent(out) :: x 
!
      end subroutine get_record_storer
!
      subroutine set_record_storer(storer, x, n)
!!
!!       Set 
!!       Written by Eirik F. Kjønstad, 2019 
!!
!!       Stores x as record n
!!
         import :: record_storer, dp
!
         implicit none 
!
         class(record_storer) :: storer 
!
         integer, intent(in) :: n 
!
         real(dp), dimension(storer%record_dim), intent(in) :: x 
!
      end subroutine set_record_storer
!
   end interface 
!
!
contains
!
!
   subroutine cycle_left_record_storer(storer)
!!
!!    Cycle left 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    If we have a list of records A, B, C, D, this routine 
!!    reorders the list to B, C, D, A. Then, when overwriting 
!!    the last record (i.e. A), we remove the originally oldest 
!!    record and insert the newest record at the end. 
!!
      implicit none
!
      class(record_storer) :: storer
!
      integer :: I, first_record 
!
      first_record = storer%record_indices(1)
!
      do I = 1, storer%n_records - 1
!
         storer%record_indices(I) =  storer%record_indices(I + 1)
!
      enddo 
!
      storer%record_indices(storer%n_records) = first_record
!
   end subroutine cycle_left_record_storer
!
!
end module record_storer_class
