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
module interval_class
!!
!!    Interval class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use kinds
!
   implicit none
!
   type interval
!
      integer :: first = -1
      integer :: last  = -1
      integer :: length = -1
!
   contains
!
   end type interval
!
   interface interval
!
      procedure :: new_interval 
!
   end interface interval 
!
contains 
!
!
   function new_interval(first, last) result(intval)
!!
!!    New interval 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      integer, intent(in) :: first, last 
!
      type(interval) :: intval 
!
      intval%first = first 
      intval%last  = last 
      intval%length = last - first + 1
!
   end function new_interval
!
!
end module interval_class
