!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module named_interval_class
!!
!!    Named interval class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2020
!!
!!    Extends the 'interval' class (first, last, length) to have a 'name' as well.
!!
   use interval_class, only: interval   
!
   implicit none
!
   type, extends(interval) :: named_interval
!
      character(len=200) :: name_ = 'no_name'
!
   end type named_interval
!
   interface named_interval
!
      procedure :: new_named_interval 
!
   end interface named_interval 
!
contains 
!
!
   function new_named_interval(name_, first, last) result(this)
!!
!!    New named interval 
!!    Written by Eirik F. Kjønstad, 2020 
!!
      implicit none 
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: first, last 
!
      type(named_interval) :: this 
!
      this%name_ = trim(name_)
!
      this%first = first 
      this%last  = last 
      this%length = last - first + 1
!
   end function new_named_interval
!
!
end module named_interval_class
