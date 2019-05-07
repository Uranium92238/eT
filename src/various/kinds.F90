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
module kinds
!
!!
!!    Kinds module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Defines a set of kinds in terms of precision. Usage:
!!
!!       In declarations:   "real(dp) :: foo", "integer(i15) :: foo_int", etc.
!!       In record-lengths: "recl=dp*200" (200 double precision numbers per record)
!!
!
   implicit none
!
   integer, parameter :: dp  = selected_real_kind(15,307)
   integer, parameter :: qp  = selected_real_kind(33,4931)
   integer, parameter :: i15 = selected_int_kind(15)
   integer, parameter :: i6  = selected_int_kind(6)
!
end module kinds
