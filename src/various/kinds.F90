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
module kinds
!
!!
!!    Kinds module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Defines a set of kinds in terms of precision. Usage:
!!
!!       In declarations:   "real(dp) :: foo", "integer(i64) :: foo_int", etc.
!!       In record-lengths: "recl=dp*200" (200 double precision numbers per record)
!!
   use, intrinsic :: iso_fortran_env, only: real32, real64, real128
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
!
   implicit none
!
   integer, parameter :: sp  = real32
   integer, parameter :: dp  = real64
   integer, parameter :: qp  = real128
   integer, parameter :: i64 = int64
   integer, parameter :: i32 = int32
   integer, parameter :: i16 = int16
   integer, parameter :: i8  = int8
!
   integer, parameter :: int_size = storage_size(i64)/8
!
end module kinds
