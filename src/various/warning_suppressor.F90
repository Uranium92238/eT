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
module warning_suppressor
!!
!!    Warning suppressor module  
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Can be used to suppress unused variable warning for variable x by "call do_nothing(x)".
!!
!
   use kinds
!
   implicit none 
!
   interface do_nothing
!!
!!    Does nothing to x!
!!
!!    NEVER use this routine just to avoid valid compiler warnings!
!!
!!    Should be used ONLY when we MUST suppress the "unused variable" compiler warning.
!! 
!!    A valid use-case is in null-object pattern routines.
!!
      procedure :: do_nothing_rank_0,           &
                   do_nothing_rank_1_real,      &
                   do_nothing_rank_1_integer,   &
                   do_nothing_rank_1_complex,   &
                   do_nothing_rank_2_real,      &
                   do_nothing_rank_2_integer,   &
                   do_nothing_rank_2_complex,   &
                   do_nothing_rank_3_real,      &
                   do_nothing_rank_3_integer,   &
                   do_nothing_rank_3_complex,   &
                   do_nothing_rank_4_real,      &
                   do_nothing_rank_4_integer,   &
                   do_nothing_rank_4_complex
!
   end interface do_nothing
!
contains
!
   subroutine do_nothing_rank_0(x)
!!
!!    Do nothing rank 0
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      class(*), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_0
!
   subroutine do_nothing_rank_1_real(x)
!!
!!    Do nothing rank 1
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      real(dp), dimension(:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_1_real
!
   subroutine do_nothing_rank_1_integer(x)
!!
!!    Do nothing rank 1
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      integer, dimension(:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_1_integer
!
   subroutine do_nothing_rank_1_complex(x)
!!
!!    Do nothing rank 1
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      complex(dp), dimension(:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_1_complex
!
   subroutine do_nothing_rank_2_real(x)
!!
!!    Do nothing rank 2
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      real(dp), dimension(:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_2_real
!
   subroutine do_nothing_rank_2_integer(x)
!!
!!    Do nothing rank 2
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      integer, dimension(:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_2_integer
!
   subroutine do_nothing_rank_2_complex(x)
!!
!!    Do nothing rank 2
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      complex(dp), dimension(:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_2_complex
!
   subroutine do_nothing_rank_3_real(x)
!!
!!    Do nothing rank 3
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      real(dp), dimension(:,:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_3_real
!
   subroutine do_nothing_rank_3_integer(x)
!!
!!    Do nothing rank 3
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      integer, dimension(:,:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_3_integer
!
   subroutine do_nothing_rank_3_complex(x)
!!
!!    Do nothing rank 3
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      complex(dp), dimension(:,:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_3_complex
!
   subroutine do_nothing_rank_4_real(x)
!!
!!    Do nothing rank 4
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      real(dp), dimension(:,:,:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_4_real
!
   subroutine do_nothing_rank_4_integer(x)
!!
!!    Do nothing rank 4
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      integer, dimension(:,:,:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_4_integer
!
   subroutine do_nothing_rank_4_complex(x)
!!
!!    Do nothing rank 4
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      complex(dp), dimension(:,:,:,:), intent(in) :: x 
!
      associate(y => x); end associate
!
   end subroutine do_nothing_rank_4_complex
!
!
end module warning_suppressor
!
