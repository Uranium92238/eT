!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module array_initialization
!
!!
!!    Array initialization module
!!
!!    Routines that initialize different types of arrays to certain values.
!!
!
   use parameters
!
   implicit none
!
contains
!
!
   subroutine zero_array(x, n)
!!
!!    Zero array
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Sets the array x of length n to zero.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(out) :: x
!
      integer :: I
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n
!
         x(I) = zero
!
      enddo
!$omp end parallel do
!
   end subroutine zero_array
!
!
   subroutine zero_array_complex(X, n)
!!
!!    Zero array
!!    Written by Eirik F. Kjønstad, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
      implicit none
!
      integer, intent(in) :: n
!
      complex(dp), dimension(n), intent(out) :: X
!
      integer :: I
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n
!
         X(I) = zero_complex
!
      enddo
!$omp end parallel do
!
   end subroutine zero_array_complex
!
!
   subroutine zero_array_int(x, n)
!!
!!    Zero array int
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Sets the integer array x of length n to zero.
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer, dimension(n), intent(out) :: x
!
      integer :: I
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n
!
         x(I) = 0
!
      enddo
!$omp end parallel do
!
   end subroutine zero_array_int
!
!
   subroutine constant_array(x, n, const)
!!
!!    Constant array
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Sets the array x of length n to const.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(out) :: x
!
      real(dp), intent(in) :: const
!
      integer :: I
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n
!
         x(I) = const
!
      enddo
!$omp end parallel do
!
   end subroutine constant_array
!
!
   subroutine identity_array(x, n)
!!
!!    Identity array
!!    Written by Ida-Marie Hoyvik, Oct 2019
!!
!!    Sets the array x of dimension nxn to be the identity matrix.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n,n), intent(out) :: x
!
      integer :: i
!
      call zero_array(x,n**2)
!
!$omp parallel do private(i) schedule(static)
      do i = 1, n
!
         x(i,i) = one
!
      enddo
!$omp end parallel do
!
   end subroutine identity_array
!
!
   subroutine copy_and_scale(alpha, X, Y, n)
!!
!!    Copy and scale
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Sets Y as:
!!
!!       Y = alpha*X
!!
!!    X and Y are vectors of length n, and alpha is a real number.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(out) :: Y
      real(dp), dimension(n), intent(in) :: X
!
      real(dp), intent(in) :: alpha
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, n
!
         Y(i) = alpha*X(i)
!
      enddo
!$omp end parallel do
!
   end subroutine copy_and_scale
!
!
   subroutine copy_and_scale_complex(alpha, X, Y, n)
!!
!!    Copy and scale
!!    Written by Sarai D. Folkestad, May 2019
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Y = alpha*X
!!
      implicit none
!
      integer, intent(in) :: n
!
      complex(dp), dimension(n), intent(out) :: Y
      complex(dp), dimension(n), intent(in) :: X
!
      complex(dp), intent(in) :: alpha
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, n
!
         Y(i) = alpha*X(i)
!
      enddo
!$omp end parallel do
!
   end subroutine copy_and_scale_complex
!
!
   subroutine copy_integer(X, Y, n, alpha)
!!
!!    Copy integer
!!    Written by Alexander C. Paul, Feb 2021
!!
!!    Sets Y as:
!!
!!       Y = X
!!
!!    X and Y are vectors of length n
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer, dimension(n), intent(out) :: Y
      integer, dimension(n), intent(in) :: X
!
      integer, optional, intent(in) :: alpha
      integer :: alpha_
!
      integer :: i
!
      alpha_ = 1
      if (present(alpha)) alpha_ = alpha
!
!$omp parallel do private(i)
         do i = 1, n
!
            Y(i) = alpha_*X(i)
!
         enddo
!$omp end parallel do
!
   end subroutine copy_integer
!
!
   subroutine set_logicals(x, n, value_)
!!
!!    Set logical
!!    Written by Alexander C. Paul, Oct 2022
!!
      implicit none
!
      integer, intent(in) :: n
!
      logical, dimension(n), intent(out) :: x
      logical, intent(in) :: value_
!
      integer :: i
!
!$omp parallel do private(i) schedule(static)
      do i = 1, n
!
         x(i) = value_
!
      enddo
!$omp end parallel do
!
   end subroutine set_logicals
!
!
end module array_initialization
