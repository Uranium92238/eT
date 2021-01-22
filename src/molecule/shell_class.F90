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
module shell_class
!
!!
!!    Shell class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad and Andreas Skeidsvoll, 2019
!!
!
   use kinds
!
   use interval_class, only: interval
   use global_out, only : output
   use memory_manager_class, only : mem
!
   implicit none
!
   type, extends(interval) :: shell ! interval: AO index range of the shell 
!
      integer    :: size_cart = -1 ! The number of basis functions in cartesian
      integer    :: l         = -1 ! The angular momentum
      integer    :: number_   = -1 ! The shell number (according to the ordering given by Libint)
!
      integer, private :: n_primitives = 0 ! Number of primitive Gaussians
!
      real(dp), dimension(:), allocatable, private :: exponents ! Exponents for primitives
      real(dp), dimension(:), allocatable, private :: coefficients ! Coefficients for primitives
!
   contains
!
      procedure :: determine_angular_momentum         => determine_angular_momentum_shell
      procedure :: determine_last_ao_index            => determine_last_ao_index_shell
!
      procedure :: initialize_exponents               => initialize_exponents_shell
      procedure :: initialize_coefficients            => initialize_coefficients_shell
!
      procedure :: destruct_exponents                 => destruct_exponents_shell
      procedure :: destruct_coefficients              => destruct_coefficients_shell
!
      procedure :: set_exponent_i                     => set_exponent_i_shell
      procedure :: get_exponent_i                     => get_exponent_i_shell
!
      procedure :: set_coefficient_i                  => set_coefficient_i_shell
      procedure :: get_coefficient_i                  => get_coefficient_i_shell
!
      procedure :: set_n_primitives                   => set_n_primitives_shell
      procedure :: get_n_primitives                   => get_n_primitives_shell
!
      procedure, nopass :: get_angular_momentum_label => get_angular_momentum_label_shell
!
      procedure :: cleanup                            => cleanup_shell 
!
   end type shell
!
!
   interface shell 
!
      procedure :: new_shell 
!
   end interface shell 
!
!
contains
!
!
   function new_shell(first, length, number_) result(sh)
!!
!!    New shell 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    first:   the first AO index of the shell 
!!    length:  the number of AOs in the shell 
!!    number_: the shell number in the full list of shells (according to Libint)
!!
      implicit none 
!
      integer, intent(in) :: first, length 
!
      integer, intent(in) :: number_ 
!
      type(shell) :: sh 
!
      sh%first = first 
      sh%length = length 
!
      call sh%determine_last_ao_index()
      call sh%determine_angular_momentum()
!
      sh%number_ = number_
!
   end function new_shell
!
!
   subroutine determine_angular_momentum_shell(sh)
!!
!!    Determine angular momentum
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the angular momentum by counting the number of basis
!!    functions in the shell (l < 10)
!!
      implicit none
!
      class(shell) :: sh
!
      integer :: i
!
      i = 0
      sh%l = -1
!
      do while (i .lt. 10)
!
         if ((2*i + 1) .eq. sh%length .or. (((i+1)*(i+2))/2) .eq. sh%length) then
!
            sh%l = i
!
         endif
!
         i = i + 1
!
      enddo
!
      if (sh%l .eq. -1) then
!
         call output%error_msg('could not determine angular momentum of shell.')
!
      endif
!
   end subroutine determine_angular_momentum_shell
!
!
   subroutine determine_last_ao_index_shell(sh)
!!
!!    Determine last AO index
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(shell) :: sh
!
      sh%last = sh%first + sh%length - 1
!
   end subroutine determine_last_ao_index_shell
!
!
   subroutine initialize_exponents_shell(sh)
!!
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      if (sh%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before alloc of exponents')
!
      if (.not. allocated(sh%exponents)) &
         call mem%alloc(sh%exponents, sh%n_primitives)
!
   end subroutine initialize_exponents_shell
!
!
   subroutine initialize_coefficients_shell(sh)
!!
!!    Initialize coefficients
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      if (sh%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before alloc of coefficients')
!
      if (.not. allocated(sh%coefficients)) &
            call mem%alloc(sh%coefficients, sh%n_primitives)
!
   end subroutine initialize_coefficients_shell
!
!
   subroutine destruct_exponents_shell(sh)
!!
!!    Destruct primitive exponents
!!    Written by Andreas Skeidsvoll, Aug 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      if (sh%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before dealloc of exponents')
!
      if (allocated(sh%exponents)) &
         call mem%dealloc(sh%exponents, sh%n_primitives)
!
   end subroutine destruct_exponents_shell
!
!
   subroutine destruct_coefficients_shell(sh)
!!
!!    Destruct coefficients
!!    Written by Andreas Skeidsvoll, Aug 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      if (sh%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before dealloc of coefficients')
!
      if (allocated(sh%coefficients)) &
         call mem%dealloc(sh%coefficients, sh%n_primitives)
!
   end subroutine destruct_coefficients_shell
!
!
   subroutine set_exponent_i_shell(sh, i, exponent)
!!
!!    Set exponent i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      integer, intent(in)   :: i
      real(dp), intent(in)       :: exponent
!
      if (i .gt. sh%n_primitives) &
         call output%error_msg('Tried to set exponent for non-exisiting primitive Gaussian')
!
      sh%exponents(i) = exponent
!
   end subroutine set_exponent_i_shell
!
!
   real(dp) function get_exponent_i_shell(sh, i)
!!
!!    Get exponent i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(in) :: sh
!
      integer, intent(in)   :: i
!
      if (i .gt. sh%n_primitives) &
         call output%error_msg('Tried to get exponent for non-exisiting primitive Gaussian')
!
      get_exponent_i_shell = sh%exponents(i)
!
   end function get_exponent_i_shell
!
!
   subroutine set_coefficient_i_shell(sh, i, coefficient)
!!
!!    Set coefficient i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      integer, intent(in)   :: i
      real(dp), intent(in)  :: coefficient
!
      if (i .gt. sh%n_primitives) &
         call output%error_msg('Tried to set coefficient for non-exisiting primitive Gaussian')
!
      sh%coefficients(i) = coefficient
!
   end subroutine set_coefficient_i_shell
!
!
   real(dp) function get_coefficient_i_shell(sh, i)
!!
!!    Get coefficient i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(in) :: sh
!
      integer, intent(in)   :: i
!
      if (i .gt. sh%n_primitives) &
         call output%error_msg('Tried to get coefficient for non-exisiting primitive Gaussian')
!
      get_coefficient_i_shell = sh%coefficients(i)
!
   end function get_coefficient_i_shell
!
!
   subroutine set_n_primitives_shell(sh, n)
!!
!!    Set number of primitives
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(inout) :: sh
!
      integer, intent(in) :: n
!
      sh%n_primitives = n
!
   end subroutine set_n_primitives_shell
!
!
   integer function get_n_primitives_shell(sh)
!!
!!    Get number of primitives
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell), intent(in) :: sh
!
      get_n_primitives_shell = sh%n_primitives
!
   end function get_n_primitives_shell
!
!
   subroutine get_angular_momentum_label_shell(l, ao, label, cartesian)
!!
!!    Get angular momentum label
!!    written by Alexander C. Paul, Dec 2019
!!
!!    Returns string containing the angular momentum and it's spatial component
!!
!!    l:          angular momentum quantum number
!!    ao:         atomic orbital of the shell
!!                NB: assuming default ordering of libint.
!!                cartesian basis functions: {xx, xy, xz, yy, yz, zz}
!!                spherical/pure functions:  m_l: {2, 1, 0, -1, -2}
!!    label:      string that is returned e.g. d_xx
!!    cartesian:  logical determining if a cartesian or "pure" basis set is used
!!
      use angular_momentum
!
      implicit none
!
      integer, intent(in) :: l
      integer, intent(in) :: ao
!
      character(len=*), intent(inout) :: label
!
      logical, intent(in) :: cartesian
!
      if (cartesian) then
!
         select case (l)
            case(0)
               label = s
            case(1)
               label = p_cart(ao)
            case(2)
               label = d_cart(ao)
            case(3)
               label = f_cart(ao)
            case(4)
               label = g
            case(5)
               label = h
            case default
               call output%error_msg('Angular momentum of atomic orbital not recognized.')
         end select
!
      else
!
         select case (l)
            case(0)
               label = s
            case(1)
               label = p(ao)
            case(2)
               label = d(ao)
            case(3)
               label = f(ao)
            case(4)
               label = g
            case(5)
               label = h
            case default
               call output%error_msg('Angular momentum of atomic orbital not recognized.')
         end select
!
      end if
!
   end subroutine get_angular_momentum_label_shell
!
!
   subroutine cleanup_shell(sh)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Dec 2019
!!
      implicit none
!
      class(shell) :: sh
!
      call sh%destruct_exponents()   
      call sh%destruct_coefficients()
!
   end subroutine cleanup_shell
!
!
end module shell_class
