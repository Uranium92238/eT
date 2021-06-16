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
   use range_class
   use global_out, only : output
   use memory_manager_class, only : mem
!
   implicit none
!
   type, extends(range_) :: shell ! AO index range of the shell
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
      procedure :: determine_angular_momentum => determine_angular_momentum_shell
!
      procedure :: initialize_exponents       => initialize_exponents_shell
      procedure :: initialize_coefficients    => initialize_coefficients_shell
!
      procedure :: destruct_exponents         => destruct_exponents_shell
      procedure :: destruct_coefficients      => destruct_coefficients_shell
!
      procedure :: set_exponent_i             => set_exponent_i_shell
      procedure :: get_exponent_i             => get_exponent_i_shell
!
      procedure :: set_coefficient_i          => set_coefficient_i_shell
      procedure :: get_coefficient_i          => get_coefficient_i_shell
!
      procedure :: set_n_primitives           => set_n_primitives_shell
      procedure :: get_n_primitives           => get_n_primitives_shell
!
      procedure :: get_angular_momentum       => get_angular_momentum_shell
      procedure :: get_angular_momentum_label => get_angular_momentum_label_shell
      procedure :: get_molden_order        => get_molden_order_shell
!
      procedure :: get_cartesian_normalization_factor &
                => get_cartesian_normalization_factor_shell
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
   function new_shell(first, length, number_) result(this)
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
      type(shell) :: this
!
      this%range_ = range_(first, length)
!
      call this%determine_angular_momentum()
!
      this%number_ = number_
!
   end function new_shell
!
!
   subroutine determine_angular_momentum_shell(this)
!!
!!    Determine angular momentum
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the angular momentum by counting the number of basis
!!    functions in the shell (l < 10)
!!
      implicit none
!
      class(shell) :: this
!
      integer :: i
!
      i = 0
      this%l = -1
!
      do while (i .lt. 10)
!
         if ((2*i + 1) .eq. this%length .or. &
            (((i+1)*(i+2))/2) .eq. this%length) then
!
            this%l = i
!
         endif
!
         i = i + 1
!
      enddo
!
      if (this%l .eq. -1) then
!
         call output%error_msg('could not determine angular momentum of shell.')
!
      endif
!
   end subroutine determine_angular_momentum_shell
!
!
   subroutine initialize_exponents_shell(this)
!!
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      if (this%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before alloc of exponents')
!
      if (.not. allocated(this%exponents)) &
         call mem%alloc(this%exponents, this%n_primitives)
!
   end subroutine initialize_exponents_shell
!
!
   subroutine initialize_coefficients_shell(this)
!!
!!    Initialize coefficients
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      if (this%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before alloc of coefficients')
!
      if (.not. allocated(this%coefficients)) &
            call mem%alloc(this%coefficients, this%n_primitives)
!
   end subroutine initialize_coefficients_shell
!
!
   subroutine destruct_exponents_shell(this)
!!
!!    Destruct primitive exponents
!!    Written by Andreas Skeidsvoll, Aug 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      if (this%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before dealloc of exponents')
!
      if (allocated(this%exponents)) &
         call mem%dealloc(this%exponents, this%n_primitives)
!
   end subroutine destruct_exponents_shell
!
!
   subroutine destruct_coefficients_shell(this)
!!
!!    Destruct coefficients
!!    Written by Andreas Skeidsvoll, Aug 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      if (this%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before dealloc of coefficients')
!
      if (allocated(this%coefficients)) &
         call mem%dealloc(this%coefficients, this%n_primitives)
!
   end subroutine destruct_coefficients_shell
!
!
   subroutine set_exponent_i_shell(this, i, exponent)
!!
!!    Set exponent i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      integer, intent(in)   :: i
      real(dp), intent(in)       :: exponent
!
      if (i .gt. this%n_primitives) &
         call output%error_msg('Tried to set exponent for non-exisiting primitive Gaussian')
!
      this%exponents(i) = exponent
!
   end subroutine set_exponent_i_shell
!
!
   function get_exponent_i_shell(this, i) result(exponent_)
!!
!!    Get exponent i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer, intent(in)   :: i
!
      real(dp) :: exponent_
!
      if (i .gt. this%n_primitives) &
         call output%error_msg('Tried to get exponent for non-exisiting primitive Gaussian')
!
      exponent_ = this%exponents(i)
!
   end function get_exponent_i_shell
!
!
   subroutine set_coefficient_i_shell(this, i, coefficient)
!!
!!    Set coefficient i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      integer, intent(in)   :: i
      real(dp), intent(in)  :: coefficient
!
      if (i .gt. this%n_primitives) &
         call output%error_msg('Tried to set coefficient for non-exisiting primitive Gaussian')
!
      this%coefficients(i) = coefficient
!
   end subroutine set_coefficient_i_shell
!
!
   real(dp) function get_coefficient_i_shell(this, i)
!!
!!    Get coefficient i
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer, intent(in) :: i
!
      if (i .gt. this%n_primitives) &
         call output%error_msg('Tried to get coefficient for non-exisiting primitive Gaussian')
!
      get_coefficient_i_shell = this%coefficients(i)
!
   end function get_coefficient_i_shell
!
!
   subroutine set_n_primitives_shell(this, n)
!!
!!    Set number of primitives
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(inout) :: this
!
      integer, intent(in) :: n
!
      this%n_primitives = n
!
   end subroutine set_n_primitives_shell
!
!
   pure function get_n_primitives_shell(this) result(n)
!!
!!    Get number of primitives
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer :: n
!
      n = this%n_primitives
!
   end function get_n_primitives_shell
!
!
   function get_angular_momentum_label_shell(this, n, cartesian) result(label)
!!
!!    Get angular momentum label
!!    written by Alexander C. Paul, Dec 2019
!!
!!    Returns string containing the angular momentum and it's spatial component
!!
!!    n:          number of the orbital of the shell
!!                NB: assuming default ordering of libint.
!!                cartesian basis functions: {xx, xy, xz, yy, yz, zz}
!!                spherical/pure functions:  m_l: {2, 1, 0, -1, -2}
!!    cartesian:  logical determining if a cartesian or "pure" basis set is used
!!    label:      string that is returned e.g. d_xx
!!
      use angular_momentum
!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer, intent(in) :: n
!
      character(len=8) :: label
!
      logical, intent(in) :: cartesian
!
      if (cartesian) then
!
         select case (this%l)
            case(0)
               label = this%get_angular_momentum()
            case(1)
               label = this%get_angular_momentum() // ' '  // p_cart(n)
            case(2)
               label = this%get_angular_momentum() // ' '  // d_cart(n)
            case(3)
               label = this%get_angular_momentum() // ' '  // f_cart(n)
            case(4)
               label = this%get_angular_momentum() // ' '  // g_cart(n)
            case(5)
               label = this%get_angular_momentum()
            case default
               call output%error_msg('Angular momentum of atomic orbital not recognized.')
         end select
!
      else
!
         select case (this%l)
            case(0)
               label = this%get_angular_momentum()
            case(1)
               label = this%get_angular_momentum() // ' ' // p(n)
            case(2)
               label = this%get_angular_momentum() // ' '  // d(n)
            case(3)
               label = this%get_angular_momentum() // ' '  // f(n)
            case(4)
               label = this%get_angular_momentum() // ' '  // g(n)
            case(5)
               label = this%get_angular_momentum()
            case default
               call output%error_msg('Angular momentum of atomic orbital not recognized.')
         end select
!
      end if
!
   end function get_angular_momentum_label_shell
!
!
   function get_molden_order_shell(this, cartesian) result(map)
!!
!!    Get Molden order
!!    Written by Alexander C. Paul, May 2021
!!
!!    Return index list mapping AOs to the order molden expects
!!
      use angular_momentum
!
      implicit none
!
      class(shell), intent(in) :: this
!
      logical, intent(in) :: cartesian
!
      integer, dimension(this%length) :: map
!
      integer :: i
!
      do i = 1, this%length
!
         if (cartesian) then
!
            select case (this%l)
               case(0)
                  map(i) = i + this%first - 1
               case(1)
                  map(i) = i + this%first - 1
               case(2)
                  map(i) = this%first - 1 + d_offsets_cart(i)
               case(3)
                  map(i) = this%first - 1 + f_offsets_cart(i)
               case(4)
                  map(i) = this%first - 1 + g_offsets_cart(i)
               case default
                  call output%error_msg('Molden cannot handle angular momentum beyond l=4.')
            end select
!
         else
!
            select case (this%l)
               case(0)
                  map(i) = i + this%first - 1
               case(1)
                  map(i) = i + this%first - 1
               case(2)
                  map(i) = this%first - 1 + d_offsets(i)
               case(3)
                  map(i) = this%first - 1 + f_offsets(i)
               case(4)
                  map(i) = this%first - 1 + g_offsets(i)
               case default
                  call output%error_msg('Molden cannot handle angular momentum beyond l=4.')
            end select
!
         end if
      end do
!
   end function get_molden_order_shell
!
!
   function get_cartesian_normalization_factor_shell(this) result(factors)
!!
!!    Get cartesian normalization factor
!!    Written by Alexander C. Paul, May 2021
!!
!!    Return factor to normalize the cartesian function determined
!!    by the angular momentum l and the basis function i
!!
      use angular_momentum
!
      implicit none
!
      class(shell), intent(in) :: this
!
      real(dp), dimension(this%length) :: factors
!
      integer:: i
!
      do i = 1, this%length
!
         select case (this%l)
            case(0)
               factors(i) = one
            case(1)
               factors(i) = one
            case(2)
               factors(i) = d_cart_normalization(i)
            case(3)
               factors(i) = f_cart_normalization(i)
            case(4)
               factors(i) = g_cart_normalization(i)
            case default
               call output%error_msg('Angular momentum of atomic orbital not recognized.')
         end select
!
      end do
!
   end function get_cartesian_normalization_factor_shell
!
!
   pure function get_angular_momentum_shell(this) result(l_letter)
!!
!!    Get angular momentum
!!    written by Alexander C. Paul, Dec 2019
!!
!!    Convert angular momentum quantum number into the letter s,p,d,f,g,h
!!
      use angular_momentum
!
      implicit none
!
      class(shell), intent(in) :: this
!
      character(len=:), allocatable :: l_letter
!
      select case (this%l)
         case(0)
            l_letter = 's'
         case(1)
            l_letter = 'p'
         case(2)
            l_letter = 'd'
         case(3)
            l_letter = 'f'
         case(4)
            l_letter = 'g'
         case(5)
            l_letter = 'h'
      end select
!
   end function get_angular_momentum_shell
!
!
   subroutine cleanup_shell(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Dec 2019
!!
      implicit none
!
      class(shell) :: this
!
      call this%destruct_exponents()
      call this%destruct_coefficients()
!
   end subroutine cleanup_shell
!
!
end module shell_class
