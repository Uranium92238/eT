!
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
module shell_class
!
!!
!!    Shell class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad and Andreas Skeidsvoll, 2019
!!
!
   use parameters
!
   use range_class, only: range_
   use global_out, only: output
   use memory_manager_class, only: mem
   use abstract_angular_momentum_class, only: abstract_angular_momentum
!
   implicit none
!
   type, extends(range_) :: shell ! AO index range of the shell
!
      integer :: number_ = -1 ! The shell number (according to the ordering given by Libint)
!
      integer, private :: n_primitives = 0 ! Number of primitive Gaussians
!
      real(dp), dimension(:), allocatable, private :: exponents ! Exponents for primitives
      real(dp), dimension(:), allocatable, private :: coefficients ! Coefficients for primitives
!
      class(abstract_angular_momentum), allocatable :: angular_momentum
!
   contains
!
      procedure :: initialize_exponents &
                => initialize_exponents_shell
      procedure :: initialize_coefficients &
                => initialize_coefficients_shell
!
      procedure :: set_exponent_i &
                => set_exponent_i_shell
      procedure :: get_exponent_i &
                => get_exponent_i_shell
!
      procedure :: set_coefficient_i &
                => set_coefficient_i_shell
      procedure :: get_coefficient_i &
                => get_coefficient_i_shell
!
      procedure :: set_n_primitives &
                => set_n_primitives_shell
      procedure :: get_n_primitives &
                => get_n_primitives_shell
!
      procedure :: get_angular_momentum &
                => get_angular_momentum_shell
!
      procedure :: get_angular_momentum_letter &
                => get_angular_momentum_letter_shell
!
      procedure :: get_angular_momentum_label &
                => get_angular_momentum_label_shell
!
      procedure :: get_molden_order &
                => get_molden_order_shell
!
      procedure :: get_normalization_factor &
                => get_normalization_factor_shell
!
      procedure :: get_overlap_of_primitives &
                => get_overlap_of_primitives_shell
!
      procedure :: get_radial_part &
                => get_radial_part_shell
!
      procedure :: get_aos_at_point &
                => get_aos_at_point_shell
!
      procedure, private :: destruct_exponents &
                         => destruct_exponents_shell
      procedure, private :: destruct_coefficients &
                         => destruct_coefficients_shell
!
      procedure :: cleanup &
                => cleanup_shell
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
   function new_shell(first, length, number_, cartesian) result(this)
!!
!!    New shell
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    first:   the first AO index of the shell
!!    length:  the number of AOs in the shell
!!    number_: the shell number in the full list of shells (according to Libint)
!!
      use angular_momentum_factory_class, only: angular_momentum_factory
!
      implicit none
!
      integer, intent(in) :: first, length
!
      integer, intent(in) :: number_
!
      logical, intent(in) :: cartesian
!
      type(shell) :: this
      type(angular_momentum_factory) :: factory
!
      this%range_ = range_(first, length)
!
      this%number_ = number_
!
      factory = angular_momentum_factory(cartesian)
      call factory%create(this%angular_momentum, length)
!
   end function new_shell
!
!
   subroutine initialize_exponents_shell(this)
!!
!!    Initialize exponents
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
   function get_angular_momentum_label_shell(this, n) result(label)
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
!!    label:      string that is returned e.g. d_xx
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer, intent(in) :: n
!
      character(len=:), allocatable :: label
!
      label = this%angular_momentum%get_label(n)
!
   end function get_angular_momentum_label_shell
!
!
   function get_molden_order_shell(this) result(map)
!!
!!    Get Molden order
!!    Written by Alexander C. Paul, May 2021
!!
!!    Return index list mapping AOs to the order molden expects
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer, dimension(this%length) :: map
!
      integer :: i
!
      do i = 1, this%length
         map(i) = this%first-1 + this%angular_momentum%get_molden_offset(i)
      end do
!
   end function get_molden_order_shell
!
!
   function get_normalization_factor_shell(this) result(factors)
!!
!!    Get normalization factor
!!    Written by Alexander C. Paul, May 2021
!!
!!    Return factor to normalize the AO determined
!!    by the angular momentum l and the basis function i
!!
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
         factors(i) = this%angular_momentum%get_normalization_factor(i)
!
      end do
!
   end function get_normalization_factor_shell
!
!
   pure function get_angular_momentum_letter_shell(this) result(l_letter)
!!
!!    Get angular momentum
!!    written by Alexander C. Paul, Dec 2019
!!
!!    Convert angular momentum quantum number into the letter s,p,d,f,g,h
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      character(len=:), allocatable :: l_letter
!
      l_letter = this%angular_momentum%get_l_string()
!
   end function get_angular_momentum_letter_shell
!
!
   pure function get_angular_momentum_shell(this) result(l)
!!
!!    Get angular momentum
!!    written by Alexander C. Paul, Dec 2019
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      integer :: l
!
      l = this%angular_momentum%get_l_integer()
!
   end function get_angular_momentum_shell
!
!
   function get_overlap_of_primitives_shell(this) result(overlap)
!!
!!    Get overlap of primitives
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
      implicit none
!
      class(shell), intent(in) :: this
!
      real(dp) :: overlap, power, base
!
      integer :: i, j
!
      overlap = zero
!
      power = three*half + real(this%get_angular_momentum(), kind=dp)
!
      do i = 1, this%n_primitives
         do j = 1, this%n_primitives
!
            base = sqrt(this%exponents(i)*this%exponents(j)) &
                  / (this%exponents(i) + this%exponents(j))
!
            overlap = overlap + this%coefficients(i)*this%coefficients(j) * (base**power)
!
         enddo
      enddo
!
      overlap = overlap*(two**power)
!
   end function get_overlap_of_primitives_shell
!
!
   function get_radial_part_shell(this, r_squared) result(radial_part)
!!
!!    Get radial part
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
      use math_utilities, only: double_factorial
!
      implicit none
!
      class(shell), intent(in) :: this
!
      real(dp), intent(in) :: r_squared
!
      real(dp) :: radial_part
      real(dp) :: overlap_primitives, normalization_constant, exponent, coefficient
      integer  :: l, i
!
      radial_part = zero
!
      overlap_primitives = this%get_overlap_of_primitives()
!
      l = this%get_angular_momentum()
!
      do i = 1, this%n_primitives
!
         exponent    = this%exponents(i)
         coefficient = this%coefficients(i)
!
!        Normalization constant for primitive containing a Racah's normalized angular part
!        (from eqn. (94) without Racah's normalization constant and (95)
!        in Giesea, T. J. HSERILib: Gaussian integral evaluation)
!
         normalization_constant = (four*exponent)**(real(l, dp)*half + three*quarter)
!
         radial_part = radial_part &
                       + normalization_constant*coefficient*exp(-exponent*r_squared)
!
      enddo
!
      normalization_constant = sqrt((two*pi)**(three*half) * overlap_primitives &
                                    * real(double_factorial(2*l-1), dp))
!
      radial_part = radial_part/normalization_constant
!
   end function get_radial_part_shell
!
!
   function get_aos_at_point_shell(this, x, y, z) result(aos_at_point)
!!
!!    Get AOs at point
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
!!    Construct Racah's normalized orbitals for a point
!!    x,y,z relative to the nucleus
!!
      use array_initialization, only: copy_and_scale
!
      implicit none
!
      class(shell), intent(in) :: this
!
      real(dp), intent(in) :: x, y, z
!
      real(dp), dimension(this%length) :: aos_at_point
!
      real(dp) :: r_squared, radial_part
!
      r_squared = x**2 + y**2 + z**2
!
      radial_part = this%get_radial_part(r_squared)
!
      call copy_and_scale(radial_part, &
                          this%angular_momentum%get_angular_part(x, y, z), &
                          aos_at_point, &
                          this%length)
!
   end function get_aos_at_point_shell
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
