   module shell_details_class
!
!!
!!    AO details class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019
!!
!
   use kinds
   use output_file_class
   use disk_manager_class
   use memory_manager_class
!
   implicit none
!
   type :: shell_details
!
      integer :: n_primitives = 0 ! Number of primitive Gaussians
!
      real(dp), dimension(:), allocatable :: exponents
      real(dp), dimension(:), allocatable :: coefficients
!
   contains
!
      procedure :: initialize_exponents      => initialize_exponents_shell_details
      procedure :: initialize_coefficients   => initialize_coefficients_shell_details
!
      procedure :: set_exponent_i            => set_exponent_i_shell_details
      procedure :: get_exponent_i            => get_exponent_i_shell_details
!
      procedure :: set_coefficient_i         => set_coefficient_i_shell_details
      procedure :: get_coefficient_i         => get_coefficient_i_shell_details
!
      procedure :: set_n_primitives          => set_n_primitives_shell_details
      procedure :: get_n_primitives          => get_n_primitives_shell_details
!
   end type shell_details
!
contains
!
!
   subroutine initialize_exponents_shell_details(ao)
!!
!!    Initialize exponents
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      if (ao%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before alloc of exponents')
!
      if (.not. allocated(ao%exponents)) allocate(ao%exponents(ao%n_primitives))
!
   end subroutine initialize_exponents_shell_details
!
!
   subroutine initialize_coefficients_shell_details(ao)
!!
!!    Initialize coefficients
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      if (ao%n_primitives == 0) &
      call output%error_msg('Number of primitive Gaussians not set before alloc of coefficients')
!
      if (.not. allocated(ao%coefficients)) allocate(ao%coefficients(ao%n_primitives))
!
   end subroutine initialize_coefficients_shell_details
!
!
   subroutine set_exponent_i_shell_details(ao, i, exponent_)
!!
!!    Set exponent
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      integer, intent(in)   :: i
      real(dp), intent(in)       :: exponent_
!
      if (i .gt. ao%n_primitives) call output%error_msg('Tried to set exponent for non-exisiting primitive Gaussian')
!
      ao%exponents(i) = exponent_
!
   end subroutine set_exponent_i_shell_details
!
!
   real(dp) function get_exponent_i_shell_details(ao, i)
!!
!!    Get exponent
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      integer, intent(in)   :: i
!
      if (i .gt. ao%n_primitives) call output%error_msg('Tried to get exponent for non-exisiting primitive Gaussian')
!
      get_exponent_i_shell_details = ao%exponents(i)
!
   end function get_exponent_i_shell_details
!
!
   subroutine set_coefficient_i_shell_details(ao, i, coefficient)
!!
!!    Set coefficient
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      integer, intent(in)   :: i
      real(dp), intent(in)       :: coefficient
!
      if (i .gt. ao%n_primitives) call output%error_msg('Tried to set coefficient for non-exisiting primitive Gaussian')
!
      ao%coefficients(i) = coefficient
!
   end subroutine set_coefficient_i_shell_details
!
!
   real(dp) function get_coefficient_i_shell_details(ao, i)
!!
!!    Get coefficient
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      integer, intent(in)   :: i
!
      if (i .gt. ao%n_primitives) call output%error_msg('Tried to get coefficient for non-exisiting primitive Gaussian')
!
      get_coefficient_i_shell_details = ao%coefficients(i)
!
   end function get_coefficient_i_shell_details
!
!
   subroutine set_n_primitives_shell_details(ao, n)
!!
!!    Set number of primitives
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(inout) :: ao
!
      integer, intent(in) :: n
!
      ao%n_primitives = n
!
   end subroutine set_n_primitives_shell_details
!
!
   integer function get_n_primitives_shell_details(ao)
!!
!!    Get number of primitives
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(shell_details), intent(in) :: ao
!
      get_n_primitives_shell_details = ao%n_primitives
!
   end function get_n_primitives_shell_details
!
!
end module shell_details_class
