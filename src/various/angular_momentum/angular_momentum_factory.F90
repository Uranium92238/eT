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
module angular_momentum_factory_class
!
!!
!!    Angular momentum factory class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Creates angular momentum object based on the length of a shell
!!    and the logical cartesian
!!
!
   use parameters
!
   implicit none
!
   type :: angular_momentum_factory
!
      logical, private :: cartesian
!
   contains
!
      procedure, public :: create => create_angular_momentum_factory
!
      procedure, nopass, private :: determine_cartesian_angular_momentum
      procedure, nopass, private :: determine_spherical_angular_momentum
!
   end type angular_momentum_factory
!
   interface angular_momentum_factory
!
      procedure :: new_angular_momentum_factory
!
   end interface angular_momentum_factory
!
!
contains
!
!
   function new_angular_momentum_factory(cartesian) result(this)
!!
!!    New angular momentum factory
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(angular_momentum_factory) :: this
!
      logical, intent(in) :: cartesian
!
      this%cartesian = cartesian
!
   end function new_angular_momentum_factory
!
!
   subroutine create_angular_momentum_factory(this, angular_momentum, n_ao_in_shell)
!!
!!    Create
!!    Written by Alexander C. Paul, June 2021
!!
!!    Constructs the angular momentum object
!!
      use global_out, only: output
!
      use abstract_angular_momentum_class
      use s_angular_momentum_class
      use cartesian_p_angular_momentum_class
      use cartesian_d_angular_momentum_class
      use cartesian_f_angular_momentum_class
      use cartesian_g_angular_momentum_class
      use cartesian_h_angular_momentum_class
      use cartesian_i_angular_momentum_class
      use spherical_p_angular_momentum_class
      use spherical_d_angular_momentum_class
      use spherical_f_angular_momentum_class
      use spherical_g_angular_momentum_class
      use spherical_h_angular_momentum_class
      use spherical_i_angular_momentum_class
!
      implicit none
!
      class(angular_momentum_factory), intent(in) :: this
      class(abstract_angular_momentum), allocatable :: angular_momentum
!
      integer, intent(in) :: n_ao_in_shell
      integer :: l
!
      if (this%cartesian) then
!
         l = this%determine_cartesian_angular_momentum(n_ao_in_shell)
!
         select case (l)
            case(0)
               angular_momentum = s_angular_momentum()
            case(1)
               angular_momentum = cartesian_p_angular_momentum()
            case(2)
               angular_momentum = cartesian_d_angular_momentum()
            case(3)
               angular_momentum = cartesian_f_angular_momentum()
            case(4)
               angular_momentum = cartesian_g_angular_momentum()
            case(5)
               angular_momentum = cartesian_h_angular_momentum()
            case(6)
               angular_momentum = cartesian_i_angular_momentum()
            case default
               call output%error_msg('Angular momentum beyond l = 6 not supported.')
         end select
!
      else
!
         l = this%determine_spherical_angular_momentum(n_ao_in_shell)
!
         select case (l)
            case(0)
               angular_momentum = s_angular_momentum()
            case(1)
               angular_momentum = spherical_p_angular_momentum()
            case(2)
               angular_momentum = spherical_d_angular_momentum()
            case(3)
               angular_momentum = spherical_f_angular_momentum()
            case(4)
               angular_momentum = spherical_g_angular_momentum()
            case(5)
               angular_momentum = spherical_h_angular_momentum()
            case(6)
               angular_momentum = spherical_i_angular_momentum()
            case default
               call output%error_msg('Angular momentum beyond l = 6 not supported.')
         end select
!
      end if
!
   end subroutine create_angular_momentum_factory
!
!
   pure function determine_cartesian_angular_momentum(n) result(l)
!!
!!    Determine cartesian angular momentum
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the angular momentum from the number of basis
!!    functions in the shell: n = (l+1)*(l+2)/2
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer :: l
!
      l = int(-three*half + sqrt(nine*quarter + two*(n-1)) )
!
   end function determine_cartesian_angular_momentum
!
!
   pure function determine_spherical_angular_momentum(n) result(l)
!!
!!    Determine spherical angular momentum
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the angular momentum from the number of basis
!!    functions in the shell: n = 2l + 1
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer :: l
!
      l = int((n - 1)/2)
!
   end function determine_spherical_angular_momentum
!
!
end module angular_momentum_factory_class
