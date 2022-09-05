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
module abstract_spherical_angular_momentum_class
!
!!
!!    Abstract spherical angular momentum class
!!    Written by Alexander C. Paul, Feb 2022
!!
!!    Stores information about spherical angular momentum
!!
!!    Libint ordering and normalization described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use abstract_angular_momentum_class
!
!
   type, abstract, extends(abstract_angular_momentum) :: abstract_spherical_angular_momentum
!
   contains
!
      procedure :: get_angular_part &
                => get_angular_part_abstract_spherical_angular_momentum
!
   end type
!
!
contains
!
!
   function get_angular_part_abstract_spherical_angular_momentum(this, x, y ,z) &
                                                                 result(angular_part)
!!
!!    Get angular part
!!    Written by Alexander C. Paul, June 2021
!!
!!    Evaluate spherical harmonic of the AOs at x, y, z
!!
!!    Solid harmonic from Molecular electronic structure theory eqn. (6.4.47)
!!
      use array_initialization, only: zero_array
      use math_utilities, only: factorial, binomial
!
      implicit none
!
      class(abstract_spherical_angular_momentum), intent(in) :: this
!
      real(dp), intent(in) :: x, y, z
!
      real(dp), dimension(this%n_functions) :: angular_part
!
      real(dp) :: C_lm_tuv, N_S_lm
      integer :: i, ml, abs_ml, floored_term
      integer :: two_v, two_v_m, two_v_max, u, t, t_max
!
      call zero_array(angular_part, this%n_functions)
!
      do i = 1, this%n_functions
!
         ml = -this%l - 1 + i
!
         abs_ml = abs(ml)
!
!        floors (l - abs_ml)/2, since abs_ml <= the non-negative l
         t_max = int((this%l - abs_ml)/2)
!
         if (ml .ge. 0) then
            two_v_m = 0
         else
            two_v_m = 1
         endif
!
!        floors (abs_ml - two_v_m)/2, since 2*v_m <= the non-negative abs_ml
         floored_term = int((abs_ml - two_v_m)/2)
!
         two_v_max = 2*floored_term + two_v_m
!
         do t = 0, t_max
            do u = 0, t
               do two_v = two_v_m, two_v_max, 2
!
!                 C_lm_tuv from Molecular electronic structure theory eqn.
!                 (6.4.48), where (two_v - two_v_m)/2 is an integer
!
                  C_lm_tuv = (-one)**(t + (two_v - two_v_m)/2) * quarter**t   &
                             * real(binomial(this%l, t) *  &
                                    binomial(this%l-t, abs_ml+t) * &
                                    binomial(t, u) * binomial(abs_ml, two_v), kind=dp)
!
                  angular_part(i) = angular_part(i)                  &
                        + C_lm_tuv * x**(2*t + abs_ml - 2*u -two_v)  &
                        * y**(2*u + two_v) * z**(this%l - 2*t - abs_ml)

!
               enddo
            enddo
         enddo
!
!        Multiply by N_S_lm from Molecular electronic structure theory eqn. (6.4.49)
!
         N_S_lm = half**abs_ml * real(factorial(this%l), kind=dp) &
                * sqrt(real(2*factorial(this%l + abs_ml) &
                            * factorial(this%l - abs_ml), kind=dp))
!
         if (ml == 0) N_S_lm = N_S_lm/sqrt(two)
!
         angular_part(i) = N_S_lm * angular_part(i)
!
      enddo
!
   end function get_angular_part_abstract_spherical_angular_momentum
!
!
end module abstract_spherical_angular_momentum_class
