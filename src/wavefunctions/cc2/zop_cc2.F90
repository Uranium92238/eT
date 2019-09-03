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
submodule (cc2_class) zop_cc2
!
!!
!!    Zeroth order properties submodule 
!!
!!    Contains routines related to the mean values, i.e. 
!!    the construction of density matrices as well as expectation 
!!    value calculation.
!!
!
   implicit none 
!
!
contains
!
!
   module subroutine prepare_for_density_cc2(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Sarai D. Folekstad, May 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      call wf%initialize_t2()
      call wf%initialize_t2bar()
!
      call wf%construct_t2()
      call wf%construct_t2bar()
!
   end subroutine prepare_for_density_cc2
!
!
end submodule zop_cc2