!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
submodule (cc3_class) initialize_destruct_cc3
!
!!
!!    Initialize destruct submodule
!!
!!    Gathers routines that initialize and destruct the CC3 type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_gs_density_cc3(wf)
!!
!!    Initialize density and CC3 corrections to the GS-density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      call mem%alloc(wf%density, wf%n_mo, wf%n_mo)
!
!     CC3 corrections to the GS-density are needed for the right transition density
      call mem%alloc(wf%GS_cc3_density_oo, wf%n_o, wf%n_o)
      call mem%alloc(wf%GS_cc3_density_vv, wf%n_v, wf%n_v)
!
   end subroutine initialize_gs_density_cc3
!
!
   module subroutine destruct_gs_density_cc3(wf)
!!
!!    Destruct density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      call mem%dealloc(wf%density, wf%n_mo, wf%n_mo)
!
!     CC3 corrections to the GS-density are needed for the right transition density
      call mem%dealloc(wf%GS_cc3_density_oo, wf%n_o, wf%n_o)
      call mem%dealloc(wf%GS_cc3_density_vv, wf%n_v, wf%n_v)
!
   end subroutine destruct_gs_density_cc3
!
!
end submodule initialize_destruct_cc3
