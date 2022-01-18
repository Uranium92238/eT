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
      class(cc3), intent(inout) :: wf
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
      class(cc3), intent(inout) :: wf
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
   module subroutine initialize_density_intermediates_cc3(wf)
!!
!!    Initialize density intermediates
!!    Written by Alexander C. Paul and Sarai D. Folkestad, Apr 2020
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      call mem%alloc(wf%r0, wf%n_singlet_states)
!
!     CC3 corrections to the left transition density are needed
!     for the excited state transition densities
      call mem%alloc(wf%L_cc3_density_oo, wf%n_o, wf%n_o)
      call mem%alloc(wf%L_cc3_density_vv, wf%n_v, wf%n_v)
!
   end subroutine initialize_density_intermediates_cc3
!
!
   module subroutine destruct_density_intermediates_cc3(wf)
!!
!!    Destruct density intermediates
!!    Written by Alexander C. Paul and Sarai D. Folkestad, Apr 2020
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      call mem%dealloc(wf%r0, wf%n_singlet_states)
!
!     CC3 corrections to the left transition density are needed
!     for the excited state transition densities
      call mem%dealloc(wf%L_cc3_density_oo, wf%n_o, wf%n_o)
      call mem%dealloc(wf%L_cc3_density_vv, wf%n_v, wf%n_v)
!
   end subroutine destruct_density_intermediates_cc3
!
!
   module subroutine initialize_eri_c1_cc3(wf)
!!
!!    Initialize eri c1
!!    Written by Alexander C. Paul, Jan 2022
!!
      use eri_1idx_transformed_tool_class, only: eri_1idx_transformed_tool
      use array_utilities, only: zero_array
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: c1
!
      if (.not. allocated(wf%eri_c1)) then
!
         wf%L_c1 = eri_cholesky_disk('C1')
         call wf%L_c1%initialize(wf%L_t1%n_J, 2, [wf%n_o, wf%n_v])
!
         call mem%alloc(c1, wf%n_v, wf%n_o)
         call zero_array(c1, wf%n_t1)
!
         call wf%construct_c1_cholesky(c1, wf%L_t1, wf%L_c1)
         wf%eri_c1 = eri_adapter(eri_1idx_transformed_tool(wf%L_t1, wf%L_c1), wf%n_o, wf%n_v)
!
         call wf%L_c1%add_observer('eri C1', wf%eri_c1%eri)
!
         call mem%dealloc(c1, wf%n_v, wf%n_o)
!
      end if
!
   end subroutine initialize_eri_c1_cc3
!
!
   module subroutine destruct_eri_c1_cc3(wf)
!!
!!    Destruct eri c1
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      if (allocated(wf%eri_c1)) then
!
         call wf%L_c1%remove_observer('eri C1')
         deallocate(wf%L_c1)
         deallocate(wf%eri_c1)
!
      end if
!
   end subroutine destruct_eri_c1_cc3
!
!
end submodule initialize_destruct_cc3
