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
submodule (ccs_class) zop_ccs
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
   module subroutine initialize_gs_density_ccs(wf)
!!
!!    Initialize density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%density)) call mem%alloc(wf%density, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_gs_density_ccs
!
!
   module subroutine destruct_gs_density_ccs(wf)
!!
!!    Destruct density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%density)) call mem%dealloc(wf%density, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_gs_density_ccs
!
!
   module subroutine prepare_for_density_ccs(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
!     For now, do nothing.
!
      write(output%unit,'(/t3,a,a,a)') 'No preparations for the density for ', &
                                       trim(wf%name_), ' wavefunction.'
!
   end subroutine prepare_for_density_ccs
!
!
   module subroutine construct_gs_density_ccs(wf)
!!
!!    Construct one-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    Constructs the one-electron density 
!!    matrix in the T1 basis
!!
!!    D_pq = < Λ | E_pq | CC >
!!
      implicit none
!
      class(ccs) :: wf
!
      call zero_array(wf%density, (wf%n_mo)**2)
!
      call wf%gs_one_el_density_ccs_oo(wf%density)
      call wf%gs_one_el_density_ccs_vo(wf%density, wf%t1bar)
!
   end subroutine construct_gs_density_ccs
!
!
   module subroutine gs_one_el_density_ccs_oo_ccs(wf, density)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ii = 2  
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, wf%n_o
!
         density(i,i) = density(i,i) + two  
!
      enddo
!$omp end parallel do
!
   end subroutine gs_one_el_density_ccs_oo_ccs
!
!
   module subroutine gs_one_el_density_ccs_vo_ccs(wf, density, tbar_ai)
!!
!!    One electron density vo
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ai = tbar_ai 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      real(dp), dimension(wf%n_v, wf%n_o) :: tbar_ai
!
      integer :: i, a
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!        
            density(wf%n_o + a, i) = density(wf%n_o + a, i) + tbar_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine gs_one_el_density_ccs_vo_ccs
!
!
   module function calculate_expectation_value_ccs(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Calculate the expectation value of a one-electron
!!    operator Â
!!
!!       < A > = < Λ | Â | CC > = sum_pq A_pq D_pq
!!
!!    where A_pq are the T1-transformed integrals
!!    and D_pq is the a one-electron density matrix
!!    in the T1-basis
!!
      implicit none
!  
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
!
      real(dp) :: expectation_value
!
      real(dp) :: ddot
!
      expectation_value = ddot(wf%n_mo**2, A, 1, density, 1)
!
   end function calculate_expectation_value_ccs
!
!
end submodule zop_ccs