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
submodule (doubles_class) zop_doubles
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

   module subroutine construct_gs_density_doubles(wf)
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
      class(doubles) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_aibj
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj
!
      call zero_array(wf%density, (wf%n_mo)**2)
!
      call wf%density_ccs_ref_ref_oo(wf%density)
      call wf%density_ccs_mu_ref_vo(wf%density, wf%t1bar)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%gs_one_el_density_doubles_ov(wf%density, wf%t1bar, t_aibj)
!
      call mem%alloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tbar_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%gs_one_el_density_doubles_oo(wf%density, tbar_aibj, t_aibj)
      call wf%gs_one_el_density_doubles_vv(wf%density, tbar_aibj, t_aibj)
!
      call mem%dealloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_gs_density_doubles
!
!
   module subroutine gs_one_el_density_doubles_oo_doubles(wf, density, tbar_akbj, t_akbi)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ij -= sum_abk t_akb,i tbar_akb,j 
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_akbj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_akbi
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_v**2)*(wf%n_o),   &
                  -one,                   &
                  t_akbi,                 & ! t_akb_i
                  (wf%n_v**2)*(wf%n_o),   &
                  tbar_akbj,              & ! tbar_akb_j
                  (wf%n_v**2)*(wf%n_o),   &
                  one,                    &
                  density,                &
                  wf%n_mo)
!
   end subroutine gs_one_el_density_doubles_oo_doubles
!
!
   module subroutine gs_one_el_density_doubles_vv_doubles(wf, density, tbar_ajci, t_bjci)
!!
!!    One electron density vv
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ab += sum_jci tbar_a,jci t_b,jci
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_ajci
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_bjci
!
      call dgemm('N', 'T',                         &
                  wf%n_v,                          &
                  wf%n_v,                          &
                  (wf%n_o**2)*(wf%n_v),            &
                  one,                             &
                  tbar_ajci,                       & ! tbar_a_jci
                  wf%n_v,                          &
                  t_bjci,                          & ! t_b_jci
                  wf%n_v,                          &
                  one,                             &
                  density(wf%n_o + 1, wf%n_o + 1), &
                  wf%n_mo)
!
   end subroutine gs_one_el_density_doubles_vv_doubles
!
!
   module subroutine gs_one_el_density_doubles_ov_doubles(wf, density, tbar_ai, t_aibj)
!!
!!    One electron density ov
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ia += sum_bj u^{ab}_ij tbar_bj = sum_bj u_ia,bj tbar_bj 
!!
!!    u^{ab}_ij = 2t_aibj - t_ajbi
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: u_aibj
!
      real(dp), dimension(:,:), allocatable :: D_ov
!
      integer :: i, a
!
      call mem%alloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      call dcopy((wf%n_v)**2*(wf%n_o)**2, t_aibj, 1, u_aibj, 1)
      call dscal((wf%n_v)**2*(wf%n_o)**2, two, u_aibj, 1)
!
      call add_1432_to_1234(-one, t_aibj, u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(D_ov, wf%n_v, wf%n_o) ! ordered v,o
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  u_aibj,        & ! u_bj_ai
                  wf%n_o*wf%n_v, &
                  tbar_ai,       & ! tbar_ai
                  1,             &
                  zero,          &
                  D_ov,          & ! D_bj
                  1)
!
      call mem%dealloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            density(i, wf%n_o + a) = density(i, wf%n_o + a) + D_ov(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_ov, wf%n_v, wf%n_o)
!
   end subroutine gs_one_el_density_doubles_ov_doubles
!
!
end submodule zop_doubles
