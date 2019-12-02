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
submodule (ccs_class) jacobian_ccs
!
!!
!!    Jacobian submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_ccs(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
!     For now, do nothing.
!
      call output%printf('- No preparations for the ' // trim(wf%name_) // &
                         ' excited state equation.', pl='v', fs='(/t3,a)')
!
   end subroutine prepare_for_jacobian_ccs
!
!
   module subroutine jacobian_transformation_ccs(wf, c)
!!
!!    Jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
!!    On exit, c is overwritten by rho.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: rho
!
!     Allocate the transformed vector & add the terms to it
!
      call mem%alloc(rho, wf%n_v, wf%n_o)
      call zero_array(rho, wf%n_t1)
!
      call wf%jacobian_ccs_a1(rho, c)
      call wf%jacobian_ccs_b1(rho, c)
!
!     Then overwrite the c vector with the transformed vector
!
      call dcopy((wf%n_o)*(wf%n_v), rho, 1, c, 1)
      call mem%dealloc(rho, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transformation_ccs
!
!
   module subroutine jacobian_ccs_a1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_b F_ab c_bi - sum_j F_ji c_aj,
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1
!
!     sum_b F_a_b c_b_i 
!
      call dgemm('N', 'N',     &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  wf%fock_ab,  &
                  wf%n_v,      &
                  c1,          &
                  wf%n_v,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
!     - sum_j c_a_j F_j_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_o,      &
                  -one,        &
                  c1,          &
                  wf%n_v,      &
                  wf%fock_ij,  &
                  wf%n_o,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
   end subroutine jacobian_ccs_a1_ccs
!
!
   module subroutine jacobian_ccs_b1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the B1 term,
!!
!!       sum_bj L_aijb c_bj = sum_bj (2 g_aijb - g_abji) c_bj,
!!
!!    and adds it to the rho1 vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_abji
      real(dp), dimension(:,:,:,:), allocatable :: L_aijb
!
      real(dp), dimension(:,:), allocatable :: c_jb
!
      type(batching_index) :: batch_b
!
      integer :: req0, req1, j, b, b_red, current_b_batch
!
      batch_b = batching_index(wf%n_v)
!
      req0 = wf%n_o*wf%n_v*wf%integrals%n_J ! L_ai^J
      req1 = 2*wf%n_v*wf%n_o**2 + & ! L_aijb and g_abji
               wf%n_v*wf%integrals%n_J ! L_ab^J
!
      call mem%batch_setup(batch_b, req0, req1)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
!        Construct L_aijb = 2 g_aijb - g_abji
!
         call mem%alloc(L_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
!
         call wf%get_voov(L_aijb,     &
                           1, wf%n_v,  &
                           1, wf%n_o,  &
                           1, wf%n_o,  &
                           batch_b%first, batch_b%last)
!
         call dscal(((wf%n_o)**2)*(wf%n_v)*(batch_b%length), two, L_aijb, 1)
!
!        Construct L_aijb = 2 g_aijb - g_abji
!
         call mem%alloc(g_abji, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_abji,                        &
                           1, wf%n_v,                    &
                           batch_b%first, batch_b%last,  &
                           1, wf%n_o,                    &
                           1, wf%n_o)
!
         call add_1432_to_1234(-one, g_abji, L_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
!
         call mem%dealloc(g_abji, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
!        Reorder c1 to do multiply with L_aijb
!
         call mem%alloc(c_jb, (wf%n_o), batch_b%length)
!
         do b = batch_b%first, batch_b%last
            do j = 1, wf%n_o
!
               b_red = b - batch_b%first + 1
!
               c_jb(j, b_red) = c1(b, j)
!
            enddo
         enddo

!
         call dgemm('N', 'N',                   &
                     (wf%n_v)*(wf%n_o),         &
                     1,                         &
                     (wf%n_o)*(batch_b%length), &
                     one,                       &
                     L_aijb,                    &
                     (wf%n_v)*(wf%n_o),         &
                     c_jb,                      &
                     (wf%n_o)*batch_b%length,   &
                     one,                       &
                     rho1,                      &
                     (wf%n_v)*(wf%n_o))
!
         call mem%dealloc(L_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
         call mem%dealloc(c_jb, (wf%n_o), (batch_b%length))
!
      enddo ! batch_b
!
   end subroutine jacobian_ccs_b1_ccs
!
!
end submodule jacobian_ccs