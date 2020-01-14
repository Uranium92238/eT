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
submodule (cc2_class) jacobian_cc2
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix 
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!   
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!  
! 
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_cc2(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      type(timings), allocatable :: timer
!
      timer = timings('Prepare for Jacobian CC2', pl='normal')
      call timer%turn_on()
!
      call wf%initialize_t2()
      call wf%construct_t2()
      call wf%save_jacobian_a1_intermediates()
!
      call timer%turn_off()
!
   end subroutine prepare_for_jacobian_cc2
!
!
   module subroutine jacobian_transformation_cc2(wf, c)
!!
!!    Jacobian transformation
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_ai = rho_ai,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable     :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj
!
      real(dp), dimension(:,:), allocatable     :: rho_ai   
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj
!
      type(timings) :: timer
!
      timer = timings('Jacobian transformation CC2', pl='normal')
      call timer%turn_on()
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      call zero_array(rho_ai, wf%n_t1)
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, c, 1, c_ai, 1)
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call wf%jacobian_doubles_a1(rho_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
      call squareup(c(wf%n_t1 + 1 : wf%n_es_amplitudes), c_aibj, wf%n_t1)
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
      call scale_diagonal(two, c_aibj, wf%n_t1)
!
      call wf%jacobian_doubles_b1(rho_ai, c_aibj)
      call wf%jacobian_doubles_c1(rho_ai, c_aibj)
      call wf%jacobian_doubles_d1(rho_ai, c_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
      call dcopy(wf%n_t1, rho_ai, 1, c, 1)
!
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     :: CC2 contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(rho_aibj, wf%n_t1**2)
!
!     Contributions from singles vector c
!
      call wf%jacobian_doubles_a2(rho_aibj, c_ai)
!
      call scale_diagonal(half, rho_aibj, wf%n_t1)
!
      call symmetric_sum(rho_aibj, wf%n_t1)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c      
!
      call wf%jacobian_cc2_b2(rho_aibj, c_aibj)
!
      call mem%dealloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!     Overwrite the incoming doubles c vector & pack in
!
      call packin(c(wf%n_t1 + 1 : wf%n_es_amplitudes), rho_aibj, wf%n_t1)
!
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_cc2
!
!
   module subroutine jacobian_cc2_b2_cc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^B2 = ε_aibj c_aibj/(1/Δ_aibj) 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)     :: rho_aibj   
!
      integer :: i, j, a, b
!
      type(timings) :: timer
!
      timer = timings('Jacobian CC2 B2', pl='verbose')
      call timer%turn_on()
!
!     c_aibj/(1/Δ_aibj) 
!
      call scale_diagonal(half, c_aibj, wf%n_o*wf%n_v)
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + c_aibj(a,i,b,j)*&
                                          (- wf%orbital_energies(i) &
                                           - wf%orbital_energies(j) &
                                           + wf%orbital_energies(wf%n_o + a) &
                                           + wf%orbital_energies(wf%n_o + b) )
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_b2_cc2
!
!
end submodule jacobian_cc2
