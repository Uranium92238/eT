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
submodule (cc2_class) jacobian_cc2
!
!!
!!    Jacobian submodule (CC2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
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
      call wf%initialize_t2()
      call wf%construct_t2()
      call wf%save_jacobian_a1_intermediates()
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
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
      type(timings) :: timer
!
      timer = new_timer('Jacobian transformation CC2')
      call timer%turn_on()
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      call zero_array(rho_ai, wf%n_t1)
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, ai) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai)
!
         enddo
      enddo
!$omp end parallel do
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
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!  
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c_aibj(a, i, b, j) = c(wf%n_o*wf%n_v + aibj)
                     c_aibj(b, j, a, i) = c(wf%n_o*wf%n_v + aibj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!$omp parallel do schedule(static) private(a, i) 
      do i =1, wf%n_o
         do a = 1, wf%n_v
!
            c_aibj(a, i, a, i) = two*c_aibj(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_doubles_b1(rho_ai, c_aibj)
      call wf%jacobian_doubles_c1(rho_ai, c_aibj)
      call wf%jacobian_doubles_d1(rho_ai, c_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
!$omp parallel do schedule(static) private(a, i, ai) 
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai) = rho_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do
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
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            rho_aibj(a, i, a, i) = half*rho_aibj(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do 
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
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!  
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj) = rho_aibj(a, i, b, j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
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
      timer = new_timer('jacobian cc2 b2')
      call timer%turn_on()
!
!     c_aibj/(1/Δ_aibj) 
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            c_aibj(a, i, a, i) = half*c_aibj(a, i, a, i) 
!
         enddo
      enddo
!$omp end parallel do
!
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
