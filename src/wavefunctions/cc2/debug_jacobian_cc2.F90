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
submodule (cc2_class) debug_jacobian_cc2
!
!!
!!    Debug Jacobian (CCS)
!!
!!    Routines debug analytical Jacobian by comparing to 
!!    Jacobian computed by numerical differentiation of 
!!    Ω 
!! 
!!       A_μν = dΩ_μ/dt_ν =  d(< μ | exp(-T) H exp(T) | R >)/dt_ν
!!
!!    The numerical differentiation is done using
!!    central differences.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine omega_for_jacobian_debug_cc2(wf, omega, t)
!!
!!    Omega for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: t
!
      real(dp), dimension(:), allocatable :: t_copy
      real(dp), dimension(:,:,:,:), allocatable :: u
      real(dp), dimension(:,:,:,:), allocatable :: t_unpacked
!
      call mem%alloc(t_copy, wf%n_gs_amplitudes)
!
      call wf%get_amplitudes(t_copy)
      call wf%set_amplitudes(t(1:wf%n_gs_amplitudes))
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      if (wf%need_g_abcd()) call wf%integrals%can_we_keep_g_pqrs_t1()
!
      call wf%construct_fock()
!
      call zero_array(omega, wf%n_es_amplitudes)
!
      call wf%omega_ccs_a1(omega(1:wf%n_gs_amplitudes))
!
     ! call mem%alloc(u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      !call squareup(t(wf%n_t1+1:wf%n_es_amplitudes), u, wf%n_t1)
!
      !call dscal(wf%n_t1**2, two, u, 1)
      !call add_packed_1432_to_unpacked_1234(-one, t(wf%n_gs_amplitudes+1:wf%n_es_amplitudes), u, wf%n_v, wf%n_o)
!
   !   call mem%alloc(t_unpacked, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!!
   !   call squareup(t(wf%n_t1+1:wf%n_es_amplitudes), t_unpacked, wf%n_t1)
   !   call copy_and_scale(two, t_unpacked, u, wf%n_t1**2)
   !   call add_1432_to_1234(-one,  t_unpacked, u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
   !   call mem%dealloc(t_unpacked, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%construct_u()
      call wf%omega_doubles_a1(omega(1:wf%n_gs_amplitudes), wf%u)
      call wf%omega_doubles_b1(omega(1:wf%n_gs_amplitudes), wf%u)
      call wf%omega_doubles_c1(omega(1:wf%n_gs_amplitudes), wf%u)
!
   !   call mem%dealloc(u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   !   call wf%construct_omega2(omega(wf%n_gs_amplitudes+1:wf%n_es_amplitudes), t(wf%n_gs_amplitudes+1:wf%n_es_amplitudes))
!
      call wf%set_amplitudes(t_copy)
      call mem%dealloc(t_copy, wf%n_gs_amplitudes)
!
   end subroutine omega_for_jacobian_debug_cc2
!
!
   module subroutine amplitudes_for_jacobian_debug_cc2(wf, t)
!!
!!    Amplitudes for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the amplitudes 
!!    with dimension n_es_amplitudes
!!
!!    Gets the t1-amplitudes 
!!    Constructs the t2-amplitudes
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: t
!
      call wf%get_amplitudes(t(1:wf%n_gs_amplitudes))
!
   !   call wf%initialize_t2()
   !   call wf%construct_t2()
!
   !   call dcopy(wf%n_t2, wf%t2, 1, t(wf%n_t1+1), 1)
   !   call wf%destruct_t2()
!
   end subroutine amplitudes_for_jacobian_debug_cc2
!
!
   module subroutine normalization_for_jacobian_debug_cc2(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Differentiation wrt. doubles amplitudes (nu > wf%n_t1)
!!    will yield factor one half on the diagonal nu = aiai
!!    and factor two on the off-diagonal aibj, ai .ne. bj
!!
!!    A_numerical_μ_aiai = dΩ_μ/dt_aiai = 1/2 ε_mu,aiai δ_μ,aiai
!!    A_numerical_μ_aibj = dΩ_μ/dt_aibj = 2 ε_mu,aibj δ_μ,aibj
!!
!!    For CC2 the doubles-doubles block of A is defined as
!!
!!    A_μ2,ν2 = delta_μ2,ν2 ε_ν2
!!
!!    This routine scales the diagonal elements by factor two
!!    and the off-diagonal elements by factor one half
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
!
      integer, intent(in) :: nu
!
      logical :: diagonal_
!
      integer :: ai, aibj, aiai
!
      if (nu .le. wf%n_t1) return
!
      diagonal_ = .false. 
!
      aibj = nu - wf%n_t1
!
      do ai = 1, wf%n_t1
!
         aiai = ai*(ai-3)/2 + 2*ai
         if (aiai == aibj) diagonal_ = .true.
!
      enddo
!
      if (diagonal_) then
!
         call dscal(wf%n_es_amplitudes, two, A_numerical_mu_nu, 1)   
!
      else     
!
         call dscal(wf%n_es_amplitudes, half, A_numerical_mu_nu, 1)  
! 
      endif
!
   end subroutine normalization_for_jacobian_debug_cc2
!
!
end submodule debug_jacobian_cc2