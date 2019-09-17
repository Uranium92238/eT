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
!!    Constructs Ω for CC2, both singles and
!!    doubles part. 
!!
!!    The construct omega routine in CC2
!!    calculates u from the integrals, hence 
!!    it is not called directly (as is done 
!!    for CCS and CCSD), but u is constructed and 
!!    the individual terms are called (ccs_a1, 
!!    doubles_a1-doubles_c1). Additionally, the 
!!    doubles part of Ω is constructed. This routine 
!!    is only used in this debug procedure.
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
!
!     Set the perturbed amplitudes
!
      call mem%alloc(t_copy, wf%n_gs_amplitudes)
!
      call wf%get_amplitudes(t_copy)
      call wf%set_amplitudes(t(1:wf%n_gs_amplitudes))
!
!     Prepare T1 integrals and Fock 
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      if (wf%integrals%get_eri_t1_mem()) call wf%integrals%can_we_keep_g_pqrs_t1()
!
      call wf%construct_fock()
!
!     Ω construction
!
      call zero_array(omega, wf%n_es_amplitudes)
!
      call wf%omega_ccs_a1(omega(1:wf%n_gs_amplitudes))
!
      call mem%alloc(u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(t(wf%n_t1+1:wf%n_es_amplitudes), u, wf%n_t1)
      call dscal(wf%n_t1**2, two, u, 1)
      call add_packed_1432_to_unpacked_1234(-one, t(wf%n_t1+1:wf%n_es_amplitudes), u, wf%n_v, wf%n_o)
!
      call wf%omega_doubles_a1(omega(1:wf%n_gs_amplitudes), u)
      call wf%omega_doubles_b1(omega(1:wf%n_gs_amplitudes), u)
      call wf%omega_doubles_c1(omega(1:wf%n_gs_amplitudes), u)
!
      call mem%dealloc(u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%construct_omega2(omega(wf%n_gs_amplitudes+1:wf%n_es_amplitudes), t(wf%n_gs_amplitudes+1:wf%n_es_amplitudes))
!
!     Reset amplitudes
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
      call wf%get_amplitudes(t(1:wf%n_t1))
!
      call wf%initialize_t2()
      call wf%construct_t2()
!
      call dcopy(wf%n_t2, wf%t2, 1, t(wf%n_t1+1), 1)
      call wf%destruct_t2()
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
!!
!!    This routine scales these diagonal elements 
!!    by factor two
!!
!!    Differentiation of Ω_μ2 wrt. doubles amplitudes
!!    gives factor two on the off-diagonal aibj, ai .ne. bj
!!
!!    A_numerical_μ2_aibj = dΩ_μ2/dt_aibj = 2 ε_μ2,aibj δ_μ2,aibj
!!
!!    For CC2, the doubles-doubles block of A is defined as
!!
!!    A_μ2,ν2 = delta_μ2,ν2 ε_ν2
!!
!!    This routine scales these off-diagonal elements 
!!    by factor half
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
         call dscal(wf%n_t2, half, A_numerical_mu_nu(wf%n_t1+1), 1)  
! 
      endif
!
   end subroutine normalization_for_jacobian_debug_cc2
!
!
end submodule debug_jacobian_cc2