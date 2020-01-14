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
submodule (mlcc2_class) debug_jacobian_mlcc2
!
!!
!!    Debug Jacobian
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
   module subroutine omega_for_jacobian_debug_mlcc2(wf, omega, t)
!!
!!    Omega for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
!!    Calculates omega 
!!    with dimension n_es_amplitudes for
!!    t given on input. For methods where 
!!
!!       n_es_amplitudes = n_gs_amplitudes
!!
!!    this simply entails calling construct omega.
!!    For e.g. CC2 this routine must be overwritten
!!    to obtain Ω_μ2
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: t
!
      real(dp), dimension(:), allocatable :: t_copy
!
      call mem%alloc(t_copy, wf%n_gs_amplitudes)
!
      call wf%get_amplitudes(t_copy)
      call wf%set_amplitudes(t(1:wf%n_gs_amplitudes))
!
      call wf%integrals%update_t1_integrals(wf%t1)
!
      wf%x2 = t(wf%n_gs_amplitudes+1:wf%n_es_amplitudes)
!
      call wf%construct_fock()
      call wf%construct_u_aibj()
!
      call zero_array(omega, wf%n_es_amplitudes)
!
      call wf%omega_ccs_a1(omega(1:wf%n_gs_amplitudes))
!
      call wf%omega_cc2_a1(omega(1:wf%n_gs_amplitudes), wf%n_cc2_o, wf%n_cc2_v, wf%first_cc2_o, wf%first_cc2_v, &
                           wf%last_cc2_o, wf%last_cc2_v)
      call wf%omega_cc2_b1(omega(1:wf%n_gs_amplitudes), wf%n_cc2_o, wf%n_cc2_v, wf%first_cc2_o, wf%first_cc2_v, &
                           wf%last_cc2_o, wf%last_cc2_v)
      call wf%omega_cc2_c1(omega(1:wf%n_gs_amplitudes), wf%n_cc2_o, wf%n_cc2_v, wf%first_cc2_o, wf%first_cc2_v)
!
      call wf%construct_omega_doubles(omega(wf%n_gs_amplitudes+1:wf%n_es_amplitudes))
!
      call wf%set_amplitudes(t_copy)
      call wf%construct_x2()
      call mem%dealloc(t_copy, wf%n_gs_amplitudes)
!
   end subroutine omega_for_jacobian_debug_mlcc2
!
!
   module subroutine amplitudes_for_jacobian_debug_mlcc2(wf, t)
!!
!!    Amplitudes for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the amplitudes 
!!    with dimension n_es_amplitudes
!!
!!    For methods where 
!!
!!       n_es_amplitudes = n_gs_amplitudes
!!
!!    this simply entails calling get_amplitudes.
!!    For e.g. CC2 this routine must be overwritten
!!    to obtain t_μ2
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: t
!
      call wf%get_amplitudes(t(1:wf%n_t1))
!
      call wf%construct_x2()
!
      call dcopy(wf%n_x2, wf%x2, 1, t(wf%n_t1 + 1), 1)
!
   end subroutine amplitudes_for_jacobian_debug_mlcc2
!
!
   module subroutine normalization_for_jacobian_debug_mlcc2(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Routine does nothing for CCS
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
!
      integer, intent(in) :: nu
!
      logical :: scale_
!
      integer :: ai, aibj, aiai
!
      if (nu .le. wf%n_t1) return
!
      scale_ = .false. 
!
      aibj = nu - wf%n_t1
!
      do ai = 1, wf%n_cc2_v*wf%n_cc2_o
!
         aiai = ai*(ai-3)/2 + 2*ai
         if (aiai == aibj) scale_ = .true.
!
      enddo
!
      if (scale_) then
!
         call dscal(wf%n_es_amplitudes, two, A_numerical_mu_nu, 1) 
!
      else
!
         call dscal(wf%n_x2, half, A_numerical_mu_nu(wf%n_t1+1), 1) 
!
      endif 
!
   end subroutine normalization_for_jacobian_debug_mlcc2
!
!
end submodule debug_jacobian_mlcc2