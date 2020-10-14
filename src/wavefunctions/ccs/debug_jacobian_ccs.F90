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
submodule (ccs_class) debug_jacobian_ccs
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
   module subroutine omega_for_jacobian_debug_ccs(wf, omega, t)
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
      class(ccs), intent(inout) :: wf
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
      call wf%set_amplitudes(t)
!
      call wf%eri%update_t1_integrals(wf%t1)
!
      call wf%construct_fock()
      call wf%construct_omega(omega)
!
      call wf%set_amplitudes(t_copy)
      call mem%dealloc(t_copy, wf%n_gs_amplitudes)
!
   end subroutine omega_for_jacobian_debug_ccs
!
!
   module subroutine amplitudes_for_jacobian_debug_ccs(wf, t)
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
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: t
!
      call wf%get_amplitudes(t)
!
   end subroutine amplitudes_for_jacobian_debug_ccs
!
!
   module subroutine normalization_for_jacobian_debug_ccs(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Routine does nothing for CCS
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
!
      integer, intent(in) :: nu
!
      call output%printf('v', 'No normalization needed for CCS (f20.5)', &
                         reals=[A_numerical_mu_nu(nu)])     
!
   end subroutine normalization_for_jacobian_debug_ccs
!
!
   module subroutine numerical_test_jacobian_ccs(wf, store_on_file)
!!
!!    Numerical test for Jacobian matrix
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the Jacobian matrix numerically
!!    by differentiating Ω wrt the ground state 
!!    amplitudes. Compares to the elements of
!!    the Jacobian matrix obtained by transforming
!!    unit vectors.
!!
!!    The optional logical parameter "store_on_file"
!!    can be used to write the two Jacobian matrices
!!    on file.
!!
!!    NOTE: This test scales as N^9 for CCS !! 
!!
!!    For now, we stop the test with error message if wf%name_ = 'cc2'
!!
      implicit none
!
      class(ccs) :: wf
!
      logical, optional :: store_on_file
!
      logical :: store_on_file_local
!
      real(dp), parameter :: dt = 1.0d-5
      real(dp), parameter :: tolerance = 1.0d-9
!
      real(dp), dimension(:), allocatable :: A_numerical_mu_nu, A_mu_nu_diff
      real(dp), dimension(:), allocatable :: omega_mu_p, omega_mu_m
      real(dp), dimension(:), allocatable :: t, t_dt, e
!
      integer :: nu, abs_max_index
!
      real(dp) :: abs_max
!
      type(sequential_file) :: A_file, A_numerical_file
!
      store_on_file_local = .false.
!
      if (present(store_on_file)) store_on_file_local = store_on_file
!
      call mem%alloc(A_numerical_mu_nu, wf%n_es_amplitudes)
      call mem%alloc(A_mu_nu_diff, wf%n_es_amplitudes)
!
      call mem%alloc(omega_mu_p, wf%n_es_amplitudes)
      call mem%alloc(omega_mu_m, wf%n_es_amplitudes)
!
      call mem%alloc(t, wf%n_es_amplitudes)
      call mem%alloc(t_dt, wf%n_es_amplitudes)
      call mem%alloc(e, wf%n_es_amplitudes)
!
      if (store_on_file_local) then
!
!        Prepare sequential file
!
         A_file = sequential_file('Jacobian_matrix', 'unformatted')
         A_numerical_file = sequential_file('Jacobian_matrix_numerical', 'unformatted')
!
         call A_file%open_('write', 'rewind')
         call A_numerical_file%open_('write', 'rewind')
!
      endif
!
      call wf%amplitudes_for_jacobian_debug(t)
!
      call wf%prepare_for_jacobian()
!
      do nu = 1, wf%n_es_amplitudes
!
!        Numerical A_μ_ν
!
         call dcopy(wf%n_es_amplitudes, t, 1, t_dt, 1)
         t_dt(nu) = t_dt(nu) - dt
!
         call wf%omega_for_jacobian_debug(omega_mu_m, t_dt)
!
         t_dt(nu) = t_dt(nu) + two*dt
!
         call wf%omega_for_jacobian_debug(omega_mu_p, t_dt)
!
!        Central difference
!
         call dcopy(wf%n_es_amplitudes, omega_mu_p, 1, A_numerical_mu_nu, 1)
         call daxpy(wf%n_es_amplitudes, -one, omega_mu_m, 1, A_numerical_mu_nu, 1)
         call dscal(wf%n_es_amplitudes, one/(two*dt), A_numerical_mu_nu, 1)
!
         call wf%normalization_for_jacobian_debug(A_numerical_mu_nu, nu)
!
!        Analytical A_μ_ν
!
         call zero_array(e, wf%n_es_amplitudes)
         e(nu) = one
!
         call wf%eri%update_t1_integrals(wf%t1)
!
         call wf%construct_fock()
!
         call wf%jacobian_transformation(e)
!
         call dcopy(wf%n_es_amplitudes, e, 1, A_mu_nu_diff, 1)
         call daxpy(wf%n_es_amplitudes, -one, A_numerical_mu_nu, 1, A_mu_nu_diff, 1)
!
         call get_abs_max_w_index(A_mu_nu_diff, wf%n_es_amplitudes, abs_max, abs_max_index)
!
!        Store on file, if requested
!
         if (store_on_file_local) then
!
            call A_file%write_(e, wf%n_es_amplitudes)
            call A_numerical_file%write_(A_numerical_mu_nu, wf%n_es_amplitudes)
!
         endif
!
         if (abs_max .gt. tolerance) then
            call output%warning_msg('Numerical and analytical Jacobian matrices do not match!')
         endif
!
      enddo
!
      call mem%dealloc(A_numerical_mu_nu, wf%n_es_amplitudes)
      call mem%dealloc(A_mu_nu_diff, wf%n_es_amplitudes)
!
      call mem%dealloc(omega_mu_p, wf%n_es_amplitudes)
      call mem%dealloc(omega_mu_m, wf%n_es_amplitudes)
!
      call mem%dealloc(t, wf%n_es_amplitudes)
      call mem%dealloc(t_dt, wf%n_es_amplitudes)
      call mem%dealloc(e, wf%n_es_amplitudes)
!
      if (store_on_file_local) then
!
         call A_file%close_('keep')
         call A_numerical_file%close_('keep')
!
      endif
!
      call output%printf('n', 'Exited numerical test for jacobian without errors.', &
                         fs='(/t3, a)')
!
   end subroutine numerical_test_jacobian_ccs
!
!
end submodule debug_jacobian_ccs
