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
submodule (ccs_class) complex_ccs
!
!!
!!    Complex submodule
!!
!!    Gathers routines that makes the CCS wavefunction complex, and that are otherwise related to
!!    the complex wavefunction.
!
   implicit none
!
!
contains
!
!
   module subroutine make_complex_ccs(wf)
!!
!!    Make complex (CCS)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call wf%make_ccs_complex()
!
   end subroutine make_complex_ccs
!
!
   module subroutine make_ccs_complex_ccs(wf)
!!
!!    Make CCS complex (CCS)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Allocates complex CCS variables, puts the real variables into the real part of the complex
!!    variables and zero into the imaginary part, and deallocates the real variables.
!!
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      wf%hf_energy_complex = cmplx(wf%hf_energy, zero, dp)
!
      if (allocated(wf%t1)) then
         call wf%initialize_t1_complex()
         wf%t1_complex = cmplx(wf%t1, zero, dp)
         call wf%destruct_t1()
      endif
!
      if (allocated(wf%t1bar)) then
         call wf%initialize_t1bar_complex()
         wf%t1bar_complex = cmplx(wf%t1bar, zero, dp)
         call wf%destruct_t1bar()
      endif
!
      if (allocated(wf%fock_ij)) then
         call wf%initialize_fock_ij_complex()
         wf%fock_ij_complex = cmplx(wf%fock_ij, zero, dp)
         call wf%destruct_fock_ij()
      endif
!
      if (allocated(wf%fock_ia)) then
         call wf%initialize_fock_ia_complex()
         wf%fock_ia_complex = cmplx(wf%fock_ia, zero, dp)
         call wf%destruct_fock_ia()
      endif
!
      if (allocated(wf%fock_ai)) then
         call wf%initialize_fock_ai_complex()
         wf%fock_ai_complex = cmplx(wf%fock_ai, zero, dp)
         call wf%destruct_fock_ai()
      endif
!
      if (allocated(wf%fock_ab)) then
         call wf%initialize_fock_ab_complex()
         wf%fock_ab_complex = cmplx(wf%fock_ab, zero, dp)
         call wf%destruct_fock_ab()
      endif
!
      if (allocated(wf%density)) then
         call wf%initialize_gs_density_complex()
         wf%density_complex = cmplx(wf%density, zero, dp)
         call wf%destruct_gs_density()
      endif
!
      call wf%integral_preparations_complex(wf%L_mo%n_J)
!
      call wf%complexify_cholesky_vectors(wf%L_mo, wf%L_mo_c)
      call wf%L_t1_c%set_equal_to(wf%L_mo_c)
!
      call wf%construct_t1_cholesky_complex(wf%t1_complex, wf%L_mo_c, wf%L_t1_c)
!
   end subroutine make_ccs_complex_ccs
!
!
   module subroutine cleanup_complex_ccs(wf)
!!
!!    Cleanup complex (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call wf%cleanup_ccs_complex()
!
   end subroutine cleanup_complex_ccs
!
!
   module subroutine cleanup_ccs_complex_ccs(wf)
!!
!!    Cleanup CCS complex (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
!!    Deallocates complex CCS variables.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      if (allocated(wf%t1_complex)) then
         call wf%destruct_t1_complex()
      endif
!
      if (allocated(wf%t1bar_complex)) then
         call wf%destruct_t1bar_complex()
      endif
!
      if (allocated(wf%fock_ij_complex)) then
         call wf%destruct_fock_ij_complex()
      endif
!
      if (allocated(wf%fock_ia_complex)) then
         call wf%destruct_fock_ia_complex()
      endif
!
      if (allocated(wf%fock_ai_complex)) then
         call wf%destruct_fock_ai_complex()
      endif
!
      if (allocated(wf%fock_ab_complex)) then
         call wf%destruct_fock_ab_complex()
      endif
!
      if (allocated(wf%density_complex)) then
         call wf%destruct_gs_density_complex()
      endif
!
      if (allocated(wf%L_t1_c)) call wf%L_t1_c%remove_observer('eri t1 complex')
!
      if (allocated(wf%eri_t1_c)) deallocate(wf%eri_t1_c)
      if (allocated(wf%L_mo_c)) deallocate(wf%L_mo_c)
      if (allocated(wf%L_t1_c)) deallocate(wf%L_t1_c)
!
   end subroutine cleanup_ccs_complex_ccs
!
!
   module subroutine construct_complex_time_derivative_ccs(wf, ddt_amplitudes_multipliers)
!!
!!    Construct complex time derivative (CCS)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the complex time derivative of the CCS amplitudes and multipliers. Does not have to
!!    be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), intent(out), dimension(2*wf%n_gs_amplitudes) :: ddt_amplitudes_multipliers
!
      call wf%construct_complex_time_derivative_amplitudes( &
         ddt_amplitudes_multipliers(1:wf%n_gs_amplitudes))
      call wf%construct_complex_time_derivative_multipliers( &
         ddt_amplitudes_multipliers(wf%n_gs_amplitudes+1:2*wf%n_gs_amplitudes))
!
   end subroutine construct_complex_time_derivative_ccs
!
!
   module subroutine construct_complex_time_derivative_amplitudes_ccs(wf, ddt_amplitudes)
!!
!!    Construct complex time derivative of amplitudes (CCS)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the complex time derivative of the CCS amplitudes (complex). The time derivative
!!    of the amplitudes is, according to Koch and Jørgensen (1990), given by
!!
!!       ddt_amplitudes = -i*omega
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_amplitudes
!
      if (allocated(wf%t1_complex)) then
!
         call wf%construct_omega_complex(ddt_amplitudes)
         call zscal(wf%n_gs_amplitudes, cmplx(zero, -one, dp), ddt_amplitudes, 1)
!
      else
!
         call output%error_msg( &
            "Need complex amplitudes to construct the time derivative of the amplitudes.")
!
      endif
!
   end subroutine construct_complex_time_derivative_amplitudes_ccs
!
!
   module subroutine construct_complex_time_derivative_multipliers_ccs(wf, ddt_multipliers)
!!
!!    Construct complex time derivative multipliers (CCS)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the complex time derivative of the CCS multipliers (complex). The time derivative
!!    of the multipliers is, according to Koch and Jørgensen (1990), given by
!!
!!       ddt_multipliers = i*multiplier_equation
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_multipliers
!
      if (allocated(wf%t1_complex) .and. allocated(wf%t1bar_complex)) then
!
         call wf%prepare_for_multiplier_equation_complex()
         call wf%construct_multiplier_equation_complex(ddt_multipliers)
         call zscal(wf%n_gs_amplitudes, cmplx(zero, one, dp), ddt_multipliers, 1)
!
      else
!
         call output%error_msg( &
            "Need complex amplitudes and multipliers to construct the time derivative.")
!
      endif
!
   end subroutine construct_complex_time_derivative_multipliers_ccs
!
!
   module subroutine complexify_cholesky_vectors(wf, L_mo, L_mo_c)
!!
!!    Complexify Cholesky vectors
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      class(abstract_eri_cholesky) :: L_mo
      class(abstract_eri_cholesky_c) :: L_mo_c
!
      type(batching_index) :: batch_q
!
      integer :: req_0, req_1, batch
!
      real(dp), dimension(:,:,:), allocatable :: L_J_pq
!
      batch_q = batching_index(wf%n_mo)
!
      req_0 = 0
      req_1 = L_mo_c%n_J*wf%n_mo &
            + L_mo%get_memory_estimate(1, wf%n_mo, 1, 1) &
            + L_mo_c%get_memory_estimate(1, wf%n_mo, 1, 1)
!
      call mem%batch_setup(batch_q, req_0, req_1, tag='Complexify Cholesky')
!
      call mem%alloc(L_J_pq, L_mo_c%n_J, wf%n_mo, batch_q%max_length)
!
      do batch = 1, batch_q%num_batches
!
         call batch_q%determine_limits(batch)
!
         call L_mo%get(L_J_pq, 1, wf%n_mo, batch_q%first, batch_q%get_last())
         call L_mo_c%set(cmplx(L_J_pq, 0.0d0, kind=dp), 1, wf%n_mo, batch_q%first, batch_q%get_last())
!
      enddo
!
      call mem%dealloc(L_J_pq, L_mo_c%n_J, wf%n_mo, batch_q%max_length)
!
      call mem%batch_finalize()
!
   end subroutine complexify_cholesky_vectors
!
!
   module subroutine integral_preparations_ccs_complex(wf, n_J)
!!
!!    Integral preparations
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_J
!
      call output%printf('m', '- Settings for integral handling:', fs='(/t3,a)')
!
      call wf%cholesky_factory_complex(n_J)
      call wf%integral_factory_complex()
!
   end subroutine integral_preparations_ccs_complex
!
!
   module subroutine integral_factory_complex(wf)
!!
!!    Integral factory complex
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=200) :: integral_storage
!
      logical :: has_room_for_eri
!
      logical :: try_to_place_eri_in_memory, eri_in_memory
!
      integral_storage = 'memory'
      call input%get_keyword('eri storage', 'integrals', integral_storage)
!
      if (trim(integral_storage) /= 'memory' .and. trim(integral_storage) /= 'none') &
         call output%error_msg('illegal value for eri storage')
!
      try_to_place_eri_in_memory = (integral_storage == 'memory' .and. wf%need_g_abcd)
      has_room_for_eri = (wf%n_mo**4)*(2*dp) < mem%get_available()/5
!
      eri_in_memory = try_to_place_eri_in_memory .and. has_room_for_eri
!
      if (eri_in_memory) then
         wf%eri_t1_c = eri_adapter_c(eri_memory_tool_c(wf%L_t1_c), wf%n_o, wf%n_v)
      else
         wf%eri_t1_c = eri_adapter_c(eri_tool_c(wf%L_t1_c), wf%n_o, wf%n_v)
      endif
!
      call wf%L_t1_c%add_observer('eri t1 complex', wf%eri_t1_c%eri)
!
      call output%printf('m', &
                         'ERI matrix in memory:       (l0)', &
                         logs=[eri_in_memory],               &
                         fs='(t6,a)')
!
   end subroutine integral_factory_complex
!
!
   module subroutine calculate_energy_omega_term_ccs_complex(wf)
!!
!!    Calculate energy_complex omega term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds multipliers dot omega to the energy_complex,
!!
!!       energy_complex += sum_mu tbar_mu Omega_mu,
!!
!!    which appears in the energy_complex expression:
!!
!!          < Lambda|H|CC > when Omega != 0.
!!
!!    This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(:), allocatable :: multipliers, omega
!
!
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
      call mem%alloc(omega, wf%n_gs_amplitudes)
!
      call wf%get_multipliers_complex(multipliers)
      call wf%construct_omega_complex(omega)
!
      wf%energy_complex = wf%energy_complex + zdot(wf%n_gs_amplitudes, multipliers, 1, omega, 1)
!
      call mem%dealloc(multipliers, wf%n_gs_amplitudes)
      call mem%dealloc(omega, wf%n_gs_amplitudes)
!
   end subroutine calculate_energy_omega_term_ccs_complex
!
!
   module subroutine calculate_energy_length_dipole_term_ccs_complex(wf, electric_field)
!!
!!    Calculate energy_complex length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds dipole part of the length gauge electromagnetic potential to the energy_complex,
!!
!!       energy_complex += 2 sum_ii (-mu·E)_ii,
!!
!!    where mu is the vector of electric dipole integral matrices
!!    and E is a uniform classical electric
!!    vector field. This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(3), intent(in) :: electric_field
!
      complex(dp), dimension(:,:,:), allocatable :: mu
!
      integer :: i
!
!     Construct t1_complex transformed dipole moment
!
      call mem%alloc(mu, wf%n_mo, wf%n_mo, 3)
      call wf%get_t1_oei_complex('dipole', mu)
!
!     Add one_complex-electron electric field contribution to the diagonal
!     of Fock and one_complex-electron integral terms
!
      do i = 1, wf%n_o
!
         wf%energy_complex = wf%energy_complex - two_complex*(mu(i, i, 1)*electric_field(1)   &
                                      + mu(i, i, 2)*electric_field(2) &
                                      + mu(i, i, 3)*electric_field(3))
!
      enddo
!
      call mem%dealloc(mu, wf%n_mo, wf%n_mo, 3)
!
   end subroutine calculate_energy_length_dipole_term_ccs_complex
!
!
end submodule complex_ccs
