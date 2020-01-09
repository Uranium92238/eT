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
      call wf%integrals%make_eri_complex()
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
end submodule complex_ccs