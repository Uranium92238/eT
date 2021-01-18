!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
submodule (ccsd_class) complex_ccsd
!
!!
!!    Complex submodule (CCSD)
!!
!!    Gathers routines that makes the CCSD wavefunction complex, and that are otherwise related to
!!    the complex wavefunction.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine make_complex_ccsd(wf)
!!
!!    Make complex (CCSD)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      call wf%make_ccs_complex()
      call wf%make_doubles_complex()
!
   end subroutine make_complex_ccsd
!
!
   module subroutine cleanup_complex_ccsd(wf)
!!
!!    Make complex (CCSD)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      call wf%cleanup_ccs_complex()
      call wf%cleanup_doubles_complex()
!
   end subroutine cleanup_complex_ccsd
!
!
   module subroutine construct_complex_time_derivative_amplitudes_ccsd(wf, ddt_amplitudes)
!!
!!    Construct complex time derivative amplitudes (CCSD)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the complex time derivative of the CCSD amplitudes (complex). The time derivative
!!    of the amplitudes is, according to Koch and Jørgensen (1990), given by
!!
!!       ddt_amplitudes = -i*omega
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_amplitudes
!
!     Multiply doubles diagonal by two to remove ai .ge. bj restriction and to incorporate
!     factor 1/2 in doubles amplitudes to match form of T2 operator in Koch and Jørgensen (1990).
!
      if (allocated(wf%t1_complex) .and. allocated(wf%t2_complex)) then
!
         call wf%construct_omega_complex(ddt_amplitudes)
         call scale_diagonal(two, ddt_amplitudes(wf%n_t1+1:wf%n_gs_amplitudes), wf%n_t1)
         call zscal(wf%n_gs_amplitudes, cmplx(zero, -one, dp), ddt_amplitudes, 1)
!
      else
!
         call output%error_msg( &
            "Need complex amplitudes to construct the time derivative of the amplitudes.")
!
      endif
!
   end subroutine construct_complex_time_derivative_amplitudes_ccsd
!
!
   module subroutine construct_complex_time_derivative_multipliers_ccsd(wf, ddt_multipliers)
!!
!!    Construct complex time derivative multipliers (CCSD)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the complex time derivative of the CCSD multipliers (complex). The time derivative
!!    of the multipliers is, according to Koch and Jørgensen (1990), given by
!!
!!       ddt_multipliers = i*multiplier_equation
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_multipliers
!
      if (allocated(wf%t1_complex) .and. allocated(wf%t1bar_complex) &
          .and. allocated(wf%t2_complex) .and. allocated(wf%t2bar_complex)) then
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
   end subroutine construct_complex_time_derivative_multipliers_ccsd
!
!
end submodule complex_ccsd