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
submodule (ccsd_class) set_get_ccsd
!
!!
!!    Set get submodule
!!
!!    Gathers routines that set and get the CCSD type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine set_amplitudes_ccsd(wf, amplitudes)
!!
!!    Set amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: amplitudes
!
      call dcopy(wf%n_t1, amplitudes, 1, wf%t1, 1)
      call dcopy(wf%n_t2, amplitudes(wf%n_t1 + 1), 1, wf%t2, 1)
!
   end subroutine set_amplitudes_ccsd
!
!
   module subroutine get_amplitudes_ccsd(wf, amplitudes)
!!
!!    Get amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes) :: amplitudes
!
      call dcopy(wf%n_t1, wf%t1, 1, amplitudes, 1)
      call dcopy(wf%n_t2, wf%t2, 1, amplitudes(wf%n_t1 + 1), 1)
!
   end subroutine get_amplitudes_ccsd
!
!
   module subroutine set_multipliers_ccsd(wf, multipliers)
!!
!!    Set multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: multipliers
!
      call dcopy(wf%n_t1, multipliers, 1, wf%t1bar, 1)
      call dcopy(wf%n_t2, multipliers(wf%n_t1 + 1), 1, wf%t2bar, 1)
!
   end subroutine set_multipliers_ccsd
!
!
   module subroutine get_multipliers_ccsd(wf, multipliers)
!!
!!    Get multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes) :: multipliers
!
      call dcopy(wf%n_t1, wf%t1bar, 1, multipliers, 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, multipliers(wf%n_t1 + 1), 1)
!
   end subroutine get_multipliers_ccsd
!
!
   module subroutine get_amplitude_block_sizes_ccsd(wf, amplitude_block_sizes)
!!
!!    Get amplitude block sizes
!!    Written by Alexander C. Paul, June 2022
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      integer, dimension(:), allocatable, intent(out) :: amplitude_block_sizes
!
      amplitude_block_sizes = [wf%n_t1, wf%n_t2]
!
   end subroutine get_amplitude_block_sizes_ccsd
!
!
end submodule set_get_ccsd