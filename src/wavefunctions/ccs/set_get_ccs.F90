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
submodule (ccs_class) set_get_ccs
!
!!
!!    Set get submodule
!!
!!    Gathers routines that set and get the CCS type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine set_amplitudes_ccs(wf, amplitudes)
!!
!!    Set amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: amplitudes
!
      call dcopy(wf%n_gs_amplitudes, amplitudes, 1, wf%t1, 1)
!
   end subroutine set_amplitudes_ccs
!
!
   module subroutine get_amplitudes_ccs(wf, amplitudes)
!!
!!    Get amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes) :: amplitudes
!
      call dcopy(wf%n_gs_amplitudes, wf%t1, 1, amplitudes, 1)
!
   end subroutine get_amplitudes_ccs
!
!
   module subroutine set_multipliers_ccs(wf, multipliers)
!!
!!    Set multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: multipliers
!
      call dcopy(wf%n_gs_amplitudes, multipliers, 1, wf%t1bar, 1)
!
   end subroutine set_multipliers_ccs
!
!
   module subroutine get_multipliers_ccs(wf, multipliers)
!!
!!    Get multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes) :: multipliers
!
      call dcopy(wf%n_gs_amplitudes, wf%t1bar, 1, multipliers, 1)
!
   end subroutine get_multipliers_ccs
!
!
   module subroutine set_fock_ccs(wf, F_pq)
!!
!!    Set Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Sets the different blocks of the Fock matrix based on the full
!!    matrix sent to the routine.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_pq
!
      integer :: i, j, a, b
!
!$omp parallel do private(i,j)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            wf%fock_ij(i,j) = F_pq(i,j)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ia(i,a) = F_pq(i, wf%n_o + a)
            wf%fock_ai(a,i) = F_pq(wf%n_o + a, i)
!
         enddo
      enddo
!$omp end parallel do 
!
!$omp parallel do private(a,b)
      do b = 1, wf%n_v
         do a = 1, wf%n_v
!
            wf%fock_ab(a,b) = F_pq(wf%n_o + a, wf%n_o + b)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine set_fock_ccs
!
!
   module subroutine set_excitation_energies_ccs(wf, energies, side)
!!
!!    Set excitation energies
!!    Written by Alexander C. Paul, Sep 2020
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_singlet_states), intent(in) :: energies
      character(len=*), intent(in) :: side
!
      if (trim(side) .eq. 'right') then 
!
         call dcopy(wf%n_singlet_states, energies, 1, &
                    wf%right_excitation_energies, 1)
!
      elseif (trim(side) .eq. 'left') then 
!
         call dcopy(wf%n_singlet_states, energies, 1, &
                    wf%left_excitation_energies, 1)
!
      elseif (trim(side) .eq. 'both') then 
!
         call dcopy(wf%n_singlet_states, energies, 1, &
                    wf%left_excitation_energies, 1)
!
         call dcopy(wf%n_singlet_states, energies, 1, &
                    wf%right_excitation_energies, 1)
!
      endif
!
   end subroutine set_excitation_energies_ccs
!
!
end submodule set_get_ccs
