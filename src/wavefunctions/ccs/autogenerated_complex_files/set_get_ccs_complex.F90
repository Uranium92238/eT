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
submodule (ccs_class) set_get_ccs_complex
!
!!
!!    Set get submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
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
   module subroutine set_amplitudes_ccs_complex(wf, amplitudes)
!!
!!    Set amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes), intent(in) :: amplitudes
!
      call zcopy(wf%n_gs_amplitudes, amplitudes, 1, wf%t1_complex, 1)
!
   end subroutine set_amplitudes_ccs_complex
!
!
   module subroutine get_amplitudes_ccs_complex(wf, amplitudes)
!!
!!    Get amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes) :: amplitudes
!
      call zcopy(wf%n_gs_amplitudes, wf%t1_complex, 1, amplitudes, 1)
!
   end subroutine get_amplitudes_ccs_complex
!
!
   module subroutine set_multipliers_ccs_complex(wf, multipliers)
!!
!!    Set multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes), intent(in) :: multipliers
!
      call zcopy(wf%n_gs_amplitudes, multipliers, 1, wf%t1bar_complex, 1)
!
   end subroutine set_multipliers_ccs_complex
!
!
   module subroutine get_multipliers_ccs_complex(wf, multipliers)
!!
!!    Get multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes) :: multipliers
!
      call zcopy(wf%n_gs_amplitudes, wf%t1bar_complex, 1, multipliers, 1)
!
   end subroutine get_multipliers_ccs_complex
!
!
   module subroutine set_fock_ccs_complex(wf, F_pq)
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
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_pq
!
      integer :: i, j, a, b
!
!$omp parallel do private(i,j)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            wf%fock_ij_complex(i,j) = F_pq(i,j)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ia_complex(i,a) = F_pq(i, wf%n_o + a)
            wf%fock_ai_complex(a,i) = F_pq(wf%n_o + a, i)
!
         enddo
      enddo
!$omp end parallel do 
!
!$omp parallel do private(a,b)
      do b = 1, wf%n_v
         do a = 1, wf%n_v
!
            wf%fock_ab_complex(a,b) = F_pq(wf%n_o + a, wf%n_o + b)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine set_fock_ccs_complex
!
!
end submodule set_get_ccs_complex
