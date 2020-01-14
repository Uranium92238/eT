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
   module subroutine t1_transform_ccs_complex(wf, Z_pq)
!!
!!    T1 transform
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Assumes that Z is in the MO basis and performs the T1 transformation,
!!
!!       Z_pq <- sum_rs X_ps Z_sr Y_qr,    i.e.    Z <- X Z Y^T
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    Here, t1 is a full MO matrix whose only non-zero block is the vir-occ
!!    part, where it is equal to t_i^a.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
   end subroutine t1_transform_ccs_complex
!
!
   module subroutine t1_transform_4_ccs_complex(wf, Z_tuvw, Z_pqrs, t1)
!!
!!    T1 transform 4 index arrays
!!    Written by Andreas Skeidsvoll, Apr 2019
!!
!!    Assumes that Z is in the MO basis and performs the T1 transformation,
!!
!!       Z_pqrs = sum_tuvw X_pt Y_qu X_rm Y_sn Z_tuvw,
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    Here, t1 is a full MO matrix whose only non-zero block is the vir-occ
!!    part, where it is equal to t_i^a.
!!    NB: needs place for an additional 2*wf%n_mo**4 + wf%n_t1 in memory.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      complex(dp), dimension(wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo), intent(in) :: Z_tuvw
      complex(dp), dimension(wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo), intent(out) :: Z_pqrs
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in) :: t1
!
   end subroutine t1_transform_4_ccs_complex
