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
   module subroutine add_t1_terms_ccs_complex(wf, Z_pq)
!!
!!    Add t1 terms
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
!!    Here, Z is assumed to be the density matrix with no T1 contributions on input
!!    - the so-called T1-transformed density matrix.
!!    The routine adds the missing T1 contributions to Z:
!!
!!       Z_pq <- sum_rs X_sp Z_sr Y_rq,    i.e.    Z <- X^T Z Y
!!
!!    where
!!
!!       X = I - t1
!!       Y = I + t1^T
!!
!!    t1 is a full MO matrix whose only non-zero block is the vir-occ part,
!!    where it is equal to t_i^a.
!!
!!    Based on t1_transform_ccs_complex by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
   end subroutine add_t1_terms_ccs_complex
!
!
   module subroutine add_t1_terms_and_transform_ccs_complex(wf, Z_pq, Z_out)
!!
!!    Add t1 terms and transform
!!    Written by Tor S. Haugland, Nov 2019 (as do_visualization)
!!
!!    Here Z, on input, is assumed to be the density matrix with no T1 contributions 
!!    - the so-called T1-transformed density matrix.
!!    The routine adds the missing T1 contributions to Z and transforms it
!!    with the AO coefficients to obtain a density as needed by the visualization tool.
!!
!!       Z_alpha,beta = (sum_pq  Z_pq C_alpha,p C_beta,q) 
!!
!!    Renamed and moved here, by Alexander C. Paul, May 2020
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(in)  :: Z_pq
      complex(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: Z_out
!
   end  subroutine add_t1_terms_and_transform_ccs_complex
