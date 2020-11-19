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
   module subroutine construct_fock_ccs_complex(wf)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_pq = h_pq + sum_k (2*g_pqkk - g_pkkq) + (effective Fock contributions)
!!
!!    Effective Fock contributions:
!!
!!       Frozen core by Sarai D. Folkestad, 2019
!!       QM/MM by Tommaso Giovannini, 2019
!!       QM/PCM by Tommaso Giovannini, 2019
!!
!!       Modified by Sarai D. Folkestad, Nov 2019
!!
!!       Added batching for N^2 memory requirement.
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine construct_fock_ccs_complex
!
!
   module subroutine add_frozen_fock_terms_ccs_complex(wf, F_pq)
!!
!!    Add frozen Fock terms
!!    Written by Sarai D. Folkestad, 2019 
!!
!!    Adds the frozen core contributions to
!!    the effective T1-transformed Fock matrix.
!!
!!    Isolated into subroutine by Eirik F. Kjønstad, 2019    
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      complex(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq 
!
   end subroutine add_frozen_fock_terms_ccs_complex
!
!
   module subroutine construct_t1_frozen_fock_terms_ccs_complex(wf, F_pq)
!!
!!    Calculate T1 Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
   end subroutine construct_t1_frozen_fock_terms_ccs_complex
!
!
   module subroutine add_t1_fock_length_dipole_term_ccs_complex(wf, electric_field)
!!
!!    Add t1 Fock length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds dipole part of the length gauge electromagnetic potential to the Fock matrix,
!!
!!       Fock matrix += -μ·E,
!!
!!    where μ is the vector of electric dipole integral matrices and E is a uniform classical electric
!!    vector field. This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(3), intent(in) :: electric_field
!
   end subroutine add_t1_fock_length_dipole_term_ccs_complex
