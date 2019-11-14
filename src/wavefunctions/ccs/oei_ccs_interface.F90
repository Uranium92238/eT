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
   module subroutine ao_to_t1_transformation_ccs(wf, x_wx, y_pq)
!!
!!    AO to T1 transformation 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Takes in an AO array and returns the T1-transformed array: 
!!
!!    x_wx (in)   array in the AO basis (w and x are AO indices)
!!    y_pq (out)  array in the T1-transformed basis (p and q are MO indices) 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)  :: x_wx 
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: y_pq 
!
   end subroutine ao_to_t1_transformation_ccs
!
!
   module subroutine ao_to_t1_transformation_ccs_complex(wf, x_wx, y_pq)
!!
!!    AO to T1 transformation 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Takes in an AO array and returns the T1-transformed array: 
!!
!!    x_wx (in)   array in the AO basis (w and x are AO indices)
!!    y_pq (out)  array in the T1-transformed basis (p and q are MO indices) 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)  :: x_wx
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: y_pq
!
   end subroutine ao_to_t1_transformation_ccs_complex
!
!
   module subroutine construct_mu_ccs(wf, mu_pqk)
!!
!!    Construct mu
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Constructs 
!!
!!       mu_pqk, k = 1, 2, 3 (x, y, z)
!!
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      real(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(out) :: mu_pqk 
!
   end subroutine construct_mu_ccs
!
!
   module subroutine construct_mu_ccs_complex(wf, mu_pqk)
!!
!!    Construct mu
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Constructs 
!!
!!       mu_pqk, k = 1, 2, 3 (x, y, z)
!!
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      complex(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(out) :: mu_pqk 
!
   end subroutine construct_mu_ccs_complex
!
!
   module subroutine construct_h_ccs(wf, h_pq)
!!
!!    Construct h
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    Constructs the one-electron Hamiltonian integrals h_pq 
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: h_pq 
!
   end subroutine construct_h_ccs
!
!
   module subroutine construct_h_ccs_complex(wf, h_pq)
!!
!!    Construct h
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    Constructs the one-electron Hamiltonian integrals h_pq 
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: h_pq
!
   end subroutine construct_h_ccs_complex
!
!
   module subroutine construct_q_ccs(wf, q_pqk)
!!
!!    Construct q
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    Constructs 
!!
!!       q_pqk, k = 1, 2, 3, 4, 5, 6 (xx, xy, xz, yy, yz, and zz)
!!
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      real(dp), dimension(wf%n_mo, wf%n_mo, 6), intent(out) :: q_pqk 
!
   end subroutine construct_q_ccs
!
!
   module subroutine construct_q_ccs_complex(wf, q_pqk)
!!
!!    Construct q
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    Constructs 
!!
!!       q_pqk, k = 1, 2, 3, 4, 5, 6 (xx, xy, xz, yy, yz, and zz)
!!
!!    in the T1 transformed basis.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
      complex(dp), dimension(wf%n_mo, wf%n_mo, 6), intent(out) :: q_pqk 
!
   end subroutine construct_q_ccs_complex
