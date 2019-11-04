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
submodule (ccs_class) oei_ccs
!
!!
!!    One-electron integrals submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!!    Submodule containing routines that can be used to construct t1-transformed one-electron integrals.
!!
!
      implicit none
!
!
contains
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
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)  :: x_wx 
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: y_pq 
!
      call wf%mo_transform(x_wx, y_pq)
      call wf%t1_transform(y_pq)
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
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)  :: x_wx
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: y_pq
!
      real(dp), dimension(:,:), allocatable :: y_pq_real
!
      call mem%alloc(y_pq_real, wf%n_mo, wf%n_mo)
!
      call wf%mo_transform(x_wx, y_pq_real)
!
      y_pq = cmplx(y_pq_real, zero, dp)
!
      call mem%dealloc(y_pq_real, wf%n_mo, wf%n_mo)
!
      call wf%t1_transform_complex(y_pq)
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
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(out) :: mu_pqk 
!
      call wf%get_mo_mu(mu_pqk)
!
      call wf%t1_transform(mu_pqk(:,:,1))
      call wf%t1_transform(mu_pqk(:,:,2))
      call wf%t1_transform(mu_pqk(:,:,3))
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
!
      complex(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(out) :: mu_pqk 
!
      real(dp), dimension(:,:,:), allocatable :: mu_pqk_real
!
      call mem%alloc(mu_pqk_real, wf%n_mo, wf%n_mo, 3)
!
      call wf%get_mo_mu(mu_pqk_real)
!
      mu_pqk = cmplx(mu_pqk_real, zero, dp)
!
      call mem%dealloc(mu_pqk_real, wf%n_mo, wf%n_mo, 3)
!
      call wf%t1_transform_complex(mu_pqk(:,:,1))
      call wf%t1_transform_complex(mu_pqk(:,:,2))
      call wf%t1_transform_complex(mu_pqk(:,:,3))
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
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: h_pq 
!
      call wf%get_mo_h(h_pq)
      call wf%t1_transform(h_pq)
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
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: h_pq
!
      real(dp), dimension(:,:), allocatable :: h_pq_real
!
      call mem%alloc(h_pq_real, wf%n_mo, wf%n_mo)
!
      call wf%get_mo_h(h_pq_real)
!
      h_pq = cmplx(h_pq_real, zero, dp)
!
      call mem%dealloc(h_pq_real, wf%n_mo, wf%n_mo)
!
      call wf%t1_transform_complex(h_pq)
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
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 6), intent(out) :: q_pqk 
!
      call wf%get_mo_q(q_pqk)
!
      call wf%t1_transform(q_pqk(:,:,1))
      call wf%t1_transform(q_pqk(:,:,2))
      call wf%t1_transform(q_pqk(:,:,3))
      call wf%t1_transform(q_pqk(:,:,4))
      call wf%t1_transform(q_pqk(:,:,5))
      call wf%t1_transform(q_pqk(:,:,6))
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
!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      complex(dp), dimension(wf%n_mo, wf%n_mo, 6), intent(out) :: q_pqk 
!
      real(dp), dimension(:,:,:), allocatable :: q_pqk_real
!
      call mem%alloc(q_pqk_real, wf%n_mo, wf%n_mo, 6)
!
      call wf%get_mo_q(q_pqk_real)
!
      q_pqk = cmplx(q_pqk_real, zero, dp)
!
      call mem%dealloc(q_pqk_real, wf%n_mo, wf%n_mo, 6)
!
      call wf%t1_transform_complex(q_pqk(:,:,1))
      call wf%t1_transform_complex(q_pqk(:,:,2))
      call wf%t1_transform_complex(q_pqk(:,:,3))
      call wf%t1_transform_complex(q_pqk(:,:,4))
      call wf%t1_transform_complex(q_pqk(:,:,5))
      call wf%t1_transform_complex(q_pqk(:,:,6))
!
   end subroutine construct_q_ccs_complex
!
!
end submodule oei_ccs
