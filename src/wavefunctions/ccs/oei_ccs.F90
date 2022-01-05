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
submodule (ccs_class) oei_ccs
!
!!
!!    One-electron integrals submodule
!!
!!    Submodule containing routines that can be used to 
!!    construct t1-transformed one-electron integrals.
!!
!
      implicit none
!
!
contains
!
!
   module subroutine get_t1_oei_ccs(wf,         &
                                    oei_type,   &
                                    oei, screening)
!!
!!    Get T1 OEI (one-electron integral)
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Calculates and saves the T1-transformed one-electron integrals specified 
!!    by 'oei_type' in the array 'oei'. 
!!
!!    Arguments:
!!
!!       oei:      (n_mo x n_mo x n_components)-array where integrals are placed 
!!       oei_type: (string) which integral to calculate 
!!
!!    Valid values of 'oei_type' are:
!!
!!       - 'hamiltonian'  One-electron Hamiltonian (h)  n_components = 1
!!       - 'overlap'      AO overlap (S)                n_components = 1
!!       - 'dipole'       Dipole moment (mu)            n_components = 3 (mu_x, mu_y, mu_z)
!!       - 'quadrupole'   Quadrupole moment (q)         n_components = 6 (q_xx, q_xy, q_xz, 
!!                                                                        q_yy, q_yz, q_zz)
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      character(len=*), intent(in) :: oei_type 
!
      real(dp), dimension(*), target, intent(out) :: oei  
!
      logical,  intent(in), optional :: screening
!
      real(dp), dimension(:,:,:), pointer :: oei_p 
!
      real(dp), dimension(:,:,:), allocatable :: oei_ao  
!
      integer :: n_components, k
!
!     Determine number of components of the OEI 
!
      n_components = wf%ao%get_n_oei_components(oei_type)
!
!     Allocate and calculate the AO integrals 
!
      call mem%alloc(oei_ao, wf%ao%n, wf%ao%n, n_components)
!
      call wf%ao%get_oei(oei_type, oei_ao, screening)
!
!     Transform the AO integrals to the T1-transformed basis
!
      oei_p(1 : wf%n_mo, 1 : wf%n_mo, 1 : n_components) => oei(1 : wf%n_mo**2 * n_components)
!
      do k = 1, n_components
!
         call wf%mo_transform(oei_ao(:,:,k), oei_p(:,:,k))
         call wf%t1_transform(oei_p(:,:,k))
!
      enddo
!
      call mem%dealloc(oei_ao, wf%ao%n, wf%ao%n, n_components)
!
   end subroutine get_t1_oei_ccs
!
!
   module subroutine get_t1_oei_ccs_complex(wf,         &
                                            oei_type,   &
                                            oei, screening)
!!
!!    Get T1 OEI (one-electron integral)
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Calculates and saves the T1-transformed one-electron integrals specified 
!!    by 'oei_type' in the array 'oei'. 
!!
!!    Arguments:
!!
!!       oei:      (n_mo x n_mo x n_components)-array where integrals are placed 
!!       oei_type: (string) which integral to calculate 
!!
!!    Valid values of 'oei_type' are:
!!
!!       - 'hamiltonian'  One-electron Hamiltonian (h)  n_components = 1
!!       - 'overlap'      AO overlap (S)                n_components = 1
!!       - 'dipole'       Dipole moment (mu)            n_components = 3 (mu_x, mu_y, mu_z)
!!       - 'quadrupole'   Quadrupole moment (q)         n_components = 6 (q_xx, q_xy, q_xz, 
!!                                                                        q_yy, q_yz, q_zz)
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      character(len=*), intent(in) :: oei_type 
!
      complex(dp), dimension(*), target, intent(out) :: oei  
!
      logical,  intent(in), optional :: screening
!
      complex(dp), dimension(:,:,:), pointer :: oei_p 
!
      real(dp), dimension(:,:), allocatable :: oei_mo ! MO integrals for a given component
!
      real(dp), dimension(:,:,:), allocatable :: oei_ao 
!
      integer :: n_components, k
!
!     Determine number of components of the OEI 
!
      n_components = wf%ao%get_n_oei_components(oei_type)
!
!     Allocate and calculate the AO integrals 
!
      call mem%alloc(oei_ao, wf%ao%n, wf%ao%n, n_components)
!
      call wf%ao%get_oei(oei_type, oei_ao, screening)
!
!     Transform the AO integrals to the T1-transformed basis
!
      oei_p(1 : wf%n_mo, 1 : wf%n_mo, 1 : n_components) => oei(1 : wf%n_mo**2 * n_components)
!
      call mem%alloc(oei_mo, wf%n_mo, wf%n_mo)
!
      do k = 1, n_components
!
         call wf%mo_transform(oei_ao(:,:,k), oei_mo)
!
         oei_p(:,:,k) = cmplx(oei_mo, zero, dp)
!
         call wf%t1_transform_complex(oei_p(:,:,k))
!
      enddo
!
      call mem%dealloc(oei_mo, wf%n_mo, wf%n_mo)
      call mem%dealloc(oei_ao, wf%ao%n, wf%ao%n, n_components)
!
   end subroutine get_t1_oei_ccs_complex
!
!
end submodule oei_ccs
