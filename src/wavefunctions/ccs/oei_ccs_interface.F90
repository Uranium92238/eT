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
   module subroutine get_t1_oei_ccs(wf,         &
                                    oei_type,   &
                                    oei)
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
      character(len=*), intent(in) :: oei_type 
      real(dp), dimension(*), target, intent(out) :: oei  
!
   end subroutine get_t1_oei_ccs
!
!
   module subroutine get_t1_oei_ccs_complex(wf,         &
                                            oei_type,   &
                                            oei)
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
      character(len=*), intent(in) :: oei_type 
      complex(dp), dimension(*), target, intent(out) :: oei  
!
   end subroutine get_t1_oei_ccs_complex
