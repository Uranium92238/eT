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
   module subroutine t1_transform_ccs(wf, Z_pq)
!!
!!    T1 transform
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq
!
   end subroutine t1_transform_ccs
!
!
   module subroutine construct_mu_ccs(wf, mu_pqk)
!!
!!    Construct mu
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
      use libint_initialization
!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(inout) :: mu_pqk 
!
   end subroutine construct_mu_ccs
!
!
   module subroutine construct_h_ccs(wf, h_pq)
!!
!!    Construct h
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: h_pq 
!
   end subroutine construct_h_ccs
!
!
   module subroutine construct_q_ccs(wf, q_pqk)
!!
!!    Construct q
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!    
!!    xx, xy, xz, yy, yz, and zz.
!!
      use libint_initialization
!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 6), intent(inout) :: q_pqk 
!
   end subroutine construct_q_ccs