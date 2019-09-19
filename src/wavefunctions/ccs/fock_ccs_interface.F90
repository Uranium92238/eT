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
   module subroutine construct_fock_ccs(wf)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad,
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine construct_fock_ccs
!
!
   module subroutine coulomb_contribution_fock_fc_ccs(wf)
!!
!!    Coulomb contribution to frozen core Fock
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
   end subroutine coulomb_contribution_fock_fc_ccs
!
!
   module subroutine exchange_contribution_fock_fc_ccs(wf)
!!
!!    Exchange contribution to frozen core fock
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
   end subroutine exchange_contribution_fock_fc_ccs
!
!
   module subroutine construct_mo_fock_fc_contribution_ccs(wf)
!!
!!    Calculate Fock frozen core
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
   end subroutine construct_mo_fock_fc_contribution_ccs
!
!
   module subroutine construct_t1_fock_fc_contribution_ccs(wf, F_pq)
!!
!!    Calculate T1 Fock frozen core
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
   end subroutine construct_t1_fock_fc_contribution_ccs
!
!