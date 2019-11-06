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
   module subroutine add_molecular_mechanics_fock_term_ccs(wf, F_pq)
!!
!!    Add molecular mechanics Fock contribution 
!!    Written by Tommaso Giovannini, 2019 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: F_pq 
!
   end subroutine add_molecular_mechanics_fock_term_ccs
!
!

   module subroutine add_pcm_fock_contribution_ccs(wf, F_pq)
!!
!!    Add PCM Fock contribution 
!!    Written by Tommaso Giovannini, 2019 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: F_pq 
!
   end subroutine add_pcm_fock_contribution_ccs
!
!
   module subroutine add_frozen_core_fock_term_ccs(wf, F_pq)
!!
!!    Add frozen core Fock contribution 
!!    Written by Sarai D. Folkestad, 2019 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq    
!
   end subroutine add_frozen_core_fock_term_ccs
!
!
   module subroutine add_frozen_hf_fock_term_ccs(wf, F_pq)
!!
!!    Add frozen HF Fock contribution 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019 
!!  
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq 
!
   end subroutine add_frozen_hf_fock_term_ccs
!
!
   module subroutine construct_t1_fock_fc_term_ccs(wf, F_pq)
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
   end subroutine construct_t1_fock_fc_term_ccs
!
!
   module subroutine construct_t1_fock_frozen_hf_term_ccs(wf, F_pq)
!!
!!    Calculate T1 Fock frozen fock contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
   end subroutine construct_t1_fock_frozen_hf_term_ccs
!
!
   module subroutine add_t1_fock_length_dipole_term_ccs(wf, electric_field)
!!
!!    Add t1 Fock length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      real(dp), dimension(3), intent(in) :: electric_field
!
   end subroutine add_t1_fock_length_dipole_term_ccs