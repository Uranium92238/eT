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
   module subroutine prepare_mos_hf(wf)
!!
!!    Prepare MOs
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine prepare_mos_hf
!
!
   module subroutine prepare_frozen_fock_contributions_hf(wf)
!!
!!    Prepare frozen Fock conttributions
!!    Written by Sarai D. Folkestad, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine prepare_frozen_fock_contributions_hf
!
!
   module subroutine remove_core_orbitals_hf(wf)
!!
!!    Remove core orbitals
!!    Written by Sarai D. Folkestad, Sep 2018 
!!
!
      implicit none 
!
      class(hf) :: wf  
!
   end subroutine remove_core_orbitals_hf
!
!
   module subroutine remove_frozen_hf_orbitals_hf(wf)
!!
!!    Remove frozen hf orbitals
!!    Written by Sarai D. Folkestad, Feb 2019
!!    Added and modified for HF by Ida-Marie Hoyvik, Oct 2019
!!
!
      use array_utilities, only : copy_and_scale
!
      implicit none
!
      class(hf), intent(inout) :: wf
!
   end subroutine remove_frozen_hf_orbitals_hf
!
!
   module subroutine construct_mo_fock_fc_contribution_hf(wf)
!!
!!    Calculate MO Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine construct_mo_fock_fc_contribution_hf
!
!
   module subroutine construct_mo_fock_frozen_hf_contribution_hf(wf)
!!
!!    Construct MO fock frozen hf  contribution
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine construct_mo_fock_frozen_hf_contribution_hf
!
!
   module subroutine initialize_orbital_coefficients_frozen_hf_hf(wf)
!!
!!    Initialize orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_orbital_coefficients_frozen_hf_hf
!
!
   module subroutine destruct_orbital_coefficients_frozen_hf_hf(wf)
!!
!!    Destruct orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_orbital_coefficients_frozen_hf_hf
!
!
   module subroutine initialize_orbital_coefficients_fc_hf(wf)
!!
!!    Initialize orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_orbital_coefficients_fc_hf
!
!
   module subroutine destruct_orbital_coefficients_fc_hf(wf)
!!
!!    Destruct orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_orbital_coefficients_fc_hf
!
!
