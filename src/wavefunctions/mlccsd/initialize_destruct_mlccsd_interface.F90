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
   module subroutine initialize_u_aibj_mlccsd(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine initialize_u_aibj_mlccsd
!
!
   module subroutine destruct_u_aibj_mlccsd(wf)
!!
!!    Destruct u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine destruct_u_aibj_mlccsd
!
!
   module subroutine initialize_t2_mlccsd(wf)
!!
!!    Initialize t2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine initialize_t2_mlccsd
!
!
   module subroutine destruct_t2_mlccsd(wf)
!!
!!    Destruct t2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine destruct_t2_mlccsd
!
!
   module subroutine initialize_amplitudes_mlccsd(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine initialize_amplitudes_mlccsd
!
!
   module subroutine destruct_amplitudes_mlccsd(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
   end subroutine destruct_amplitudes_mlccsd
!
!
    module subroutine initialize_orbital_coefficients_cc2_mlccsd(wf)
!!
!!    Initialize orbital coefficients cc2
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine initialize_orbital_coefficients_cc2_mlccsd
!
!
    module subroutine destruct_orbital_coefficients_cc2_mlccsd(wf)
!!
!!    Destruct orbital coefficients cc2
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine destruct_orbital_coefficients_cc2_mlccsd
!
!
    module subroutine initialize_orbital_energies_cc2_mlccsd(wf)
!!
!!    Initialize orbital energies cc2
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine initialize_orbital_energies_cc2_mlccsd
!
!
    module subroutine destruct_orbital_energies_cc2_mlccsd(wf)
!!
!!    Destruct orbital energies cc2
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine destruct_orbital_energies_cc2_mlccsd
!
!
    module subroutine initialize_O_o_mlccsd(wf)
!!
!!    Initialize O_o
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine initialize_O_o_mlccsd
!
!
    module subroutine destruct_O_o_mlccsd(wf)
!!
!!    Destruct O_o
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine destruct_O_o_mlccsd
!
!
    module subroutine initialize_O_v_mlccsd(wf)
!!
!!    Initialize O_v 
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine initialize_O_v_mlccsd
!
!
    module subroutine destruct_O_v_mlccsd(wf)
!!
!!    Destruct O_v
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(mlccsd) :: wf
!
    end subroutine destruct_O_v_mlccsd
