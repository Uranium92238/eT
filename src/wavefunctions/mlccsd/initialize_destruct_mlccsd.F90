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
submodule (mlccsd_class) initialize_desctruct_mlccsd
!
!!
!!    Initialize destruct submodule (MLCC2)
!!
!!    Gathers routines that initialize and destruct the MLCC2 allocatable variables.
!!
!
   implicit none
!
!
contains
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
      call mem%alloc(wf%u_aibj, wf%n_cc2_v + wf%n_ccsd_v, &
                                wf%n_cc2_o + wf%n_ccsd_o, &
                                wf%n_cc2_v + wf%n_ccsd_v, &
                                wf%n_cc2_o + wf%n_ccsd_o)
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
      if (allocated(wf%u_aibj)) call mem%dealloc(wf%u_aibj, wf%n_cc2_v + wf%n_ccsd_v, &
                                                            wf%n_cc2_o + wf%n_ccsd_o, &
                                                            wf%n_cc2_v + wf%n_ccsd_v, &
                                                            wf%n_cc2_o + wf%n_ccsd_o)
!
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
      call mem%alloc(wf%t2, wf%n_t2)
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
      if (allocated(wf%t2)) call mem%dealloc(wf%t2, wf%n_t2)
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
      call wf%initialize_t1()
      call wf%initialize_t2()
      call wf%initialize_x2()
      call wf%initialize_u_aibj()
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
      call wf%destruct_t1()
      call wf%destruct_t2()
      call wf%destruct_x2()
      call wf%destruct_u_aibj()
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
      call mem%alloc(wf%orbital_coefficients_cc2, wf%n_ao, wf%n_mo)
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
      if (allocated (wf%orbital_coefficients_cc2)) &
        call mem%dealloc(wf%orbital_coefficients_cc2, wf%n_ao, wf%n_mo)
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
      call mem%alloc(wf%orbital_energies_cc2, wf%n_mo)
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
      if (allocated (wf%orbital_energies_cc2)) &
        call mem%dealloc(wf%orbital_energies_cc2, wf%n_mo)
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
      call mem%alloc(wf%O_o,  wf%n_ccsd_o + wf%n_cc2_o,  wf%n_ccsd_o + wf%n_cc2_o)
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
      if (allocated (wf%O_o)) &
        call mem%dealloc(wf%O_o,  wf%n_ccsd_o + wf%n_cc2_o,  wf%n_ccsd_o + wf%n_cc2_o)
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
      call mem%alloc(wf%O_v,  wf%n_ccsd_v + wf%n_cc2_v,  wf%n_ccsd_v + wf%n_cc2_v)
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
      if (allocated (wf%O_v)) &
        call mem%dealloc(wf%O_v,  wf%n_ccsd_v + wf%n_cc2_v,  wf%n_ccsd_v + wf%n_cc2_v)
!
    end subroutine destruct_O_v_mlccsd
!
!
end submodule initialize_desctruct_mlccsd