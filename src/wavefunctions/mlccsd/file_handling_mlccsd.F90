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
submodule (mlccsd_class) file_handling_mlccsd
!
!!
!!    File handling submodule
!!
!!    Gathers routines that save wavefunction parameters to file,
!!    and reads them from file, plus other routines related to the
!!    handling of the files that belong to the wavefunction.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine get_restart_vector_mlccsd(wf, storer, vector, energy, &
                                              restart_from, restart_to)
!!
!!    Get restart vector
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Gets start vector and energy from file and
!!    handles the basis transformations according to:
!!
!!    restart from "right" to "left"
!!    L^a_i = 2R^a_i
!!    L^ab_ij = 4R^ab_ij - 2R^ba_ij
!!
!!    restart from "left" to "right"
!!    R^a_i = 1/2 L^a_i
!!    R^ab_ij = 1/6 (2L^ab_ij + L^ba_ij)
!!
      use array_initialization, only: zero_array, copy_and_scale
      use reordering, only: construct_packed_contravariant
      use reordering, only: construct_packed_covariant
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      type(amplitude_file_storer), intent(inout) :: storer
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
!
      real(dp), intent(out) :: energy
!
      character(len=*), intent(in) :: restart_from, restart_to
!
      real(dp), dimension(:), allocatable :: temp
!
      real(dp) :: alpha
      integer  :: n_amplitudes_read, n_v, n_o
!
      n_o = wf%n_cc2_o + wf%n_ccsd_o
      n_v = wf%n_cc2_v + wf%n_ccsd_v
!
      call storer%read_(vector, energy, n_amplitudes_read)
!
      if (n_amplitudes_read == wf%n_t1) then
!
         call zero_array(vector(wf%n_t1+1:wf%n_es_amplitudes), wf%n_x2)
!
      else if (n_amplitudes_read == wf%n_es_amplitudes) then
!
!        Handle restart "left" from "right" and "right" from "left"
         if (restart_from /= restart_to) then
!
            call mem%alloc(temp, wf%n_x2)
            alpha = one
!
            if (restart_from == 'left' .and. restart_to == 'right') then
!
!              2R^ab_ij - R^ba_ij
               alpha = two
               call construct_packed_contravariant(vector(wf%n_t1+1:wf%n_es_amplitudes), &
                                                   temp, n_v, n_o)
!
            else if (restart_from == 'right' .and. restart_to == 'left') then
!
!              1/3 (2L^ab_ij + L^ba_ij)
               alpha = half
               call construct_packed_covariant(vector(wf%n_t1+1:wf%n_es_amplitudes), &
                                               temp, n_v, n_o)
!
            end if
!
            call dscal(wf%n_t1, alpha, vector, 1)
            call copy_and_scale(alpha, temp, vector(wf%n_t1+1:wf%n_es_amplitudes), wf%n_x2)
            call mem%dealloc(temp, wf%n_x2)
!
         end if
!
      else
!
         call output%error_msg('Did not recognize number of amplitudes in (a0) &
                               &expected (i0) or (i0) found (i0) amplitudes.', &
                               ints=[wf%n_es_amplitudes, wf%n_t1, n_amplitudes_read], &
                               chars=[storer%get_filename()])
!
      end if
!
   end subroutine get_restart_vector_mlccsd
!
!
   module subroutine save_mlcc_orbitals_mlccsd(wf)
!!
!!    Save MLCC orbitals
!!    Written by Sarai D. Folkestad
!!
!!    File format:
!!
!!    1: n_ccs_o, n_ccs_v, n_cc2_o, n_cc2_v, n_ccsd_o, n_ccsd_v, orbitals
!!    2: orbital energies
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      call wf%orbital_coefficients_mlcc_file%open_('write','rewind')
!
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_ccs_o)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_ccs_v)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_cc2_o)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_cc2_v)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_ccsd_o)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_ccsd_v)
!
      call wf%orbital_coefficients_mlcc_file%write_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
      call wf%orbital_coefficients_mlcc_file%close_('keep')
!
!     Print MLCC orbital energies to file
!
      call wf%orbital_energies_mlcc_file%open_('write', 'rewind')
      call wf%orbital_energies_mlcc_file%write_(wf%orbital_energies, wf%n_mo)
      call wf%orbital_energies_mlcc_file%close_('keep')
!
   end subroutine save_mlcc_orbitals_mlccsd
!
!
   module subroutine read_mlcc_orbitals_mlccsd(wf)
!!
!!    Read MLCC orbitals
!!    Written by Sarai D. Folkestad
!!
!!    File format:
!!
!!    1: n_ccs_o, n_ccs_v, n_cc2_o, n_cc2_v, n_ccsd_o, n_ccsd_v, orbitals
!!    2: orbital energies
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      call wf%orbital_coefficients_mlcc_file%open_('read','rewind')
!
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_ccs_o)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_ccs_v)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_cc2_o)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_cc2_v)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_ccsd_o)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_ccsd_v)
!
      call wf%orbital_coefficients_mlcc_file%read_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
      call wf%orbital_coefficients_mlcc_file%close_('keep')
!
!     Print MLCC orbital energies to file
!
      call wf%orbital_energies_mlcc_file%open_('read', 'rewind')
      call wf%orbital_energies_mlcc_file%read_(wf%orbital_energies, wf%n_mo)
      call wf%orbital_energies_mlcc_file%close_('keep')
!
   end subroutine read_mlcc_orbitals_mlccsd
!
!
end submodule file_handling_mlccsd
