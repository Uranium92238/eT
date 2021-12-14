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
submodule (mlcc2_class) file_handling_mlcc2
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
   module subroutine read_doubles_vector_mlcc2(wf, file_, vector, expected_n, read_n)
!!
!!    Read doubles vector
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    File format:
!!    n_t1, t1, n_t2, t2
!!
!!    read_n: optionally adds the number of amplitudes read to read_n
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      integer, intent(in) :: expected_n
      real(dp), dimension(expected_n), intent(out) :: vector
!
      integer, intent(inout) :: read_n
!
      integer(i64) :: n
      integer      :: n_t2, iostat
!
      call file_%read_(n, i64+1 + (wf%n_t1+1)*dp, iostat)
      n_t2 = int(n)
!
      if (.not. is_iostat_end(iostat)) then
!
         if (n_t2 .ne. expected_n) then
            call output%error_msg('Wrong number of doubles amplitudes in (a0): &
                                  &expected (i0), found (i0)', &
                                  chars=[file_%get_name()], ints=[expected_n, n_t2])
         end if
!
         call file_%read_(vector, expected_n)
!
         read_n = read_n + n_t2
!
      end if
!
   end subroutine read_doubles_vector_mlcc2
!
!
   module subroutine save_doubles_vector_mlcc2(wf, file_, vector, n)
!!
!!    Save doubles vector
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    File format: energy, n_t1, t1, n_x2, t2
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      integer, intent(in) :: n
      real(dp), dimension(n), intent(in) :: vector
!
      integer(i64) :: m
!
      m = int(n, i64)
!
      call file_%write_(m, dp*(wf%n_t1+1) + i64+1)
      call file_%write_(vector, n)
!
   end subroutine save_doubles_vector_mlcc2
!
!
   module subroutine save_mlcc_orbitals_mlcc2(wf)
!!
!!    Save MLCC orbitals
!!    Written by Sarai D. Folkestad
!!
!!    File format:
!!
!!    1: n_ccs_o, n_ccs_v, n_cc2_o, n_cc2_v, orbitals
!!    2: orbital energies
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      call wf%orbital_coefficients_mlcc_file%open_('write','rewind')
!
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_ccs_o)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_ccs_v)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_cc2_o)
      call wf%orbital_coefficients_mlcc_file%write_(wf%n_cc2_v)
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
   end subroutine save_mlcc_orbitals_mlcc2
!
!
   module subroutine read_mlcc_orbitals_mlcc2(wf)
!!
!!    Read MLCC orbitals
!!    Written by Sarai D. Folkestad
!!
!!    File format:
!!
!!    1: n_ccs_o, n_ccs_v, n_cc2_o, n_cc2_v, orbitals
!!    2: orbital energies
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      call wf%orbital_coefficients_mlcc_file%open_('read','rewind')
!
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_ccs_o)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_ccs_v)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_cc2_o)
      call wf%orbital_coefficients_mlcc_file%read_(wf%n_cc2_v)
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
   end subroutine read_mlcc_orbitals_mlcc2
!
!
   module subroutine read_excitation_vector_file_mlcc2(wf, file_, vector, energy, read_n)
!!
!!    Read excitation vector file
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Reads excitation vector from file structured as follows:
!!    excitation_energy, n_t1, X1
!!
!!    read_n: optionally adds the number of amplitudes read to read_n
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: vector
!
      real(dp), intent(out) :: energy
!
      integer, intent(inout), optional :: read_n
      integer :: n
!
      n = 0
!
      call file_%open_('read', 'rewind')
!
      call file_%read_(energy)
      call wf%read_singles_vector(file_, vector, n)
      call wf%read_doubles_vector(file_, vector(wf%n_t1+1:), wf%n_x2, n)
!
      call file_%close_
!
      if (present(read_n)) read_n = n
!
   end subroutine read_excitation_vector_file_mlcc2
!
!
   module subroutine save_excitation_vector_on_file_mlcc2(wf, file_, vector, energy)
!!
!!    Save excitation vector on file
!!    Written by Alexander C. Paul, Sep 2020
!!
!!    Writes excitation vector o file structured as follows:
!!    excitation_energy, n_t1, X1
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: vector
!
      real(dp), intent(in) :: energy
!
      call file_%open_('write', 'rewind')
!
      call file_%write_(energy)
      call wf%save_singles_vector(file_, vector)
      call wf%save_doubles_vector(file_, vector(wf%n_t1+1:), wf%n_x2)
!
      call file_%close_
!
   end subroutine save_excitation_vector_on_file_mlcc2
!
!
   module subroutine get_restart_vector_mlcc2(wf, file_, vector, energy, &
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
      use array_utilities, only: copy_and_scale, zero_array
      use reordering, only: construct_packed_contravariant
      use reordering, only: construct_packed_covariant
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      type(stream_file), intent(inout) :: file_
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
      integer  :: n_amplitudes_read
!
      call wf%read_excitation_vector_file(file_, vector, energy, n_amplitudes_read)
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
            if (restart_from == 'right' .and. restart_to == 'left') then
!
!              2R^ab_ij - R^ba_ij
               alpha = two
               call construct_packed_contravariant(vector(wf%n_t1+1:wf%n_es_amplitudes), &
                                                   temp, wf%n_cc2_v, wf%n_cc2_o)
!
            else if (restart_from == 'left' .and. restart_to == 'right') then
!
!              1/3 (2L^ab_ij + L^ba_ij)
               alpha = half
               call construct_packed_covariant(vector(wf%n_t1+1:wf%n_es_amplitudes), &
                                               temp, wf%n_cc2_v, wf%n_cc2_o)
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
                               chars=[file_%get_name()])
!
      end if
!
   end subroutine get_restart_vector_mlcc2
!
!
end submodule file_handling_mlcc2
