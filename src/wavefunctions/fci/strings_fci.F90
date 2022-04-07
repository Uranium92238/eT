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
submodule (fci_class) strings_fci
!
!!
!! FCI strings management submodule
!!
!! Gathers routines that generate and make operation on FCI strings
!!
!
   use math_utilities, only: binomial
!
   implicit none
!
!
contains
!
!
   module subroutine generate_strings(wf)
!!
!!    Generate strings
!!    Written by Enrico Ronca, 2020
!!
!!    Generate strings
!!
!!    Bit strings with lowest lying orbital to the right.
!!    Occupied orbital -> 1 and unoccupied orbital -> 0
!!
!!    Example - 2 electrons:
!!
!!       Orbitals 1 and 2 occupied -> 11
!!       Orbitals 1 and 3 occupied -> 101
!!       Orbitals 2 and 4 occupied -> 1010
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      type(timings), allocatable :: timer
!
      timer = timings('Generate strings', pl='n')
      call timer%turn_on
!
      call generate_strings_recursive(wf%alpha_strings,     &
                                      wf%n_alpha_strings,   &
                                      wf%n_mo,              &
                                      wf%n_alpha)
!
      call generate_strings_recursive(wf%beta_strings,   &
                                      wf%n_beta_strings, &
                                      wf%n_mo,           &
                                      wf%n_beta)
!
      call timer%turn_off
!
   end subroutine generate_strings
!
!
   recursive module subroutine generate_strings_recursive(strings,       &
                                                          n_strings,     &
                                                          n_orbitals,    &
                                                          n_electrons)
!!
!!    Generate strings recursive
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    Bit strings with lowest lying orbital to the right.
!!    Occupied orbital -> 1 and unoccupied orbital -> 0
!!
!!    Example - 2 electrons:
!!
!!       Orbitals 1 and 2 occupied -> 11
!!       Orbitals 1 and 3 occupied -> 101
!!       Orbitals 2 and 4 occupied -> 1010
!!
!!    Generate strings recursively:
!!
!!    Termination conditions:
!!
!!       1. Fill one electron in n orbitals
!!
!!       2. Fill n electrons into n orbitals
!!
!!    Recursive calls (two cases)
!!
!!       1. Fill (n_electrons) electrons into (n_orbitals - 1) orbitals
!!
!!       2. Fill (n_electrons - 1) electrons into (n_orbitals - 1) orbitals, and add an
!!          electron to the n'th orbital
!!
      implicit none
!
      integer, intent(in) :: n_electrons
      integer, intent(in) :: n_orbitals
      integer, intent(in) :: n_strings
!
      integer, dimension(n_strings), intent(inout):: strings
!
      integer :: n_strings_minus_o
      integer :: n_strings_minus_o_and_e
      integer :: i
!
      integer, dimension(:), allocatable :: strings_minus_o
      integer, dimension(:), allocatable :: strings_minus_o_and_e
!
      if (n_electrons == 0) then
         strings(1) = 0
         return
      end if
!
      if (n_electrons == 1) then
         strings = fill_electron_in_n_orbitals(n_orbitals)
         return
      end if
!
      if (n_electrons == n_orbitals) then
         strings(1) = fill_n_electrons_in_n_orbitals(n_orbitals)
         return
      end if
!
!     Fill all electrons in first (n_orbitals - 1) orbitals.
!     Last orbital is vacant.
!
      n_strings_minus_o = binomial(n_orbitals - 1, n_electrons)
!
      call mem%alloc(strings_minus_o, n_strings_minus_o)
!
      call generate_strings_recursive(strings_minus_o,   &
                                      n_strings_minus_o, &
                                      n_orbitals - 1,    &
                                      n_electrons)
!
      strings(1 : n_strings_minus_o) = strings_minus_o
!
      call mem%dealloc(strings_minus_o, n_strings_minus_o)
!
!     Fill (n_electrons - 1) electrons in first (n_orbitals - 1) orbitals and add electron
!     to the last orbital.
!
      n_strings_minus_o_and_e = binomial(n_orbitals - 1, n_electrons - 1)
      call mem%alloc(strings_minus_o_and_e, n_strings_minus_o_and_e)
!
      call generate_strings_recursive(strings_minus_o_and_e,     &
                                      n_strings_minus_o_and_e,   &
                                      n_orbitals - 1,            &
                                      n_electrons - 1)
!
      do i = 1, n_strings_minus_o_and_e
!
         strings_minus_o_and_e(i) = create_p(strings_minus_o_and_e(i), p=n_orbitals)
!
      enddo
!
      strings(n_strings_minus_o + 1 :) = strings_minus_o_and_e(:)
!
      call mem%dealloc(strings_minus_o_and_e, n_strings_minus_o_and_e)
!
   end subroutine generate_strings_recursive
!
!
   pure function create_p(string_in, p) result(string_out)
!!
!!    Create p
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      integer, intent(in) :: string_in
      integer, intent(in) :: p
!
      integer :: string_out
!
      string_out = ior(string_in, ishft(1, p - 1))
!
   end function create_p
!
!
   pure function annihilate_p(string_in, p) result(string_out)
!!
!!    Annihilate p
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      integer, intent(in) :: string_in
      integer, intent(in) :: p
!
      integer :: string_out
!
      string_out = ieor(string_in, ishft(1, p - 1))
!
   end function annihilate_p
!
!
   pure function fill_electron_in_n_orbitals(n) result(strings)
!!
!!    Fill electron in n orbitals
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    Strings for filling a single electron in n orbitals
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer, dimension(n) :: strings
!
      integer :: i
!
      do i = 1, n
         strings(i) = ishft(1, i - 1)
      end do
!
   end function fill_electron_in_n_orbitals
!
!
   pure function fill_n_electrons_in_n_orbitals(n) result(string)
!!
!!    Fill n electrons in n orbitals
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    String for n electrons in n orbitals
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer :: string
!
      integer :: i
!
      string = 0
!
      do i = 1, n
!
         string = ior(string, ishft(1, i - 1))
!
      end do
!
   end function fill_n_electrons_in_n_orbitals
!
!
   pure function orbital_occupied_in_string(string, orbital_number) result(occupied)
!!
!!    Orbital occupied in string
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      integer, intent(in) :: string, orbital_number
!
      logical :: occupied
!
      occupied = iand(ishft(1, orbital_number - 1), string) /= 0
!
   end function orbital_occupied_in_string
!
!
   module function generate_address_for_string(string, n_orbitals, n_electrons) result(address)
!!
!!    Generate address for string
!!    Written by Enrico Ronca, 2020
!!
!!    Generate the address of the string, i.e., the corresponding
!!    index in the FCI wavefunction.
!!
!!    The strings are ordered according to reverse lexical ordering,
!!    see Helgaker et al., Molecular electronic structure theory, p. 556-558
!!    and specifically equation (11.8.3).
!!
!!    The calculation is performed with the following steps:
!!
!!    - Loops through the bits of string and checks whether the corresponding orbital is occupied.
!!    - Compute how many strings have been produced by in order to fill this orbital (offset)
!!    - Accumulates the result to obtain the address for the specific string
!!
!!
      implicit none
!
      integer, intent(in)  :: n_orbitals, n_electrons
      integer, intent(in)  :: string
      integer :: address
!
      integer :: offset, orbital, n_remaining_electrons
!
      address               = 1
      n_remaining_electrons = n_electrons
!
      do orbital = n_orbitals, 1, -1
!
         if (n_remaining_electrons == 0 .or. &
             orbital == n_remaining_electrons) exit
!
         if (orbital_occupied_in_string(string, orbital)) then
!
            offset = get_address_offset(orbital, n_remaining_electrons)
!
            address = address + offset
            n_remaining_electrons = n_remaining_electrons - 1
!
         end if
!
      end do
!
   end function generate_address_for_string
!
!
   function get_address_offset(orbital, n_electrons) result(offset)
!!
!!    Get address offset
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    See Helgaker et al., Molecular electronic structure theory, p. 556-558
!!
!!    Returns the (diagonal) arc weight ending at vertex Y^1_o,n_e
!!
      implicit none
!
      integer, intent(in) :: orbital, n_electrons
      integer :: offset
!
      if (n_electrons > orbital) &
         call output%error_msg('too many electrons for the given number of orbitals')
!
      if (orbital == 0 .or. n_electrons == 0) &
         call output%error_msg('either 0 orbitals or electrons in ' &
                               //'call to get_address_offset')
!
      if (orbital > 1) then
!
         offset = binomial(orbital - 1, n_electrons)
!
      else
!
         offset = 0
!
      endif
!
   end function get_address_offset
!
!
   pure function sign_create_p_annihilate_q(p, q, string) result(sign_)
!!
!!    Sign create p annihilate q
!!    Written by Enrico Ronca, 2020
!!
!!    Calculates sign of
!!
!!       a_p^dagger a_q | k >,
!!
!!    where k is given by the string.
!!    p is assumed to be unoccupied in the string, and q to be occupied.
!!
!!    The sign is given by counting the number of occupied orbitals
!!    between p and q. Also see Helgaker et al., Molecular electronic structure theory,
!!    p. 8, equation (1.3.17)
!!
      implicit none
!
      integer, intent(in)  :: p, q
      integer, intent(in) :: string
!
      integer :: sign_
      integer :: exponent
!
      exponent = mod(n_occupations_between_p_and_q(p, q, string), 2)
      sign_ = (-1)**exponent
!
   end function sign_create_p_annihilate_q
!
!
   pure function sign_create_or_annihilate_p(p, string) result(sign_)
!!
!!    Sign create or annihilate p
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2022
!!
!!    Calculates sign of
!!
!!       a_p | k > or a_p^dagger | k >
!!
!!    where k is given by the string.
!!    p is assumed to be occupied in the string.
!!
      implicit none
!
      integer, intent(in)  :: p
      integer, intent(in) :: string
!
      integer :: sign_
      integer :: exponent_
!
      exponent_ = poppar(ishft(string, -p))
      sign_ = (-1)**exponent_
!
   end function sign_create_or_annihilate_p
!
!
   pure function n_occupations_between_p_and_q(p, q, string) result(n_occupations)
!!
!!    N occupations between p and q
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      integer, intent(in)  :: p, q
      integer, intent(in) :: string
!
      integer :: n_occupations
!
      integer :: temporary_string
!
!     Make string with ones at positions between, not including, p and q
      if (p > q) then
!
         temporary_string = ishft(1, p - 1) - ishft(1, q)
!
      else
!
         temporary_string = ishft(1, q - 1) - ishft(1, p)
!
      end if
!
!     Count the ones appearing in both string and temporary_string
      n_occupations = popcnt(iand(string, temporary_string))
!
   end function n_occupations_between_p_and_q
!
!
   module subroutine construct_excitation_maps(wf, n_o, n_v, n_strings, strings, excitation_maps)
!!
!!    Construct excitation maps
!!    Written by Enrico Ronca, 2020
!!
!!    Generates a table for the strings relationship in terms of a
!!    creation-annihilating operator pair.
!!
!!    Constructs the excitation_maps array for a set of strings;
!!
!!    The array is n_strings-by-(n_o + n_o * n_v)-by-4 and stores information
!!    about the new string generated by acting with a_p^dagger a_q on an existing string.
!!
!!    For a given string, we have an (n_o + n_o * n_v)-by-4 array,
!!    where the first n_o rows contain info about occupied-occupied excitations
!!
!!       a_p^dagger a_p | k >  (does not change the string)
!!
!!    and the next n_o * n_v rows contain info occupied-virtual exciations
!!
!!       a_p^dagger a_q | k >  (p virtual and q occupied, does generate a new string)
!!
!!    For a given string (i.e., a given | k >) and a given operation (a_p^dagger a_p or a_p^dagger a_q),
!!    there are four columns containing:
!!
!!       1. The virtual orbital of the generated string
!!       2. The occupied orbital of the generated string
!!       3. The address of the generated string
!!       4. The prefactor/sign of the a_p^dagger a_q operation (+ or -)
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      integer, intent(in) :: n_o, n_v
!
      integer, intent(in) :: n_strings
      integer, dimension(n_strings), intent(in) :: strings
!
      integer, dimension(n_strings, n_o*n_v + n_o, 4), intent(out) :: excitation_maps
!
      integer :: s, i, a, ai
      integer :: address
      integer :: string
      integer :: new_string
!
      integer, dimension(wf%n_mo) :: occupied, virtual
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct excitation maps', pl='n')
      call timer%turn_on
!
!$omp parallel do private(s, i, a, ai, occupied, virtual, string, new_string, address) &
!$omp shared(n_o, n_v, n_strings, strings, excitation_maps)
      do s = 1, n_strings
!
         string = strings(s)
!
         call wf%get_occupied_and_virtuals_in_string(string, occupied, virtual)
!
         do i = 1, n_o
!
            excitation_maps(s, i, 1) = occupied(i)
            excitation_maps(s, i, 2) = occupied(i)
            excitation_maps(s, i, 3) = s
            excitation_maps(s, i, 4) = 1
!
         end do
!
         do i = 1, n_o
            do a = 1, n_v
!
               new_string = create_p(annihilate_p(string, p=occupied(i)), p=virtual(a))
!
               address = generate_address_for_string(new_string, wf%n_mo, n_o)
!
               ai = n_v*(i - 1) + a
!
               excitation_maps(s, ai + n_o, 1) = virtual(a)
               excitation_maps(s, ai + n_o, 2) = occupied(i)
               excitation_maps(s, ai + n_o, 3) = address
               excitation_maps(s, ai + n_o, 4) = sign_create_p_annihilate_q(virtual(a), occupied(i), string)
!
            end do
         end do
!
      end do
!$omp end parallel do
!
      call timer%turn_off
!
   end subroutine construct_excitation_maps
!
!
   module subroutine get_occupied_and_virtuals_in_string(wf, string, occupied, virtual)
!!
!!    Get occupied and virtuals in string
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    On exit, "occupied" is a list of occupied orbitals in the string (with zeros at the end),
!!    and similarly for "virtual".
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      integer, intent(in) :: string
!
      integer, dimension(wf%n_mo), intent(out) :: occupied
      integer, dimension(wf%n_mo), intent(out) :: virtual
!
      integer :: p, count_o, count_v
!
      occupied = 0
      virtual = 0
!
      count_o = 0
      count_v = 0
!
      do p = 1, wf%n_mo
!
         if (orbital_occupied_in_string(string, p)) then
!
            count_o = count_o + 1
            occupied(count_o) = p
!
         else
!
            count_v = count_v + 1
            virtual(count_v) = p
!
         end if
!
      end do
!
   end subroutine get_occupied_and_virtuals_in_string
!
!
   module subroutine calculate_occupied_indices_in_strings(n_strings, n_orbitals, n_o, strings, occupied_indices)
!!
!!    Calculate occupied indices in strings
!!    Written by Enrico Ronca, 2020
!!
!!    Constructs occupied_indices (n_strings-by-n_o)
!!
!!    For each string, occupied_list contains the orbital indices that are occupied in the string
!!
      implicit none
!
      integer,                            intent(in)  :: n_strings
      integer,                            intent(in)  :: n_orbitals
      integer, dimension(n_strings),      intent(in)  :: strings
!
      integer, intent(in) :: n_o
      integer, dimension(n_strings, n_o), intent(out) :: occupied_indices
!
      integer :: s
!
      integer :: string, p, count_o
!
      occupied_indices = 0
!
      do s = 1, n_strings
!
         string = strings(s)
!
         count_o = 0
!
         do p = 1, n_orbitals
!
            if (orbital_occupied_in_string(string, p)) then
!
               count_o = count_o + 1
!
               occupied_indices(s, count_o) = p
!
            endif
         enddo
!
      enddo
!
   end subroutine calculate_occupied_indices_in_strings
!
!
   module subroutine construct_destruction_strings_and_signs(wf, &
                                                             destruction_strings, destruction_signs, &
                                                             n_strings, strings, n_electrons)
!!
!!    Construct destruction strings and signs
!!    Written by Alexander C. Paul, Eirik F. Kjønstad, Sarai D. Folkestad,
!!    and Enrico Ronca, 2022
!!
!!
!!    Constructs destruction_strings and destruction_signs, which contains, for each string (|k>)
!!    and each MO (p), information about the string generated from creation:
!!
!!        a_p | k >.
!!
!!    For each string and MO:
!!
!!       destruction_strings: address of the generated determinant (zero if p is unoccupied)
!!       destruction_signs:   sign of generated determinant (zero if p is unoccupied), see eq. 1.2.3 MEST
!!
!
      use array_utilities, only: zero_array_int, zero_array
!
      implicit none
!
      class(fci), intent(in) :: wf
!
      integer, intent(in) :: n_strings, n_electrons
!
      integer,  dimension(n_strings), intent(in) :: strings

      real(dp),  dimension(n_strings, wf%n_mo), intent(out) :: destruction_signs
      integer,  dimension(n_strings, wf%n_mo), intent(out)  :: destruction_strings
!
      integer :: p, s, string
      integer :: new_string, new_address
!
      call zero_array_int(destruction_strings, n_strings*wf%n_mo)
      call zero_array(destruction_signs, n_strings*wf%n_mo)
!
      do s = 1, n_strings
!
         string = strings(s)
!
         do p = 1, wf%n_mo
!
            if (orbital_occupied_in_string(string, p)) then
!
               new_string = annihilate_p(string, p)
!
               new_address = generate_address_for_string(new_string, wf%n_mo, n_electrons - 1)
!
               destruction_strings(s, p) = new_address
               destruction_signs(s, p) = real(sign_create_or_annihilate_p(p, string), kind=dp)
!
            end if
!
         end do
      end do
!
   end subroutine construct_destruction_strings_and_signs
!
!
   module subroutine construct_creation_strings_and_signs(wf, &
                                                          creation_strings, creation_signs, &
                                                          n_strings, strings, n_electrons)
!!
!!    Construct creation strings and signs
!!    Written by Alexander C. Paul, Eirik F. Kjønstad, Sarai D. Folkestad,
!!    and Enrico Ronca, 2022
!!
!!    Constructs creation_strings and creation_signs, which contains, for each string (|k>)
!!    and each MO (p), information about the string generated from creation:
!!
!!        a^dagger_p | k >.
!!
!!    For each string and MO:
!!
!!       creation_strings: address of the generated determinant (zero if p is occupied)
!!       creation_signs:   sign of generated determinant (zero if p is occupied), see eq. 1.2.3 MEST
!!
      use array_utilities, only: zero_array_int, zero_array
!
      implicit none
!
      class(fci), intent(in) :: wf
!
      integer, intent(in) :: n_strings, n_electrons
!
      integer, dimension(n_strings), intent(in) :: strings
!
      real(dp), dimension(n_strings, wf%n_mo), intent(out) :: creation_signs
      integer, dimension(n_strings, wf%n_mo), intent(out)  :: creation_strings
!
      integer :: p, s, string
      integer :: new_string, new_address
!
      call zero_array_int(creation_strings, n_strings*wf%n_mo)
      call zero_array(creation_signs, n_strings*wf%n_mo)
!
      do s = 1, n_strings
!
         string = strings(s)
!
         do p = 1, wf%n_mo
!
            if (.not. orbital_occupied_in_string(string, p)) then
!
               new_string = create_p(string, p)
!
               new_address = generate_address_for_string(new_string, wf%n_mo, n_electrons + 1)
!
               creation_strings(s, p) = new_address
               creation_signs(s, p)   = real(sign_create_or_annihilate_p(p, string), kind=dp)
!
            end if
!
         end do
      end do
!
   end subroutine construct_creation_strings_and_signs
!
!
end submodule strings_fci
