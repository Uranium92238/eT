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
module rohf_class
!
!!
!!    Restricted open-shell Hartree-Fock (ROHF) class module
!!    Written by Sarai D. Folkestad 2020
!!
!!    See Tsuchimochi and Scuseria, J. Chem. Phys. 133, 141102 (2010)
!!
!

   use parameters
   use cuhf_class, only: cuhf
!
   use memory_manager_class, only: mem
   use global_out,           only: output
   use global_in,            only: input
!
   implicit none
!
   type, extends(cuhf) :: rohf
!
      character(len=200), private :: coupling_parameters = 'guest-saunders'
!
      real(dp), dimension(3), private :: A = [half, half, half]
      real(dp), dimension(3), private :: B = [half, half, half]
!
   contains
!
      procedure, private :: make_rohf_orbitals        => make_rohf_orbitals_rohf
      procedure :: finalize_gs                        => finalize_gs_rohf
      procedure :: read_settings                      => read_settings_rohf
      procedure :: read_rohf_settings                 => read_rohf_settings_rohf
      procedure :: write_orbital_info                 => write_orbital_info_rohf
!
   end type rohf
!
   interface rohf
!
      procedure :: new_rohf
!
   end interface rohf
!
contains
!
!
   function new_rohf() result(wf)
!!
!!    New rohf
!!    Written by Sarai D. Folkestad 2020
!!
      implicit none
!
      type(rohf) :: wf
!
      wf%name_ = 'rohf'
      wf%cumulative_fock = .false.
      wf%fractional_uniform_valence = .false.
!
      call wf%read_settings()
!
      call wf%print_banner()
!
   end function new_rohf
!
!
   subroutine read_settings_rohf(wf)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Designed to be overwritten by descendants.
!!
      implicit none
!
      class(rohf) :: wf
!
      call wf%read_hf_settings()
      call wf%read_uhf_settings()
      call wf%read_rohf_settings()
!
   end subroutine read_settings_rohf
!
!
   subroutine read_rohf_settings_rohf(wf)
!!
!!    Read rohf settings
!!    Written by Sarai D. Folkestad, Nov 2021
!!
!!    Reads settings specific to the wavefunction.
!!
      implicit none
!
      class(rohf) :: wf
!
      call input%get_keyword('rohf coupling parameters', 'solver scf', wf%coupling_parameters)
!
      if (trim(wf%coupling_parameters) == 'guest-saunders') then
!
         wf%A = [half, half, half]
         wf%B = [half, half, half]
!
      else if (trim(wf%coupling_parameters) == 'mcweeny-diercksen') then
!
         wf%A = [one/three, one/three, two/three]
         wf%B = [two/three, one/three, one/three]
!
      else if (trim(wf%coupling_parameters) == 'faegri-manne') then
!
         wf%A = [half, one, half]
         wf%B = [half, zero, half]
!
      else
!
         call output%error_msg('did not recognize the ROHF coupling parameters')
!
      endif
!
      call output%warning_msg('ROHF orbital energies are not uniquely defined, &
                           & and may not fulfill Koopmans theorem')
!
   end subroutine read_rohf_settings_rohf
!
!
   subroutine make_rohf_orbitals_rohf(wf)
!!
!!    Make ROHF orbitals
!!    Written by Sarai D. Folkestad, Nov 2021
!!
!!    Makes ROHF orbitals from converged CUHF alpha and beta orbitals
!!    Given the parameters A and B (not uniquely defined)
!!
!!    Diagonalizas the block diagonal matrix
!!
!!          R_xx = A_x * F_alpha + B_x * F_beta
!!
!!    x = c - core
!!      = o - open
!!      = v - virtual
!!
!!    See Tsuchimochi and Scuseria, J. Chem. Phys. 133, 141102 (2010) or
!!    B. N. Plakhutin, E. V. Gorelik and N. N. Breslavskaya, J. Chem. Phys. 125, 204110 (2006)
!!
!!    NOTE: the off-diagonal blocks of R are assumed zero, since the
!!    CUHF equations have been solved prior to the call to this routine
!!
      use array_utilities, only: symmetric_sandwich, symmetric_sandwich_right_transposition
      use array_utilities, only: generalized_diagonalization_symmetric
      implicit none
!
      class(rohf), intent(inout) :: wf
!
      real(dp), dimension(:), allocatable :: occupation_numbers
      real(dp), dimension(:,:), allocatable :: C_no_basis, F_no_a, F_no_b, R, SC
      real(dp), dimension(:,:), allocatable :: F, S
!
      integer:: i, j
!
      call mem%alloc(occupation_numbers, wf%n_mo)
      call mem%alloc(C_no_basis, wf%ao%n, wf%n_mo)
!
      call wf%construct_natural_occupation_number_basis(C_no_basis, occupation_numbers)
!
      call mem%alloc(F_no_a, wf%n_mo, wf%n_mo)
      call symmetric_sandwich(F_no_a, &
                              wf%ao_fock_a,               &
                              C_no_basis, &
                              wf%ao%n, wf%n_mo)
!
      call mem%alloc(F_no_b, wf%n_mo, wf%n_mo)
      call symmetric_sandwich(F_no_b, &
                              wf%ao_fock_b,               &
                              C_no_basis, &
                              wf%ao%n, wf%n_mo)
!
      call mem%alloc(R, wf%n_mo, wf%n_mo, set_zero=.true.)
!
      do i = 1, wf%n_mo
         do j = 1, wf%n_mo
!
            if (abs(occupation_numbers(i) - one) .lt. 1.0d-8 .and. &
                abs(occupation_numbers(j) - one) .lt. 1.0d-8) then ! xx = cc
!
               R(i, j) = wf%A(1)*F_no_a(i, j) + wf%B(1)*F_no_b(i, j)
!
            else if (abs(occupation_numbers(i) - half) .lt. 1.0d-8 .and. &
                     abs(occupation_numbers(j) - half) .lt. 1.0d-8) then ! xx = oo
!
               R(i, j) = wf%A(2)*F_no_a(i, j) + wf%B(2)*F_no_b(i, j)
!
            else if (abs(occupation_numbers(i)) .lt. 1.0d-8 .and. &
                     abs(occupation_numbers(j)) .lt. 1.0d-8) then  ! xx = vv
!
               R(i, j) = wf%A(3)*F_no_a(i, j) + wf%B(3)*F_no_b(i, j)
!
            endif
!
         enddo
      enddo
!
      call mem%dealloc(F_no_a, wf%n_mo, wf%n_mo)
      call mem%dealloc(F_no_b, wf%n_mo, wf%n_mo)
!
!     Transform back to AO basis
      call mem%alloc(SC, wf%ao%n, wf%n_mo)
!
      call dgemm('n', 'n', &
                  wf%ao%n, &
                  wf%n_mo, &
                  wf%ao%n, &
                  one,     &
                  wf%ao%s, &
                  wf%ao%n, &
                  C_no_basis, &
                  wf%ao%n, &
                  zero,    &
                  SC,      &
                  wf%ao%n)
!
      call symmetric_sandwich_right_transposition(wf%ao_fock, &
                                                  R, &
                                                  SC, wf%n_mo, wf%ao%n)
!
      call mem%dealloc(SC, wf%ao%n, wf%n_mo)
      call mem%dealloc(R, wf%n_mo, wf%n_mo)!
      call mem%dealloc(occupation_numbers, wf%n_mo)
      call mem%dealloc(C_no_basis, wf%ao%n, wf%n_mo)
!
!     Diagonalize to obtain orbital coefficients
!
      call mem%alloc(F, wf%n_mo, wf%n_mo)
      call mem%alloc(S, wf%n_mo, wf%n_mo)
!
      call wf%ao_to_reduced_ao_transformation(F, wf%ao_fock)
      call wf%ao%get_reduced_ao_metric(S) ! RAO metric
!
      call generalized_diagonalization_symmetric(F, S, wf%n_mo, wf%orbital_energies)
      call wf%set_orbital_coefficients_from_reduced_ao_C(F, wf%orbital_coefficients)
!
      call mem%dealloc(F, wf%n_mo, wf%n_mo)
      call mem%dealloc(S, wf%n_mo, wf%n_mo)
!
   end subroutine make_rohf_orbitals_rohf
!
!
   subroutine finalize_gs_rohf(wf)
!!
!!    Finalize gs
!!    Written by Sarai D. Folkestad, Nov 2021
!!
!!    Prepares a single set of orbitals,
!!    common for alpha and beta electrons
!!
      implicit none
!
      class(rohf), intent(inout) :: wf
!
      call wf%make_rohf_orbitals()
!
   end subroutine finalize_gs_rohf
!
!
   subroutine write_orbital_info_rohf(wf)
!!
!!    Write orbital info
!!    Written by Alexander C. Paul Nov 2020
!!
      implicit none
!
      class(rohf), intent(inout) :: wf
!
      call wf%hf%write_orbital_info()
!
   end subroutine write_orbital_info_rohf
!
!
end module rohf_class
