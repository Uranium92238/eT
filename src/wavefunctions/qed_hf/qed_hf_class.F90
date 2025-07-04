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
module qed_hf_class
!
!!
!!    QED Hartree-Fock (QED-HF) class module
!!    Written by Tor S. Haugland, May 2021
!!
!!    Solves Hamiltonians of the form
!!       H = H_e + w b^+ b + g_b (b + b^+) + g_s^2/w
!!    using the mean-field approximation. The solution is
!!    written in terms of a HF state and a coherent state.
!!
!!    One example is g_b = g_s = lambda * sqrt(w/2) * (d dot e)
!!    where lambda=sqrt(1/epsilon V), d is the dipole and e
!!    is the transversal polarization vector.
!!
!!    Check out the QED tool for more details.
!!
!
   use hf_class, only: hf
   use qed_tool_class, only: qed_tool
!
   use parameters
   use memory_manager_class, only: mem
   use global_out, only: output
   use qed_ao_eri_getter_class, only: qed_ao_eri_getter
!
   implicit none
!
   type, extends(hf) :: qed_hf
!
      type(qed_tool), allocatable :: qed
!
   contains
!
!     Overwritten routines
!
      procedure :: prepare                         => prepare_qed_hf
      procedure :: cleanup                         => cleanup_qed_hf
      procedure :: print_summary                   => print_summary_qed_hf
!
      procedure :: construct_mo_fock_fc_term       => construct_mo_fock_fc_term_qed_hf
      procedure :: update_fock_and_energy          => update_fock_and_energy_qed_hf
!
      procedure :: get_ao_h                        => get_ao_h_qed_hf
!
   end type qed_hf
!
!
   interface qed_hf
!
      procedure :: new_qed_hf
!
   end interface qed_hf
!
!
contains
!
!
   function new_qed_hf() result(wf)
!!
!!    New qed hf
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    Adapted by Tor S. Haugland, Jan 2020
!!
      use citation_printer_class,   only : eT_citations
!
      implicit none
!
      type(qed_hf) :: wf
!
      wf%name_ = 'qed-rhf'
!
      call wf%read_settings()
      call wf%print_banner()
!
      call eT_citations%add('QED-HF')
!
   end function new_qed_hf
!
!
   subroutine prepare_qed_hf(wf, centers, embedding, charge)
!!
!!    Prepare
!!    Written by Tor S. Haugland, Jan 2020
!!
!!    Initializes files, writes the restart file used for consistency checks
!!    and constructs screening vectors
!!
      use atomic_center_class, only: atomic_center
!
      implicit none
!
      class(qed_hf) :: wf
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers
      integer, intent(in), optional :: charge
!
      logical, intent(in), optional :: embedding
!
      if (present(embedding)) then
         if (embedding) call output%error_msg("Embedding not supported in QED-HF")
      end if
!
      call wf%general_preparations(centers, embedding, charge)
!
      wf%qed = qed_tool()
      call wf%qed%initialize(wf%ao)
!
      wf%eri_getter = qed_ao_eri_getter(wf%ao, wf%qed)
!
   end subroutine prepare_qed_hf
!
!
   subroutine cleanup_qed_hf(wf)
!!
!!    Cleanup
!!    Written by Tor S. Haugland, Jan 2020
!!
      implicit none
!
      class(qed_hf) :: wf
!
      call wf%qed%cleanup()
      call wf%hf%cleanup()
!
   end subroutine cleanup_qed_hf
!
!
   subroutine update_fock_and_energy_qed_hf(wf, cumulative)
!!
!!    Update Fock and energy
!!    Written by Tommaso Giovannini, July 2019
!!    Modified by Tor S. Haugland, June 2021
!!       Added qed%update_coherent_state
!!
!!    Call either cumulative or non-cumulative updating depending on options.
!!    In QED-HF, the photon state coherent state is also updated.
!!
      implicit none
!
      class(qed_hf), intent(inout) :: wf
      logical, intent(in) :: cumulative
!
      real(dp), dimension(:,:), allocatable :: h
!
      if (wf%qed%is_optimizing_photons) &
         call wf%qed%update_coherent_state(wf%ao_density)
!
      if (.not. cumulative) then
!
         call wf%update_G_non_cumulative()
!
      else
!
         call wf%update_G_cumulative()
!
      endif
!
!     Set the two-electron part
!
      call dcopy(wf%ao%n**2, wf%ao_G, 1, wf%ao_fock, 1)
      call dscal(wf%ao%n**2, half, wf%ao_fock, 1)
!
!     Add the one-electron part
!
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
      call daxpy(wf%ao%n**2, one, h, 1, wf%ao_fock, 1)
!
!     Energy
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, h)
!
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!
   end subroutine update_fock_and_energy_qed_hf
!
!
   subroutine print_summary_qed_hf(wf)
!!
!!    Print Summary
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(qed_hf), intent(inout) :: wf
!
      call wf%hf%print_summary()
      call wf%qed%print_summary()
!
   end subroutine print_summary_qed_hf
!
!
   subroutine construct_mo_fock_fc_term_qed_hf(wf, mo_fc_fock)
!!
!!    Calculate MO Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the frozen core contribution to
!!    the fock matrix
!!
!!       F_pq = (2 g_wxyz D^FC_yz - g_wyzx D^FC_xy) C_pw C_qx
!!
!!    in preparation of FC-CC
!!
      implicit none
!
      class(qed_hf) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: mo_fc_fock
!
      real(dp), dimension(:,:), allocatable :: D
      real(dp), dimension(:,:), allocatable :: ao_F_fc
!
      if (.not. wf%qed%is_basis_complete) &
         call output%error_msg("Frozen core not suppported (use quadrupole oei)")
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
!     Add frozen core contribution
!
      call dgemm('N', 'T',                      &
                  wf%ao%n,                      &
                  wf%ao%n,                      &
                  wf%n_frozen_core_orbitals,    &
                  two,                          &
                  wf%orbital_coefficients_fc,   &
                  wf%ao%n,                      &
                  wf%orbital_coefficients_fc,   &
                  wf%ao%n,                      &
                  zero,                         &
                  D,                            &
                  wf%ao%n)
!
      call daxpy(wf%ao%n**2, half, D, 1, wf%frozen_CCT, 1)
!
!     Construct the frozen core contribution to the active Fock matrix
!
      call mem%alloc(ao_F_fc, wf%ao%n, wf%ao%n)
!
      call wf%get_G(D, ao_F_fc)
!
      call wf%mo_transform(ao_F_fc, mo_fc_fock)
!
      call mem%dealloc(ao_F_fc, wf%ao%n, wf%ao%n)
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      call wf%destruct_orbital_coefficients_fc()
!
   end subroutine construct_mo_fock_fc_term_qed_hf
!
!
   subroutine get_ao_h_qed_hf(wf, h)
!!
!!    Get AO h
!!    Written by Tor S. Haugland, Oct 2021
!!
!!    Get the Hamiltonian one-electron integrals h in
!!    the AO basis
!!       h_wx = h^e_wx + h^QED_wx
!!
      implicit none
!
      class(qed_hf), intent(in) :: wf
      real(dp), dimension(wf%ao%n, wf%ao%n) :: h
!
      call dcopy(wf%ao%n**2, wf%ao%h, 1, h, 1)
      call wf%qed%add_ao_h_contribution(h)
!
   end subroutine get_ao_h_qed_hf
!
!
end module qed_hf_class
