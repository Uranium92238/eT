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
submodule (hf_class) mo_hf
!
!!
!!    MO submodule 
!!
!!    Collects the routines to construct MO quantities in HF theory.
!!    These are used e.g. in the MO-based based MO-SCF-DIIS.
!!
!
   implicit none
!
!
contains
!
!
!
   module function get_max_roothan_hall_mo_gradient_hf(wf) result(max_gradient)
!!
!!    Get max of Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       E^(1) = - 4 F_ai,
!!
!!    and returns the maximum absolute value of E^(1).
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp) :: max_gradient
!
      real(dp), dimension(:,:), allocatable :: gradient
!
      call mem%alloc(gradient, wf%n_v, wf%n_o)
!
      call wf%get_roothan_hall_mo_gradient(gradient)
!
      max_gradient = get_abs_max(gradient, wf%n_o*wf%n_v)
!
      call mem%dealloc(gradient, wf%n_v, wf%n_o)
!
   end function get_max_roothan_hall_mo_gradient_hf
!
!
   module subroutine get_roothan_hall_mo_gradient_hf(wf, G)
!!
!!    Get packed Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       G_ai = E^(1)_ai = - 4 F_ai
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o), intent(out) :: G
!
      integer :: a, i
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            G(a, i) = - four*wf%mo_fock(i, wf%n_o + a)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_roothan_hall_mo_gradient_hf
!
!
   module subroutine roothan_hall_update_orbitals_mo_hf(wf)
!!
!!    Roothan-Hall update of orbitals
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Update orbitals after a Roothan-Hall step in the MO basis
!!
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: C_old
!
      call wf%do_roothan_hall_mo(wf%mo_fock, wf%orbital_energies)
!
!     Update of the orbital coefficients in the new mo basis
!
      call mem%alloc(C_old, wf%n_ao, wf%n_mo)
      call dcopy(wf%n_ao*wf%n_mo, wf%orbital_coefficients, 1, C_old, 1)
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  one,                       &
                  C_old,                     &
                  wf%n_ao,                   &
                  wf%W_mo_update,            &
                  wf%n_mo,                   &
                  zero,                      &
                  wf%orbital_coefficients,   &
                  wf%n_ao)
!
      call mem%dealloc(C_old, wf%n_ao, wf%n_mo)
!
   end subroutine roothan_hall_update_orbitals_mo_hf
!
!
   module subroutine do_roothan_hall_mo_hf(wf, F, e)
!!
!!    Do Roothan-Hall
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Solves the equation 
!!
!!       F W = W e, 
!!
!!    where W (called mo_basis_update) is the 
!!    transformation matrix from old to new MO basis
!!
!!       C_new = C_old W
!!
!!    'F' : Fock matrix in MO basis
!!
!!    'e' : orbital energies 
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)    :: F
      real(dp), dimension(wf%n_mo),          intent(inout) :: e
!
      real(dp), dimension(:), allocatable :: work
!
      integer :: info
!
!     Solve F W = W e
!
      info = 0
!
      call mem%alloc(work, 4*wf%n_mo)
!
      call dcopy(wf%n_mo**2, F, 1, wf%W_mo_update, 1)
!
      call dsyev('V', 'L',              &
                  wf%n_mo,              &
                  wf%W_mo_update,       &
                  wf%n_mo,              &
                  e,                    &
                  work,                 &
                  4*(wf%n_mo),          &
                  info)
!
      call mem%dealloc(work, 4*wf%n_mo)
!
      if (info .ne. 0) then
!
         call output%error_msg('Could not solve Roothan-Hall equations.' // &
                              ' "Dsyev" finished with info: (i0)', ints=[info])
!
      endif
!
   end subroutine do_roothan_hall_mo_hf
!
!
   module subroutine update_fock_and_energy_mo_hf(wf, prev_ao_density)
!!
!!    Update Fock and energy
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Constructs Fock matrix in the AO basis
!!    and calculates the energy.
!!
!!    Transforms the active Fock matrix to the new MO basis.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
      real(dp), dimension(:,:), allocatable :: G
!
      type(timings) :: timer
!
      timer = timings('AO Fock construction', pl='normal')
      call timer%turn_on()
!
      if (present(prev_ao_density)) then 
!
!        Construct the two electron part of the Fock matrix (G),
!        and add the contribution to the Fock matrix
!
         call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
!
         call mem%alloc(G, wf%n_ao, wf%n_ao)
!
         call wf%construct_ao_G(wf%ao_density,                    &
                                G,                                &
                                C_screening=(wf%name_ == 'mlhf'))
!
         call daxpy(wf%n_ao**2, one, G, 1, wf%ao_fock, 1)
!
         call mem%dealloc(G, wf%n_ao, wf%n_ao)
!
         call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)         
!
      else
!
!        AO fock construction and energy calculation
!
!        Construct the two electron part of the Fock matrix (G),
!        and add the contribution to the Fock matrix
!
         call wf%construct_ao_G(wf%ao_density,                    &
                                wf%ao_fock,                       &
                                C_screening=(wf%name_ == 'mlhf'))
!
!        Add the one-electron part
!
         call daxpy(wf%n_ao**2, one, wf%ao_h, 1, wf%ao_fock, 1)
!
      endif
!
      call timer%turn_off()
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, wf%ao_h)
!
!     Transformation of the AO fock in the MO basis
!
      call wf%mo_transform(wf%ao_fock, wf%mo_fock)
!
   end subroutine update_fock_and_energy_mo_hf
!
!
   module subroutine prepare_for_roothan_hall_mo_hf(wf)
!!
!!    Construct idempotent density and prepare for Roothan-Hall
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Constructs the ao Fock matrix and
!!    performs a Roothan-Hall step to get the
!!    initial idempotent density. 
!!
!!    The routine also prints the number of 
!!    electrons and the energy of the initial guess.
!!
!!    These steps are performed in 
!!
!!       prepare_for_roothan_hall_hf,
!!
!!    this routine also allocates the MO 
!!    quantities mo_fock and W_mo_update and
!!    initializew W_mo_update to the identity matrix
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%prepare_for_roothan_hall()
!
!     Allocate active mo specific arrays
!     and construct them
!
      call wf%initialize_W_mo_update()
      call wf%initialize_mo_fock()
!
      call identity_array(wf%W_mo_update, wf%n_mo)
!
   end subroutine prepare_for_roothan_hall_mo_hf
!
!
   module subroutine read_for_scf_restart_mo_hf(wf)
!!
!!    Read for SCF restart
!!    Written by Sarai D. Folkestad and Linda Goletto, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%read_orbital_coefficients()
      call wf%update_ao_density()
      call wf%read_orbital_energies()
!
!     Allocate active mo specific arrays
!     and construct them
!
      call wf%initialize_W_mo_update()
      call wf%initialize_mo_fock()
!
      call identity_array(wf%W_mo_update, wf%n_mo)
!
   end subroutine read_for_scf_restart_mo_hf
!
!
   module subroutine construct_mo_fock_hf(wf)
!!
!!    Construct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the MO Fock matrix F_pq using the current AO
!!    Fock and the orbital coefficients.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%ao_fock,                &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  X,                         & ! X = F^ao C
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  wf%mo_fock,                & ! F = C^T F^ao C
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_ao, wf%n_mo)
!
   end subroutine construct_mo_fock_hf
!
!
   module subroutine get_fock_ov_hf(wf, F)
!!
!!    Get HF equations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the occupied-virtual block of the Fock MO matrix,
!!    and returns the result in the array F.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v)   :: F ! F_ia
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_v)
!
      call dgemm('N', 'N',                                  &
                  wf%n_ao,                                  &
                  wf%n_v,                                   &
                  wf%n_ao,                                  &
                  one,                                      &
                  wf%ao_fock,                               &
                  wf%n_ao,                                  &
                  wf%orbital_coefficients(1, wf%n_o + 1),   &
                  wf%n_ao,                                  &
                  zero,                                     &
                  X,                                        &
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_o,                    &
                  wf%n_v,                    &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  F,                         &
                  wf%n_o)
!
      call mem%dealloc(X, wf%n_ao, wf%n_v)
!
   end subroutine get_fock_ov_hf
!
!
end submodule mo_hf
