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
submodule (fci_class) contract_fci
!
!!
!! FCI contraction submodule
!!
!! Performs the transformation of a vector c by the Hamiltonian
!!
!!    sigma = H c
!!
!! The transformation is performed in line with the N-resolution method of
!! Ch. 11.8 in Molecular Electronic Structure Theory (MEST), by Helgaker et al.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine hamiltonian_transformation_no_spin_symmetry_fci(wf, c, sigma)
!!
!!    Hamiltonian transformation no spin symmetry
!!    Written by Enrico Ronca, 2020
!!
!!    Performs the transformation of a vector c by the Hamiltonian
!!
!!       sigma = H c
!!
!!    On exit, the sigma vector is stored in c.
!!
!!    Constructs the matrices:
!!
!!       D_rsKaKb = sum_JaJb <KaKb|E_rs|JaJb> C_JaJb
!!
!!       G_pqKaKb = 1/2 sum_rs h2e_pqrs D_rsKaKb
!!
!!    And the sigma vector is then given as
!!
!!       sigma_IaIb = sum_KaKb sum_pq <IaIb|E_pq|KaKb> G_pqKaKb
!!
!!    See MEST, eqs. (11.8.13) - (11.8.15)
!!
!!    Note that the one-electron contribution to sigma (eq. 11.8.10)
!!    is included through the modifications of the integrals (h2e = effective_2e_hamiltonian).
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: c
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(out) :: sigma
!
      real(dp), dimension(:,:,:,:), allocatable :: D_rsKaKb
      real(dp), dimension(:,:,:,:), allocatable :: G_pqKaKb
!
      type(timings), allocatable :: timer
!
      timer = timings('Hamiltonian transformation FCI', pl='m')
      call timer%turn_on
!
      call mem%alloc(D_rsKaKb, wf%n_mo, wf%n_mo, wf%n_alpha_strings, wf%n_beta_strings)
      call wf%construct_D(c, D_rsKaKb)
!
      call mem%alloc(G_pqKaKb, wf%n_mo, wf%n_mo, wf%n_alpha_strings, wf%n_beta_strings)
      call dgemm('N', 'N',                       &
                  wf%n_mo**2,                    &
                  wf%n_determinants,             &
                  wf%n_mo**2,                    &
                  half,                          &
                  wf%effective_2e_hamiltonian,   &
                  wf%n_mo**2,                    &
                  D_rsKaKb,                      &
                  wf%n_mo**2,                    &
                  zero,                          &
                  G_pqKaKb,                      &
                  wf%n_mo**2)
!
      call mem%dealloc(D_rsKaKb, wf%n_mo, wf%n_mo, wf%n_alpha_strings, wf%n_beta_strings)
!
      call wf%construct_sigma(G_pqKaKb, sigma)
!
      call daxpy(wf%n_determinants, wf%get_nuclear_repulsion(), c, 1, sigma, 1)
!
      call mem%dealloc(G_pqKaKb, wf%n_mo, wf%n_mo, wf%n_alpha_strings, wf%n_beta_strings)
!
      call timer%turn_off
!
   end subroutine hamiltonian_transformation_no_spin_symmetry_fci
!
!
   module subroutine construct_D(wf, c, D_rsKaKb)
!!
!!    Construct D
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    Constructs the D matrix defined in eqs. (11.8.3) and (11.8.16) of MEST
!!
!!       D_rsKaKb = sum_JaJb <KaKb|E_rs|JaJb> C_JaJb
!!                = sum_Ja <Ka|E_pq^a|Ja> C_JaKb + sum_Jb <Kb|E_pq^b|Jb> C_KaJb
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(fci), intent(in) :: wf
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in)  :: c
      real(dp), dimension(wf%n_mo, wf%n_mo, &
                          wf%n_alpha_strings, wf%n_beta_strings), intent(out) :: D_rsKaKb
!
      integer :: Ja, Jb, rs, Ka, r, s, Kb
!
      real(dp) :: sign_
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct D matrix', pl='v')
      call timer%turn_on
!
      call zero_array(D_rsKaKb, wf%n_mo**2 * wf%n_determinants)
!
!      D_rsKaKb = sum_Ja <Ka|E_rs^a|Ja> C_JaKb
!
!$omp parallel do private(Kb, Ja, rs, Ka, r, s, sign_)
      do Kb = 1, wf%n_beta_strings
         do Ja = 1, wf%n_alpha_strings
            do rs = 1, wf%n_alpha_excitations
!
               r  = wf%excitation_maps_alpha(Ja, rs, 1)
               s  = wf%excitation_maps_alpha(Ja, rs, 2)
               Ka = wf%excitation_maps_alpha(Ja, rs, 3)
!
               sign_ = real(wf%excitation_maps_alpha(Ja, rs, 4), kind=dp)
!
               D_rsKaKb(r, s, Ka, Kb) = D_rsKaKb(r, s, Ka, Kb) + sign_ * c(Ja, Kb)
!
            enddo
         end do
      end do
!$omp end parallel do
!
!      D_rsKaKb += sum_Jb <Kb|E_rs^b|Jb> C_KaJb
!
!$omp parallel do private(Ka, Jb, rs, Kb, r, s, sign_)
      do Ka = 1, wf%n_alpha_strings
         do Jb = 1, wf%n_beta_strings
            do rs = 1, wf%n_beta_excitations

               r = wf%excitation_maps_beta(Jb, rs, 1)
               s = wf%excitation_maps_beta(Jb, rs, 2)
               Kb = wf%excitation_maps_beta(Jb, rs, 3)

               sign_ = real(wf%excitation_maps_beta(Jb, rs, 4), kind=dp)

               D_rsKaKb(r, s, Ka, Kb) = D_rsKaKb(r, s, Ka, Kb) + sign_ * c(Ka, Jb)
!
            end do
         end do
      end do
!$omp end parallel do
!
      call timer%turn_off
!
   end subroutine construct_D
!
!
   module subroutine construct_sigma(wf, G_pqKaKb, sigma)
!!
!!    Construct sigma
!!    Written by Enrico Ronca and Sarai D. Folkestad, 2020-2021
!!
!!    Constructs the transformed vector sigma, defined in eq. (11.8.15) of
!!    MEST.
!!
!!       sigma_IaIb = sum_KaKb sum_pq <IaIb|E_pq|KaKb> G_pqKaKb
!!                  = sum_Ka sum_pq <Ia|E_pq^a|Ka> G_pqKaIb
!!                  + sum_Kb sum_pq <Ib|E_pq^b|Kb> G_pqIaKb
!!
!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(fci), intent(in)  :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, &
                          wf%n_alpha_strings, wf%n_beta_strings), intent(in)  :: G_pqKaKb
!
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(out) :: sigma
!
      integer :: Ia, Ib, Ka, Kb, p, q, pq
      real(dp) :: sign_
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct sigma', pl='v')
      call timer%turn_on
!
      call zero_array(sigma, wf%n_determinants)
!
!     sigma_IaIb = sum_KaKb sum_pq <Ia|E_pq^a|Ka> G_pqKaIb
!
!$omp parallel do private(Ia, Ib, pq, p, q, sign_, Ka)
      do Ib = 1, wf%n_beta_strings
         do Ka = 1, wf%n_alpha_strings
            do pq = 1, wf%n_alpha_excitations
!
               p  = wf%excitation_maps_alpha(Ka, pq, 1)
               q  = wf%excitation_maps_alpha(Ka, pq, 2)
               Ia = wf%excitation_maps_alpha(Ka, pq, 3)
!
               sign_ = real(wf%excitation_maps_alpha(Ka, pq, 4), kind=dp)
!
               sigma(Ia,Ib) = sigma(Ia,Ib) + sign_ * G_pqKaKb(p, q, Ka, Ib)
!
            enddo
         end do
      end do
!$omp end parallel do
!
!     sigma_IaIb += sum_Kb sum_pq <Ib|E_pq^b|Kb> G_pqIaKb
!
!$omp parallel do private(Ia, Ib, pq, p, q, sign_, Kb)
      do Ia = 1, wf%n_alpha_strings
         do Kb = 1, wf%n_beta_strings
            do pq = 1, wf%n_beta_excitations

               p  = wf%excitation_maps_beta(Kb, pq, 1)
               q  = wf%excitation_maps_beta(Kb, pq, 2)
               Ib = wf%excitation_maps_beta(Kb, pq, 3)
!
               sign_ = real(wf%excitation_maps_beta(Kb, pq, 4), kind=dp)
!
               sigma(Ia, Ib) = sigma(Ia, Ib) + sign_ * G_pqKaKb(p, q, Ia, Kb)
!
            end do
!
         end do
      end do
!$omp end parallel do
!
      call timer%turn_off
!
   end subroutine construct_sigma
!
!
end submodule contract_fci
