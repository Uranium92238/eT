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
submodule (fci_class) properties_fci
!
!!
!! Properties FCI submodule
!!
!! Contains routines to calculate properties in FCI
!!
!
   implicit none
!
!
contains
!
!
   module function get_electronic_dipole_fci(wf) result(mu_electronic)
!!
!!    Get electronic dipole
!!    Written by Alexander C. Paul, May 2022
!!
      use array_utilities, only: symmetric_sandwich_right_transposition
!
      implicit none
!
      class(fci), intent(in) :: wf
!
      real(dp), dimension(3) :: mu_electronic
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: ao_mu_pqk, mu_pqk
!
      call mem%alloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
      call mem%alloc(ao_mu_pqk, wf%ao%n, wf%ao%n, 3)
!
      call wf%ao%get_oei('dipole', ao_mu_pqk)
!
      do k = 1, 3
!
         call wf%mo_transform(ao_mu_pqk(:,:,k), mu_pqk(:,:,k))
         mu_electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k), wf%density)
!
      enddo
!
      call mem%dealloc(ao_mu_pqk, wf%ao%n, wf%ao%n, 3)
      call mem%dealloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
!
   end function get_electronic_dipole_fci
!
!
   module function get_electronic_quadrupole_fci(wf) result(q_electronic)
!!
!!    Get electronic quadrupole
!!    Written by Alexander C. Paul, May 2022
!!
      use array_utilities, only: symmetric_sandwich_right_transposition
!
      implicit none
!
      class(fci), intent(in) :: wf
!
      real(dp), dimension(6) :: q_electronic
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: ao_q_pqk, q_pqk
!
      call mem%alloc(q_pqk, wf%n_mo, wf%n_mo, 6)
      call mem%alloc(ao_q_pqk, wf%ao%n, wf%ao%n, 6)
!
      call wf%ao%get_oei('quadrupole', ao_q_pqk)
!
      do k = 1, 6
!
         call wf%mo_transform(ao_q_pqk(:,:,k), q_pqk(:,:,k))
         q_electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k), wf%density)
!
      enddo
!
      call mem%dealloc(q_pqk, wf%n_mo, wf%n_mo, 6)
      call mem%dealloc(ao_q_pqk, wf%ao%n, wf%ao%n, 6)
!
   end function get_electronic_quadrupole_fci
!
!
   module function calculate_expectation_value_fci(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Calculate the expectation value of a one-electron operator A
!!
!!       < A > = < FCI| A |FCI > = sum_pq A_pq D_pq
!!
      implicit none
!
      class(fci), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
      real(dp) :: expectation_value, ddot
!
      expectation_value = ddot(wf%n_mo**2, A, 1, density, 1)
!
   end function calculate_expectation_value_fci
!
!
   module subroutine construct_gs_density_fci(wf)
!!
!!    Construct GS Density
!!    Written by Alexander C. Paul and Sarai D. Folkestad, May 2022
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call wf%construct_density(wf%ci_coefficients(:,:,1), &
                                wf%ci_coefficients(:,:,1), &
                                wf%density)
!
   end subroutine construct_gs_density_fci
!
!
   module subroutine construct_density(wf, fci_vector_1, fci_vector_2, density)
!!
!!    Construct Density
!!    Written by Alexander C. Paul and Sarai D. Folkestad, May 2022
!!
!!    Constructs the density matrix defined in eqs. (11.8.3) and (11.8.16) of MEST
!!
!!       density = <fci_1|E_pq|fci_2>
!!               = sum_JaKaKb C_KaKb <Ka|E_pq^a|Ja> C_JaKb
!!               + sum_JbKbKa C_KaKb <Kb|E_pq^b|Jb> C_KaJb
!!
!!    where the alpha determinants Ja and Ka are related by a single excitation
!!    and analogously for the beta determinants Jb and Kb
!!
      use array_utilities, only:zero_array
!
      implicit none
!
      class(fci), intent(in) :: wf
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: fci_vector_1
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: fci_vector_2
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct FCI density matrix', pl='v')
      call timer%turn_on
!
      call zero_array(density, wf%n_mo**2)
      call wf%add_alpha_density(density, fci_vector_1, fci_vector_2)
      call wf%add_beta_density(density, fci_vector_1, fci_vector_2)
!
      call timer%turn_off
!
   end subroutine construct_density
!
!
   module subroutine add_alpha_density(wf, density, fci_vector_1, fci_vector_2)
!!
!!    Add alpha Density
!!    Written by Alexander C. Paul and Sarai D. Folkestad, May 2022
!!
!!    Adds the density matrix for alpha spin to 'density'
!!
!!       density += sum_JaKaKb C_KaKb <Ka|E_pq^a|Ja> C_JaKb
!!
!!    where the alpha determinants Ja and Ka are related by a single excitation
!!
      use array_utilities, only:zero_array
      use omp_lib
!
      implicit none
!
      class(fci), intent(in) :: wf
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: fci_vector_1
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: fci_vector_2
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(:,:,:), allocatable :: D
!
      integer :: Ja, Ka, Kb, pq, p, q
      integer :: n_threads = 1, thread_n = 1
      real(dp) :: sign_
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(D, wf%n_mo, wf%n_mo, n_threads)
      call zero_array(D, wf%n_mo**2*n_threads)
!
!$omp parallel do private(Ka, Kb, Ja, pq, p, q, sign_, thread_n) shared(D, wf)
      do Ja = 1, wf%n_alpha_strings
!
!$       thread_n = omp_get_thread_num() + 1
!
         do pq = 1, wf%n_alpha_excitations
!
               p  = wf%excitation_maps_alpha(Ja, pq, 1)
               q  = wf%excitation_maps_alpha(Ja, pq, 2)
               Ka = wf%excitation_maps_alpha(Ja, pq, 3)
!
               sign_ = real(wf%excitation_maps_alpha(Ja, pq, 4), kind=dp)
!
            do Kb = 1, wf%n_beta_strings
!
               D(p, q, thread_n) = D(p, q, thread_n) &
                                 + sign_ * fci_vector_1(Ka, Kb) * fci_vector_2(Ja, Kb)
!
            enddo
         end do
      end do
!$omp end parallel do
!
      do thread_n = 1, n_threads
         call daxpy(wf%n_mo**2, one, D(:,:,thread_n), 1, density, 1)
      end do
!
      call mem%dealloc(D, wf%n_mo, wf%n_mo, n_threads)
!
   end subroutine add_alpha_density
!
!
   module subroutine add_beta_density(wf, density, fci_vector_1, fci_vector_2)
!!
!!    Add beta Density
!!    Written by Alexander C. Paul and Sarai D. Folkestad, May 2022
!!
!!    Adds the density matrix for beta spin to 'density'
!!
!!       density += sum_JbKbKa C_KaKb <Kb|E_pq^b|Jb> C_KaJb
!!
!!    where the beta determinants Jb and Kb are related by a single excitation
!!
      use array_utilities, only:zero_array
      use omp_lib
!
      implicit none
!
      class(fci), intent(in) :: wf
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: fci_vector_1
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: fci_vector_2
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(:,:,:), allocatable :: D
!
      integer :: Jb, Kb, Ka, pq, p, q
      integer :: n_threads = 1, thread_n = 1
      real(dp) :: sign_
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(D, wf%n_mo, wf%n_mo, n_threads)
      call zero_array(D, wf%n_mo**2*n_threads)
!
!$omp parallel do private(Ka, Kb, Jb, pq, p, q, sign_, thread_n) shared(D, wf)
      do Jb = 1, wf%n_beta_strings
!
!$       thread_n = omp_get_thread_num() + 1
!
         do pq = 1, wf%n_beta_excitations
!
            p = wf%excitation_maps_beta(Jb, pq, 1)
            q = wf%excitation_maps_beta(Jb, pq, 2)
            Kb = wf%excitation_maps_beta(Jb, pq, 3)
            sign_ = real(wf%excitation_maps_beta(Jb, pq, 4), kind=dp)
!
            do Ka = 1, wf%n_alpha_strings
!
               D(p, q, thread_n) = D(p, q, thread_n) &
                                 + sign_ * fci_vector_1(Ka, Kb) * fci_vector_2(Ka, Jb)
!
            end do
         end do
      end do
!$omp end parallel do
!
      do thread_n = 1, n_threads
         call daxpy(wf%n_mo**2, one, D(:,:,thread_n), 1, density, 1)
      end do
!
      call mem%dealloc(D, wf%n_mo, wf%n_mo, n_threads)
!
   end subroutine add_beta_density
!
!
end submodule properties_fci
