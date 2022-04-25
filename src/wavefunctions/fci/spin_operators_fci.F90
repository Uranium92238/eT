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
submodule (fci_class) spin_operators_fci
!
!!
!! FCI Spin Operators submodule
!!
!! Gathers routines to calculate Spin expectation values at the FCI level
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
   module function sm_sp_expectation_value(wf, ci_coefficients) result(sm_sp)
!!
!!    S- S+ expectation value
!!    Written by Alexander C. Paul, Eirik F. Kjønstad, Sarai D. Folkestad,
!!    and Enrico Ronca, 2022
!!
!!    Computes:
!!
!!       <FCI|S- S+ |FCI> = sum_IaIb sum_JaJb C_IaIb C_JaJb <IaIb|S- S+ |JaJb>
!!
!!    We apply S+:
!!
!!       S+|JaJb> = sum_p c^dagger_pa c_pb A_Ja A_Jb |vac> ,
!!                = sum_p A_Ka A_Kb |vac> P^(Jpa)P^(Jpb)
!!
!!    where
!!
!!       c^dagger_pa - creates an alpha electron in orbital p
!!       c_pb        - annihilates a beta electron in orbital p
!!       A_Ja/A_Jb   - Are strings of ordered alpha/beta creation operators defining |JaJb>
!!       A_Ka        - Is the string of ordered alpha creation operators that in addition to
!!                     the operators in A_Ja includes c^dagger_pa
!!       A_Kb        - Is the string of ordered beta creation operators that includes the
!!                     operators in A_Jb EXCEPT c^dagger_pb
!!       P^(Jpa)     - Is the change in parity resulting in reordering c^dagger_pa A_Ja -> A_Ka
!!       P^(Jpb)     - Is the change in parity resulting in reordering c_pb A_Jb -> A_Kb
!!
!!    Using
!!
!!       (S+|FCI>)^dagger = <FCI|S-
!!
!!    we evaluate the expectation value as a dot-product of S+|FCI> vectors
!!
      use omp_lib, only: omp_get_max_threads, omp_get_thread_num
      use array_utilities, only: zero_array
!
      implicit none
!
      class(fci), intent(in) :: wf
      real(dp) :: sm_sp
!
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(in) :: ci_coefficients
!
      real(dp), dimension(:,:), allocatable :: creation_alpha_signs, destruction_beta_signs
      integer,  dimension(:,:), allocatable :: creation_alpha_strings, destruction_beta_strings
!
      real(dp),  dimension(:,:,:), allocatable :: sp_ci_coefficients
!
      real(dp) :: parity_JaJb, ddot
!
      integer :: Ja, Jb, p, Ka, Kb, n_alpha_plus_1_strings, n_beta_minus_1_strings
      integer :: n_threads = 1, thread_n = 1
!
!$    n_threads = omp_get_max_threads()
!
      n_alpha_plus_1_strings = binomial(wf%n_mo, wf%n_alpha + 1)
      n_beta_minus_1_strings = binomial(wf%n_mo, wf%n_beta  - 1)
!
      call mem%alloc(sp_ci_coefficients, n_alpha_plus_1_strings, &
                                         n_beta_minus_1_strings, &
                                         n_threads)
!
      call mem%alloc(creation_alpha_strings, wf%n_alpha_strings, wf%n_mo)
      call mem%alloc(creation_alpha_signs,   wf%n_alpha_strings, wf%n_mo)
!
      call mem%alloc(destruction_beta_strings, wf%n_beta_strings, wf%n_mo)
      call mem%alloc(destruction_beta_signs,   wf%n_beta_strings, wf%n_mo)
!
!     Strings (and corresponding signs) obtained by creating an alpha electron
      call wf%construct_creation_strings_and_signs(creation_alpha_strings, &
                                                   creation_alpha_signs,   &
                                                   wf%n_alpha_strings,     &
                                                   wf%alpha_strings,       &
                                                   wf%n_alpha)
!
!     Strings (and corresponding signs) obtained by destructing a beta electron
      call wf%construct_destruction_strings_and_signs(destruction_beta_strings,  &
                                                      destruction_beta_signs,    &
                                                      wf%n_beta_strings,         &
                                                      wf%beta_strings,           &
                                                      wf%n_beta)
!
      call zero_array(sp_ci_coefficients, n_alpha_plus_1_strings * &
                                          n_beta_minus_1_strings * &
                                          n_threads)
!
!$omp parallel do &
!$omp shared(sp_ci_coefficients, ci_coefficients, creation_alpha_strings, &
!$omp destruction_beta_strings , creation_alpha_signs, destruction_beta_signs) &
!$omp private(p, Jb, Ja, Ka, Kb, parity_JaJb, thread_n)
      do Jb = 1, wf%n_beta_strings
!
!$       thread_n = omp_get_thread_num() + 1
!
         do p = 1, wf%n_mo
!
            if (destruction_beta_signs(Jb, p) == 0) cycle
!
            do Ja = 1, wf%n_alpha_strings
!
               if (creation_alpha_signs(Ja, p) == 0) cycle
!
               Ka = creation_alpha_strings(Ja, p)
               Kb = destruction_beta_strings(Jb, p)
!
               parity_JaJb = destruction_beta_signs(Jb, p)*creation_alpha_signs(Ja, p)
!
               sp_ci_coefficients(Ka, Kb, thread_n) = sp_ci_coefficients(Ka, Kb, thread_n) &
                                                    + ci_coefficients(Ja, Jb)*parity_JaJb
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      do thread_n = 2, n_threads
!
         call daxpy(n_alpha_plus_1_strings*n_beta_minus_1_strings, one, &
                     sp_ci_coefficients(:,:,thread_n), 1, sp_ci_coefficients(:,:,1), 1)
!
      enddo 
!
      sm_sp = ddot(n_alpha_plus_1_strings*n_beta_minus_1_strings, &
                   sp_ci_coefficients, 1, sp_ci_coefficients, 1)
!
      call mem%dealloc(sp_ci_coefficients, n_alpha_plus_1_strings, &
                                           n_beta_minus_1_strings, &
                                           n_threads)
!
      call mem%dealloc(creation_alpha_strings, wf%n_alpha_strings, wf%n_mo)
      call mem%dealloc(creation_alpha_signs, wf%n_alpha_strings, wf%n_mo)
!
      call mem%dealloc(destruction_beta_strings, wf%n_beta_strings, wf%n_mo)
      call mem%dealloc(destruction_beta_signs, wf%n_beta_strings, wf%n_mo)
!
   end function sm_sp_expectation_value
!
!
   module function calculate_spin_squared_fci(wf, state) result(s2)
!!
!!    Calculate S^2
!!    Written by Alexander C. Paul, Eirik F. Kjønstad, Sarai D. Folkestad,
!!    and Enrico Ronca, 2022
!!
!!    Calculates <fci|S^2|fci> via MEST eqs. (2.4.36) and (2.4.37):
!!
!!    <fci|S^2|fci> = <fci| (S-S+ + S_z (S_z + 1)) |fci>
!!                  = <fci| S-S+ |fci> + 1/4[(N_alpha-N_beta)^2 + 2(N_alpha-N_beta)]
!!
      implicit none
!
      class(fci), intent(in) :: wf
!
      integer, intent(in) :: state
!
      real(dp) :: s2, sz_contribution, sm_sp_contribution
!
      sm_sp_contribution = wf%sm_sp_expectation_value(wf%ci_coefficients(:,:,state))

      sz_contribution = quarter*(two*real(wf%n_alpha - wf%n_beta, dp)  &
                                   + real(wf%n_alpha - wf%n_beta, dp)**2)
!
      s2 = sz_contribution + sm_sp_contribution
!
   end function calculate_spin_squared_fci
!
!
   module subroutine calculate_spin_multiplicities_fci(wf, multiplicities)
!!
!!    Calculate spin multiplicities
!!    Written by Enrico Ronca, 2022
!!
!!    Calculates multiplicity = 2S + 1, where
!!
!!       S = sqrt(S^2 + 1/4) - 1/2
!!
      implicit none
!
      class(fci), intent(in) :: wf
!
      real(dp), dimension(wf%n_states), intent(out) :: multiplicities
!
      integer :: state
!
      real(dp) :: S, S2
!
      type(timings), allocatable :: timer
!
      if (wf%n_alpha == 0 .or. wf%n_beta == 0) then
         multiplicities = wf%ao%get_n_electrons() + one
         return
      end if
!
      if (wf%n_v == 0) then
         multiplicities = one
         return
      end if
!
      timer = timings('Calculate spin multiplicities', pl='n')
      call timer%turn_on
!
      do state = 1, wf%n_states
!
         S2 = wf%calculate_spin_squared(state)
!
         S = sqrt(S2 + quarter) - half
!
         multiplicities(state) = two * S + one
!
      end do
!
      call timer%turn_off
!
   end subroutine calculate_spin_multiplicities_fci
!
!
end submodule spin_operators_fci
