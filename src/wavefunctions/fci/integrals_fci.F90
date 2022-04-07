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
submodule (fci_class) fci_integrals
!
!!
!! FCI Integrals submodule
!!
!! Gathers routines to prepare integrals in MO basis for a FCI calculation
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_g_pqrs(wf)
!!
!!    Construct g_pqrs
!!    Written by Enrico Ronca and Tor S. Haugland, 2022
!!
!!    Performs Cholesky decomposition of the electron repulsion integral matrix.
!!    and tranform the integrals in the mo basis.
!!
      use eri_cd_class, only: eri_cd
      use eri_tool_class, only: eri_tool
      use eri_cholesky_memory_class,   only: eri_cholesky_memory
!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      type(eri_cd), allocatable :: eri_cholesky_solver
      type(eri_cholesky_memory), allocatable :: L_mo
      type(eri_tool), allocatable :: eri
!
      integer :: n_J
!
      call output%printf('v', 'Doing Cholesky decomposition of the ERIs.')
!
      eri_cholesky_solver = eri_cd(wf%ao)
      call eri_cholesky_solver%run(wf%ao)
      call eri_cholesky_solver%diagonal_test(wf%ao)
!
      L_mo = eri_cholesky_memory()
      n_J = eri_cholesky_solver%get_n_cholesky()
      call L_mo%initialize(n_J, 1, [wf%n_mo])
!
      eri = eri_tool(L_mo)
!
      call eri_cholesky_solver%construct_cholesky_mo_vectors(wf%ao, wf%ao%n, wf%n_mo, &
                                                             wf%orbital_coefficients, L_mo)
!
      call eri_cholesky_solver%cleanup()
!
      call eri%get(wf%g_pqrs,    &
                   1, wf%n_mo,   &
                   1, wf%n_mo,   &
                   1, wf%n_mo,   &
                   1, wf%n_mo)
!
      deallocate(L_mo)
      deallocate(eri)
!
   end subroutine construct_g_pqrs
!
!
   module subroutine construct_h_pq(wf)
!!
!!    Construct h_pq
!!    Written by Enrico Ronca and Eirik F. Kjønstad, 2022
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: h_ao
!
      call mem%alloc(h_ao, wf%ao%n, wf%ao%n)
!
      call wf%ao%get_oei('hamiltonian', h_ao)
!
      call wf%mo_transform(h_ao, wf%h_pq)
!
      call mem%dealloc(h_ao, wf%ao%n, wf%ao%n)
!
   end subroutine construct_h_pq
!
!
   module subroutine construct_hamiltonian_integrals(wf)
!!
!!    Construct Hamiltonian integrals
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct Hamiltonian integrals', pl='n')
      call timer%turn_on
!
      call wf%construct_h_pq()
      call wf%construct_g_pqrs()
      call wf%construct_effective_2e_hamiltonian()
!
      call timer%turn_off
!
   end subroutine construct_hamiltonian_integrals
!
!
   module subroutine construct_effective_2e_hamiltonian(wf)
!!
!!    Construct effective 2-electron Hamiltonian
!!    Written by Enrico Ronca, 2020
!!
!!    Modify two-electron integrals to include the effective one-electron integrals
!!
!!       k_pq = (h_pq - 1/2 sum_r g_prrq)
!!
!!    defined in eq. (11.8.8) in MEST
!!
!!    The one-electron integrals modify the ERIs according to
!!
!!       h2e_rrpq =  (g_pqrr + 1/N_e k_pq)
!!       h2e_pqrr =  (g_rrpq + 1/N_e k_pq)
!!
!!    this is to avoid separate computation of the one-electron contribution to sigma.
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      real(dp) :: n_electrons
!
      integer :: p,q, r
!
      real(dp), dimension(:,:), allocatable :: X_pq
!
      n_electrons = real(wf%n_alpha + wf%n_beta, kind=dp)
!
      call dcopy(wf%n_mo**4, wf%g_pqrs, 1, wf%effective_2e_hamiltonian, 1)
!
      call mem%alloc(X_pq, wf%n_mo, wf%n_mo)
!
      call dcopy(wf%n_mo**2, wf%h_pq, 1, X_pq, 1)
!
!$omp parallel do private(p, q, r)
      do q = 1, wf%n_mo
         do p = 1, wf%n_mo
            do r = 1, wf%n_mo
!
               X_pq(p, q) = X_pq(p, q) - half*wf%g_pqrs(p, r, r, q)
!
            end do
         end do
      end do
!$omp end parallel do
!
      call dscal(wf%n_mo**2, one/n_electrons, X_pq, 1)
!
!$omp parallel do private(p, q, r)
      do q = 1, wf%n_mo
         do p = 1, wf%n_mo
            do r = 1, wf%n_mo
!
               wf%effective_2e_hamiltonian(r,r,p,q) = &
                     wf%effective_2e_hamiltonian(r,r,p,q) + X_pq(p,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(p, q, r)
      do r = 1, wf%n_mo
         do q = 1, wf%n_mo
            do p = 1, wf%n_mo
!
               wf%effective_2e_hamiltonian(p,q,r,r) = &
                  wf%effective_2e_hamiltonian(p,q,r,r) + X_pq(p,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_pq, wf%n_mo, wf%n_mo)
!
   end subroutine construct_effective_2e_hamiltonian
!
!
end submodule fci_integrals
