!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
submodule (mlcc2_class) omega_mlcc2
!
!!
!!    Omega submodule
!!
!!    Based on the CC2 routines written by Sarai D. Folkestad 
!!    and Eirik F. Kjønstad
!!
!!    Routines to construct
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
   module subroutine construct_omega_mlcc2(wf, omega)
!!
!!    Construct omega
!!    Written by Sarai D. Folkestad, Jan 2019
!!
!!    Directs the construction of the omega vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      type(timings) :: timer
!
      timer = timings('omega MLCC2')
      call timer%turn_on()
!
      call zero_array(omega,wf%n_gs_amplitudes)
!
      call wf%omega_ccs_a1(omega)
!
      call wf%construct_x2()
      call wf%construct_u_aibj()
!
      call wf%omega_cc2_a1(omega, wf%n_cc2_o, wf%n_cc2_v, wf%first_cc2_o, wf%first_cc2_v, &
                           wf%last_cc2_o, wf%last_cc2_v)
      call wf%omega_cc2_b1(omega, wf%n_cc2_o, wf%n_cc2_v, wf%first_cc2_o, wf%first_cc2_v, &
                           wf%last_cc2_o, wf%last_cc2_v)
      call wf%omega_cc2_c1(omega, wf%n_cc2_o, wf%n_cc2_v, wf%first_cc2_o, wf%first_cc2_v)
!
      call timer%turn_off()
!
   end subroutine construct_omega_mlcc2
!
!
   module subroutine omega_cc2_a1_mlcc2(wf, omega, n_cc2_o, n_cc2_v, first_cc2_o, &
                                       first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Omega MLCC2 A1 term
!!    Adapted by Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_bcj u_bicj g_abjc
!!
!!    and adds it to the projection vector omega.
!!
!!    The term is calculated while batching over index a.
!!
!!    Index restrictions:
!!
!!       b, i, c, j : CC2 orbitals 
!!
!!       a : unrestricted
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
!
!     Local variables
!
      real(dp), dimension(:,:,:), allocatable :: L_Jcj, L_Jab, L_aJb, X_Jbi
!
      type(batching_index) :: batch_a
!
      integer :: req0, req1
!
      integer :: current_a_batch
!
      type(timings) :: timer
!  
      timer = timings('omega mlcc2 a1')
      call timer%turn_on()
!
!     Using L_Jjc_t1 =  L_Jjc_mo = L_Jcj_mo
!
      call mem%alloc(L_Jcj, wf%integrals%n_J, n_cc2_v, n_cc2_o)
      call wf%integrals%get_cholesky_mo(L_Jcj,              &
                                       wf%n_o + first_cc2_v,&
                                       wf%n_o + last_cc2_v, &
                                       first_cc2_o,         &
                                       last_cc2_o)
!
!     X_Jbi = u_bicj L_Jcj
!
      call mem%alloc(X_Jbi, wf%integrals%n_J, n_cc2_v, n_cc2_o)
!
      call dgemm('N', 'N',          &
                  wf%integrals%n_J, &
                  n_cc2_o*n_cc2_v,  &
                  n_cc2_o*n_cc2_v,  &
                  one,              &
                  L_Jcj,            &
                  wf%integrals%n_J, &
                  wf%u_aibj,        & ! u_cjbi
                  n_cc2_o*n_cc2_v,  &
                  zero,             &
                  X_Jbi,            &
                  wf%integrals%n_J)
!
      call mem%dealloc(L_Jcj, wf%integrals%n_J, n_cc2_v, n_cc2_o)
!
      req0 = 0
      req1 = 2*(n_cc2_v)*(wf%integrals%n_J)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(L_Jab, wf%integrals%n_J, batch_a%length, n_cc2_v)
         call wf%integrals%get_cholesky_t1(L_Jab,                 &
                                          wf%n_o + batch_a%first, &
                                          wf%n_o + batch_a%last,  & 
                                          wf%n_o + first_cc2_v,   &
                                          wf%n_o + last_cc2_v)
!
         call mem%alloc(L_aJb, batch_a%length, wf%integrals%n_J, n_cc2_v)
         call sort_123_to_213(L_Jab, L_aJb, wf%integrals%n_J, batch_a%length, n_cc2_v)
         call mem%dealloc(L_Jab, wf%integrals%n_J, batch_a%length, n_cc2_v)
!
         call dgemm('N','N',                          &
                     batch_a%length,                  &
                     n_cc2_o,                         &
                     wf%integrals%n_J*n_cc2_v,        &
                     one,                             &
                     L_aJb,                           & 
                     batch_a%length,                  &
                     X_Jbi,                           & 
                     wf%integrals%n_J*n_cc2_v,        &
                     one,                             &
                     omega(batch_a%first,first_cc2_o),& 
                     wf%n_v)
!
         call mem%dealloc(L_aJb, batch_a%length, wf%integrals%n_J, n_cc2_v)
!
      enddo ! batch_a
!
      call mem%dealloc(X_Jbi, wf%integrals%n_J, n_cc2_v, n_cc2_o)
!
      call timer%turn_off()
!
   end subroutine omega_cc2_a1_mlcc2
!
!
   module subroutine omega_cc2_b1_mlcc2(wf, omega, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Omega MLCC2 B1 term
!!    Adapted by Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_bkj g_kbji * u_ajbk,
!!
!!    and adds it to the omega vector
!!
!!    Index restrictions:
!!
!!       a, j, b, k : CC2 orbitals
!!
!!       i : unrestricted
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf

      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kbji
      real(dp), dimension(:,:,:,:), allocatable :: g_jbki
!
      type(timings) :: timer 
!  
      timer = timings('omega mlcc2 b1')
      call timer%turn_on()
!
!     g_kbji ordered as g_jbki
!
      call mem%alloc(g_kbji, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call wf%get_ovoo(g_kbji, &
                       first_cc2_o, last_cc2_o, &
                       first_cc2_v, last_cc2_v, &
                       first_cc2_o, last_cc2_o, &
                       1, wf%n_o)
!
      call mem%alloc(g_jbki, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call sort_1234_to_3214(g_kbji, g_jbki, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call mem%dealloc(g_kbji, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
!     omega_Ai += - sum_BKJ g_KB,Ji * u_AJ,BK
!
      call dgemm('N', 'N',                &
                  n_cc2_v,                &
                  wf%n_o,                 &
                  (n_cc2_o**2)*n_cc2_v,   &
                  -one,                   &
                  wf%u_aibj,              &! u_a_jbk
                  n_cc2_v,                &
                  g_jbki,                 &! g_jbk_i
                  (n_cc2_o**2)*n_cc2_v,   &
                  one,                    &
                  omega(first_cc2_v, 1),  & ! omega_a_i
                  wf%n_v)
!
      call mem%dealloc(g_jbki, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine omega_cc2_b1_mlcc2
!
!
   module subroutine omega_cc2_c1_mlcc2(wf, omega, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v)
!!
!!    Omega MLCC2 C1 term
!!    Adapted by Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the C1 term,
!!
!!       C1: sum_bj u_ai,bj * F_jb,
!!
!!    and adds it to the omega vector.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
      real(dp), dimension(:,:), allocatable :: F_bj, omega_ai
!
      integer :: b, j, a, i
!
      type(timings) :: timer
!  
      timer = timings('omega mlcc2 c1')
      call timer%turn_on()
!
!     Collect correct fock matrix elements
!
      call mem%alloc(F_bj, n_cc2_v, n_cc2_o)
!
!$omp parallel do private(b, j)
      do b = 1, n_cc2_v
         do j = 1, n_cc2_o
!
            F_bj(b,j) = wf%fock_ia(j + first_cc2_o - 1, b + first_cc2_v - 1)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(omega_ai, n_cc2_v, n_cc2_o)
!
      call dgemm('N','N',              &
                  (n_cc2_v)*(n_cc2_o), &
                  1,                   &
                  (n_cc2_v)*(n_cc2_o), &
                  one,                 &
                  wf%u_aibj,           & ! u_ai_bj
                  (n_cc2_v)*(n_cc2_o), &
                  F_bj,                &
                  (n_cc2_v)*(n_cc2_o), &
                  zero,                &
                  omega_ai,            &
                  (n_cc2_v)*(n_cc2_o))
!
      call mem%dealloc(F_bj, n_cc2_v, n_cc2_o)
!
      do a = 1, n_cc2_v
         do i = 1, n_cc2_o
!
            omega(a + first_cc2_v - 1, i + first_cc2_o - 1) = &
                     omega(a + first_cc2_v - 1, i + first_cc2_o - 1) + omega_ai(a, i)
                           
!
         enddo
      enddo
!
      call mem%dealloc(omega_ai, n_cc2_v, n_cc2_o)
!
      call timer%turn_off()
!
    end subroutine omega_cc2_c1_mlcc2
!
!
   module subroutine construct_omega_doubles_mlcc2(wf, omega2)
!!
!!    Construct Omega doubles
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the doubles part of omega for MLCC2
!!
!!    Note that this is not used to solve MLCC2 equations, 
!!    but it is only used for debug of Jacobian matrix
!!
!!    omega_aibj = 1/Δ_aibj g_aibj + e_aibj s_aibj 
!!
!!    Note that it is calculated in the biorthonormal basis!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_x2), intent(out) :: omega2
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer :: a, i, b, j, ai, bj, aibj, aiai
!
      call mem%alloc(g_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      call wf%get_vovo(g_aibj, &
                        wf%first_cc2_v, wf%last_cc2_v, &
                        wf%first_cc2_o, wf%last_cc2_o, &
                        wf%first_cc2_v, wf%last_cc2_v, &
                        wf%first_cc2_o, wf%last_cc2_o)
!
      call packin(omega2, g_aibj, wf%n_cc2_v*wf%n_cc2_o)
!
      call mem%dealloc(g_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      do a = 1, wf%n_cc2_v
         do i = 1, wf%n_cc2_o
            do b = 1, wf%n_cc2_v
               do j = 1, wf%n_cc2_o
!
                  ai = wf%n_cc2_v*(i - 1) + a
                  bj = wf%n_cc2_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai,bj) - 3)/2 + ai + bj
!
                  omega2(aibj) = omega2(aibj) + wf%x2(aibj)*&
                        (wf%orbital_energies(a + wf%first_cc2_v - 1 + wf%n_o)&
                        +wf%orbital_energies(b + wf%first_cc2_v - 1 + wf%n_o)&
                        -wf%orbital_energies(i + wf%first_cc2_o - 1)&
                        -wf%orbital_energies(j + wf%first_cc2_o - 1))
!
               enddo
            enddo
         enddo
      enddo
!
      do ai = 1, wf%n_cc2_v*wf%n_cc2_o
!
         aiai = ai*(ai - 3)/2 + ai*2
!
         omega2(aiai) = half*omega2(aiai)
!
      enddo
!
   end subroutine construct_omega_doubles_mlcc2
!
!
end submodule omega_mlcc2
