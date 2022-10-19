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
submodule (hf_class) tdhf_transformation
!
!!
!!    TDHF transformation submodule
!!
!!    Contains the RPA and Tamm-Dancoff (CIS/CCS) transformations
!!    And relevant routines for the calculation of these transformations.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine ao_G_with_non_symmetric_D_hf(wf, D, G)
!!
!!    AO G with non symmetric D
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Computes
!!
!!       G = (2g_wxyz - g_wzyx) D_zy
!!
!!    Where D is non-symmetric
!!
!
      use omp_lib
      use array_initialization, only: zero_array
!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: G
!
      integer :: s1s2, s1s2_packed, s1, s2, w, x, w_red, x_red
      integer :: s3s4, s3s4_packed, s3, s4, y, z, y_red, z_red
!
      real(dp) :: deg_12, deg_34, deg_
!
      real(dp), dimension(wf%ao%max_sh_size**4), target :: g_ABCD
      real(dp), dimension(:,:,:,:), pointer :: g_ABCD_p
!
      integer :: n_threads, thread
!
      real(dp), dimension(:,:,:), allocatable :: G_thread
!
      call zero_array(G, wf%ao%n**2)
!
      n_threads = 1
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(G_thread, wf%ao%n, wf%ao%n, n_threads, set_zero=.true.)
!
!$omp parallel do                                                                 &
!$omp private(s1, s2, s3, s4, s1s2, s1s2_packed, s3s4, deg_12, deg_34, deg_,&
!$omp w, x, y, z, w_red, x_red, y_red, z_red, g_ABCD, g_ABCD_p, thread, s3s4_packed) schedule(dynamic)
      do s1s2 = 1, wf%ao%n_sig_eri_shp
!
         thread = 1
!$       thread = thread + omp_get_thread_num()
!
         s1s2_packed = wf%ao%cs_eri_max_indices(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = wf%ao%cs_eri_max_indices(s1s2_packed, 2)
!
         if (s1 == s2) then
            deg_12 = half
         else
            deg_12 = one
         endif
!
         do s3s4 = 1, wf%ao%n_sig_eri_shp
!
            s3s4_packed = wf%ao%cs_eri_max_indices(s3s4, 3)
!
            s3 = wf%ao%cs_eri_max_indices(s3s4_packed, 1)
            s4 = wf%ao%cs_eri_max_indices(s3s4_packed, 2)
!
            if (s3 .gt. s1) cycle
!
            if (s3 == s4) then
               deg_34 = half
            else
               deg_34 = one
            endif
!
            if (s1 == s3) then
               deg_ = half
            else
               deg_ = one
            endif
!
            deg_ = deg_*deg_12*deg_34
!
            call wf%get_ao_g(g_ABCD, s1, s2, s3, s4)
!
            g_ABCD_p(1 : wf%ao%shells(s1)%length, 1 : wf%ao%shells(s2)%length, &
                  1 : wf%ao%shells(s3)%length, 1 : wf%ao%shells(s4)%length) &
                  => g_ABCD(1 : (wf%ao%shells(s1)%length)   &
                                 *(wf%ao%shells(s2)%length) &
                                 *(wf%ao%shells(s3)%length) &
                                 *(wf%ao%shells(s4)%length))
!
            do z = wf%ao%shells(s4)%first, wf%ao%shells(s4)%get_last()
               z_red = z - wf%ao%shells(s4)%first + 1
!
               do y = wf%ao%shells(s3)%first, wf%ao%shells(s3)%get_last()
                  y_red = y - wf%ao%shells(s3)%first + 1
!
                  do x = wf%ao%shells(s2)%first, wf%ao%shells(s2)%get_last()
                     x_red = x - wf%ao%shells(s2)%first + 1
!
                     do w = wf%ao%shells(s1)%first, wf%ao%shells(s1)%get_last()
                        w_red = w - wf%ao%shells(s1)%first + 1
!
!                       J-terms
!
                        G_thread(w,x,thread) = G_thread(w,x,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(y,z)
                        G_thread(x,w,thread) = G_thread(x,w,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(y,z)
                        G_thread(w,x,thread) = G_thread(w,x,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(z,y)
                        G_thread(x,w,thread) = G_thread(x,w,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(z,y)
!
                        G_thread(y,z,thread) = G_thread(y,z,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(w,x)
                        G_thread(z,y,thread) = G_thread(z,y,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(w,x)
                        G_thread(y,z,thread) = G_thread(y,z,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(x,w)
                        G_thread(z,y,thread) = G_thread(z,y,thread) &
                                             + deg_*two*g_ABCD_p(w_red, x_red, y_red, z_red)*D(x,w)
!
!                       K-terms
!
                        G_thread(w,z,thread) = G_thread(w,z,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(y, x)
                        G_thread(x,z,thread) = G_thread(x,z,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(y, w)
                        G_thread(w,y,thread) = G_thread(w,y,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(z, x)
                        G_thread(x,y,thread) = G_thread(x,y,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(z, w)
!
                        G_thread(z,w,thread) = G_thread(z,w,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(x, y)
                        G_thread(z,x,thread) = G_thread(z,x,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(w, y)
                        G_thread(y,w,thread) = G_thread(y,w,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(x, z)
                        G_thread(y,x,thread) = G_thread(y,x,thread) &
                                             - deg_*g_ABCD_p(w_red, x_red, y_red, z_red)*D(w, z)
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      do thread = 1, n_threads
!
         call daxpy(wf%ao%n**2, one, G_thread(1, 1, thread), 1, G, 1)
!
      enddo
!
      call mem%dealloc(G_thread, wf%ao%n, wf%ao%n, n_threads)
!
   end subroutine ao_G_with_non_symmetric_D_hf
!
!
   module subroutine construct_tdhf_D_hf(wf, D, b)
!!
!!    Construct TDHF D
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Constructs the D matrix
!!
!!       D_wx = b_ai C_wa C_xi
!!
!!    Where C is the orbital coefficient matrix,
!!    and the indices a and i denote virtual and occupied
!!    orbitals respectively
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: D
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: b
!
      real(dp), dimension(:,:), allocatable :: X_wi
!
      call mem%alloc(X_wi, wf%ao%n, wf%n_o)
!
      call dgemm('N', 'N',                                &
                  wf%ao%n,                                &
                  wf%n_o,                                 &
                  wf%n_v,                                 &
                  one,                                    &
                  wf%orbital_coefficients(:, wf%n_o + 1), &
                  wf%ao%n,                                &
                  b,                                      &
                  wf%n_v,                                 &
                  zero,                                   &
                  X_wi,                                   &
                  wf%ao%n)
!
      call dgemm('N', 'T',                 &
                  wf%ao%n,                 &
                  wf%ao%n,                 &
                  wf%n_o,                  &
                  one,                     &
                  X_wi,                    &
                  wf%ao%n,                 &
                  wf%orbital_coefficients, &
                  wf%ao%n,                 &
                  zero,                    &
                  D,                       &
                  wf%ao%n)
!
      call mem%dealloc(X_wi, wf%ao%n, wf%n_o)
!
   end subroutine construct_tdhf_D_hf
!
!
   module subroutine get_vo_block_hf(wf, X, X_ai)
!!
!!    Get vo block
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Constructs the virtual-occupied block of
!!    the matrix X^MO
!!
!!       X^MO_ai = X_wx C_wa C_xi
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: X
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: X_ai
!
      real(dp), dimension(:,:), allocatable :: X_wi
!
      call mem%alloc(X_wi, wf%ao%n, wf%n_o)
!
      call dgemm('N', 'N',                 &
                  wf%ao%n,                 &
                  wf%n_o,                  &
                  wf%ao%n,                 &
                  one,                     &
                  X,                       &
                  wf%ao%n,                 &
                  wf%orbital_coefficients, &
                  wf%ao%n,                 &
                  zero,                    &
                  X_wi,                    &
                  wf%ao%n)
!
      call dgemm('T', 'N',                                &
                  wf%n_v,                                 &
                  wf%n_o,                                 &
                  wf%ao%n,                                &
                  one,                                    &
                  wf%orbital_coefficients(:, wf%n_o + 1), &
                  wf%ao%n,                                &
                  X_wi,                                   &
                  wf%ao%n,                                &
                  zero,                                   &
                  X_ai,                                   &
                  wf%n_v)
!
      call mem%dealloc(X_wi, wf%ao%n, wf%n_o)
!
   end subroutine get_vo_block_hf
!
!
   module subroutine A_transformation_hf(wf, b, sigma)
!!
!!    A transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    sigma_ai = (2g_aijb - g_abji)b_bj + (e_a - e_i)b_ai
!!             = C_wa Cxi (2g_wxyz - g_wyzx) D^T_yz  + (e_a - e_i)b_ai
!!
!!    1. Construct D_wx = b_ai C_wa C_xi
!!
!!    2. Make G(D^T) for non-symmetric D
!!
!!    3. Construct vo-block of G in the MO basis
!!
!!    4. Add one electron terms
!!
      use array_utilities, only: transpose_
!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: b
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma
!
      real(dp), dimension(:,:), allocatable :: D, D_transpose, G, G_ai
!
      integer :: a, i
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
      call wf%construct_tdhf_D(D, b)
!
      call mem%alloc(D_transpose, wf%ao%n, wf%ao%n)
      call transpose_(D, D_transpose, wf%ao%n)
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      call mem%alloc(G, wf%ao%n, wf%ao%n)
      call wf%ao_G_with_non_symmetric_D(D_transpose, G)
      call mem%dealloc(D_transpose, wf%ao%n, wf%ao%n)
!
      call mem%alloc(G_ai, wf%n_v, wf%n_o)
      call wf%get_vo_block(G, G_ai)
      call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            sigma(a,i) = sigma(a,i) + G_ai(a,i) &
                        + (wf%orbital_energies(a + wf%n_o) - wf%orbital_energies(i))*b(a,i)
!
         enddo
      enddo
!
      call mem%dealloc(G_ai, wf%n_v, wf%n_o)
!
   end subroutine A_transformation_hf
!
!
   module subroutine B_transformation_hf(wf, b, sigma)
!!
!!    B transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    sigma_ai = (2g_aibj - g_ajbi)b_bj
!!             = C_wa Cxi (2g_wxyz - g_wyzx) D_yz
!!
!!    1. Construct D_wx = b_ai C_wa C_xi
!!
!!    2. Make G(D) for non-symmetric D
!!
!!    3. Construct vo-block of G in the MO basis
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: b
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma
!
      real(dp), dimension(:,:), allocatable :: D, G, G_ai
!
      integer :: a, i
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
      call wf%construct_tdhf_D(D, b)
!
      call mem%alloc(G, wf%ao%n, wf%ao%n)
      call wf%ao_G_with_non_symmetric_D(D, G)
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      call mem%alloc(G_ai, wf%n_v, wf%n_o)
      call wf%get_vo_block(G, G_ai)
      call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            sigma(a,i) = sigma(a,i) + G_ai(a,i)
!
         enddo
      enddo
!
      call mem%dealloc(G_ai, wf%n_v, wf%n_o)
!
   end subroutine B_transformation_hf
!
!
   module subroutine tamm_dancoff_transformation_hf(wf, b, sigma)
!!
!!    Tamm-Dancoff transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    sigma = Ab
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: b
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: sigma
!
      call zero_array(sigma, wf%n_v*wf%n_o)
      call wf%A_transformation(b, sigma)
!
   end subroutine tamm_dancoff_transformation_hf
!
!
   module subroutine rpa_transformation_hf(wf, b, sigma)
!!
!!    RPA transformation
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    RPA eigenvalue problem:
!!
!!       [  A  B ][X] = omega [X]
!!       [ -B -A ][Y]         [Y]
!!
!!    But equations are currently solved as:
!!
!!       (A + B)(A - B) [X + Y] = omega^2[X + Y]
!!
!!    Therefore we construct,
!!
!!       sigma = (A + B)(A - B) b = (A + B) c, for c = (A - B) b
!!
!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), target, intent(in) :: b
      real(dp), dimension(wf%n_v, wf%n_o), target, intent(out) :: sigma
!
      real(dp), dimension(:,:), allocatable :: c
!
      call zero_array(sigma, wf%n_v*wf%n_o)

      call mem%alloc(c, wf%n_v, wf%n_o, set_zero=.true.)
!
      call wf%A_transformation(b, c)
      call wf%B_transformation(b, c)
!
      call wf%A_transformation(c, sigma)

      call dscal(wf%n_v*wf%n_o, -one, c, 1)
      call wf%B_transformation(c, sigma)
!
      call mem%dealloc(c, wf%n_v, wf%n_o)
!
   end subroutine rpa_transformation_hf
!
!
   module subroutine tdhf_summary_hf(wf)
!!
!!    TDHF summary
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Prints the TDHF excitation energies.
!!
      implicit none
!
      class(hf) :: wf
!
      integer :: I
!
      call output%printf('m', '- TDHF excitation energies:', fs='(/t3,a)')
!
      call output%printf('m', 'Excitation energy', fs='(/t39,a)')
      call output%print_separator('m', 42, '-', fs='(t27,a)')
      call output%printf('m', ' State                (Hartree)             (eV)', fs='(t6,a)')
      call output%print_separator('m', 63, '-', fs='(t6,a)')
!
      do I = 1, wf%n_tdhf_vectors
!
         call output%printf('m', '(i4)             (f19.12)   (f19.12)',                  &
                            ints=[I], reals=[wf%tdhf_excitation_energies(I),              &
                            wf%tdhf_excitation_energies(I)*Hartree_to_eV], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('m', 63, '-', fs='(t6,a)')
!
      call output%printf('m', 'eV/Hartree (CODATA 2014): (f11.8)', &
                         reals=[Hartree_to_eV], fs='(t6,a)')
!
   end subroutine tdhf_summary_hf
!
!
   module subroutine get_rpa_preconditioner_hf(wf, preconditioner)
!!
!!    Get RPA preconditioner
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Preconditioner is square of the orbital energy differences
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_v*wf%n_o) :: preconditioner
!
      integer :: I
!
      call wf%get_orbital_differences(preconditioner, wf%n_v*wf%n_o)
!
      do I = 1, wf%n_v*wf%n_o
!
         preconditioner(I) = preconditioner(I)**2
!
      enddo
!
   end subroutine get_rpa_preconditioner_hf
!
!
   module subroutine get_rpa_response_preconditioner_hf(wf, preconditioner, frequency)
!!
!!    Get RPA preconditioner
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Preconditioner is square of the orbital energy differences
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(2*wf%n_v*wf%n_o), intent(out) :: preconditioner
!
      real(dp), intent(in) :: frequency
!
      integer :: I
!
      real(dp), dimension(:), allocatable :: eps
!
      call mem%alloc(eps, wf%n_v*wf%n_o)
!
      call wf%get_orbital_differences(eps, wf%n_v*wf%n_o)
!
      do I = 1, wf%n_v*wf%n_o
!
         preconditioner(I) = eps(I) - frequency
         preconditioner(I + wf%n_v*wf%n_o) = eps(I) + frequency
!
      enddo
      call mem%dealloc(eps, wf%n_v*wf%n_o)

!
   end subroutine get_rpa_response_preconditioner_hf
!
!
   module subroutine get_tdhf_start_vector_hf(wf, trial, c, restart)
!!
!!    Get TDHF start vectors
!!    Written by Sarai D. Folkestad, May 2021
!!
      use array_initialization, only: zero_array
      use array_utilities, only: get_n_lowest
!
      implicit none
!
      class(hf) :: wf
!
      integer :: trial
!
      real(dp), dimension(wf%n_v*wf%n_o) :: c
!
      logical :: restart
!
      real(dp), dimension(:), allocatable :: eps, lowest_eps
      integer, dimension(:),  allocatable :: start_indices
!
      if (restart) then
!
         if (wf%tdhf_files(trial)%exists()) then
!
            call output%printf('m', 'Requested restart. Reading in solution from file.', &
                               fs='(/t3,a)')
!
            call wf%read_tdhf_vector(c, trial)
            return
!
         endif
!
         call output%warning_msg('asked for restart of TDHF excitation energies, &
                                 &but no vectors found on disk')
!
      endif
!
      call mem%alloc(eps, wf%n_v*wf%n_o)
      call wf%get_orbital_differences(eps, wf%n_v*wf%n_o)
!
      call mem%alloc(start_indices, trial)
!
      call mem%alloc(lowest_eps, trial)
      call get_n_lowest(trial, wf%n_v*wf%n_o, eps, lowest_eps, start_indices)
      call mem%dealloc(lowest_eps, trial)
      call mem%dealloc(eps, wf%n_v*wf%n_o)
!
      call zero_array(c,  wf%n_v*wf%n_o)
!
      c(start_indices(trial)) = one
!
      call mem%dealloc(start_indices, trial)
!
   end subroutine get_tdhf_start_vector_hf
!
!
   module subroutine get_tamm_dancoff_preconditioner_hf(wf, preconditioner)
!!
!!    Get Tamm-Dancoff preconditioner
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Preconditioner is orbital energy differences
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_v*wf%n_o) :: preconditioner
!
      call wf%get_orbital_differences(preconditioner, wf%n_v*wf%n_o)
!
   end subroutine get_tamm_dancoff_preconditioner_hf
!
!
   module subroutine get_dipole_gradient(wf, g, k)
!!
!!    Get dipole gradient
!!    Written by Sarai D. Folkestadm Mar 2022
!!
!!    NOTE: In order to keep the definitions of A and B as in the
!!    eigenvalue equation, we define our orbital rotation operator
!!    as
!!
!!       exp(iK), where
!!       K = 1/sqrt(2) sum_n K_n E_n, for positive and negative n
!!
!!    This results in the prefactor sqrt(2).
!!    The index n is a super-index denoting an excitation
!!
!!       E_n = E_ai
!!       E_-n = E_ia = E_ai^dagger
!!
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(out) :: g
!
      integer, intent(in) :: k
!
      real(dp), dimension(:,:,:), allocatable :: mu_wxk
      real(dp), dimension(:,:),   allocatable :: mu_pqk
!
      integer :: a, i, ai
!
      real(dp) :: prefactor
!
      call mem%alloc(mu_wxk, wf%ao%n, wf%ao%n, 3)
      call wf%ao%get_oei('dipole', mu_wxk)
!
      call mem%alloc(mu_pqk, wf%n_mo, wf%n_mo)
      call wf%mo_transform(mu_wxk(:,:,k), mu_pqk)
!
      prefactor = sqrt(two)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai =  wf%n_v*(i - 1) + a
!
            g(ai) = prefactor*mu_pqk(wf%n_o + a, i)
            g(ai + wf%n_v*wf%n_o) = -prefactor*mu_pqk(i, wf%n_o + a)
!
         enddo
      enddo
!
      call mem%dealloc(mu_wxk, wf%ao%n, wf%ao%n, 3)
      call mem%dealloc(mu_pqk, wf%n_mo, wf%n_mo)
!
   end subroutine get_dipole_gradient
!
!
   module subroutine get_S2_transformation(wf, trial, transform)
!!
!!    Get S2 transformation
!!    Written by Sarai D. Folkestad, Mar 2022
!!
!!    Computes the metric of the response equation
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(in) :: trial
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(out) :: transform
!
      integer :: a, i, ai
!
      transform = zero
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai =  wf%n_v*(i - 1) + a
!
            transform(ai) = trial(ai)
            transform(ai + wf%n_v*wf%n_o) = - trial(ai + wf%n_v*wf%n_o)
!
         enddo
      enddo
!
   end subroutine get_S2_transformation
!
!
   module subroutine get_E2_transformation(wf, trial, transform)
!!
!!    Get E2 transformation
!!    Written by Sarai D. Folkestad, Mar 2022
!!
!!    Computes the Hessian
!!
!!       E^[2] = (  A -B )
!!               ( -B  A )
!!
!!    NOTE: the sign of B is related to
!!    its definition and the antisymmetry of the
!!    property (e.g. dipole) gradient.
!!
!!    Optionally the property gradient can be symmetric and
!!    -B -> B, see e.q. McWheeny, Methods of molecular quantum mechanics
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(in) :: trial
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(out) :: transform
!
      transform = zero
!
      call wf%A_transformation(trial(1:wf%n_v*wf%n_o), transform(1:wf%n_v*wf%n_o))
      call wf%B_transformation(-trial(wf%n_v*wf%n_o+1:), transform(1:wf%n_v*wf%n_o))
!
      call wf%B_transformation(-trial(1:wf%n_v*wf%n_o), transform(wf%n_v*wf%n_o+1:))
      call wf%A_transformation(trial(wf%n_v*wf%n_o+1:), transform(wf%n_v*wf%n_o+1:))
!
   end subroutine get_E2_transformation
!
!
   module subroutine print_polarizability_hf(wf, polarizability, frequency)
!!
!!    Print polarizability
!!    Written by Regina Matveeva, Feb 2022
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(3,3), intent(in) :: polarizability
      real(dp), intent(in) :: frequency
!
      character(len=4), dimension(3) :: operator_ = ['mu_x', 'mu_y', 'mu_z']
!
      integer :: k, l
!
      call output%printf('m', '(a0) polarizability for frequency (e11.4):', &
         chars=[convert_to_uppercase(trim(wf%name_))], reals=[frequency], fs='(/t3,a)')
!
      call output%print_separator('m', 44, '-', fs='(t3,a)')
!
      do k = 1, 3
         do l = 1,3
!
            call output%printf('m', '<< ' // operator_(k) // ', ' //  &
                              operator_(l) // ' >>' // ': (f19.12)', &
                              reals=[polarizability(k,l)], fs='(t6,a)')
!
         end do
      end do
!
      call output%print_separator('m', 44, '-', fs='(t3,a)')

!
   end subroutine print_polarizability_hf
!
!
   module subroutine get_response_transformation(wf, trial, transform, frequency)
!!
!!    Get response transformation
!!    Written by Sarai D. Folkestad, Mar 2022
!!
!!    Computes
!!
!!       Y = (E^[2] - omega S^[2]) X
!!
!!    where X is the trial and Y is the transform
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp), intent(in) :: frequency
!
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(in) :: trial
      real(dp), dimension(2*wf%n_o*wf%n_v), intent(out) :: transform
      real(dp), dimension(:), allocatable :: S2_transform
!
      call mem%alloc(S2_transform, 2*wf%n_v*wf%n_o)
!
      call wf%get_E2_transformation(trial, transform)
      call wf%get_S2_transformation(trial, S2_transform)
!
      call daxpy(2*wf%n_o*wf%n_v, -frequency, S2_transform, 1, transform, 1)
!
      call mem%dealloc(S2_transform, 2*wf%n_v*wf%n_o)
!
   end subroutine get_response_transformation
!
!
end submodule tdhf_transformation
