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
submodule (hf_class) ao_fock 
!
!!
!!    AO Fock submodule 
!!
!!    Collects the routines used in the construction of the AO Fock matrix.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine update_fock_and_energy_hf(wf, prev_ao_density)
!!
!!    Wrapper for Update Fock and energy
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Call either cumulative or non-cumulative updating depending on options
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
      if (.not. present(prev_ao_density)) then
!
          call wf%update_fock_and_energy_non_cumulative()
!
      else
!
         if(.not.wf%system%mm_calculation.and..not.wf%system%pcm_calculation) then
!
            call wf%update_fock_and_energy_cumulative(prev_ao_density)
!
         else
!
            if(wf%system%mm%forcefield.eq.'non-polarizable') then
!
               call wf%update_fock_and_energy_cumulative(prev_ao_density)
!
            else
!
               call wf%update_fock_and_energy_non_cumulative()
!
            endif
!
         endif
!
      endif
!
!
   end subroutine update_fock_and_energy_hf
!
!
   module subroutine update_fock_and_energy_cumulative_hf(wf, prev_ao_density)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed. By cumulatively
!!    we mean using the density change to build the Fock matrix
!!    in the iterative loop.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in) :: prev_ao_density
!
      real(dp), dimension(:, :), allocatable :: G
!
      type(timings) :: timer
!
      timer = timings('AO Fock construction', pl='normal')
      call timer%turn_on()
!
!     Fock
!
!     Construct the two electron part of the Fock matrix (G),
!     and add the contribution to the Fock matrix
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
!
      call mem%alloc(G, wf%n_ao, wf%n_ao)
      call wf%construct_ao_G(wf%ao_density, G)
!
      call daxpy(wf%n_ao**2, one, G, 1, wf%ao_fock, 1)
!
      call mem%dealloc(G, wf%n_ao, wf%n_ao)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)
!
      call timer%turn_off()
!
!     Energy
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, wf%ao_h)
!
!     Note! Must be called after wf%calculate_hf_energy_from_fock !!
!
      if(wf%system%mm_calculation) call wf%calculate_mm_energy_terms()
      if(wf%system%pcm_calculation) call wf%calculate_pcm_energy_terms()
!
   end subroutine update_fock_and_energy_cumulative_hf
!
!
   module subroutine update_fock_and_energy_non_cumulative_hf(wf)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified by Tommaso Giovannini, May 2019 for QM/MM
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed.
!!
      implicit none
!
      class(hf) :: wf
!
      type(timings) :: timer
!
      timer = timings('AO Fock construction', pl='normal')
      call timer%turn_on()
!
!     Construct the two electron part of the Fock matrix (G),
!     and add the contribution to the Fock matrix
!
      call wf%construct_ao_G(wf%ao_density, wf%ao_fock)
!
!     Add the one-electron part
!
      call daxpy(wf%n_ao**2, one, wf%ao_h, 1, wf%ao_fock, 1)
!
!     QM/MM and PCM specific work
!
      if (wf%system%mm_calculation) call wf%add_mm_fock_terms()
      if (wf%system%pcm_calculation) call wf%add_pcm_fock_term()
!
      call timer%turn_off()
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, wf%ao_h)
!
!     QM/MM and PCM energy terms (NOTE! these must be called after calculate_hf_energy_from_fock !!)
!
      if(wf%system%mm_calculation) call wf%calculate_mm_energy_terms()
      if(wf%system%pcm_calculation) call wf%calculate_pcm_energy_terms()
!
   end subroutine update_fock_and_energy_non_cumulative_hf
!
!
   module subroutine construct_ao_G_hf(wf, D, G, C_screening)
!!
!!    Construct AO G(D)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Modified by Sarai D. Folkestad, Oct 2019 for G construction
!!    only
!!
!!    Calculates
!!
!!       G(D)_αβ = sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where D is the inactive AO Density.
!!
!!    Modified by Eirik F. Kjønstad, Jan 2020. Added C_screening optional.
!!
!!    C_screening: If true, G(D) will be constructed in the AO basis, as usual, but the 
!!                 Coulomb and exchange screening will target the MO basis G(D). Used in 
!!                 MLHF when constructing G(Da). Default: false.
!!
      class(hf)   :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: D
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: G
!
      integer :: thread = 0, n_threads = 1
!
      real(dp), dimension(:,:), allocatable :: sp_Density_schwarz
      real(dp), dimension(:,:,:), allocatable :: G_thread
!
      integer :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      logical, optional, intent(in) :: C_screening
      logical :: C_screening_local 
!
      C_screening_local = .false.
      if (present(C_screening)) C_screening_local = C_screening
!
      if (C_screening_local) &
         call output%printf('v', 'Will perform MO screening in G(D) construction')
!
!     Construct the density screening vector and the maximum element in the density
!
      call mem%alloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
!
      call wf%construct_sp_Density_schwarz(sp_Density_schwarz, D)
      max_D_schwarz = get_abs_max(sp_Density_schwarz, wf%system%n_s**2)
!
!     Compute number of significant ERI shell pairs (the G construction
!     only loops over these shell pairs) and the maximum element
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
      max_eri_schwarz = get_abs_max(wf%sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
!     Construct the two electron part of the Fock matrix, using the screening vectors
!     and parallelizing over available threads (each gets its own copy of the Fock matrix)
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(G_thread, wf%n_ao, wf%n_ao, n_threads) ! [G(thread 1) G(thread 2) ...]
      call zero_array(G_thread, (wf%n_ao**2)*n_threads)
!
      if (.not. C_screening_local) then 
!
         call wf%construct_ao_G_thread_terms(G_thread, D, n_threads, max_D_schwarz, max_eri_schwarz, &
                                         sp_density_schwarz, n_sig_sp,                 &
                                         wf%coulomb_threshold, wf%exchange_threshold,  &
                                         wf%libint_epsilon, wf%system%shell_limits)
!
      else 
!
         call wf%construct_ao_G_thread_terms_mo_screened(G_thread, D, n_threads, max_D_schwarz, &
                                         max_eri_schwarz, &
                                         sp_density_schwarz, n_sig_sp,                 &
                                         wf%coulomb_threshold, wf%exchange_threshold,  &
                                         wf%libint_epsilon, wf%system%shell_limits)         
!
      endif
!
      call mem%dealloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix,
!     and symmetrize the result
!
      call zero_array(G, (wf%n_ao**2))
!
      do thread = 1, n_threads
!
         call daxpy(wf%n_ao**2, one, G_thread(1, 1, thread), 1, G, 1)
!
      enddo
!
      call mem%dealloc(G_thread, wf%n_ao, wf%n_ao, n_threads)
!
      call symmetric_sum(G, wf%n_ao)
      call dscal(wf%n_ao**2, half, G, 1)
!
   end subroutine construct_ao_G_hf
!
!
   module subroutine construct_ao_G_thread_terms_hf(wf, F, D, n_threads, max_D_schwarz,   &
                                          max_eri_schwarz, sp_density_schwarz, n_sig_sp,  &
                                          coulomb_thr, exchange_thr, precision_thr, shells)
!!
!!    Construct AO G
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ =+ sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ (= G(D)_αβ),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get G(D)_αβ.
!!
!!    Routine partly based on Hartree-Fock implementation shipped with 
!!    the Libint 2 integral package by E. Valeev. 
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads, n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz 
      real(dp), intent(in) :: coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
      real(dp) :: deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                             &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34,     &
!$omp w, x, y, z, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset,         &
!$omp temp4, temp5, temp6, temp7, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,        &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                                &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4), &
                           sp_density_schwarz_s1s2)
!
               temp8 = max(sp_density_schwarz_s3s2, &
                           sp_density_schwarz_s3s1, &
                           sp_density_schwarz(s4,s2), &
                           sp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr .and. temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4,         &
                  precision_thr/max(temp7,temp8), thread, skip, shells(s1)%length, shells(s2)%length, &
                  shells(s3)%length, shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*(shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*(shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           temp3 = eighth*temp*d3
                           temp4 = eighth*temp*d4
                           temp5 = eighth*temp*d5
                           temp6 = eighth*temp*d6
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
!
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
!
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call set_coulomb_precision_c(wf%libint_epsilon)
!
   end subroutine construct_ao_G_thread_terms_hf
!
!
   module subroutine construct_ao_G_thread_terms_mo_screened_hf(wf, F, D, n_threads, max_D_schwarz,   &
                                          max_eri_schwarz, sp_density_schwarz, n_sig_sp,  &
                                          coulomb_thr, exchange_thr, precision_thr, shells)
!!
!!    Construct AO G MO screened
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ =+ sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ (= G(D)_αβ),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get G(D)_αβ.
!!
!!    Routine partly based on Hartree-Fock implementation shipped with 
!!    the Libint 2 integral package by E. Valeev. 
!!
!!    This routine modifies the existing construct AO G routine by:
!!
!!       - Constructing the maximum MO coefficients list:
!!
!!          C_max(s) = max_p | C_wp | for AOs w in the shell s
!!
!!       - Use this MO coefficients list to screen (Coulomb, exchange) for the precision of 
!!
!!          G_pq = G_alpha,beta C_alpha,p C_beta,q        (not G_alpha,beta)
!!
!!    Used to construct G(Da), where Da is the active density, in MLHF. In that case 
!!    the MOs are local, which means that screening for G(Da) in the MO basis is more efficient  
!!    than screening for G(Da) in the AO basis.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads, n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz 
      real(dp), intent(in) :: coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
      real(dp) :: deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip, q
!
      type(timings), allocatable :: C_max_timer
!
      real(dp), dimension(:), allocatable :: C_max
!
      call mem%alloc(C_max, wf%system%n_s) 
!
      C_max_timer = timings('MO coefficients maximums screening vector', 'n')
      call C_max_timer%turn_on()
!
      do s1 = 1, wf%system%n_s 
!
         C_max(s1) = zero 
!
         do q = 1, wf%n_mo 
            do w = shells(s1)%first, shells(s1)%last
!
               if (abs(wf%orbital_coefficients(w, q)) .gt. C_max(s1)) then 
!
                  C_max(s1) = abs(wf%orbital_coefficients(w, q))
!
               endif
!
            enddo
         enddo
      enddo
!
      call C_max_timer%turn_off()
!
!$omp parallel do                                                                             &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34,     &
!$omp w, x, y, z, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset,         &
!$omp temp4, temp5, temp6, temp7, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,        &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                                &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)
!
               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4)*C_max(s1)*C_max(s2),  &
                           sp_density_schwarz_s1s2*C_max(s3)*C_max(s4))
!
               temp8 = max(sp_density_schwarz_s3s2*C_max(s1)*C_max(s4),    &
                           sp_density_schwarz_s3s1*C_max(s2)*C_max(s4),    &
                           sp_density_schwarz(s4,s2)*C_max(s1)*C_max(s3),  &
                           sp_density_schwarz(s1,s4)*C_max(s2)*C_max(s3))
!
               if (temp8*temp .lt. exchange_thr .and. temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4,              &
                                                          precision_thr/max(temp7,temp8), &
                                                          thread,                         &
                                                          skip,                           &
                                                          shells(s1)%length,              &
                                                          shells(s2)%length,              &
                                                          shells(s3)%length,              &
                                                          shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*&
                         (shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*&
                                    (shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           temp3 = eighth*temp*d3
                           temp4 = eighth*temp*d4
                           temp5 = eighth*temp*d5
                           temp6 = eighth*temp*d6
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
!
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
!
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call set_coulomb_precision_c(wf%libint_epsilon)
      call mem%dealloc(C_max, wf%system%n_s) 
!
   end subroutine construct_ao_G_thread_terms_mo_screened_hf
!
!
   module subroutine construct_coulomb_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,     &
                                                   sp_density_schwarz, &
                                                   n_sig_sp, coulomb_thr, precision_thr, shells)
!!
!!    AO Fock Coulomb construction loop
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the Coulomb two-electron part of the Fock matrix,
!!
!!       F_αβ = F_αβ + sum_γδ g_αβγδ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get the Coulomb part of G(D)_αβ.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads,  n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: d1, d2, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp7, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp1, temp2, d1, d2, thread, thread_offset,                            &
!$omp temp7, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,                                &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                            &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4), sp_density_schwarz_s1s2)
!
               if (temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp7, thread, skip, shells(s1)%length, shells(s2)%length,    &
                  shells(s3)%length, shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*(shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*(shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_coulomb_ao_G_hf
!
!
   module subroutine construct_exchange_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz, &
                                          sp_density_schwarz, &
                                           n_sig_sp, exchange_thr, precision_thr, shells)
!!
!!    AO Fock exchange construction loop
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ = F_αβ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!!    Note: the contributions from each thread need to be added to a single
!!    n_ao x n_ao matrix & symmetrized to get the exchange part of G(D)_αβ.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads,  n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp3, temp4, temp5, temp6, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp3, d3, d4, d5, d6, thread, thread_offset,                           &
!$omp temp4, temp5, temp6, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,           &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                            &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. exchange_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. exchange_thr) cycle ! Screened out shell pair
!
               temp8 = max(sp_density_schwarz_s3s2,   &
                           sp_density_schwarz_s3s1,   &
                           sp_density_schwarz(s4,s2), &
                           sp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp8, thread, skip, shells(s1)%length, shells(s2)%length,    &
                  shells(s3)%length, shells(s4)%length)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%length)*(shells(s2)%length)*(shells(s3)%length)*(shells(s4)%length)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%length*(shells(s2)%length*(shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp3 = eighth*temp*d3
                           temp4 = eighth*temp*d4
                           temp5 = eighth*temp*d5
                           temp6 = eighth*temp*d6
!
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_exchange_ao_G_hf
!
!
   module subroutine construct_sp_eri_schwarz_hf(wf)
!!
!!    Construct shell-pair electronic-repulsion-integral Schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of g_wxwx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively.
!!
!
      use array_utilities, only: quicksort_with_index_descending
!
      implicit none
!
      class(hf) :: wf
!
      integer, dimension(:),  allocatable :: sp_eri_schwarz_index_list
      real(dp), dimension(:), allocatable :: sorted_sp_eri_schwarz
!
!     Local variables
!
      integer :: s1, s2, s1s2, sp
!
      real(dp) :: maximum
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      type(interval) :: A_interval, B_interval
!
      type(timings) :: timer, timer_sort
!
      timer = timings('Construct ERI screening list', 'v')
      call timer%turn_on()
!
!     Set the maximum element in each shell pair
!
      call set_coulomb_precision_c(1.0d-50)
!
!$omp parallel do private(s1, s2, s1s2, A_interval, B_interval, g, maximum) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
            s1s2 = (max(s1,s2)*(max(s1,s2)-3)/2) + s1 + s2
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call wf%system%construct_ao_g_wxyz(g, s1, s2, s1, s2)
!
            maximum = get_abs_max(g, ((A_interval%length)*(B_interval%length))**2)
!
            wf%sp_eri_schwarz(s1s2, 1) = sqrt(maximum)
!
            wf%sp_eri_schwarz_list(s1s2, 1) = s1
            wf%sp_eri_schwarz_list(s1s2, 2) = s2
!
         enddo
      enddo
!$omp end parallel do
!
!     Sort the sp_eri_schwarz vector and use the resulting index list
!     to resort the sp_eri_schwarz_list matrix
!
      call mem%alloc(sp_eri_schwarz_index_list, wf%system%n_s*(wf%system%n_s + 1)/2)
      call mem%alloc(sorted_sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
      call dcopy(wf%system%n_s*(wf%system%n_s + 1)/2, &
                  wf%sp_eri_schwarz(:, 1), 1, sorted_sp_eri_schwarz, 1)
!
      timer_sort = timings('Construct ERI screening list (time to sort)','v')
      call timer_sort%turn_on()
!
      call quicksort_with_index_descending(sorted_sp_eri_schwarz, sp_eri_schwarz_index_list, &
                                                wf%system%n_s*(wf%system%n_s + 1)/2)
!
      call timer_sort%turn_off()
!
      call dcopy(wf%system%n_s*(wf%system%n_s + 1)/2, &
                  wf%sp_eri_schwarz(:, 1), 1, wf%sp_eri_schwarz(:, 2), 1)
!
      call dcopy(wf%system%n_s*(wf%system%n_s + 1)/2, sorted_sp_eri_schwarz, 1, &
                  wf%sp_eri_schwarz(:, 1), 1)

      call mem%dealloc(sorted_sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
!$omp parallel do private(sp)
      do sp = 1, wf%system%n_s*(wf%system%n_s + 1)/2
!
         wf%sp_eri_schwarz_list(sp, 3) = sp_eri_schwarz_index_list(sp)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(sp_eri_schwarz_index_list, wf%system%n_s*(wf%system%n_s + 1)/2)
!
      call timer%turn_off()
!
   end subroutine construct_sp_eri_schwarz_hf
!
!
   module subroutine construct_sp_density_schwarz_hf(wf, sp_density_schwarz, D)
!!
!!    Construct shell-pair density schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of D_wx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s) :: sp_density_schwarz
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), dimension(:,:), allocatable, target :: D_red
      real(dp), dimension(:,:), contiguous, pointer :: D_red_p => null()
!
      type(interval) :: A_interval, B_interval
!
      integer :: s1, s2
      integer :: n_threads = 1, thread = 0
!
      real(dp) :: maximum
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(D_red,wf%system%max_shell_size**2,n_threads)
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D_red_p, maximum, thread) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
!$          thread = omp_get_thread_num()
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            D_red_p(1:A_interval%length,1:B_interval%length) => D_red(1:A_interval%length*B_interval%length,thread+1)
!
            D_red_p = D(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            maximum = get_abs_max(D_red_p, (A_interval%length)*(B_interval%length))
!
            nullify(D_red_p)
!
            sp_density_schwarz(s1, s2) = maximum
            sp_density_schwarz(s2, s1) = maximum
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_red,wf%system%max_shell_size**2,n_threads)
!
   end subroutine construct_sp_density_schwarz_hf
!
!
   module subroutine get_n_sig_eri_sp_hf(wf, n_sig_sp)
!!
!!    Get number of significant ERI shell-pairs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the number of significant shell pairs. The threshold
!!    determines how small the largest element of g_wxwx in a shell
!!    pair AB (w in A, x in B) to be ignored completely in the Fock
!!    construction loop.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(inout) :: n_sig_sp
!
      integer :: s1s2
!
      n_sig_sp = 0
!
      do s1s2 = 1, wf%system%n_s*(wf%system%n_s + 1)/2
!
         if (wf%sp_eri_schwarz(s1s2, 1)*wf%sp_eri_schwarz(1, 1) .lt. sqrt(wf%libint_epsilon)) then
!
            exit
!
         else
!
            n_sig_sp = n_sig_sp + 1
!
         endif
!
      enddo
!
   end subroutine get_n_sig_eri_sp_hf
!
!
   module subroutine construct_ao_G_1der_hf(wf, G_ao, D_ao)
!!
!!    Construct AO G 1der
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβqk = sum_γδ g_αβγδqk D_γδ - 1/2 * sum_γδ g_αδγβqk D_γδ (= G(D)_αβqk),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: G_ao
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_ao
!
      real(dp), dimension((wf%system%max_shell_size**4)*3*4), target :: g_ABCDqk
      real(dp), dimension(:,:,:,:,:,:), pointer, contiguous :: g_ABCDqk_p
!
      integer :: A, B, C, D, D_max, w, x, y, z, n_sig_sp, AB, AB_packed
      integer :: w_red, x_red, y_red, z_red, tot_dim, k, q, n_threads, thread
!
      real(dp) :: d1, d2, d3, d4, d5, d6
!
      integer, dimension(4) :: atoms
!
      real(dp) :: deg, deg_AB, deg_CD, deg_AB_CD
!
      real(dp), dimension(3,4) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
!
      real(dp), dimension(:,:,:,:,:), allocatable :: G_ao_t
!
      thread = 0
      n_threads = 1
!$    n_threads = omp_get_max_threads()
      call mem%alloc(G_ao_t, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms, n_threads)
      G_ao_t = zero
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
!
!$omp parallel do private(A, B, C, D, D_max, atoms, deg, deg_CD, deg_AB, deg_AB_CD, g_ABCDqk, g_ABCDqk_p, &
!$omp w, x, y, z, w_red, x_red, y_red, z_red, temp, temp1, temp2, temp3, temp4, temp5, temp6, &
!$omp d1, d2, d3, d4, d5, d6, thread, q, k, tot_dim, AB, AB_packed) schedule(dynamic)
      do AB = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
!
         if (wf%sp_eri_schwarz(AB, 1)*wf%sp_eri_schwarz(1, 1) < wf%coulomb_threshold) cycle
         AB_packed = wf%sp_eri_schwarz_list(AB, 3)
!
         A = wf%sp_eri_schwarz_list(AB_packed, 1)
         B = wf%sp_eri_schwarz_list(AB_packed, 2)
!
         atoms(1) = wf%system%shell2atom(A)
         atoms(2) = wf%system%shell2atom(B)
!
         deg_AB = real(2-B/A, kind=dp)
!
            do C = 1, A
!
               D_max = (C/A)*B + (1-C/A)*C
               atoms(3) = wf%system%shell2atom(C)
!
               do D = 1, D_max
!
                  deg_CD    = real(2-D/C, kind=dp)
                  deg_AB_CD = min(1-C/A+2-min(D/B,B/D), 2)

                  deg = deg_AB*deg_CD*deg_AB_CD ! Shell degeneracy
!
                  atoms(4) = wf%system%shell2atom(D)
!
                  call wf%system%construct_ao_g_wxyz_1der(g_ABCDqk, A, B, C, D)
!
                  tot_dim = (wf%system%shell_limits(A)%length)*(wf%system%shell_limits(B)%length)&
                              *(wf%system%shell_limits(C)%length)*(wf%system%shell_limits(D)%length)&
                              *3*4
!
                  g_ABCDqk(1:tot_dim) = deg*g_ABCDqk(1:tot_dim)
!
                  g_ABCDqk_p(1 : wf%system%shell_limits(A)%length, 1 : wf%system%shell_limits(B)%length, &
                             1 : wf%system%shell_limits(C)%length, 1 : wf%system%shell_limits(D)%length, &
                             1 : 3, 1 : 4) => g_ABCDqk(1 : tot_dim)
!
                  do z = wf%system%shell_limits(D)%first, wf%system%shell_limits(D)%last
!
                     z_red = z - wf%system%shell_limits(D)%first + 1
!
                     do y = wf%system%shell_limits(C)%first, wf%system%shell_limits(C)%last
!
                        y_red = y - wf%system%shell_limits(C)%first + 1
!
                        d1 = D_ao(y, z)
!
                        do x = wf%system%shell_limits(B)%first, wf%system%shell_limits(B)%last
!
                           x_red = x - wf%system%shell_limits(B)%first + 1
!
                           d3 = D_ao(x, y)
                           d5 = D_ao(x, z)
!
                           do w = wf%system%shell_limits(A)%first, wf%system%shell_limits(A)%last
!
                              d2 = D_ao(w, x)
                              d4 = D_ao(w, y)
                              d6 = D_ao(w, z)
!
                              w_red = w - wf%system%shell_limits(A)%first + 1
!
                              do k = 1, 4
                                 do q = 1, 3
                                    temp(q, k) = g_ABCDqk_p(w_red, x_red, y_red, z_red, q, k)
                                 enddo
                              enddo
!
                              temp1 = half*temp*d1
                              temp2 = half*temp*d2
!
                              temp3 = eighth*temp*d3
                              temp4 = eighth*temp*d4
                              temp5 = eighth*temp*d5
                              temp6 = eighth*temp*d6
!
                              do k = 1, 4
                                 do q = 1, 3
!
                                    G_ao_t(w, x, q, atoms(k), thread+1) = G_ao_t(w, x, q, atoms(k), thread+1) + temp1(q, k)
                                    G_ao_t(y, x, q, atoms(k), thread+1) = G_ao_t(y, x, q, atoms(k), thread+1) - temp6(q, k)
                                    G_ao_t(y, z, q, atoms(k), thread+1) = G_ao_t(y, z, q, atoms(k), thread+1) + temp2(q, k)
                                    G_ao_t(w, z, q, atoms(k), thread+1) = G_ao_t(w, z, q, atoms(k), thread+1) - temp3(q, k)
                                    G_ao_t(x, z, q, atoms(k), thread+1) = G_ao_t(x, z, q, atoms(k), thread+1) - temp4(q, k)
                                    G_ao_t(w, y, q, atoms(k), thread+1) = G_ao_t(w, y, q, atoms(k), thread+1) - temp5(q, k)
!
                                 enddo
                              enddo
!
                           enddo
                        enddo
                     enddo
                  enddo
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
      G_ao = zero
!
      do thread = 1, n_threads
!
         call daxpy(3*wf%system%n_atoms*wf%n_ao**2, one, G_ao_t(1,1,1,1,thread), 1, G_ao, 1)
!
      enddo
!
      call mem%dealloc(G_ao_t, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms, n_threads)
!
   end subroutine construct_ao_G_1der_hf
!
!
   module subroutine add_pcm_fock_term_hf(wf)
!!
!!    Update Fock PCM 
!!    Written by Tommaso Giovannini, Oct 2019 
!!
!!    The QM Fock is updated with the contributions coming 
!!    from the PCM:
!!       q*V_wx
!!
!!    Done by interfacing to PCMSolver
!!
      implicit none
!
      class(hf) :: wf
!
      integer :: i
!     
! 
      call zero_array(wf%pcm_fock, wf%n_ao*wf%n_ao)
      call zero_array(wf%system%pcm%pcm_rhs, wf%system%pcm%n_tesserae)
!      
!     electrostatic potential contracted with density : \sum_i V_mu(D_mu)(r_i)
!
      call wf%construct_ao_electrostatics(0, 1, 'prop', &
                                          wf%system%pcm%n_tesserae, &
                                          wf%system%pcm%grid_coord*bohr_to_angstrom, &
                                          property_points=wf%system%pcm%pcm_rhs, &
                                          ao_density=wf%ao_density) 
!      
!     solve q=D^-1 (V(D)) 
! 
      call wf%system%pcm%set_surface_function('NucMEP')
!                                          
      call wf%system%pcm%compute_asc('NucMEP', 'NucASC')
!                                          
      call wf%system%pcm%get_surface_function('NucASC')
!                                          
      call output%print_separator('verbose', 67, fs='(/t3,a)')
!
      call output%printf('v', 'Atom         PCM ASC            PCM RHS', fs='(t6,a)')
!
      do i = 1, wf%system%pcm%n_tesserae
!
         call output%printf('v', '(i4)      (e13.6)      (e13.6)', ints=[i], &
                            reals=[wf%system%pcm%charges(i), &
                            wf%system%pcm%pcm_rhs(i)], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('verbose', 67)
!
!     Fock creation: F_munu = \sum_i q_i V_munu(r_i)
!
      call wf%construct_ao_electrostatics(0, 0, 'fock', &
                                          wf%system%pcm%n_tesserae, &
                                          wf%system%pcm%grid_coord*bohr_to_angstrom, & 
                                          elec_fock=wf%pcm_fock, &
                                          charges=wf%system%pcm%charges) 
!
      wf%ao_fock = wf%ao_fock + half*wf%pcm_fock
!
      call output%print_matrix('debug', 'QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'PCM Fock', wf%pcm_fock, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM/PCM Fock', wf%ao_fock, wf%n_ao, wf%n_ao)
!
!
   end subroutine add_pcm_fock_term_hf
!
!
   module subroutine add_mm_fock_terms_hf(wf)
!!
!!    Update Fock with polarizable QM/MM terms
!!    For now: QM/FQ model (see mm_class and output file)
!!    Written by Tommaso Giovannini, July 2019 for QM/MM
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:), allocatable                   :: potential_points
      integer :: i
!
      if (wf%system%mm%forcefield .eq. 'non-polarizable') then
!
         call daxpy(wf%n_ao**2, half, wf%nopol_h_wx, 1, wf%ao_fock, 1)
!
      elseif (wf%system%mm%forcefield .eq. 'fq') then
!
         if(.not.allocated(potential_points)) call mem%alloc(potential_points, wf%system%mm%n_atoms)
!
         call zero_array(wf%pol_emb_fock,wf%n_ao*wf%n_ao)
         call zero_array(wf%system%mm%pol_emb_rhs,wf%system%mm%n_variables)
!
!        electrostatic potential contracted with density : \sum_i V_mu(D_mu)(r_i)
!
         call wf%construct_ao_electrostatics(0,1,'prop',wf%system%mm%n_atoms,wf%system%mm%coordinates, &
                                             property_points=potential_points,ao_density=wf%ao_density) 
!
!        rhs for fq: -chi - V(D)
!
         wf%system%mm%pol_emb_rhs(1:wf%system%mm%n_atoms) = -wf%system%mm%chi + potential_points
!
!        solve q=D^-1 (-chi - V(D))
!
         call dgemm('N', 'N',                   &
                     wf%system%mm%n_variables,  &
                     1,                         &
                     wf%system%mm%n_variables,  &
                     one,                       &
                     wf%system%mm%fq_matrix,    &
                     wf%system%mm%n_variables,  &
                     wf%system%mm%pol_emb_rhs,  &
                     wf%system%mm%n_variables,  &
                     zero,                      &
                     wf%system%mm%pol_emb_lhs,  &
                     wf%system%mm%n_variables)
!
!
         call output%print_separator('verbose', 67, fs='(/t3,a)')
!
         call output%printf('v', 'Atom          FQ LHS             FQ RHS       &
                            & QM Potential@FQs', fs='(t6,a)')
!
         do i = 1, wf%system%mm%n_atoms
!
            call output%printf('v', '(i4)      (e13.6)      (e13.6)      (e13.6)', &
                               ints=[i], reals=[wf%system%mm%pol_emb_lhs(i), &
                               wf%system%mm%pol_emb_rhs(i), &
                               potential_points(i)], fs='(t6,a)')
!
         enddo
!
         call output%print_separator('verbose', 67)
!
!        Put FQ charges into charge (I am discrading langrangian multipliers)
!
         wf%system%mm%charge = wf%system%mm%pol_emb_lhs(1:wf%system%mm%n_atoms)
!
!        Fock creation: F_munu = \sum_i q_i V_munu(r_i)
!
         call wf%construct_ao_electrostatics(0,0,'fock',wf%system%mm%n_atoms,wf%system%mm%coordinates, & 
                                             elec_fock=wf%pol_emb_fock,charges=wf%system%mm%charge) 
!
         call daxpy(wf%n_ao**2, half, wf%pol_emb_fock, 1, wf%ao_fock, 1)
!
         call output%print_matrix('debug', 'QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
         call output%print_matrix('debug', 'FQ Fock', wf%pol_emb_fock, wf%n_ao, wf%n_ao)
         call output%print_matrix('debug', 'QM/FQ Fock', wf%ao_fock, wf%n_ao, wf%n_ao)
!
         call mem%dealloc(potential_points, wf%system%mm%n_atoms)
!
      else
!
         call output%error_msg('did not recognize the force field')
!
      endif
!
   end subroutine add_mm_fock_terms_hf
!
!
   module subroutine set_screening_and_precision_thresholds_hf(wf, gradient_threshold)
!!
!!    Set screening and precision thresholds
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the screening thresholds for Coulomb and exchange
!!    integrals given the convergence threshold for the gradient
!!
!!       coulomb_threshold  = gradient_threshold * 1.0d-3
!!       exchange_threshold = gradient_threshold * 1.0d-3
!!
!!       libint_epsilon = (gradient_threshold * 1.0d-3)**2
!!
!!    unless stricter thresholds are already set on input or by default.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), intent(in) :: gradient_threshold
!
      if (wf%coulomb_threshold .gt. gradient_threshold*1.0d-3)  &
                  wf%coulomb_threshold  = gradient_threshold*1.0d-3
!
      if (wf%exchange_threshold .gt. gradient_threshold*1.0d-3) &
                  wf%exchange_threshold = gradient_threshold*1.0d-3
!
      wf%coulomb_threshold  = min(wf%coulomb_threshold,  1.0d-12)
      wf%exchange_threshold = min(wf%exchange_threshold, 1.0d-10)
!
      if (wf%libint_epsilon .gt. (wf%coulomb_threshold)**2) &
                  wf%libint_epsilon = (wf%coulomb_threshold)**2
!
   end subroutine set_screening_and_precision_thresholds_hf
!
!
end submodule ao_fock
