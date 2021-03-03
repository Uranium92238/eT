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
   module subroutine update_fock_and_energy_hf(wf, cumulative)
!!
!!    Wrapper for Update Fock and energy
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Call either cumulative or non-cumulative updating depending on options
!!
      implicit none
!
      class(hf), intent(inout) :: wf
      logical, intent(in) :: cumulative
!
      if (.not. cumulative) then
!
         call wf%update_G_non_cumulative()
!
     else
!
         call wf%update_G_cumulative(wf%previous_ao_density)
!
     endif
!
!     Set the two-electron part
!
      call dcopy(wf%ao%n**2,  wf%ao_G, 1, wf%ao_fock, 1)
!
!     Add the one-electron part
!
      call daxpy(wf%ao%n**2, one, wf%ao%h, 1, wf%ao_fock, 1)
!
!     Embedding contribution
!
      if (wf%embedded) then
!
!        Update embedding if necessary
!
         call wf%embedding%update(wf%ao, wf%ao_density)
!
!        Add embedding term to fock matrix
!
         call daxpy(wf%ao%n**2, one, wf%ao%v, 1, wf%ao_fock, 1)
!
      endif
!
!     Energy
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, wf%ao%h)
!
   end subroutine update_fock_and_energy_hf
!
!
   module subroutine update_G_cumulative_hf(wf, prev_ao_density)
!!
!!    Update G cumulative
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed. By cumulatively
!!    we mean using the density change to build the Fock matrix
!!    in the iterative loop.
!!
!!    Modified by SDF, Sep 2020:
!!    
!!    Cumulative construction of two-electron part and not
!!    full Fock matrix
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in) :: prev_ao_density
!
      real(dp), dimension(:, :), allocatable :: G
!
      type(timings) :: timer
!
      timer = timings('AO G construction', pl='normal')
      call timer%turn_on()
!
      call output%printf('v', 'Fock matrix construction using density differences')
!
!     Determine density differences
!
      call daxpy(wf%ao%n**2, -one, prev_ao_density, 1, wf%ao_density, 1)
!
!     Construct the two electron part of the Fock matrix (G)
!
      call mem%alloc(G, wf%ao%n, wf%ao%n)
      call wf%construct_ao_G(wf%ao_density, G)
      call daxpy(wf%ao%n**2, one, G, 1, wf%ao_G, 1)
      call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
!     Restore the density
!
      call daxpy(wf%ao%n**2, one, prev_ao_density, 1, wf%ao_density, 1)
!
      call timer%turn_off()
!
   end subroutine update_G_cumulative_hf
!
!
   module subroutine update_G_non_cumulative_hf(wf)
!!
!!    Update G
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      type(timings) :: timer
!
      timer = timings('AO G construction', pl='normal')
      call timer%turn_on()
!
!     Construct the two electron part of the Fock matrix (G),
!     and add the contribution to the Fock matrix
!
      call wf%construct_ao_G(wf%ao_density, wf%ao_G)
!
      call timer%turn_off()
!
   end subroutine update_G_non_cumulative_hf
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
      class(hf), intent(inout):: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in)    :: D
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(inout) :: G
!
      integer :: thread = 0, n_threads = 1
!
      real(dp), dimension(:,:), allocatable :: shp_Density_schwarz
      real(dp), dimension(:,:,:), allocatable :: G_thread
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      real(dp) :: eri_precision
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
      call mem%alloc(shp_density_schwarz, wf%ao%n_sh, wf%ao%n_sh)
!
      call wf%construct_shp_Density_schwarz(shp_Density_schwarz, D)
      max_D_schwarz = get_abs_max(shp_Density_schwarz, wf%ao%n_sh**2)
!
!     Compute number of significant ERI shell pairs (the G construction
!     only loops over these shell pairs) and the maximum element
!
      max_eri_schwarz = get_abs_max(wf%ao%cs_eri_max, wf%ao%n_sh*(wf%ao%n_sh + 1)/2)
!
!     Construct the two electron part of the Fock matrix, using the screening vectors
!     and parallelizing over available threads (each gets its own copy of the Fock matrix)
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(G_thread, wf%ao%n, wf%ao%n, n_threads) ! [G(thread 1) G(thread 2) ...]
      call zero_array(G_thread, (wf%ao%n**2)*n_threads)
!
      eri_precision = wf%ao%get_libint_epsilon()
!
      if (.not. C_screening_local) then 
!
         call wf%construct_ao_G_thread_terms(G_thread, D, n_threads, max_D_schwarz,           &
                                             max_eri_schwarz, shp_density_schwarz, wf%ao%n_sig_eri_shp, &
                                             wf%coulomb_threshold, wf%exchange_threshold,     &
                                             eri_precision, wf%ao%shells)
!
      else 
!
         call wf%construct_ao_G_thread_terms_mo_screened(G_thread, D, n_threads, max_D_schwarz,   &
                                                         max_eri_schwarz, shp_density_schwarz,    &
                                                         wf%ao%n_sig_eri_shp, wf%coulomb_threshold,         &
                                                         wf%exchange_threshold, eri_precision,    &
                                                         wf%ao%shells)         
!
      endif
!
      call mem%dealloc(shp_density_schwarz, wf%ao%n_sh, wf%ao%n_sh)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix,
!     and symmetrize the result
!
      call zero_array(G, (wf%ao%n**2))
!
      do thread = 1, n_threads
!
         call daxpy(wf%ao%n**2, one, G_thread(1, 1, thread), 1, G, 1)
!
      enddo
!
      call mem%dealloc(G_thread, wf%ao%n, wf%ao%n, n_threads)
!
      call symmetric_sum(G, wf%ao%n)
      call dscal(wf%ao%n**2, half, G, 1)
!
   end subroutine construct_ao_G_hf
!
!
   module subroutine construct_ao_G_thread_terms_hf(wf, F, D, n_threads, max_D_schwarz,   &
                                          max_eri_schwarz, shp_density_schwarz, n_sig_shp,  &
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
      integer, intent(in) :: n_threads, n_sig_shp
!
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
!
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz 
      real(dp), intent(in) :: coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in) :: shp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, shp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
      real(dp) :: deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: shp_density_schwarz_s1s2, shp_density_schwarz_s3s2, shp_density_schwarz_s3s1
!
      real(dp), dimension(wf%ao%max_sh_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                             &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34,     &
!$omp w, x, y, z, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset,         &
!$omp temp4, temp5, temp6, temp7, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,        &
!$omp shp_eri_schwarz_s1s2, shp_density_schwarz_s1s2, deg_12_34,                                &
!$omp shp_density_schwarz_s3s2, shp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_shp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%ao%n ! Start column of thread's Fock matrix
!
         shp_eri_schwarz_s1s2 = wf%ao%cs_eri_max(s1s2, 1)
!
         s1s2_packed = wf%ao%cs_eri_max_indices(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = wf%ao%cs_eri_max_indices(s1s2_packed, 2)
!
         shp_density_schwarz_s1s2 = shp_density_schwarz(s1, s2)
         if (shp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            shp_density_schwarz_s3s2 = shp_density_schwarz(s3, s2)
            shp_density_schwarz_s3s1 = shp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = shp_eri_schwarz_s1s2*wf%ao%cs_eri_max(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(shp_density_schwarz(s3,s4), &
                           shp_density_schwarz_s1s2)
!
               temp8 = max(shp_density_schwarz_s3s2, &
                           shp_density_schwarz_s3s1, &
                           shp_density_schwarz(s4,s2), &
                           shp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr .and. temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%ao%get_eri(g, s1, s2, s3, s4, precision_thr/max(temp7,temp8), skip)
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
   end subroutine construct_ao_G_thread_terms_hf
!
!
   module subroutine construct_ao_G_thread_terms_mo_screened_hf(wf, F, D, n_threads, max_D_schwarz,   &
                                          max_eri_schwarz, shp_density_schwarz, n_sig_shp,  &
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
      integer, intent(in) :: n_threads, n_sig_shp
!
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
!
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz 
      real(dp), intent(in) :: coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in) :: shp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, shp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
      real(dp) :: deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: shp_density_schwarz_s1s2, shp_density_schwarz_s3s2, shp_density_schwarz_s3s1
!
      real(dp), dimension(wf%ao%max_sh_size**4) :: g
!
      integer :: thread = 0, skip, q
!
      type(timings), allocatable :: C_max_timer
!
      real(dp), dimension(:), allocatable :: C_max
!
      call mem%alloc(C_max, wf%ao%n_sh) 
!
      C_max_timer = timings('MO coefficients maximums screening vector', 'n')
      call C_max_timer%turn_on()
!
      do s1 = 1, wf%ao%n_sh 
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
!$omp shp_eri_schwarz_s1s2, shp_density_schwarz_s1s2, deg_12_34,                              &
!$omp shp_density_schwarz_s3s2, shp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_shp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%ao%n ! Start column of thread's Fock matrix
!
         shp_eri_schwarz_s1s2 = wf%ao%cs_eri_max(s1s2, 1)
!
         s1s2_packed = wf%ao%cs_eri_max_indices(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = wf%ao%cs_eri_max_indices(s1s2_packed, 2)
!
         shp_density_schwarz_s1s2 = shp_density_schwarz(s1, s2)
         if (shp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            shp_density_schwarz_s3s2 = shp_density_schwarz(s3, s2)
            shp_density_schwarz_s3s1 = shp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = shp_eri_schwarz_s1s2*wf%ao%cs_eri_max(s3s4, 2)
!
               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(shp_density_schwarz(s3,s4)*C_max(s1)*C_max(s2),  &
                           shp_density_schwarz_s1s2*C_max(s3)*C_max(s4))
!
               temp8 = max(shp_density_schwarz_s3s2*C_max(s1)*C_max(s4),    &
                           shp_density_schwarz_s3s1*C_max(s2)*C_max(s4),    &
                           shp_density_schwarz(s4,s2)*C_max(s1)*C_max(s3),  &
                           shp_density_schwarz(s1,s4)*C_max(s2)*C_max(s3))
!
               if (temp8*temp .lt. exchange_thr .and. temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%ao%get_eri(g, s1, s2, s3, s4, precision_thr/max(temp7,temp8), skip)
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
      call mem%dealloc(C_max, wf%ao%n_sh) 
!
   end subroutine construct_ao_G_thread_terms_mo_screened_hf
!
!
   module subroutine construct_coulomb_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,     &
                                                   shp_density_schwarz, &
                                                   n_sig_shp, coulomb_thr, precision_thr, shells)
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
      integer, intent(in) :: n_threads,  n_sig_shp
!
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
!
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, precision_thr
!
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in)               :: shp_density_schwarz
!
      real(dp) :: d1, d2, shp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp7, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: shp_density_schwarz_s1s2, shp_density_schwarz_s3s2, shp_density_schwarz_s3s1
!
      real(dp), dimension(wf%ao%max_sh_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp1, temp2, d1, d2, thread, thread_offset,                            &
!$omp temp7, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,                                &
!$omp shp_eri_schwarz_s1s2, shp_density_schwarz_s1s2, deg_12_34,                            &
!$omp shp_density_schwarz_s3s2, shp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_shp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%ao%n ! Start column of thread's Fock matrix
!
         shp_eri_schwarz_s1s2 = wf%ao%cs_eri_max(s1s2, 1)
!
         s1s2_packed = wf%ao%cs_eri_max_indices(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = wf%ao%cs_eri_max_indices(s1s2_packed, 2)
!
         shp_density_schwarz_s1s2 = shp_density_schwarz(s1, s2)
         if (shp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            shp_density_schwarz_s3s2 = shp_density_schwarz(s3, s2)
            shp_density_schwarz_s3s1 = shp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = shp_eri_schwarz_s1s2*wf%ao%cs_eri_max(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(shp_density_schwarz(s3,s4), shp_density_schwarz_s1s2)
!
               if (temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%ao%get_eri(g, s1, s2, s3, s4, precision_thr/temp7, skip)
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
                                          shp_density_schwarz, &
                                           n_sig_shp, exchange_thr, precision_thr, shells)
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
      integer, intent(in) :: n_threads,  n_sig_shp
!
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
!
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, exchange_thr, precision_thr
!
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in) :: shp_density_schwarz
!
      real(dp) :: d3, d4, d5, d6, shp_eri_schwarz_s1s2
      real(dp) :: temp, temp3, temp4, temp5, temp6, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: shp_density_schwarz_s1s2, shp_density_schwarz_s3s2, shp_density_schwarz_s3s1
!
      real(dp), dimension(wf%ao%max_sh_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp3, d3, d4, d5, d6, thread, thread_offset,                           &
!$omp temp4, temp5, temp6, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,           &
!$omp shp_eri_schwarz_s1s2, shp_density_schwarz_s1s2, deg_12_34,                            &
!$omp shp_density_schwarz_s3s2, shp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_shp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%ao%n ! Start column of thread's Fock matrix
!
         shp_eri_schwarz_s1s2 = wf%ao%cs_eri_max(s1s2, 1)
!
         s1s2_packed = wf%ao%cs_eri_max_indices(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = wf%ao%cs_eri_max_indices(s1s2_packed, 2)
!
         shp_density_schwarz_s1s2 = shp_density_schwarz(s1, s2)
         if (shp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. exchange_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            shp_density_schwarz_s3s2 = shp_density_schwarz(s3, s2)
            shp_density_schwarz_s3s1 = shp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = shp_eri_schwarz_s1s2*wf%ao%cs_eri_max(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. exchange_thr) cycle ! Screened out shell pair
!
               temp8 = max(shp_density_schwarz_s3s2,   &
                           shp_density_schwarz_s3s1,   &
                           shp_density_schwarz(s4,s2), &
                           shp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%ao%get_eri(g, s1, s2, s3, s4, precision_thr/temp8, skip)
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
                           wxyz = shells(s1)%length*(shells(s2)%length*&
                                  (shells(s3)%length*(z_red-1)+y_red-1)+x_red-1)+w_red
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
   module subroutine construct_shp_density_schwarz_hf(wf, shp_density_schwarz, D)
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
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh) :: shp_density_schwarz
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
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
      call mem%alloc(D_red,wf%ao%max_sh_size**2,n_threads)
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D_red_p, maximum, thread) schedule(dynamic)
      do s1 = 1, wf%ao%n_sh
         do s2 = 1, s1
!
!$          thread = omp_get_thread_num()
!
            A_interval = wf%ao%shells(s1)
            B_interval = wf%ao%shells(s2)
!
            D_red_p(1:A_interval%length,1:B_interval%length) => D_red(1:A_interval%length*B_interval%length,thread+1)
!
            D_red_p = D(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            maximum = get_abs_max(D_red_p, (A_interval%length)*(B_interval%length))
!
            nullify(D_red_p)
!
            shp_density_schwarz(s1, s2) = maximum
            shp_density_schwarz(s2, s1) = maximum
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_red,wf%ao%max_sh_size**2,n_threads)
!
   end subroutine construct_shp_density_schwarz_hf
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
      real(dp), dimension(wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers) :: G_ao
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D_ao
!
      real(dp), dimension((wf%ao%max_sh_size**4)*3*4), target :: g_ABCDqk
      real(dp), dimension(:,:,:,:,:,:), pointer, contiguous :: g_ABCDqk_p
!
      integer :: A, B, C, D, D_max, w, x, y, z, AB, AB_packed
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
      call mem%alloc(G_ao_t, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers, n_threads)
      G_ao_t = zero
!
!$omp parallel do private(A, B, C, D, D_max, atoms, deg, deg_CD, deg_AB, deg_AB_CD, g_ABCDqk, g_ABCDqk_p, &
!$omp w, x, y, z, w_red, x_red, y_red, z_red, temp, temp1, temp2, temp3, temp4, temp5, temp6, &
!$omp d1, d2, d3, d4, d5, d6, thread, q, k, tot_dim, AB, AB_packed) schedule(dynamic)
      do AB = 1, wf%ao%n_sig_eri_shp
!
!$       thread = omp_get_thread_num()
!
         if (wf%ao%cs_eri_max(AB, 1)*wf%ao%cs_eri_max(1, 1) < wf%coulomb_threshold) cycle
         AB_packed = wf%ao%cs_eri_max_indices(AB, 3)
!
         A = wf%ao%cs_eri_max_indices(AB_packed, 1)
         B = wf%ao%cs_eri_max_indices(AB_packed, 2)
!
         atoms(1) = wf%ao%shell_to_center(A)
         atoms(2) = wf%ao%shell_to_center(B)
!
         deg_AB = real(2-B/A, kind=dp)
!
            do C = 1, A
!
               D_max = (C/A)*B + (1-C/A)*C
               atoms(3) = wf%ao%shell_to_center(C)
!
               do D = 1, D_max
!
                  deg_CD    = real(2-D/C, kind=dp)
                  deg_AB_CD = min(1-C/A+2-min(D/B,B/D), 2)

                  deg = deg_AB*deg_CD*deg_AB_CD ! Shell degeneracy
!
                  atoms(4) = wf%ao%shell_to_center(D)
!
                  call wf%ao%get_eri_1der(g_ABCDqk, A, B, C, D)
!
                  tot_dim = (wf%ao%shells(A)%length)*(wf%ao%shells(B)%length)&
                              *(wf%ao%shells(C)%length)*(wf%ao%shells(D)%length)&
                              *3*4
!
                  g_ABCDqk(1:tot_dim) = deg*g_ABCDqk(1:tot_dim)
!
                  g_ABCDqk_p(1 : wf%ao%shells(A)%length, 1 : wf%ao%shells(B)%length, &
                             1 : wf%ao%shells(C)%length, 1 : wf%ao%shells(D)%length, &
                             1 : 3, 1 : 4) => g_ABCDqk(1 : tot_dim)
!
                  do z = wf%ao%shells(D)%first, wf%ao%shells(D)%last
!
                     z_red = z - wf%ao%shells(D)%first + 1
!
                     do y = wf%ao%shells(C)%first, wf%ao%shells(C)%last
!
                        y_red = y - wf%ao%shells(C)%first + 1
!
                        d1 = D_ao(y, z)
!
                        do x = wf%ao%shells(B)%first, wf%ao%shells(B)%last
!
                           x_red = x - wf%ao%shells(B)%first + 1
!
                           d3 = D_ao(x, y)
                           d5 = D_ao(x, z)
!
                           do w = wf%ao%shells(A)%first, wf%ao%shells(A)%last
!
                              d2 = D_ao(w, x)
                              d4 = D_ao(w, y)
                              d6 = D_ao(w, z)
!
                              w_red = w - wf%ao%shells(A)%first + 1
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
         call daxpy(3*wf%n_atomic_centers*wf%ao%n**2, one, G_ao_t(1,1,1,1,thread), 1, G_ao, 1)
!
      enddo
!
      call mem%dealloc(G_ao_t, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers, n_threads)
!
   end subroutine construct_ao_G_1der_hf
!
!
   module subroutine set_screening_and_precision_thresholds_hf(wf, gradient_threshold)
!!
!!    Set screening and precision thresholds
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the screening thresholds for Coulomb and exchange
!!    integrals given the convergence threshold for the gradient
!!    unless other thresholds are already set on input.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), intent(in) :: gradient_threshold
!
      real(dp) :: epsilon ! Libint ERI precision 
!
!     If not specified by user, set Coulomb and exchange thresholds 
!     (For ref., this will give 10-12 and 10-10 for gradient threshold: 10-6)
!
      if (.not. input%is_keyword_present('coulomb threshold', 'solver scf')) then 
!
         wf%coulomb_threshold  = gradient_threshold*1.0d-6
!
      endif 
!
      if (.not. input%is_keyword_present('exchange threshold', 'solver scf')) then 
!
         wf%exchange_threshold = gradient_threshold*1.0d-4
!
      endif
!
!     If not requested by user, set Libint integral accuracy according to thresholds 
!
      if (.not. input%is_keyword_present('integral precision', 'solver scf')) then 
!
!        Tighten the default threshold if it is larger than 
!        the lowest screening threshold squared
!
         epsilon = min(wf%coulomb_threshold, wf%exchange_threshold)**2
!
      else
!
!        Warn user about tampering with integral precision (Libint epsilon)
!
         call output%warning_msg("Setting a specific integral precision can cause instabilities.&
                                 & Be aware that the 'precision' given to Libint is not exact&
                                 & but approximate. Don't use unless you know what you are doing.")
!
      endif 
!
!     If not requested by user, set integral cutoff
!
      if (.not. input%is_keyword_present('integral cutoff', 'solver scf')) then
!
         wf%integral_cutoff = sqrt(epsilon)
!
      endif
!
!     Currently, exchange threshold must be >= coulomb threshold in the Fock construction 
!
      if (wf%exchange_threshold .lt. wf%coulomb_threshold) then 
!
         call output%printf('m', 'Exchange threshold is restricted to being higher than, or&
                                 & equal, to the Coulomb threshold. Setting them equal to&
                                 & each other.')
!
         wf%exchange_threshold = wf%coulomb_threshold
!
      endif
!
      call wf%ao%set_libint_epsilon(epsilon)
!
   end subroutine set_screening_and_precision_thresholds_hf
!
!
end submodule ao_fock
