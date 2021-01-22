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
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in), optional :: prev_ao_density
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
      class(hf) :: wf
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in) :: prev_ao_density
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
      class(hf) :: wf
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
      class(hf)   :: wf
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in)    :: D
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(inout) :: G
      logical, optional, intent(in) :: C_screening
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
      integer, intent(in) :: n_threads, n_sig_shp
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz 
      real(dp), intent(in) :: coulomb_thr, exchange_thr, precision_thr
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in) :: shp_density_schwarz
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
      integer, intent(in) :: n_threads, n_sig_shp
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz 
      real(dp), intent(in) :: coulomb_thr, exchange_thr, precision_thr
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in) :: shp_density_schwarz
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
      integer, intent(in) :: n_threads,  n_sig_shp
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, precision_thr
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in)               :: shp_density_schwarz
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
      integer, intent(in) :: n_threads,  n_sig_shp
      type(interval), dimension(wf%ao%n_sh), intent(in) :: shells
      real(dp), dimension(wf%ao%n, wf%ao%n*n_threads)   :: F
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, exchange_thr, precision_thr
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh), intent(in) :: shp_density_schwarz
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
      real(dp), dimension(wf%ao%n_sh, wf%ao%n_sh) :: shp_density_schwarz
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
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
      real(dp), dimension(wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers) :: G_ao
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D_ao
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
      real(dp), intent(in) :: gradient_threshold
!
   end subroutine set_screening_and_precision_thresholds_hf
