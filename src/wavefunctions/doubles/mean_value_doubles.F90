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
submodule (doubles_class) mean_value_doubles
!
!!
!!    Mean-value submodule
!!
!!    Contains routines related to the mean values, i.e.
!!    the construction of density matrices as well as expectation
!!    value calculation.
!!
!!    The ground state density is constructed as follows:
!!
!!          D_pq = < Lambda| E_pq |CC >
!!    where:
!!          < Lambda| = < HF| + sum_mu tbar_mu < mu| exp(-T)
!!
!!
!!    In general a CC density matrix can be written as:
!!
!!          D_pq = < X| e^(-T) E_pq e^T |Y >
!!
!!    where X and Y are left and right vectors with contributions from
!!    a reference determinant and excited determinants (< mu|, |nu >):
!!
!!          D_pq =             X_ref < HF| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_mu  < mu| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_ref < HF| e^(-T) E_pq e^T |mu >  Y_mu
!!                 + sum_mu,nu X_mu  < mu| e^(-T) E_pq e^T |nu >  Y_nu
!!
!!    Depending on the type of density matrix (Ground state, transition ,
!!    excited state, interstate transition) different states and thus different
!!    amplitudes X_ref, X_mu, Y_ref and Y_mu will contribute.
!!
!!    In EOM theory the states can be written as the following vectors:
!!
!!          |CC >     = R_0 = (1, 0)
!!          |Lambda > = L_0 = (1, tbar_mu)
!!          |R_k >    = R_k = (-sum_mu(tbar_mu*R_mu), R_mu)
!!          |L_k >    = L_k = (0, L_mu)
!!
!!    The routine names derive from the contribution of the vectors:
!!
!!       ref_ref: first component of the vector for the left and right state
!!
!!       mu_ref:  second component of the vector for the left and
!!                first component of the vector for the right state
!!
!!       ref_mu:  first component of the vector for the left and
!!                second component of the vector for the right state
!!
!!       mu_nu:   second component of the vector for the left and right state
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_gs_density_doubles(wf)
!!
!!    Construct GS density
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    Constructs the one-electron density matrix in the T1 basis
!!
!!    D_pq = < Lambda| E_pq |CC >
!!
!!    Contributions to the density are split up as follows:
!!    D_pq = D_pq(ref-ref) + sum_mu tbar_mu D_pq(mu-ref)
!!
!!    The second term is separated in "mu_ref_density_terms",
!!    as it is used for the left transition density as well.
!!
      implicit none
!
      class(doubles) :: wf
      type(timings)  :: timer
!
      real(dp), dimension(:), allocatable :: tbar
!
      timer = timings('Ground state density', pl='m')
      call timer%turn_on
!
      call mem%alloc(tbar, wf%n_t1 + wf%n_t2)
      call wf%get_full_multipliers(tbar)
!
      call wf%mu_ref_density_terms(wf%density, 0, tbar)
!
      call mem%dealloc(tbar, wf%n_t1 + wf%n_t2)
!
      call wf%density_ccs_ref_ref_oo(wf%density)
!
      call timer%turn_off
!
   end subroutine construct_gs_density_doubles
!
!
   module subroutine mu_ref_density_terms_doubles(wf, density, state, L)
!!
!!    Density mu ref terms
!!    Written by Alexander C. Paul, May 2021
!!
!!    Constructs terms of the form:
!!       sum_mu L_mu < mu| E_pq |HF >
!!
!!    corresponding to terms of the ground state density
!!    and the left transition density.
!!
      use array_initialization, only: zero_array
      use reordering, only: squareup
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
!
      integer, intent(in) :: state
!
!     L might only be contiguous in the ranges (1:n_t1) and (1+t1: n_t1+n_t2)
!     as L can also be a combined array of t1bar + t2bar
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in) :: L
!
      real(dp), dimension(:,:,:,:), allocatable :: L2, t2
!
      type(timings)     :: timer
      character(len=40) :: timer_name
!
      write(timer_name, '(a,i0,a)') 'Doubles contribution to <', state,'|E_pq|0>'
      timer = timings(trim(timer_name), pl='v')
!
      call zero_array(density, wf%n_mo**2)
      call wf%ccs%mu_ref_density_terms(density, state, L(1:wf%n_t1))
!
      call timer%turn_on() ! Doubles contribution
!
      call mem%alloc(t2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t2, wf%n_v*wf%n_o)
!
      call wf%density_doubles_mu_ref_ov(density, L(1:wf%n_t1), t2)
!
      call mem%alloc(L2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(L(wf%n_t1+1 : wf%n_t1+wf%n_t2), L2, wf%n_t1)
!
      call wf%density_doubles_mu_ref_oo(density, L2, t2)
      call wf%density_doubles_mu_ref_vv(density, L2, t2)
!
      call mem%dealloc(t2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine mu_ref_density_terms_doubles
!
!
   module subroutine density_doubles_mu_ref_oo_doubles(wf, density, tbar_akbj, t_akbi)
!!
!!    One electron density excited-determinant/reference oo-term
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    explicit term in this routine:
!!          D_ij -= sum_abk t_akb,i tbar_akb,j
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_akbj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_akbi
!
      call dgemm('T', 'N',          &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_v**2*wf%n_o, &
                  -one,             &
                  t_akbi,           & ! t_akb_i
                  wf%n_v**2*wf%n_o, &
                  tbar_akbj,        & ! tbar_akb_j
                  wf%n_v**2*wf%n_o, &
                  one,              &
                  density,          &
                  wf%n_mo)
!
   end subroutine density_doubles_mu_ref_oo_doubles
!
!
   module subroutine density_doubles_mu_ref_vv_doubles(wf, density, tbar_ajci, t_bjci)
!!
!!    One electron density excited-determinant/reference vv-term
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!          D_ab += sum_jci tbar_a,jci t_b,jci
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_ajci
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_bjci
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_o**2*wf%n_v, &
                  one,              &
                  tbar_ajci,        & ! tbar_a_jci
                  wf%n_v,           &
                  t_bjci,           & ! t_b_jci
                  wf%n_v,           &
                  one,              &
                  density(wf%n_o + 1, wf%n_o + 1), &
                  wf%n_mo)
!
   end subroutine density_doubles_mu_ref_vv_doubles
!
!
   module subroutine density_doubles_mu_ref_ov_doubles(wf, density, tbar_ai, t_aibj)
!!
!!    One electron density excited-determinant/reference ov-term
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!          D_ia += sum_bj u^ab_ij tbar_bj = sum_bj u_ia,bj tbar_bj
!!
!!          u^{ab}_ij = 2t_aibj - t_ajbi
!!
      use reordering, only: add_1432_to_1234
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: u_aibj
!
      real(dp), dimension(:,:), allocatable :: D_ov
!
      integer :: i, a
!
      call mem%alloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_v**2*wf%n_o**2, t_aibj, 1, u_aibj, 1)
      call dscal(wf%n_v**2*wf%n_o**2, two, u_aibj, 1)
!
      call add_1432_to_1234(-one, t_aibj, u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(D_ov, wf%n_v, wf%n_o) ! ordered v,o
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  u_aibj,        & ! u_bj_ai
                  wf%n_o*wf%n_v, &
                  tbar_ai,       & ! tbar_ai
                  1,             &
                  zero,          &
                  D_ov,          & ! D_bj
                  1)
!
      call mem%dealloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            density(i, wf%n_o + a) = density(i, wf%n_o + a) + D_ov(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_ov, wf%n_v, wf%n_o)
!
   end subroutine density_doubles_mu_ref_ov_doubles
!
!
end submodule mean_value_doubles
