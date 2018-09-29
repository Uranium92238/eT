submodule (ccsd_class) jacobian_ccsd
!
!!
!!    Jacobian submodule (ccsd)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Andreas Skeidsvoll, 2018
!!
!
   implicit none
!
!
contains
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai^B1 = sum_bj F_jb (2*c_ai_bj  -  c_aj_bi)
!!              = sum_bj F_jb v_ai_jb
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      real(dp), dimension(:,:), allocatable :: v_ai_jb
!
!     Construct v_aibj = 2*c_aibj - c_ajbi
!
!               v_ai_jb(ai,jb) = 2*c_ai_bj(ai,bj) - c_ai_bj(aj,bi)
!                        1234              1243             1342
!
      call mem%alloc(v_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      v_ai_jb = zero
!
      call add_1243_to_1234(two, c_ai_bj, v_ai_jb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one, c_ai_bj, v_ai_jb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  v_ai_jb,           &
                  (wf%n_o)*(wf%n_v), &
                  wf%fock_ia,        & ! F_jb
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_a_i,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(v_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    rho_ai^C1 = - sum_bjk L_jikb c_aj_bk
!!              = - sum_bjk (2*g_jikb - g_kijb) c_aj_bk
!!
!!    The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!    where c_ai_bj(ai,bj) = c_aibj above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      real(dp), dimension(:,:), allocatable :: g_ji_kb
      real(dp), dimension(:,:), allocatable :: L_jb_ki
      real(dp), dimension(:,:), allocatable :: c_a_jbk
!
!     Construct the integral g_ji_kb = sum_J L_ji_J * L_kb_J
!
      call mem%alloc(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call get_ooov(g_ji_kb)
!
!     Constructing L_jikb = 2*g_jikb - g_kijb
!
!                  L_jb_ki(jb,ki) = 2*g_ji_kb(ji,kb) - g_ji_kb(ki,jb)
!                           1234              1432             3412
!
!
      call mem%alloc(L_jb_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
      L_jb_ki = zero
!
      call add_1432_to_1234(two, g_ji_kb, L_jb_ki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call add_3412_to_1234(-one, g_ji_kb, L_jb_ki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ji_kb, (wf%n_o)**2, (wf%n_o)*(wf%n_v))
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  c_ai_bj,                & ! c_a_jbk
                  wf%n_v,                 &
                  L_jb_ki,                & ! L_jkb_i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  rho_a_i,                &
                  wf%n_v)
!
      call mem%dealloc(L_jb_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
   module subroutine jacobian_ccsd_a2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^A2 = sum_c g_aibc c_cj - sum_k g_aikj c_bk
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^A2,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      real(dp), dimension(:,:), allocatable :: g_ai_kj ! g_aikj
      real(dp), dimension(:,:), allocatable :: g_ai_jk ! g_aikj reordered
      real(dp), dimension(:,:), allocatable :: g_ai_bc ! g_aibc, batching over b
!
      real(dp), dimension(:,:), allocatable :: rho_ba_ij ! rho_ai_bj, term 1 (see below)
!
!     Batching variables
!
      integer(i15) :: required = 0
      integer(i15) :: current_b_batch = 0
      integer(i15) :: aib_offset = 0
!
      type(batching_index) :: batch_b
!
!     :: Term 1. - sum_k g_aikj c_bk ::
!
!     Calculate g_ai_kj
!
      call mem%alloc(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call wf%get_vooo(g_ai_kj)
!
!     Reorder to g_ai_jk = g_aikj = g_ai_kj
!
      call mem%alloc(g_ai_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call sort_1234_to_1243(g_ai_kj, g_ai_jk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ai_kj, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call mem%alloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_a_i,                &
                  wf%n_v,               &
                  g_ai_jk,              & ! g_aij_k
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  rho_ba_ij,            & ! rho_b_aij
                  wf%n_v)
!
      call mem%dealloc(g_ai_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      call add_3124_to_1234(one, rho_ba_ij, rho_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ba_ij, (wf%n_v)**2, (wf%n_o)**2)
!
!     :: Term 2. rho_ai_bj =+ sum_c g_aibc c_cj ::
!
!     We do the matrix multiplication as g_aib_c c_cj,
!     batching over the b index.
!
      required = wf%integrals%get_required_vvvo() + (wf%n_v)*(wf%n_o)*(wf%n_v)*(batch_b%length)
!
!     Initialize batching variable
!
      call batch_b%init(wf%n_v)
      call mem%num_batch(batch_b, required)
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Determine limits for the current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
         call mem%alloc(g_ai_bc, (wf%n_v)*(wf%n_o), (wf%n_v)*(batch_b%length))
!
         call get_vovv(g_ai_bc,       &
                       1,             &
                       wf%n_v,        &
                       1,             &
                       wf%n_o,        &
                       batch_b%first, &
                       batch_b%last,  &
                       1,             &
                       wf%n_v)
!
!        Calculate the contribution to rho_aib_j = sum_c g_aib_c c_cj
!
         aib_offset = index_three(1, 1, batch_b%first, wf%n_v, wf%n_o)
!
         call dgemm('N', 'N',                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_ai_bc,                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_a_i,                              &
                     wf%n_v,                             &
                     one,                                &
                     rho_ai_bj(aib_offset,1),            &
                     (wf%n_o)*(wf%n_v)**2)
!
          call mem%dealloc(g_ai_bc, (wf%n_v)*(wf%n_o), wf%n_v*(batch_b%length))
!
      enddo ! End of batches over b
!
   end subroutine jacobian_ccsd_a2_ccsd
!
!
end submodule jacobian_ccsd
