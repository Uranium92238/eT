submodule (ccsd_class) omega_ccsd
!
!!
!!    Omega submodule (ccsd)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!
   implicit none
!
!
!
!
contains
!
!

   module subroutine omega_ccsd_a1_ccsd(wf, omega1)
!!
!!    Omega A1 term
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd g_adkc * u_ki^cd,
!!
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf.
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
!     Batching variables
!
      integer(i15) :: required        = 0
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      real(dp), dimension(:,:), allocatable :: u_dk_ci, t_dk_ci, g_ad_kc
!
      integer(i15) :: ad_dim
!
!     u_ki^cd = 2*t_ki^cd - t_ik^cd (ordered as u_dk_ci)
!
      call mem%alloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_dk_ci, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      u_dk_ci = - t_dk_ci
!
      call add_1432_to_1234(two, t_dk_ci, u_dk_ci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Prepare for batching
!
!     Estimated memory required to construct g_adkc
!
    !  required = wf%integrals%get_required_vvov()
!
!     Initialization of the batching variable
!
      call batch_a%init(wf%n_v)                ! Initialize batching index a
      call mem%num_batch(batch_a, required) ! Determine batching information
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index
!
         call batch_a%determine_limits(current_a_batch)
!
         ad_dim = (batch_a%length)*(wf%n_v)
!
!        Form g_ad_kc = g_adkc
!
         call mem%alloc(g_ad_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
         call wf%get_vvov(g_ad_kc,                       &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
         call dgemm('N','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     g_ad_kc,                    & ! g_a_dkc
                     batch_a%length,             &
                     u_dk_ci,                    & ! u_dkc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     omega1(batch_a%first,1),    &
                     wf%n_v)
!
         call mem%dealloc(g_ad_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
      enddo ! End of batches of the index a
!
      call mem%dealloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_a1_ccsd
!
   module subroutine omega_ccsd_b1_ccsd(wf, omega1)
!!
!!    Omega B1
!!    Written by Sarai D. Fokestad, Eirik F. Kjønstad,
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!!    Calculates the B1 term,
!!
!!       B1:   - sum_ckl u_kl^ac * g_kilc,
!!
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
      real(dp), dimension(:,:), allocatable :: g_lc_ki ! g_kilc
      real(dp), dimension(:,:), allocatable :: t_al_ck ! t_kl^ac
      real(dp), dimension(:,:), allocatable :: u_al_ck ! u_kl^ac = 2 t_kl^ac - t_lk^ac
!
!     Get g_ki_lc = g_kilc
!
      call mem%alloc(g_lc_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
      call wf%get_ovoo(g_lc_ki,    &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_o)
!
!     Form u_al_ck = u_kl^ac = 2 * t_kl^ac - t_lk^ac
!
!     Squareup amplitudes and reorder: t_ak_cl to t_al_ck
!     u_al_ck = 2 * t_kl^ac - t_lk^ac = 2 * t_al_ck(al,ck) - t_al_ck(ak,cl)
!
      call mem%alloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*wf%n_o)
      call squareup_and_sort_1234_to_1432(wf%t2, t_al_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(u_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      u_al_ck = zero
!
      call add_1432_to_1234(-one, t_al_ck, u_al_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_al_ck, 1, u_al_ck, 1)
!
      call mem%dealloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!
      call dgemm('N','N',                 &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  u_al_ck,                & ! u_a_lck
                  wf%n_v,                 &
                  g_lc_ki,                & ! g_lck_i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  omega1,                 &
                  wf%n_v)
!
!     Deallocate remaining vectors
!
      call mem%dealloc(u_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(g_lc_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
   end subroutine omega_ccsd_b1_ccsd
!
!
   module subroutine omega_ccsd_c1_ccsd(wf, omega1)
!!
!!    Omega C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!!    Calculates the C1 term of omega,
!!
!!       C1: sum_ck F_kc*u_ai_ck,
!!
!!    and adds it to the projection vector (omega1) of
!!    the wavefunction object wf
!!
!!    u_ai_kc = 2*t_ck_ai - t_ci_ak
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
      real(dp), dimension(:,:), allocatable :: F_c_k ! F_kc
!
      real(dp), dimension(:,:), allocatable :: u_ai_ck
      real(dp), dimension(:,:), allocatable :: t_ai_ck
!
      integer(i15) :: k, c
!
!     Form u_ai_ck = u_ik^ac = 2*t_ik^ac - t_ki^ac
!
      call mem%alloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_ai_ck, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      u_ai_ck = zero
!
      call add_1432_to_1234(-one, t_ai_ck, u_ai_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ai_ck, 1, u_ai_ck, 1)
!
      call mem%dealloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder the Fock matrix F_ck = F_kc
!
      call mem%alloc(F_c_k, wf%n_v, wf%n_o)
!
      call sort_12_to_21(wf%fock_ia, F_c_k, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  F_c_k,             & ! Equal to F_kc
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega1,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(F_c_k, wf%n_v, wf%n_o)
!
   end subroutine omega_ccsd_c1_ccsd
!
   module subroutine omega_ccsd_a2_ccsd(wf, omega2)
!!
!!    Omega A2 term
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!      
!!    A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!!
!!    Structure: Batching over both a and b for A2.2.
!!                t^+_ci_dj = t_ci_dj + t_di_cj
!!                t^-_ci_dj = t_ci_dj - t_di_cj
!!                g^+_ac_bd = g_ac_bd + g_bc_ad
!!                g^-_ac_bd = g_ac_bd - g_bc_ad
!!
!!                omega_A2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!!                omega_A2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)), intent(inout):: omega2
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: g_aibj
      real(dp), dimension(:,:), allocatable :: g_ac_bd
      real(dp), dimension(:,:), allocatable :: g_p_ab_cd
      real(dp), dimension(:,:), allocatable :: g_m_ab_cd
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_p_cd_ij
      real(dp), dimension(:,:), allocatable :: t_m_cd_ij
!
!     Reordered omega 2
!
      real(dp), dimension(:,:), allocatable :: omega2_p_ab_ij
      real(dp), dimension(:,:), allocatable :: omega2_m_ab_ij
!
!     Indices
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ab = 0, ca = 0, bc = 0, cd = 0, ad = 0, db = 0, ac = 0, bd = 0
      integer(i15) :: ai = 0, aj = 0, bj = 0, bi = 0, ci = 0, cj = 0, dj = 0, di = 0
      integer(i15) :: ij = 0
!
      integer(i15) :: aibj = 0, biaj = 0, cidj = 0, dicj = 0
!
!     Batching and memory handling variables
!
      integer(i15) :: required = 0
!
      integer(i15) :: current_a_batch = 0
      integer(i15) :: current_b_batch = 0
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
      integer(i15) :: threads = 0
      integer(i15) :: omp_get_num_threads
!
!     Timing variables
!
      real(dp) :: time_non_integral_part
      real(dp) :: a2_begin_time, a2_end_time
      real(dp) :: a2_begin_time_1, a2_end_time_1
!
      time_non_integral_part = zero
!
      call cpu_time(a2_begin_time)
!
!     :: Calculate the A2.1 term of omega ::
!
!     Create g_ai_bj
!
      call mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_vovo(g_ai_bj,    &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_o)
!
!     Add A2.1 to Omega 2
!
      call add_to_packed(omega2, g_ai_bj, (wf%n_o)*(wf%n_v))

      call mem%dealloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!    ::  Calculate the A2.2 term  of omega ::
!
      !required = wf%get_vvvv_required_mem() + 4*(wf%n_o**2)*(wf%n_v**2) + 2*(wf%n_v**4)
!
!     Initialize batching variables
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%num_two_batch(batch_a, batch_b, required)
!
!     Start looping over a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Determine the limits for the a-batch
!
         call batch_a%determine_limits(current_a_batch)
!
!        Start looping over b-batches
!
         do current_b_batch = 1, batch_b%num_batches
!
!           Determine the limits for the b-batch
!
            call batch_b%determine_limits(current_b_batch)
!
!           Allocate g_ca_db
!
            call mem%alloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!           Get g_ac_bd
!
            call cpu_time(a2_begin_time_1)
!
            call wf%get_vvvv(g_ac_bd,        &
                              batch_a%first, &
                              batch_a%last,  &
                              1,             &
                              wf%n_v,        &
                              batch_b%first, &
                              batch_b%last,  &
                              1,             &
                              wf%n_v)
!
            call cpu_time(a2_end_time_1)
            time_non_integral_part = time_non_integral_part &
                                   - a2_end_time_1 + a2_begin_time_1
!
            if (current_b_batch .eq. current_a_batch) then
!
!              Allocate for +-g, +-t
!
               call mem%alloc(g_p_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
               call mem%alloc(g_m_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
               call mem%alloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
               call mem%alloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
               g_p_ab_cd = zero
               g_m_ab_cd = zero
               t_p_cd_ij = zero
               t_m_cd_ij = zero
!
!              Reorder g_ca_db to g_ab_cd and t_ci_dj to t_cd_ij
!
!$omp parallel do schedule(static) private(d,b,a,ac,cd,bd,bc,ab,ad,i,j,ij,ci,dj,cj,di,cidj,dicj)
               do c = 1, wf%n_v
!
                  do d = 1, c
!
                     cd = index_packed(c, d)
!
                     do a = 1, batch_a%length
!
                        ac = index_two(a, c, batch_a%length)
                        ad = index_two(a, d, batch_a%length)
!
                        do  b = 1, batch_b%length
                           if ((a+batch_a%first-1) .ge. (b+batch_b%first-1)) then
!
                              bd = index_two(b, d, batch_b%length)
                              bc = index_two(b, c, batch_b%length)
!
                              ab = index_packed(a, b)
!
                              g_p_ab_cd(ab, cd) = g_ac_bd(ac, bd) + g_ac_bd(ad, bc)
                              g_m_ab_cd(ab, cd) = g_ac_bd(ac, bd) - g_ac_bd(ad, bc)
!
                             if(c .ne. d) then
                               g_p_ab_cd(ab, cd) = two*g_p_ab_cd(ab, cd)
                               g_m_ab_cd(ab, cd) = two*g_m_ab_cd(ab, cd)
                             endif
!
                           endif
                        enddo
                     enddo
!
                    do i = 1, wf%n_o
                       do j = 1, i
!
                          ij = index_packed(i, j)
!
                           ci = index_two(c, i, wf%n_v)
                           dj = index_two(d, j, wf%n_v)
                           cj = index_two(c, j, wf%n_v)
                           di = index_two(d, i, wf%n_v)
!
                           cidj = index_packed(ci, dj)
                           dicj = index_packed(cj, di)
!
                           t_p_cd_ij(cd, ij) = wf%t2(cidj, 1) + wf%t2(dicj, 1)
                           t_m_cd_ij(cd, ij) = wf%t2(cidj, 1) - wf%t2(dicj, 1)
!
                       enddo
                    enddo
                 enddo
              enddo
!$omp end parallel do
!
!              Dellocate g_ac_bd
!
               call mem%dealloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!              Allocate omega +-
!
              call mem%alloc(omega2_p_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
              call mem%alloc(omega2_m_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
!
!              omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
!
              call dgemm('N','N',                      &
                          packed_size(batch_a%length), &
                          packed_size(wf%n_o),         &
                          packed_size(wf%n_v),         &
                          one/four,                    &
                          g_p_ab_cd,                   &
                          packed_size(batch_a%length), &
                          t_p_cd_ij,                   &
                          packed_size(wf%n_v),         &
                          zero,                        &
                          omega2_p_ab_ij,              &
                          packed_size(batch_a%length))
!
              call dgemm('N','N',                      &
                          packed_size(batch_a%length), &
                          packed_size(wf%n_o),         &
                          packed_size(wf%n_v),         &
                          one/four,                    &
                          g_m_ab_cd,                   &
                          packed_size(batch_a%length), &
                          t_m_cd_ij,                   &
                          packed_size(wf%n_v),         &
                          zero,                        &
                          omega2_m_ab_ij,              &
                          packed_size(batch_a%length) )
!
!             Deallocate +-g, +-t
!
              call mem%dealloc(g_p_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
              call mem%dealloc(g_m_ab_cd, packed_size(batch_a%length), packed_size(wf%n_v))
              call mem%dealloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
              call mem%dealloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
              do i = 1, wf%n_o
                 do j = 1, i
!
                    ij = index_packed(i, j)
!
                    do a = 1, batch_a%length
!
                       ai = index_two(a + batch_a%first - 1, i, wf%n_v) ! A is full-CCSD-space a index
                       aj = index_two(a + batch_a%first - 1, j, wf%n_v) ! A is full-CCSD-space a index
!
                       do b = 1, batch_b%length
!
                          if ((a+batch_a%first-1) .ge. (b+batch_b%first-1)) then
                             bj = index_two(b + batch_b%first - 1, j, wf%n_v) ! B is full-CCSD-space b index
                             bi = index_two(b + batch_b%first - 1, i, wf%n_v) ! B is full-CCSD-space b index
!
                             ab = index_packed(a, b)
!
                             aibj = index_packed(ai, bj)
                             biaj = index_packed(bi, aj)
!
!                            Reorder into omega2_aibj
!
                             omega2(aibj,1) = omega2(aibj, 1) &
                                                   + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                             if (aibj .ne. biaj) then
                                omega2(biaj,1) = omega2(biaj, 1) &
                                                   + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
                             endif
                          endif
!
                       enddo
                    enddo
                 enddo
              enddo
!
!              Deallocate omega +-
!
               call mem%dealloc(omega2_p_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
               call mem%dealloc(omega2_m_ab_ij, packed_size(batch_a%length), packed_size(wf%n_o))
            else
!
!              Allocate for +-g, +-t
!
               call mem%alloc(g_p_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call mem%alloc(g_m_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call mem%alloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
               call mem%alloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
               g_p_ab_cd = zero
               g_m_ab_cd = zero
               t_p_cd_ij = zero
               t_m_cd_ij = zero
!
!$omp parallel do schedule(static) private(d,b,a,ac,cd,bd,bc,ab,ad,i,j,ij,ci,dj,cj,di,cidj,dicj)
               do c = 1, wf%n_v
                  do d = 1, c
!
                     cd = index_packed(c, d)
!
                     do a = 1, batch_a%length
!
                        ac = index_two(a, c, batch_a%length)
                        ad = index_two(a, d, batch_a%length)
!
                        do  b = 1, batch_b%length
!
                             bd = index_two(b, d, batch_b%length)
                             bc = index_two(b, c, batch_b%length)
!
                             ab = index_two(a, b, batch_a%length)
!
                             g_p_ab_cd(ab, cd) = g_ac_bd(ac, bd) + g_ac_bd(ad, bc)
                             g_m_ab_cd(ab, cd) = g_ac_bd(ac, bd) - g_ac_bd(ad, bc)
!
                            if(c .ne. d) then
                              g_p_ab_cd(ab, cd) = two*g_p_ab_cd(ab, cd)
                              g_m_ab_cd(ab, cd) = two*g_m_ab_cd(ab, cd)
                            endif
!
                       enddo
                    enddo
!
                    do i = 1, wf%n_o
                       do j = 1, i
!
                           ij = index_packed(i, j)
!
                           ci = index_two(c, i, wf%n_v)
                           dj = index_two(d, j, wf%n_v)
                           cj = index_two(c, j, wf%n_v)
                           di = index_two(d, i, wf%n_v)
!
                           cidj = index_packed(ci, dj)
                           dicj = index_packed(cj, di)
!
                           t_p_cd_ij(cd, ij) = wf%t2(cidj, 1) + wf%t2(dicj, 1)
                           t_m_cd_ij(cd, ij) = wf%t2(cidj, 1) - wf%t2(dicj, 1)
!
                       enddo
                    enddo
                 enddo
              enddo
!$omp end parallel do
!
!              Dellocate g_ac_bd
!
               call mem%dealloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
!              Allocate omega +-
!
               call mem%alloc(omega2_p_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
               call mem%alloc(omega2_m_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
!
!               omega2_ab_ij = sum_(cd) g_ab_cd*t_cd_ij
!
              call dgemm('N','N',                            &
                          (batch_a%length)*(batch_b%length), &
                          packed_size(wf%n_o),               &
                          packed_size(wf%n_v),               &
                          one/four,                          &
                          g_p_ab_cd,                         &
                          (batch_a%length)*(batch_b%length), &
                          t_p_cd_ij,                         &
                          packed_size(wf%n_v),               &
                          zero,                              &
                          omega2_p_ab_ij,                    &
                          (batch_a%length)*(batch_b%length))
!
              call dgemm('N','N',                            &
                          (batch_a%length)*(batch_b%length), &
                          packed_size(wf%n_o),               &
                          packed_size(wf%n_v),               &
                          one/four,                          &
                          g_m_ab_cd,                         &
                          (batch_a%length)*(batch_b%length), &
                          t_m_cd_ij,                         &
                          packed_size(wf%n_v),               &
                          zero,                              &
                          omega2_m_ab_ij,                    &
                          (batch_a%length)*(batch_b%length))
!
!              Deallocate +-g, +-t
!
               call mem%dealloc(g_p_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call mem%dealloc(g_m_ab_cd, (batch_a%length)*(batch_b%length), packed_size(wf%n_v))
               call mem%dealloc(t_p_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
               call mem%dealloc(t_m_cd_ij, packed_size(wf%n_v), packed_size(wf%n_o))
!
               do i = 1, wf%n_o
                  do j = 1, i
!
                     ij = index_packed(i, j)
!
                     do a = 1, batch_a%length
!
                        ai = index_two(a + batch_a%first - 1, i, wf%n_v) ! A is full-space a index
                        aj = index_two(a + batch_a%first - 1, j, wf%n_v) ! A is full-space a index
!
                        do b = 1, batch_b%length
!
                              bj = index_two(b + batch_b%first - 1, j, wf%n_v) ! B is full-space b index
                              bi = index_two(b + batch_b%first - 1, i, wf%n_v) ! B is full-space b index
!
                              ab = index_two(a, b, batch_a%length)
!
                              aibj = index_packed(ai, bj)
                              biaj = index_packed(bi, aj)
!
!                             Reorder into omega2_aibj
!
                              omega2(aibj,1) = omega2(aibj, 1) &
                                          + omega2_p_ab_ij(ab, ij) + omega2_m_ab_ij(ab, ij)
!
                              if (aibj .ne. biaj) then
!
                                 omega2(biaj,1) = omega2(biaj, 1) &
                                          + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
!
                              endif
!
                        enddo
                     enddo
                  enddo
               enddo
!
!              Deallocate omega +-
!
               call mem%dealloc(omega2_p_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
               call mem%dealloc(omega2_m_ab_ij, (batch_a%length)*(batch_b%length), packed_size(wf%n_o))
!
            endif
!
         enddo ! End batching over b
      enddo ! End batching over a
!
   end subroutine omega_ccsd_a2_ccsd
end submodule
