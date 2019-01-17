submodule (ccsd_class) omega_ccsd
!
!!
!!    Omega submodule (ccsd)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Andreas Skeidsvoll, 2018
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
!
   module subroutine construct_omega_ccsd(wf, omega)
!!
!!    Construct omega (CCSD)
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!    and Andreas Skeidsvoll, 2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: omega1, omega2
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call mem%alloc(omega2, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2, 1)
!
!     Set the omega vector to zero
!
      omega1 = zero
      omega2 = zero
!
!     Construct singles contributions
!
      call wf%omega_ccsd_a1(omega1)
      call wf%omega_ccsd_b1(omega1)
      call wf%omega_ccsd_c1(omega1)
!
      call wf%omega_ccs_a1(omega1)
!
!     Construct doubles contributions
!
      call wf%omega_ccsd_a2(omega2)
      call wf%omega_ccsd_b2(omega2)
      call wf%omega_ccsd_c2(omega2)
      call wf%omega_ccsd_d2(omega2)
      call wf%omega_ccsd_e2(omega2)
!
      call dcopy((wf%n_o)*(wf%n_v), omega1, 1, omega, 1)
      call dcopy((wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v)+1)/2, omega2, 1, omega((wf%n_o)*(wf%n_v)+1, 1), 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
      call mem%dealloc(omega2, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2, 1)
!
   end subroutine construct_omega_ccsd
!
!
   module subroutine omega_ccsd_a1_ccsd(wf, omega1)
!!
!!    Omega A1 term
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!    and Andreas Skeidsvoll, 2018
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
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      real(dp), dimension(:,:), allocatable :: u_dk_ci, t_dk_ci, g_ad_kc
!
      integer(i15) :: ad_dim, rec0, rec1
!
      type(timings) :: ccsd_a1_timer
!
      call ccsd_a1_timer%init('omega ccsd a1')
      call ccsd_a1_timer%start()
!
!     u_ki^cd = 2*t_ki^cd - t_ik^cd (ordered as u_dk_ci)
!
      call mem%alloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_dk_ci, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      u_dk_ci = zero
      call daxpy(((wf%n_v)*(wf%n_o))**2, -one, t_dk_ci, 1, u_dk_ci, 1)
!
      call add_1432_to_1234(two, t_dk_ci, u_dk_ci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Prepare for batching
!
      rec0 = wf%n_o*wf%integrals%n_J*wf%n_v
!
      rec1 = wf%n_v*wf%integrals%n_J + wf%n_v**2*(wf%n_o)
!
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_a, rec0, rec1)
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
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
      call ccsd_a1_timer%freeze()
      call ccsd_a1_timer%switch_off()
!
   end subroutine omega_ccsd_a1_ccsd
!
!
   module subroutine omega_ccsd_b1_ccsd(wf, omega1)
!!
!!    Omega B1
!!    Written by Sarai D. Fokestad, Eirik F. Kjønstad,
!!    and Andreas Skeidsvoll, 2018
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
      type(timings) :: ccsd_b1_timer 
!  
      call ccsd_b1_timer%init('omega ccsd b1')
      call ccsd_b1_timer%start()
!
!     Form u_al_ck = u_kl^ac = 2 * t_kl^ac - t_lk^ac
!     Square up amplitudes and reorder: t_ak_cl to t_al_ck
!
      call mem%alloc(g_lc_ki, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
      call wf%get_ovoo(g_lc_ki)
!
!     u_al_ck = 2 * t_kl^ac - t_lk^ac = 2 * t_al_ck(al,ck) - t_al_ck(ak,cl)
!
      call mem%alloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call squareup(wf%t2, t_al_ck, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      u_al_ck = zero
      call daxpy(((wf%n_v)*(wf%n_o))**2, -one, t_al_ck, 1, u_al_ck, 1)
!
      call add_1432_to_1234(two, t_al_ck, u_al_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
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
      call ccsd_b1_timer%freeze()
      call ccsd_b1_timer%switch_off()
!
   end subroutine omega_ccsd_b1_ccsd
!
!
   module subroutine omega_ccsd_c1_ccsd(wf, omega1)
!!
!!    Omega C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
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
      type(timings) :: ccsd_c1_timer 
!
      call ccsd_c1_timer%init('omega ccsd c1')
      call ccsd_c1_timer%start()
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
      call ccsd_c1_timer%freeze()
      call ccsd_c1_timer%switch_off()
!
   end subroutine omega_ccsd_c1_ccsd
!
!
   module subroutine omega_ccsd_a2_ccsd(wf, omega2)
!!
!!    Omega A2 term
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    and Andreas Skeidsvoll, 2018
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
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
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
      integer(i15) :: a, b, c, d
      integer(i15) :: i, j
!
      integer(i15) :: ab, bc, cd, ad, ac, bd
      integer(i15) :: ai, aj, bj, bi, ci, cj, dj, di
      integer(i15) :: ij
!
      integer(i15) :: aibj = 0, biaj = 0, cidj = 0, dicj = 0
!
!     Batching and memory handling variables
!
      integer(i15) :: rec0, rec1_a, rec1_b, rec2
!
      integer(i15) :: current_a_batch = 0
      integer(i15) :: current_b_batch = 0    
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
      type(timings) :: ccsd_a2_timer, ccsd_a2_integral_timer
!
      call ccsd_a2_timer%init('omega ccsd a2')
      call ccsd_a2_integral_timer%init('omega ccsd a2 g_abcd')
!      
      call ccsd_a2_timer%start()
!
!     :: Calculate the A2.1 term of omega ::
!
!     Create g_ai_bj
!
      call mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_vovo(g_ai_bj)
!
!     Add A2.1 to Omega 2
!
      call add_to_packed(omega2, g_ai_bj, (wf%n_o)*(wf%n_v))

      call mem%dealloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!    ::  Calculate the A2.2 term  of omega ::
!
      rec0 = 2*(packed_size(wf%n_o))*(packed_size(wf%n_v))
!
      rec1_a = wf%integrals%n_J*wf%n_v 
      rec1_b = wf%integrals%n_J*wf%n_v 
!
      rec2 = 2*wf%n_v**2 + 2*(packed_size(wf%n_o))
!
!     Initialize batching variables
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_b, rec0, rec1_a, rec1_b, rec2)
!
!     Start looping over a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_b_batch = 1, batch_b%num_batches
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(g_ac_bd, (wf%n_v)*(batch_a%length), (wf%n_v)*(batch_b%length))
!
            call ccsd_a2_integral_timer%start()
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
            call ccsd_a2_integral_timer%freeze()
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
                  do d = 1, c
!
                     cd = (max(c,d)*(max(c,d)-3)/2) + c + d
!
                     do a = 1, batch_a%length
!
                        ac = batch_a%length*(c - 1) + a
                        ad = batch_a%length*(d - 1) + a 
!
                        do  b = 1, batch_b%length
                           if ((a+batch_a%first-1) .ge. (b+batch_b%first-1)) then
!
                              bd = batch_b%length*(d - 1) + b 
                              bc = batch_b%length*(c - 1) + b
!
                              ab = (max(a,b)*(max(a,b)-3)/2) + a + b
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
!
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
!
                                omega2(biaj,1) = omega2(biaj, 1) &
                                                   + omega2_p_ab_ij(ab, ij) - omega2_m_ab_ij(ab, ij)
!
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
!
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
                        ac =  batch_a%length*(c - 1) + a
                        ad =  batch_a%length*(d - 1) + a
!
                        do  b = 1, batch_b%length
!
                             bd = batch_b%length*(d - 1) + b
                             bc = batch_b%length*(c - 1) + b
!
                             ab = index_two(a, b, batch_a%length)
!
                             g_p_ab_cd(ab, cd) = g_ac_bd(ac, bd) + g_ac_bd(ad, bc)
                             g_m_ab_cd(ab, cd) = g_ac_bd(ac, bd) - g_ac_bd(ad, bc)
!
                            if(c .ne. d) then
!
                              g_p_ab_cd(ab, cd) = two*g_p_ab_cd(ab, cd)
                              g_m_ab_cd(ab, cd) = two*g_m_ab_cd(ab, cd)
!
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
                        ai = wf%n_v*(i - 1) + a + batch_a%first - 1 ! A is full-space a index
                        aj = wf%n_v*(j - 1) + a + batch_a%first - 1 ! A is full-space a index
!
                        do b = 1, batch_b%length
!
                           if (a + batch_a%first - 1 .ge. b + batch_b%first - 1) then
!
                              bj = wf%n_v*(j - 1) + b + batch_b%first - 1 ! B is full-space b index
                              bi = wf%n_v*(i - 1) + b + batch_b%first - 1 ! B is full-space b index
!
                              ab = batch_a%length*(b - 1) + a
!
                              aibj = max(ai, bj)*(max(ai, bj)-3)/2 + ai + bj
                              biaj = max(bi, aj)*(max(bi, aj)-3)/2 + bi + aj
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
      call ccsd_a2_timer%freeze()
!
      call ccsd_a2_timer%switch_off()
      call ccsd_a2_integral_timer%switch_off()
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine omega_ccsd_b2_ccsd(wf, omega2)
!!
!!    Omega B2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    and Andreas Skeidsvoll, 2018
!!
!!    Omega B2 = sum_(kl) t_ak_bl*(g_kilj + sum_(cd) t_ci_dj * g_kc_ld)
!!
!!    Structure: g_kilj is constructed first and reordered as g_kl_ij.
!!    Then the contraction over cd is performed, and the results added to g_kl_ij.
!!    t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kc_ld
      real(dp), dimension(:,:), allocatable :: g_kl_cd
      real(dp), dimension(:,:), allocatable :: g_kl_ij
      real(dp), dimension(:,:), allocatable :: g_ki_lj
!
!     Reordered T2 apmlitudes
!
      real(dp), dimension(:,:), allocatable :: t_cd_ij
!
!     Reordered omega
!
      real(dp), dimension(:,:), allocatable :: omega_ab_ij
      real(dp), dimension(:,:), allocatable :: omega_ai_bj
!
      type(timings) :: ccsd_b2_timer
!
      call ccsd_b2_timer%init('omega ccsd b2')
      call ccsd_b2_timer%start()
!
!     Allocate and construct g_ki_lj
!
      call mem%alloc(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
!
      call wf%get_oooo(g_ki_lj)
!
      call mem%alloc(g_kl_ij, (wf%n_o)*(wf%n_o),(wf%n_o)*(wf%n_o))
!
      call sort_1234_to_1324(g_ki_lj, g_kl_ij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ki_lj, (wf%n_o)*(wf%n_o), (wf%n_o)*(wf%n_o))
!
!     Allocate and construct g_kc_ld
!
      call mem%alloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_ld)
!
!     Reorder g_kc_ld as g_kl_cd
!
      call mem%alloc(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
      call sort_1234_to_1324(g_kc_ld, g_kl_cd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kc_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_ci_dj as t_cd_ij
!
      call mem%alloc(t_cd_ij, (wf%n_v)**2, (wf%n_o)**2)
      call squareup_and_sort_1234_to_1324(wf%t2, t_cd_ij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_cd,     &
                  (wf%n_o)**2, &
                  t_cd_ij,     &
                  (wf%n_v)**2, &
                  one,         &
                  g_kl_ij,     &
                  (wf%n_o)**2)
!
!     Deallocate t_cd_ij and g_kl_cd
!
      call mem%dealloc(g_kl_cd, (wf%n_o)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     omega_ab_ij = sum_(kl) t_ab_kl*X_kl_ij
!
      call mem%alloc(omega_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_cd_ij,     & ! t_ab_kl
                  (wf%n_v)**2, &
                  g_kl_ij,     &
                  (wf%n_o)**2, &
                  zero,        &
                  omega_ab_ij, &
                  (wf%n_v)**2)
!
      call mem%dealloc(t_cd_ij, (wf%n_v)**2, (wf%n_o)**2)
      call mem%dealloc(g_kl_ij, (wf%n_o)**2, (wf%n_o)**2)
!
!     Reorder into omega2
!
      call mem%alloc(omega_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_1324(omega_ab_ij, omega_ai_bj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(omega_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      call add_to_packed(omega2, omega_ai_bj, (wf%n_o)*(wf%n_v))
      call mem%dealloc(omega_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call ccsd_b2_timer%freeze()
      call ccsd_b2_timer%switch_off()
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   module subroutine omega_ccsd_c2_ccsd(wf, omega2)
!!
!!    Omega C2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    and Andreas Skeidsvoll, 2018
!!
!!    Omega C2 = -1/2 * sum_(ck) t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!!                    - sum_(ck) t_bk_ci*(g_kj_ac - sum_(dl)t_al_dj * g_kd_lc)
!!                    - 1/2 * sum_ck u_jk^bc g_acki
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kd_lc
      real(dp), dimension(:,:), allocatable :: g_dl_ck
      real(dp), dimension(:,:), allocatable :: g_ki_ac
      real(dp), dimension(:,:), allocatable :: g_ai_ck
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_ai_dl
      real(dp), dimension(:,:), allocatable :: t_ck_bj
      real(dp), dimension(:,:), allocatable :: t_bk_cj
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:), allocatable :: X_ai_ck
      real(dp), dimension(:,:), allocatable :: Y_ai_bj
!
!     Reordered U2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: u_ck_bj
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! Holds term temporarily
      real(dp), dimension(:,:), allocatable :: omega_a_batch
!
!     Indices
!
      integer(i15) :: a, b, c
      integer(i15) :: i, j, k
!
      integer(i15) :: ai, aj, bi, bj, ac, ki, ck
      integer(i15) :: ai_batch
!
      integer(i15) :: aibj
!
!     Batching and memory handling
!
      integer(i15) :: rec0, rec1
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      type(timings) :: ccsd_c2_timer
!
      call ccsd_c2_timer%init('omega ccsd c2')
      call ccsd_c2_timer%start()
!
!     Sort t_al_di = t_li^ad as t_ai_dl (1234 to 1432)
!
      call mem%alloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call squareup_and_sort_1234_to_1432(wf%t2, t_ai_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Allocate and construct g_kd_lc
!
      call mem%alloc(g_kd_lc,(wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kd_lc)
!
!     Sort g_kd_lc to g_dl_ck (1234 to 2341)
!
      call mem%alloc(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call sort_1234_to_2341(g_kd_lc, g_dl_ck, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kd_lc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     -1/2*sum_(dl) t_ai_dl*g_dl_ck = X_ai_ck
!
      call mem%alloc(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -half,             &
                  t_ai_dl,           &
                  (wf%n_o)*(wf%n_v), &
                  g_dl_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate L_ia_J and g_dl_ck,
!
      call mem%dealloc(g_dl_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form u_ck_bj = u_jk^bc = u_kj^cb = 2 * t_jk^bc - t_kj^bc
!
      call mem%alloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_ck_bj, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      u_ck_bj = zero
!
      call add_1432_to_1234(-one, t_ck_bj, u_ck_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ck_bj, 1, u_ck_bj, 1)
!
      call mem%dealloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Allocate a holder for - 1/2 * sum_ck u_jk^bc g_acki,
!     constructed in batches over the a index below
!
      call mem%alloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      omega2_ai_bj = zero
!
!     Constructing g_ki_ac
!
!     Prepare for batching
!
      rec0 = wf%n_o**2*wf%integrals%n_J
!
      rec1 = wf%n_v*wf%integrals%n_J + (wf%n_o)*(wf%n_v**2)
!
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_a, rec0, rec1)
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
!        Allocate and construct g_ki_ac
!
         call mem%alloc(g_ki_ac, (wf%n_o)**2, (batch_a%length)*(wf%n_v))
!
         call wf%get_oovv(g_ki_ac,                       &
                           1, wf%n_o,                    &
                           1, wf%n_o,                    &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v)
!
!        X_ai_ck = X_ai_ck + g_ki_ac
!
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = wf%n_v*(k-1) + c
!
               do i = 1, wf%n_o
!
                  ki = wf%n_o*(i - 1) + k
!
                  do a = 1, batch_a%length
!
                     ac = batch_a%length*(c-1) + a
                     ai = wf%n_v*(i-1) + a + batch_a%first - 1 
!
                     X_ai_ck(ai, ck) = X_ai_ck(ai, ck) + g_ki_ac(ki, ac)
!
                  enddo
               enddo
            enddo
         enddo
!
!        Calculate the contribution to the term
!
!        omega_ai_bj = - 1/2 * sum_ck u_jk^bc g_acki
!
         call mem%alloc(g_ai_ck, (batch_a%length)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
         call sort_1234_to_3241(g_ki_ac, g_ai_ck, wf%n_o, wf%n_o, batch_a%length, wf%n_v)
!
         call mem%dealloc(g_ki_ac, (wf%n_o)**2, (batch_a%length)*(wf%n_v))
!
!        - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_ai_ck u_ck_bj
!
         call mem%alloc(omega_a_batch, batch_a%length*wf%n_o, wf%n_v*wf%n_o)
!
         call dgemm('N','N',                 &
                  (wf%n_o)*(batch_a%length), &
                  (wf%n_o)*(wf%n_v),         &
                  (wf%n_o)*(wf%n_v),         &
                  -one/two,                  &
                  g_ai_ck,                   &
                  (wf%n_o)*(batch_a%length), &
                  u_ck_bj,                   &
                  (wf%n_o)*(wf%n_v),         &
                  zero,                      &
                  omega_a_batch, &
                  (wf%n_o)*batch_a%length)
!
         call mem%dealloc(g_ai_ck, (batch_a%length)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
         do j = 1, wf%n_o
            do b = 1, wf%n_v
!
               bj = wf%n_v*(j-1) + b
!
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     ai_batch = batch_a%length*(i-1) + a
                     ai = wf%n_v*(i-1) + a + batch_a%first - 1
!
                     omega2_ai_bj(ai, bj) = omega_a_batch(ai_batch, bj)               
!
                  enddo
               enddo
            enddo
         enddo
!
         call mem%dealloc(omega_a_batch, batch_a%length*wf%n_o, wf%n_v*wf%n_o)
!
      enddo ! End of batching
!
!     Deallocate reordered u_ck_bj vector
!
      call mem%dealloc(u_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the - 1/2 * sum_ck u_jk^bc g_acki term to omega
!
!     Omega2(aibj, 1) =+ omega2_ai_bj(ai,bj) + omega2_bj_ai(bj,ai)
!
      call symmetric_sum(omega2_ai_bj, (wf%n_o)*(wf%n_v))         ! symmetrize
      call add_to_packed(omega2, omega2_ai_bj, (wf%n_o)*(wf%n_v)) ! add to packed
!
      call mem%dealloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Reorder t_bk_cj = t_kj^bc as t_ck_bj
!
      call mem%alloc(t_bk_cj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_bk_cj, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_3214(t_bk_cj, t_ck_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_bk_cj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form intermediate Y_ai_bj = - sum_(ck) X_ai_ck*t_ck_bj
!
      call mem%alloc(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  X_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  t_ck_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ai_bj,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(t_ck_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Omega_aibj,1 = P_ai_bj ( 1/2*Y_ai_bj + Y_aj_bi )
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = wf%n_v*(i - 1) + a
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aj = wf%n_v*(j - 1) + a
                     bi = wf%n_v*(i - 1) + b
!
                     aibj = max(ai, bj)*(max(ai, bj)-3)/2 + ai + bj
!
                     omega2(aibj, 1) = omega2(aibj, 1) + half*Y_ai_bj(ai, bj) + Y_ai_bj(aj, bi) &
                                                        + half*Y_ai_bj(bj, ai) + Y_ai_bj(bi, aj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(Y_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call ccsd_c2_timer%freeze()
      call ccsd_c2_timer%switch_off()
!
   end subroutine omega_ccsd_c2_ccsd
!
!
   module subroutine omega_ccsd_d2_ccsd(wf, omega2)
!!
!!    Omega D2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    and Andreas Skeidsvoll, 2018
!!
!!    Calculates the D2 term,
!!
!!      D2: sum_ck u_jk^bc g_aikc
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!    where
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!    The first and second terms are referred to as D2.1 and D2.2.
!!    All terms are added to the omega vector of the wavefunction object wf.
!!
!!    The routine adds the terms in the following order: D2.2, D2.1
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing D2.2 & D2.1
!
      real(dp), dimension(:,:), allocatable :: g_ld_kc ! g_ldkc
      real(dp), dimension(:,:), allocatable :: L_ld_kc ! L_ldkc = 2 * g_ldkc - g_lckd
      real(dp), dimension(:,:), allocatable :: t_ai_dl ! t_il^ad
      real(dp), dimension(:,:), allocatable :: u_ai_dl ! u_il^ad = 2 * t_il^ad - t_li^ad
      real(dp), dimension(:,:), allocatable :: u_ai_ld ! u_il^ad = 2 * t_il^ad - t_li^ad
      real(dp), dimension(:,:), allocatable :: Z_ai_kc ! An intermediate, see below
      real(dp), dimension(:,:), allocatable :: g_ai_kc ! g_aikc
!
      type(timings) :: ccsd_d2_timer
!
      call ccsd_d2_timer%init('omega ccsd d2')
      call ccsd_d2_timer%start()
!
!     :: Calculate the D2.2 term of omega ::
!
!     Form L_ld_kc = L_ldkc = 2*g_ld_kc(ld,kc) - g_ld_kc(lc,kd)
!
      call mem%alloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%get_ovov(g_ld_kc)
!
      call mem%alloc(L_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ld_kc = zero
!
      call add_1432_to_1234(-one, g_ld_kc, L_ld_kc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, g_ld_kc, 1, L_ld_kc, 1)
!
      call mem%dealloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form u_ai_ld = u_il^ad = 2 * t_il^ad - t_li^ad
!
      call mem%alloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_ai_dl, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      u_ai_dl = zero
!
      call add_1432_to_1234(-one, t_ai_dl, u_ai_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ai_dl, 1, u_ai_dl, 1)
!
      call mem%dealloc(t_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_1243(u_ai_dl, u_ai_ld, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(u_ai_dl, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the intermediate Z_ai_kc = sum_dl u_ai_ld L_ld_kc and set it to zero
!
      call mem%alloc(Z_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  L_ld_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Z_ai_kc,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form the D2.2 term, 1/4 sum_kc Z_ai_kc u_kc_bj = 1/4 sum_kc Z_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call mem%alloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one/four,          &
                  Z_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
!     Some justification for the above matrix multiplication. We have
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_ai_ld(ai,ld) = u_il^ad,
!     which means that u_ai_ld(bj,kc)^T = u_ai_ld(kc,bj) = u_kj^cb = u_jk^bc.
!
      call mem%dealloc(Z_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     :: Calculate the D2.1 term of omega ::
!
!     Form g_ai_kc = g_aikc
!
      call mem%alloc(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_voov(g_ai_kc)
!
!     Calculate the D2.1 term, sum_ck u_jk^bc g_aikc = sum_ck g_ai_kc(ai,kc) u_ai_ld(bj,kc)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_ai_kc,           &
                  (wf%n_o)*(wf%n_v), &
                  u_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega2_ai_bj,      &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(g_ai_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add the D2.1 term to the omega vector
!
      call symmetrize_and_add_to_packed(omega2, omega2_ai_bj, (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call ccsd_d2_timer%freeze()
      call ccsd_d2_timer%switch_off()
!
   end subroutine omega_ccsd_d2_ccsd
!
!
   module subroutine omega_ccsd_e2_ccsd(wf, omega2)
!
!     Omega E2
!     Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!     and Andreas Skeidsvoll, 2018
!
!     Calculates the E2 term,
!
!      E2: sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd)
!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!
!     where
!
!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!
!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!
!     Both are permuted added to the projection vector element omega2(ai,bj) of
!     the wavefunction object wf.
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
!     Vectors for E2.1 term
!
      real(dp), dimension(:,:), allocatable :: omega2_bj_ai ! For storing the E2.1 term temporarily
      real(dp), dimension(:,:), allocatable :: g_ld_kc      ! g_ldkc
      real(dp), dimension(:,:), allocatable :: u_bl_dk      ! u_kl^bd
      real(dp), dimension(:,:), allocatable :: t_bk_dl      ! t_kl^bd
      real(dp), dimension(:,:), allocatable :: u_bk_dl      ! u_kl^bd
      real(dp), dimension(:,:), allocatable :: X_b_c        ! An intermediate, see below for definition
      real(dp), dimension(:,:), allocatable :: t_cj_ai      ! t_ij^ac
!
!     Vectors for E2.2 term
!
      real(dp), dimension(:,:), allocatable :: omega2_ai_bj ! For storing the E2.2 term temporarily
      real(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition
!
      type(timings) :: ccsd_e2_timer 
!
      call ccsd_e2_timer%init('omega ccsd e2')
      call ccsd_e2_timer%start()
!
!     :: Calculate the E2.1 term of omega ::
!
!     Form u_bl_dk = u_kl^bd = u_bk_dl
!                  = 2 * t_kl^bd - t_lk^bd = 2 * t_bk_dl(bk,dl) - t_bk_dl(bl, dk)
!
      call mem%alloc(t_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call squareup(wf%t2, t_bk_dl, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      u_bk_dl = zero
!
      call add_1432_to_1234(-one, t_bk_dl, u_bk_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_bk_dl, 1, u_bk_dl, 1)
!
      call mem%dealloc(t_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_bl_dk, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call sort_1234_to_1432(u_bk_dl, u_bl_dk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(u_bk_dl, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form g_ld_kc = g_ldkc
!
      call mem%alloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%get_ovov(g_ld_kc)
!
!     Make the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call mem%alloc(X_b_c, wf%n_v, wf%n_v)
!
!     First, copy the virtual-virtual Fock matrix into the intermediate
!
      call dcopy((wf%n_v)**2, wf%fock_ab, 1, X_b_c, 1) ! X_b_c = F_bc
!
!     Then, add the second contribution,
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  u_bl_dk,              & ! u_b_ldk
                  wf%n_v,               &
                  g_ld_kc,              & ! g_ldk_c
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_b_c,                &
                  wf%n_v)
!
!     Form t_cj_ai = t_ij^ac = t_ji^ca
!
      call mem%alloc(t_cj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call squareup(wf%t2, t_cj_ai, (wf%n_v)*(wf%n_o))
!
!     Form the E2.1 term
!
      call mem%alloc(omega2_bj_ai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  X_b_c,                &
                  wf%n_v,               &
                  t_cj_ai,              & ! t_c_jai
                  wf%n_v,               &
                  zero,                 &
                  omega2_bj_ai,         & ! omega2_b_jai
                  wf%n_v)
!
!     Add the E2.1 term to the omega vector and deallocations
!
      call symmetrize_and_add_to_packed(omega2, omega2_bj_ai, (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(omega2_bj_ai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
      call mem%dealloc(X_b_c, wf%n_v, wf%n_v)
!
!     :: Calculate E2.2 term of omega ::
!
!     Make the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc
!                                 = F_kj  + sum_cdl u_lj^dc g_kcld
!                                 = F_k_j + sum_cdl g_k_cld * u_cld_j
!
      call mem%alloc(Y_k_j, wf%n_o, wf%n_o)
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj
!
      call dcopy((wf%n_o)**2, wf%fock_ij, 1, Y_k_j, 1)
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that
!     Y_k_j = F_k_j + sum_cdl g_k_cld u_cld_j
!
!     Note:
!
!     g_ld_kc(kc,ld) = g_kcld              -> pretend that this is g_k_cld
!     u_b_ldk(c,ldj) = u_jl^cd (= u_lj^dc) -> pretend that this is u_cld_j
!
      call dgemm('N','N',              &
                 wf%n_o,               &
                 wf%n_o,               &
                 (wf%n_o)*(wf%n_v)**2, &
                 one,                  &
                 g_ld_kc,              & ! g_k_cld
                 wf%n_o,               &
                 u_bl_dk,              & ! u_cld_j
                 (wf%n_o)*(wf%n_v)**2, &
                 one,                  &
                 Y_k_j,                &
                 wf%n_o)
!
      call mem%dealloc(u_bl_dk, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
      call mem%dealloc(g_ld_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate the E2.2 term,
!     - sum_k t_aib_k Y_k_j = - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
!
!     Note: t_cj_ai = t_ji^ca => t_cj_ai(ai,bk) = t_ik^ab;
!     thus, we can treat t_cj_ai as t_aib_k = t_ik^ab.
!
      call mem%alloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one,                 &
                 t_cj_ai,              & ! t_aib_k
                 (wf%n_o)*(wf%n_v)**2, &
                 Y_k_j,                &
                 wf%n_o,               &
                 zero,                 &
                 omega2_ai_bj,         & ! omega2_aib_j
                 (wf%n_o)*(wf%n_v)**2)
!
!     Deallocate Y_k_j and the amplitudes
!
      call mem%dealloc(Y_k_j, wf%n_o, wf%n_o)
      call mem%dealloc(t_cj_ai, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Add the E2.2 term to the omega vector
!
      call symmetrize_and_add_to_packed(omega2, omega2_ai_bj, (wf%n_o)*(wf%n_v))
!
!     Deallocate the E2.2 term
!
      call mem%dealloc(omega2_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call ccsd_e2_timer%freeze()
      call ccsd_e2_timer%switch_off()
!
   end subroutine omega_ccsd_e2_ccsd
!
!
end submodule omega_ccsd 
