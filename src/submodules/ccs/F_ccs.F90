submodule (ccs_class) F_ccs
!
!!
!!    F-transformation submodule (ccs)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    Routines for the linear transform of
!!    vectors by the F matrix 
!!
!!    ρ = F * c,
!!
!!    where
!!   
!!    F_μ,ν = < Λ' | [[ exp(T) H_0 exp(T), τ_μ ], τ_ν ] | R >,
!!
!!    Where < Λ' | = < R | + sum_μ tbar_μ < μ |
!!  
! 
   implicit none
!
!
contains
!
!
   module subroutine F_transform_vector_ccs(wf, c)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018     
!!
!!    Directs the transformation by the F matrix.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: rho
!
      integer(i15) :: a, i, ai
!
      call mem%alloc(rho, wf%n_v, wf%n_o)
      rho = zero
!
      call wf%F_ccs_a1_0(c, rho)
!  
      call wf%F_ccs_a1_1(c, rho)
      call wf%F_ccs_b1_1(c, rho)
      call wf%F_ccs_c1_1(c, rho)
!
!$omp parallel do private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai, 1) = rho(a, i)
!
         enddo
      enddo
!$omp end parallel do      
!
      call mem%dealloc(rho, wf%n_v, wf%n_o)
!
   end subroutine F_transform_vector_ccs
!
!
   module subroutine F_ccs_a1_0_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,0 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,0_ai = 2 L_iajb c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_iajb)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     L_iajb = 2 g_iajb - g_ibja (ordered as L_aibj)
!
      L_aibj = zero
!
      call add_2143_to_1234(two, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     rho_ai += 2 L_iajb c_bj
!
      call dgemm('N', 'N', &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  two,                 &
                  L_aibj,              &
                  (wf%n_v)*(wf%n_o),   &
                  c_ai,                & ! c_bj
                  (wf%n_v)*(wf%n_o),   & 
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccs_a1_0_ccs
!
!
   module subroutine F_ccs_a1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1_ai = - (F_ib tbar_aj + F_ja tbar_bi) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable :: X_ij
      real(dp), dimension(:,:), allocatable :: X_ji
!
!     Term 1 : - F_ib tbar_aj c_bj
!
      call mem%alloc(X_ij, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',    &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ia, &
                  wf%n_o,     &
                  c_ai,       & ! c_bj
                  wf%n_v,     &
                  zero,       &
                  X_ij,       &
                  wf%n_o)
!
      call dgemm('N', 'T',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%t1bar,   & !tbar_aj
                  wf%n_v,     &
                  X_ji,       &
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)
!
      call mem%dealloc(X_ij, wf%n_o, wf%n_o)
!
!     Term 2 : - F_ja tbar_bi c_bj
!
      call mem%alloc(X_ji, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',    &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  c_ai,       & ! c_bj
                  wf%n_v,     &
                  wf%t1bar,   & ! tbar_bi
                  wf%n_v,     &
                  zero,       &
                  X_ji,       &
                  wf%n_o)
!
      call dgemm('T', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &  
                  wf%n_o,     &
                  -one,       &
                  wf%fock_ia, & ! F_ja
                  wf%n_o,     &
                  X_ji,       &
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)
!
      call mem%dealloc(X_ji, wf%n_o, wf%n_o)
!
   end subroutine F_ccs_a1_1_ccs
!
!
   module subroutine F_ccs_b1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation B1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_B1,1_ai = - (L_ikjb tbar_ak + L_jkia tbar_bk) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjb
      real(dp), dimension(:,:,:,:), allocatable :: L_kibj
!
      real(dp), dimension(:,:), allocatable :: X_ki
      real(dp), dimension(:,:), allocatable :: X_kj
!
!     Term 1: - L_ikjb tbar_ak c_bj
!
!     L_ikjb = 2 g_ikjb - g_jkib (ordered as L_kibj)
!
      call mem%alloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_ooov(g_ikjb)   
!
      call mem%alloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)  
!
      L_kibj = zero
      call add_2143_to_1234(two, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     X_ki = L_kibj c_bj
!
      call mem%alloc(X_ki, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  L_kibj,              & ! L_ki_jb
                  wf%n_o**2,           &
                  c_ai,                & ! c_bj
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_ki,                &
                  wf%n_o**2)
!
!      rho_ai += tbar_ak X_ki
!
      call dgemm('N', 'N',    &  
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%t1bar,   & ! tbar_ak
                  wf%n_v,     &
                  X_ki,       &
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v) 
!
!     Term 2: - L_jkia tbar_bk c_bj
!
!     NOTE: We will now pretend that L_kibj(k, i, b, j) = L_ikjb is
!           L_kjai(k, j, a, i) = L_jkia in order to not make integrals twice
!
!     X_kj = tbar_bk c_bj
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',    &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%t1bar,   & ! tbar_bk
                  wf%n_v,     &
                  c_ai,       & ! c_bj
                  wf%n_v,     &
                  zero,       &
                  X_kj,       &
                  wf%n_o)
!
      call dgemm('T', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_o)**2,         &
                  -one,                &
                  L_kibj,              & ! L_kj_ai
                  wf%n_o**2,           &
                  X_kj,                &
                  wf%n_o**2,           &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o) 
!
   end subroutine F_ccs_b1_1_ccs
!
!
   module subroutine F_ccs_c1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation C1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_C1,1_ai = (L_cajb tbar_ci + L_cbia tbar_cj) c_bj
!!                = (X_iajb + X_jbia) c_bj
!!
!!    In batches of c
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: X_iajb
      real(dp), dimension(:,:,:,:), allocatable :: X_aibj
      real(dp), dimension(:,:,:,:), allocatable :: g_cajb
      real(dp), dimension(:,:,:,:), allocatable :: L_cajb
!
      type(batching_index) :: batch_c
!
      integer(i15) :: current_c_batch, req0, req1
!
      call mem%alloc(X_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      X_iajb = zero
!
!     Prepare for batching over index c
!
      req0 = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
!
      req1 = max(2*(wf%n_o)*(wf%n_v**2), (wf%n_o)*(wf%n_v**2) + (wf%integrals%n_J)*(wf%n_v))
!
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1)
!
      do current_c_batch = 1, batch_c%num_batches
!
!        Determine the limits for the current c-batch
!
         call batch_c%determine_limits(current_c_batch)
!
!        L_cajb = 2 g_cajb - g_cbja
!
         call mem%alloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_cajb,                        &
                           batch_c%first, batch_c%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1,wf%n_v)
!
         call mem%alloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dcopy((batch_c%length)*(wf%n_v**2)*(wf%n_o), g_cajb, 1, L_cajb, 1)
         call dscal((batch_c%length)*(wf%n_v**2)*(wf%n_o), two, L_cajb, 1)
!
         call add_1432_to_1234(-one, g_cajb, L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_cajb, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('T', 'N',                      &
                     wf%n_o,                       &
                     (wf%n_v**2)*wf%n_o,           &
                     batch_c%length,               &
                     one,                          &
                     wf%t1bar(batch_c%first, 1),   & ! tbar_c_i
                     wf%n_v,                       &
                     L_cajb,                       & ! L_c_ajb
                     batch_c%length,               &
                     one,                          &
                     X_iajb,                       & ! X_i_ajb
                     wf%n_o)
!
         call mem%dealloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! batches of c
!
      call mem%alloc(X_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2143(X_iajb, X_aibj, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     rho_ai += X_aibj c_bj
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  X_aibj,              & ! X_ai_bj
                  (wf%n_v)*(wf%n_o),   &
                  c_ai,                & ! c_bj
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
!     rho_ai += X_bjai c_bj (Using bj <-> ai in X_aibj)
!
      call dgemm('T', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  X_aibj,              & ! X_bj_ai
                  (wf%n_v)*(wf%n_o),   &
                  c_ai,                & ! c_bj
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccs_c1_1_ccs
!
end submodule F_ccs
