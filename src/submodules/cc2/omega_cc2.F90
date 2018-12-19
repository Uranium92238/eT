submodule (cc2_class) omega_cc2
!
!!
!!    Omega submodule (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto and Alexander Paul, 2018
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
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      omega = zero
!
      call wf%omega_ccs_a1(omega)
!
      call wf%omega_cc2_a1(omega)
      call wf%omega_cc2_b1(omega)
      call wf%omega_cc2_c1(omega)
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf, omega)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: g_bi_cj
      real(dp), dimension(:,:), allocatable :: L_bj_ci
      real(dp), dimension(:,:), allocatable :: g_ab_jc
!
      integer(i15) :: b, j, c, i
      integer(i15) :: bj, ci, bi, cj
!
      integer(i15) :: rec0, rec1_b, rec1_c, rec2
!
      integer(i15) :: current_b_batch, current_c_batch
!
      type(batching_index) :: batch_b, batch_c
!
      rec0 = 0
!
      rec1_b = (wf%integrals%n_J)*(wf%n_v**2)
      rec1_c = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
!
      rec2 =  (wf%n_o**2)*(wf%n_v**2) + (wf%n_o)*(wf%n_v**3)
!
      call batch_b%init(wf%n_v)
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_b, batch_c, rec0, rec1_b, rec1_c, rec2)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
            call mem%alloc(g_bi_cj, (batch_b%length)*(wf%n_o), (batch_c%length)*(wf%n_o))
            call mem%alloc(L_bj_ci, (batch_b%length)*(wf%n_o), (batch_c%length)*(wf%n_o))
!
            call wf%get_vovo(g_bi_cj,                      &
                              batch_b%first, batch_b%last, &
                              1, wf%n_o,                   &
                              batch_c%first, batch_c%last, &
                              1, wf%n_o)
!
            do b = 1, (batch_b%length)
               do  j = 1, wf%n_o
                   do c = 1, (batch_c%length)
                      do i = 1, wf%n_o
!
                         bj = batch_b%length*(j-1) + b
                         ci = batch_c%length*(i-1) + c
                         bi = batch_b%length*(i-1) + b
                         cj = batch_c%length*(j-1) + c
!                        
                         L_bj_ci(bj,ci) = - (two*g_bi_cj(bi,cj)/( wf%fock_diagonal(b + batch_b%first - 1 + wf%n_o, 1) &
                                                                + wf%fock_diagonal(c + batch_c%first - 1 + wf%n_o, 1) &
                                                                - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))) &
                                               + g_bi_cj(bj,ci)/( wf%fock_diagonal(b + batch_b%first - 1 + wf%n_o, 1) &
                                                                + wf%fock_diagonal(c + batch_c%first - 1 + wf%n_o, 1)  &
                                                                - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))
                      enddo
                   enddo
               enddo
            enddo
!
            call mem%dealloc(g_bi_cj, (batch_b%length)*(wf%n_o), (batch_c%length)*(wf%n_o))
!
            call mem%alloc(g_ab_jc, (batch_b%length)*(wf%n_v), (batch_c%length)*(wf%n_o))
!
            call wf%get_vvov(g_ab_jc,                       &
                              1, wf%n_v,                    &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_o,                    &
                              batch_c%first, batch_c%last)
!
            call dgemm('N','N',                                   &
                        wf%n_v,                                   &
                        wf%n_o,                                   &
                        (batch_b%length)*(batch_c%length)*wf%n_o, &
                        one,                                      &
                        g_ab_jc,                                  &
                        wf%n_v,                                   &
                        L_bj_ci,                                  &
                        (batch_b%length)*(batch_c%length)*wf%n_o, &
                        one,                                      &
                        omega,                                    &
                        wf%n_v)
!
            call mem%dealloc(g_ab_jc, (batch_b%length)*(wf%n_v), (batch_c%length)*(wf%n_o))
            call mem%dealloc(L_bj_ci, (batch_b%length)*(wf%n_o), (batch_c%length)*(wf%n_o))
!
         enddo
      enddo
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, omega)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: g_aj_bk
      real(dp), dimension(:,:), allocatable :: g_jb_ki
      real(dp), dimension(:,:), allocatable :: g_kb_ji
!
      integer(i15) :: a, b, j, k
      integer(i15) :: aj, bk
!
      call mem%alloc(g_aj_bk, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%get_vovo(g_aj_bk)
!
      do a = 1, wf%n_v
         do j = 1, wf%n_o
            do b = 1, wf%n_v
               do k = 1, wf%n_o
!
                  aj = wf%n_v*(j-1) + a
                  bk = wf%n_v*(k-1) + b
!
                  g_aj_bk(aj,bk) = - g_aj_bk(aj,bk)/(wf%fock_diagonal(a + wf%n_o, 1) + wf%fock_diagonal(b + wf%n_o, 1)&
                                                   - wf%fock_diagonal(j, 1) - wf%fock_diagonal(k, 1))
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%alloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call wf%get_ovoo(g_jb_ki)
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_o**2*wf%n_v, &
                  one,              &
                  g_aj_bk,          &
                  wf%n_v,           &
                  g_jb_ki,          &
                  wf%n_o**2*wf%n_v, &
                  one,              &
                  omega,            &
                  wf%n_v)
!
      call mem%dealloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call mem%alloc(g_kb_ji, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call wf%get_ovoo(g_kb_ji)
!
      call mem%alloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call sort_1234_to_3214(g_kb_ji, g_jb_ki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kb_ji, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call dgemm('N','N',            &
                  wf%n_v,            &
                  wf%n_o,            &
                  (wf%n_o**2)*wf%n_v,&
                  -two,              &
                  g_aj_bk,           &
                  wf%n_v,            &
                  g_jb_ki,           &
                  (wf%n_o**2)*wf%n_v,&
                  one,               &
                  omega,             &
                  wf%n_v)
!
      call mem%dealloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%dealloc(g_aj_bk, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
   end subroutine omega_cc2_b1_cc2
!
!
   module subroutine omega_cc2_c1_cc2(wf, omega)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: g_aibj
      real(dp), dimension(:,:), allocatable :: L_aijb
!
      integer(i15) :: i, j, a, b
      integer(i15) :: ai, aj, bi, bj, jb
!
      call mem%alloc(g_aibj, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
      call wf%get_vovo(g_aibj)
!
      call mem%alloc(L_aijb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
!$omp parallel do schedule(static) private(i, j, a, b, ai, aj, bi, bj, jb)
      do i = 1, wf%n_o 
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = wf%n_v*(i-1) + a
               aj = wf%n_v*(j-1) + a
!
               do b = 1, wf%n_v
!
                  bi = wf%n_v*(i-1) + b
                  bj = wf%n_v*(j-1) + b
                  jb = wf%n_o*(b-1) + j
!
                  L_aijb(ai, jb) = -two*(g_aibj(ai, bj)/(wf%fock_diagonal(a + wf%n_o, 1) + wf%fock_diagonal(b + wf%n_o, 1) &
                                                      - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))) &
                                    + (g_aibj(aj, bi)/(wf%fock_diagonal(a + wf%n_o, 1) + wf%fock_diagonal(b + wf%n_o, 1) &
                                                      - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))) 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_aibj, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
      call dgemm('N', 'N', &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  one,           &
                  L_aijb,        &
                  wf%n_v*wf%n_o, &
                  wf%fock_ia,    & ! F_jb_1
                  wf%n_v*wf%n_o, &
                  one,           &
                  omega,         &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(L_aijb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
   end subroutine omega_cc2_c1_cc2
!
end submodule
