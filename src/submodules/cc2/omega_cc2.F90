submodule (cc2_class) omega_cc2
!
!!
!!    Omega submodule (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
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
!!    and Linda Goletto, 2018
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
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: g_bi_cj
      real(dp), dimension(:,:), allocatable :: l_bj_ci
      real(dp), dimension(:,:), allocatable :: g_ab_jc
!
      integer(i15) :: b, j, c, i
      integer(i15) :: bj, ci, bi, cj
!
      call mem%alloc(g_bi_cj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%alloc(l_bj_ci, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%get_vovo(g_bi_cj)
!
      do b = 1, wf%n_v
         do  j = 1, wf%n_o
             do c = 1, wf%n_v
                do i = 1, wf%n_o
                   bj = index_two(b,j,wf%n_v)
                   ci = index_two(c,i,wf%n_v)
                   bi = index_two(b,i,wf%n_v)
                   cj = index_two(c,j,wf%n_v)
!                  
                   l_bj_ci= - (two*g_bi_cj(bi,cj)/(wf%fock_ab(b,b) + wf%fock_ab(c,c) - &
                           wf%fock_ij(i,i) - wf%fock_ij(j,j))) + g_bi_cj(bj,ci)/(wf%fock_ab(b,b) &
                           + wf%fock_ab(c,c) - wf%fock_ij(i,i) - wf%fock_ij(j,j))
                enddo
             enddo
         enddo
      enddo
!
      call mem%dealloc(g_bi_cj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call mem%alloc(g_ab_jc, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_vvov(g_ab_jc)
!
              call dgemm('N','N',                            &
                          wf%n_v,                            &
                          wf%n_o,                            &
                          wf%n_v**2*wf%n_o,                  &
                          one,                               &
                          g_ab_jc,                           &
                          wf%n_v,                            &
                          l_bj_ci,                           &
                          wf%n_v**2*wf%n_o,                  &
                          one,                               &
                          omega,                             &
                          wf%n_v)
!
      call mem%dealloc(g_ab_jc, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(l_bj_ci, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, omega)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
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
                  aj = index_two(a,j,wf%n_v)
                  bk = index_two(b,k,wf%n_v)
                  g_aj_bk(aj,bk) = g_aj_bk(aj,bk)/(- wf%fock_ab(a,a) &
                     - wf%fock_ab(b,b)  + wf%fock_ij(j,j) + wf%fock_ij(k,k))
               enddo
            enddo
         enddo
      enddo
!
      call mem%alloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call wf%get_ovoo(g_jb_ki)
!
              call dgemm('N','N',                            &
                          wf%n_v,                            &
                          wf%n_o,                            &
                          wf%n_o**2*wf%n_v,                  &
                          -one,                              &
                          g_aj_bk,                           &
                          wf%n_v,                            &
                          g_jb_ki,                           &
                          wf%n_o**2*wf%n_v,                  &
                          one,                               &
                          omega,                             &
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
      call sort_1234_to_3214(g_kb_ji,g_jb_ki,wf%n_o,wf%n_v,wf%n_o,wf%n_o)
!
      call mem%dealloc(g_kb_ji, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
              call dgemm('N','N',                            &
                          wf%n_v,                            &
                          wf%n_o,                            &
                          wf%n_o**2*wf%n_v,                  &
                          two,                               &
                          g_aj_bk,                           &
                          wf%n_v,                            &
                          g_jb_ki,                           &
                          wf%n_o**2*wf%n_v,                  &
                          one,                               &
                          omega,                             &
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
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine omega_cc2_c1_cc2
!
end submodule
