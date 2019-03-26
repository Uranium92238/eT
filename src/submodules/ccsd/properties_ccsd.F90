submodule (ccsd_class) properties_ccsd
!
!!
!!    Properties submodule (CCSD)
!!    Written by Josefine H. Andersen, 2019
!!
!!    <insert description of actions>
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_etaX_ccsd(wf, Xoperator, etaX)
!!
!!    Construct etaX
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
!     etaX_ai:
!
      call wf%construct_etaX_ccs_singles(Xoperator, etaX)
      call wf%construct_etaX_singles_q1(Xoperator, etaX)
      call wf%construct_etaX_singles_q2(Xoperator, etaX)
!
!     etaX_aibj:
!
      call wf%construct_etaX_doubles_q1(Xoperator, etaX)
      call wf%construct_etaX_doubles_q2(Xoperator, etaX)
!
   end subroutine construct_etaX_ccsd
!
!
   module subroutine construct_etaX_ccs_singles_ccsd(wf, Xoperator, etaX)
!!
!!    Construct CCS term of CCSD etaX singles
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    2*X_ia
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!
      real(dp), dimension(wf%n_t1, 1), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: X_ia
!
      call mem%alloc(X_ia, wf%n_t1, 1)
      call  wf%get_operator_ov(Xoperator, X_ia)
!
      call dscal(wf%n_t1, two, X_ia, 1)
      
      call sort_12_to_21(X_ia, etaX, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_ia, wf%n_t1, 1)
!
   end subroutine construct_etaX_ccs_singles_ccsd
!
!
   module subroutine construct_etaX_singles_q1_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q1 term of etaX singles
!!    Written by Josefine H. Andersen
!!
!!    Q1 = sum_c tb_ci X_ca - sum_k tb_ak X_ik
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout)   :: etaX
!      
      real(dp), dimension(:,:), allocatable :: etaX_temp
!
      real(dp), dimension(:,:), allocatable :: X_ca
      real(dp), dimension(:,:), allocatable :: X_ik
!      
      call mem%alloc(etaX_temp, wf%n_v, wf%n_o)
!
!     :: First term  sum_c tb_ci X_ca
!
      call mem%alloc(X_ca, wf%n_v*wf%n_v, 1)
      call wf%get_operator_vv(Xoperator, X_ca)
!
      call dgemm('T','N',   &
                 wf%n_v,    &
                 wf%n_o,    &
                 wf%n_v,    &
                 one,       &
                 X_ca,      &
                 wf%n_v,    &
                 wf%t1bar,  &
                 wf%n_v,    &
                 zero,      &
                 etaX_temp, &
                 wf%n_v)
!         
      call mem%dealloc(X_ca, wf%n_v*wf%n_v, 1)
!      
!     :: Second term  - sum_k tb_ak X_ik
!
      call mem%alloc(X_ik, wf%n_o*wf%n_o, 1)
      call wf%get_operator_oo(Xoperator, X_ik)
!      
      call dgemm('N','T',   &
                 wf%n_v,    &
                 wf%n_o,    &
                 wf%n_o,    &
                 -one,      &
                 wf%t1bar,  &
                 wf%n_v,    &
                 X_ik,      &
                 wf%n_o,    &
                 one,       &
                 etaX_temp, &
                 wf%n_v)
!         
      call mem%dealloc(X_ik, wf%n_o*wf%n_o, 1)
!
!     Add eta_temp to etaX
!
      call daxpy(wf%n_t1, -one, etaX_temp, 1, etaX, 1) !OBS negative sign due to
                                                       !unsolved bug
!
      call mem%dealloc(etaX_temp, wf%n_v, wf%n_o)
!
   end subroutine construct_etaX_singles_q1_ccsd
!
!
   module subroutine construct_etaX_singles_q2_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q2 term of etaX singles
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Q2 = - sum_ckdl (tb_ckal X_id t_ckdl + tb_ckdi X_la t_ckdl)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: etaX_temp
!
      real(dp), dimension(:,:), allocatable :: X_id     ! X_la
!
      real(dp), dimension(:,:), allocatable :: tb_ck_al ! tb_ckdi
      real(dp), dimension(:,:), allocatable :: t_lck_d
!
      real(dp), dimension(:,:), allocatable :: I_a_d    ! intermediate, first term
      real(dp), dimension(:,:), allocatable :: I_l_i    ! intermediate, second term
!      
      real(dp), parameter :: one = 1.0
!
      call mem%alloc(etaX_temp, wf%n_v, wf%n_o)
!
      call mem%alloc(X_id, wf%n_o, wf%n_v)
      call wf%get_operator_ov(Xoperator, X_id)
!
!     Read multipliers and squareup
!
      call mem%alloc(tb_ck_al, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      call squareup(wf%t2bar, tb_ck_al, (wf%n_v)*(wf%n_o))
!      
!     Read amplitudes and order as t_lck_d = t_kl^cd
!
      call mem%alloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      call squareup_and_sort_1234_to_2341(wf%t2, t_lck_d, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!      
!     :: First term: - sum_ckdl tb_ckal X_id t_ckdl
!
!     Form the intermediate I_a_d = sum_ckl tb_a_lck t_lck_d = sum_ckl tb_ckal t_kl^cd
!  
      call mem%alloc(I_a_d, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tb_ck_al,             & 
                  wf%n_v,               &
                  t_lck_d,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  I_a_d,                &
                  wf%n_v)
!
!     Add   - sum_ckdl tb_ckal X_id t_kl^cd
!         = - sum_d I_a_d X_id 
!         = - sum_d I_a_d X_i_a^T(d,i)
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  I_a_d,      &
                  wf%n_v,     &
                  X_id,       &
                  wf%n_o,     &
                  zero,       &
                  etaX_temp,   &
                  wf%n_v)          
!
      call mem%dealloc(I_a_d, wf%n_v, wf%n_v)
!
!     :: Second term: sum_ckdl tb_ckdi X_la t_ckdl
!
!     Form the intermediate X_l_i = sum_ckd t_l_ckd b_ckd_i  = sum_ckd b_ckdi t_kl^cd
!
      call mem%alloc(I_l_i, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_lck_d,              & 
                  (wf%n_o),             &
                  tb_ck_al,             & 
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  I_l_i,                &
                  wf%n_o)
!
      call mem%dealloc(t_lck_d, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      call mem%dealloc(tb_ck_al, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     Add - sum_ckdl b_ckdi X_la t_kl^cd = - sum_l X_la I_l_i = - sum_l X_i_a^T(a,l) I_l_i(l,i)
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  X_id,       &
                  wf%n_o,     &
                  I_l_i,      &
                  wf%n_o,     &
                  one,        &
                  etaX_temp,  &
                  wf%n_v)
!
      call mem%dealloc(I_l_i, wf%n_o, wf%n_o)
      call mem%dealloc(X_id, wf%n_o, wf%n_v)
!
!     Add eta_temp to etaX
!
      call daxpy(wf%n_t1, -one, etaX_temp, 1, etaX, 1) !OBS: negative sign due to
                                                       !unsolved bug 
!
      call mem%dealloc(etaX_temp, wf%n_v, wf%n_o)
!
   end subroutine construct_etaX_singles_q2_ccsd
!
!
   module subroutine construct_etaX_doubles_q1_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q1 term of etaX doubles
!!    Written by Josefine H. Andersen, Feb 2019
!!    
!!    Q1 = 2 X_jb tb_ai - X_ib tb_aj
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: etaX_ai_bj
      real(dp), dimension(:,:), allocatable :: etaX_temp
!      
      real(dp), dimension(:,:), allocatable :: tb_ai
      real(dp), dimension(:,:), allocatable :: X_ai
!
      integer :: i, a, j, b, ai, bj
!
      call mem%alloc(X_ai, wf%n_o, wf%n_v)
      call wf%get_operator_ov(Xoperator, X_ai)
!      
      call mem%alloc(tb_ai, wf%n_v, wf%n_o)    
      call dcopy(wf%n_t1, wf%t1bar, 1, tb_ai, 1) 
!
      call mem%alloc(etaX_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      etaX_ai_bj = zero
!      
!     :: First and second term
!
       do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!  
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  etaX_ai_bj(ai, bj) = etaX_ai_bj(ai, bj) - (X_ai(i, b))*tb_ai(a, j) &
                                    + two*(X_ai(j, b))*tb_ai(a, i)
!
               enddo
            enddo
!
         enddo
      enddo
!
      call mem%dealloc(tb_ai, wf%n_v, wf%n_o) 
      call mem%dealloc(X_ai, wf%n_o, wf%n_v)
!
      call mem%alloc(etaX_temp, wf%n_t2, 1)
      etaX_temp = zero ! necessary
!
      call symmetrize_and_add_to_packed(etaX_temp, etaX_ai_bj, (wf%n_v*wf%n_o))
!
      call mem%dealloc(etaX_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     Add temporary etaX to etaX
!
      call daxpy(wf%n_t2, -one, etaX_temp, 1, etaX(wf%n_t1+1:wf%n_es_amplitudes, 1), 1)
!
      call mem%dealloc(etaX_temp, wf%n_t2, 1)
!
   end subroutine construct_etaX_doubles_q1_ccsd
!
!
   module subroutine construct_etaX_doubles_q2_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q2 term of etaX doubles
!!    Written by Josefine H. Andersen, Feb 2019
!!    
!!    Q2 = sum_c tb_aicj X_cb - sum_k tb_aibk X_jk
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: etaX_temp
      real(dp), dimension(:,:), allocatable :: etaX_ai_bj
      real(dp), dimension(:,:), allocatable :: etaX_aij_b
!
      real(dp), dimension(:,:), allocatable :: tb_ai_bj
      real(dp), dimension(:,:), allocatable :: tb_aij_c
!
      real(dp), dimension(:,:), allocatable :: X_cb
      real(dp), dimension(:,:), allocatable :: X_jk
!
      real(dp), parameter :: one = 1.0
!
      call mem%alloc(tb_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      tb_ai_bj = zero
      call squareup(wf%t2bar, tb_ai_bj, (wf%n_v)*(wf%n_o))
!
!     Reorder multipiers to tb_aij_c
!
      call mem%alloc(tb_aij_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      tb_aij_c = zero
      call sort_1234_to_1243(tb_ai_bj, tb_aij_c, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: First term: sum_c b_aicj X_cb 
!
      call mem%alloc(X_cb, wf%n_v, wf%n_v) 
      call wf%get_operator_vv(Xoperator, X_cb)
!
      call mem%alloc(etaX_aij_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
      etaX_aij_b = zero
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_v,               &
                  one,                  &
                  tb_aij_c,             &
                  (wf%n_v)*(wf%n_o)**2, &
                  X_cb,                 &
                  wf%n_v,               &
                  zero,                 &
                  etaX_aij_b,           &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(tb_aij_c, (wf%n_v)*(wf%n_o)**2, wf%n_v)      
      call mem%dealloc(X_cb, wf%n_v, wf%n_v)
!
!     Add etaX_aij_b to temporary etaX
!
      call mem%alloc(etaX_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      etaX_ai_bj = zero
      call add_1243_to_1234(one, etaX_aij_b, etaX_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(etaX_aij_b, (wf%n_v)*(wf%n_o)**2, wf%n_v)
!
!     :: Second term: -sum_k tb_aick X_jk
!
      call mem%alloc(X_jk, wf%n_o, wf%n_o)
      call wf%get_operator_oo(Xoperator, X_jk)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  tb_ai_bj,             &
                  (wf%n_o)*(wf%n_v)**2, &
                  X_jk,                 &
                  wf%n_o,               &
                  one,                  &
                  etaX_ai_bj,           &
                  (wf%n_o)*(wf%n_v)**2)
!          
      call mem%dealloc(tb_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      call mem%dealloc(X_jk, wf%n_o, wf%n_o)
!
      call mem%alloc(etaX_temp, wf%n_t2, 1)
      etaX_temp = zero
      call symmetrize_and_add_to_packed(etaX_temp, etaX_ai_bj, wf%n_o*wf%n_v)      
!
      call mem%dealloc(etaX_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     add temporary etaX to etaX vector
!
      call daxpy(wf%n_t2, -one, etaX_temp, 1, etaX(wf%n_t1+1:, 1), 1)
!
      call mem%dealloc(etaX_temp, wf%n_t2, 1)
!      
   end subroutine construct_etaX_doubles_q2_ccsd
!
!
   module subroutine construct_csiX_ccsd(wf, Xoperator, csiX)
!!
!!    Construct csiX
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: csiX
!
      call wf%construct_csiX_singles(Xoperator, csiX)
!      
      call wf%construct_csiX_doubles(Xoperator, csiX)
!
   end subroutine construct_csiX_ccsd
!
!
   module subroutine construct_csiX_singles_ccsd(wf, Xoperator, csiX)
!!
!!    Construct csiX singles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
!!    X_ai + sum_ck u_aick X_kc
!!    
!!    where u_aick = 2t_ck_ai - t_ci_ak
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: csiX
!      
      real(dp), dimension(:,:), allocatable   :: csiX_temp
!      
      real(dp), dimension(:,:), allocatable   :: u_ai_ck
      real(dp), dimension(:,:), allocatable   :: t_ai_ck
!      
      real(dp), dimension(:,:), allocatable   :: X_kc
      real(dp), dimension(:,:), allocatable   :: X_ck
!
      real(dp), parameter :: one = 1.0, two = 2.0
!
      call mem%alloc(csiX_temp, wf%n_t1, 1)
      csiX_temp = zero
!
!     :: First term: X_ai
!
      call wf%get_operator_vo(Xoperator, csiX_temp)     
!
!     :: Second term: u_aick X_kc
!
      call mem%alloc(X_kc, wf%n_o*wf%n_v, 1)
      call mem%alloc(X_ck, wf%n_v*wf%n_o, 1)
!      
      call wf%get_operator_ov(Xoperator, X_kc)
      call sort_12_to_21(X_kc, X_ck, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_kc, wf%n_o*wf%n_v, 1)
!
      call mem%alloc(t_ai_ck, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
      call squareup(wf%t2, t_ai_ck,  wf%n_o*wf%n_v)
!
!     Form u_aick = 2 t_ai_ck - t_ak_ci
!
      call mem%alloc(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      u_ai_ck = zero
!
      call add_1432_to_1234(-one, t_ai_ck, u_ai_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o) 
      call daxpy((wf%n_o)**2 * (wf%n_v)**2, two, t_ai_ck, 1, u_ai_ck, 1)
!
      call mem%dealloc(t_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Form sum_ck u_ai_ck X_kc
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_ai_ck,           &
                  (wf%n_o)*(wf%n_v), &
                  X_ck,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  csiX_temp,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(X_ck, wf%n_v*wf%n_o, 1)
!
!     Add temporary csiX to csiX
!
      call daxpy(wf%n_t1, one, csiX_temp, 1, csiX(1:wf%n_t1, 1), 1)
!
      call mem%dealloc(csiX_temp, wf%n_t1, 1)
!
   end subroutine construct_csiX_singles_ccsd
!
!
   module subroutine construct_csiX_doubles_ccsd(wf, Xoperator, csiX)
!!
!!    Construct csiX doubles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
!!    sum_c t_aicj X_bc - sum_k t_aibk X_kj
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: csiX
!
      real(dp), dimension(:,:), allocatable :: csiX_temp
      real(dp), dimension(:,:), allocatable :: csiX_bj_ai
      real(dp), dimension(:,:), allocatable :: csiX_ai_bj
!
      real(dp), dimension(:,:), allocatable :: t_cj_ai
!
      real(dp), dimension(:,:), allocatable :: X_bc
      real(dp), dimension(:,:), allocatable :: X_kj
!
      real(dp), parameter :: one = 1.0
!
      call mem%alloc(t_cj_ai, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      call squareup(wf%t2, t_cj_ai, wf%n_v*wf%n_o)
!
!     :: First term: sum_c t_aicj X_bc
!
      call mem%alloc(X_bc, wf%n_v, wf%n_v)
      call wf%get_operator_vv(Xoperator, X_bc)
!
      call mem%alloc(csiX_bj_ai, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      csiX_bj_ai = zero
!
      call dgemm('N','N',            &
                 wf%n_v,             &
                 wf%n_v*(wf%n_o)**2, &
                 wf%n_v,             &
                 one,                &
                 X_bc,               &
                 wf%n_v,             &
                 t_cj_ai,            &
                 wf%n_v,             &
                 zero,               &
                 csiX_bj_ai,         &
                 wf%n_v)
!         
      call mem%dealloc(X_bc, wf%n_v, wf%n_v)
!
!     Add csiX_bj_ai to temporary csiX vector
!
      call mem%alloc(csiX_temp, wf%n_t2, 1)
      csiX_temp = zero
      call symmetrize_and_add_to_packed(csiX_temp, csiX_bj_ai, wf%n_v*wf%n_o)
!
      call mem%dealloc(csiX_bj_ai, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!      
!     :: Second term: -sum_k t_ai_bk X_kj
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
      call wf%get_operator_oo(Xoperator, X_kj) 
!
      call mem%alloc(csiX_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      csiX_ai_bj = zero
!      
      call dgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one,                 &
                 t_cj_ai,              &
                 wf%n_o*(wf%n_v)**2,   &
                 X_kj,                 &
                 wf%n_o,               &
                 zero,                 &
                 csiX_ai_bj,           &
                 wf%n_o*(wf%n_v)**2)
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(t_cj_ai, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     Add csi_ai_bj to temporary csiX vector
!
      call symmetrize_and_add_to_packed(csiX_temp, csiX_ai_bj, wf%n_v*wf%n_o)
!
      call mem%dealloc(csiX_ai_bj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     Add doubles contribution to csiX
!
      call daxpy(wf%n_t2, one, csiX_temp, 1, csiX(wf%n_t1+1:wf%n_es_amplitudes, 1), 1)
!
      call mem%dealloc(csiX_temp, wf%n_t2, 1)
!
   end subroutine construct_csiX_doubles_ccsd
!
!
   module subroutine get_eom_contribution_ccsd(wf, etaX, csiX, Xoperator)
!!
!!    Add EOM contribution to csiX vector, CCSD
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(inout), optional :: Xoperator
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(in)    :: csiX
!
!
      !call wf%get_eom_xcc_contribution(etaX, csiX)
      call wf%get_eom_doubles_contribution(Xoperator, etaX)
!
      call wf%get_eom_xcc_contribution(etaX, csiX)
!
   end subroutine get_eom_contribution_ccsd
!
!
   module subroutine get_eom_doubles_contribution_ccsd(wf, Xoperator, etaX)
!!
!!    Build EOM contribution to etaX in CCSD
!!    Written by Josefine H. Andersen, February 2019
!!
!!    = sum_ck tb_aick X_ck + sum_ckdl tb_aick u_ckdl X_ld
!!
!!    where u_ckdl = 2*t_ckdl - t_cldk
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: etaX_temp
!
      real(dp), dimension(:,:), allocatable :: tb_ai_bj
      real(dp), dimension(:,:), allocatable :: t_ck_dl
      real(dp), dimension(:,:), allocatable :: u_ck_dl
!
      real(dp), dimension(:,:), allocatable :: I_ai_dl
!
      real(dp), dimension(:,:), allocatable :: X_ck
      real(dp), dimension(:,:), allocatable :: X_ld
      real(dp), dimension(:,:), allocatable :: X_dl
!
      real(dp), parameter :: one = 1.0, two = 2.0
!
!     :: First term: sum_ck tb_aick X_ck
!
      call mem%alloc(tb_ai_bj, wf%n_v*wf%n_o,  wf%n_v*wf%n_o)
      call squareup(wf%t2bar, tb_ai_bj, wf%n_v*wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v*wf%n_o, 1)
      call wf%get_operator_vo(Xoperator, X_ck)
!
      call mem%alloc(etaX_temp, wf%n_v*wf%n_o, 1)
!      
      call dgemm('N','N',           &
                 (wf%n_v)*(wf%n_o), &
                 1,                 &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 tb_ai_bj,          &
                 (wf%n_v)*(wf%n_o), &
                 X_ck,              &
                 (wf%n_v)*(wf%n_o), &
                 zero,              &
                 etaX_temp,         &
                 (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_ck, wf%n_v*wf%n_o, 1)
!
!     :: Second term: sum_ckdl tb_aick u_ckdl X_ld
!
      call mem%alloc(X_ld, wf%n_o, wf%n_v)
      call wf%get_operator_ov(Xoperator, X_ld)
!      
      call mem%alloc(X_dl, wf%n_v, wf%n_o)
      call sort_12_to_21(X_ld, X_dl, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_ld, wf%n_o, wf%n_v)
!
!     Form u_ai_ck
!
      call mem%alloc(t_ck_dl, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      call squareup(wf%t2, t_ck_dl, wf%n_v*wf%n_o)
!
      call mem%alloc(u_ck_dl, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      u_ck_dl = zero
!
      call add_1432_to_1234(-one, t_ck_dl, u_ck_dl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ck_dl, 1, u_ck_dl, 1)
!
      call mem%dealloc(t_ck_dl, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     Form the intermediate I_ai_dl = sum_ck tb_ai_ck u_ck_dl
!
      call mem%alloc(I_ai_dl, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!      
      call dgemm('N','N',           &
                 (wf%n_v)*(wf%n_o), &
                 (wf%n_v)*(wf%n_o), &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 tb_ai_bj,          &
                 (wf%n_v)*(wf%n_o), &
                 u_ck_dl,           &
                 (wf%n_v)*(wf%n_o), &
                 zero,              &
                 I_ai_dl,           &
                 (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(tb_ai_bj, wf%n_v*wf%n_o,  wf%n_v*wf%n_o)
      call mem%dealloc(u_ck_dl, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
!     Form sum_dl I_ai_dl X_ld^T
!
      call dgemm('N','N',           &
                 (wf%n_v)*(wf%n_o), &
                 1,                 &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 I_ai_dl,           &
                 (wf%n_v)*(wf%n_o), &
                 X_dl,              &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 etaX_temp,         &
                 (wf%n_v)*(wf%n_o))
!
            call mem%dealloc(I_ai_dl, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
            call mem%dealloc(X_dl, wf%n_v, wf%n_o)
!
!     Add EOM contribution to etaX
!
      call daxpy(wf%n_t1, -one, etaX_temp, 1, etaX, 1)
!
      call mem%dealloc(etaX_temp, wf%n_v*wf%n_o, 1)
!
   end subroutine get_eom_doubles_contribution_ccsd
!
!
end submodule properties_ccsd
