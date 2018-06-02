submodule (ccs_class) fock
!
!!
!!    Fock submodule
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    initialize_fock_matrix:  allocates and zeroes fock matrix blocks
!!    destruct_fock_matrix:    deallocates fock matrix blocks
!!
!!    construct_fock:          constructs T1_transformed Fock matrix.
!!    one_electron_t1:         T1-transformation of h_pq (called by construct_fock).
!!
!
   implicit none
!
   character(len=40) :: integral_type
!
contains
!
!
    module subroutine initialize_fock_matrix_ccs(wf)
!!
!!     Initialize Fock matrix
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Allocates and sets Fock matrix blocks (ij, ia, ai, ab) to zero
!!     before calling the Fock matrix constructor.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ij)) call wf%mem%alloc(wf%fock_ij, wf%n_o, wf%n_o)
      if (.not. allocated(wf%fock_ia)) call wf%mem%alloc(wf%fock_ia, wf%n_o, wf%n_v)
      if (.not. allocated(wf%fock_ai)) call wf%mem%alloc(wf%fock_ai, wf%n_v, wf%n_o)
      if (.not. allocated(wf%fock_ab)) call wf%mem%alloc(wf%fock_ab, wf%n_v, wf%n_v)
!
      wf%fock_ij = zero
      wf%fock_ia = zero
      wf%fock_ai = zero
      wf%fock_ab = zero
!
   end subroutine initialize_fock_matrix_ccs
!
!
    module subroutine destruct_fock_matrix_ccs(wf)
!!
!!     Destruct Fock matrix
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!     Deallocates the Fock matrix blocks if they are allocated.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ij)) call wf%mem%dealloc(wf%fock_ij, wf%n_o, wf%n_o)
      if (allocated(wf%fock_ia)) call wf%mem%dealloc(wf%fock_ia, wf%n_o, wf%n_v)
      if (allocated(wf%fock_ai)) call wf%mem%dealloc(wf%fock_ai, wf%n_v, wf%n_o)
      if (allocated(wf%fock_ab)) call wf%mem%dealloc(wf%fock_ab, wf%n_v, wf%n_v)
!
   end subroutine destruct_fock_matrix_ccs
!
!
   module subroutine construct_fock_ccs(wf)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Constructs the T1-transformed Fock matrix blocks (occ/vir-occ/vir),
!!    and saves the result in the class variables fock_pq.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: fock_ao
      real(dp), dimension(:,:), allocatable :: fock_matrix
!
      type(file) :: h1ao_file
!
      real(dp), dimension(:,:), allocatable :: h1ao    ! AO matrix h_αβ
      real(dp), dimension(:,:), allocatable :: h1ao_sq ! Unpacked AO matrix h_αβ
      real(dp), dimension(:,:), allocatable :: h1mo    ! MO matrix h_pq
      real(dp), dimension(:,:), allocatable :: X       ! An intermediate
!
!     Indices
!
      integer(i15) :: i = 0, j = 0, k = 0, a = 0, b = 0
      integer(i15) :: kj = 0, ii = 0, ij = 0, kk = 0, ik = 0
      integer(i15) :: jj = 0, ji = 0, ai = 0, ib = 0, bi = 0
      integer(i15) :: aj = 0, ja = 0, ab = 0, ia = 0
!
!     Useful orbital information
!
      integer(i15) :: n_ao_sq_packed = 0 ! Dimension of packed (n_ao x n_ao) matrix
!
!     Two electron integrals
!
      real(dp), dimension(:,:), allocatable :: g_ij_kl
      real(dp), dimension(:,:), allocatable :: g_ab_ij
      real(dp), dimension(:,:), allocatable :: g_ai_jb
      real(dp), dimension(:,:), allocatable :: g_ia_jk
      real(dp), dimension(:,:), allocatable :: g_ai_jk
!
!     :: One-electron contribution ::
!
!     Open file to read the one-electron AO integrals
!
      h1ao_file%name = 'MLCC_AOINT'
      call wf%disk%open_file(h1ao_file, 'formatted', 'readwrite', 'sequential')
!
      n_ao_sq_packed = packed_size(wf%n_ao)
!
      call wf%mem%alloc(h1ao, n_ao_sq_packed, 1)
      h1ao = zero
!
      read(h1ao_file%unit, *) (h1ao(i,1), i = 1, n_ao_sq_packed)
!
      call wf%disk%close_file(h1ao_file)
!
!     Squareup h_αβ
!
      call wf%mem%alloc(h1ao_sq, wf%n_ao, wf%n_ao)
      h1ao_sq = zero
!
      call squareup(h1ao, h1ao_sq, wf%n_ao)
!
      call wf%mem%dealloc(h1ao, n_ao_sq_packed, 1)
!
!     Transform to MO basis, h_αβ -> h_pq
!
      call wf%mem%alloc(X, wf%n_ao, wf%n_mo) ! Intermediate
!
      call dgemm('N','N',     &
                  wf%n_ao,    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  one,        &
                  h1ao_sq,    &
                  wf%n_ao,    &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  zero,       &
                  X,          &
                  wf%n_ao)
!
      call wf%mem%dealloc(h1ao_sq, wf%n_ao, wf%n_ao)
!
      call wf%mem%alloc(h1mo, wf%n_mo, wf%n_mo)
!
      call dgemm('T','N',     &
                  wf%n_mo,    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  one,        &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  X,          &
                  wf%n_ao,    &
                  zero,       &
                  h1mo,       &
                  wf%n_mo)
!
!     T1-transform h_pq, then copy it to fock_matrix array
!
      call wf%mem%alloc(fock_matrix, wf%n_mo, wf%n_mo)
      fock_matrix = zero
!
      call wf%one_electron_t1(h1mo, fock_matrix)
!
      call wf%mem%dealloc(h1mo, wf%n_mo, wf%n_mo)
!
!     Deallocate intermediate X and fock_ao
!
      call wf%mem%dealloc(X,wf%n_ao, wf%n_mo)
!
!     :: Two-electron occupied-occupied block: F_ij = h_ij + sum_k (2*g_ijkk - g_ikkj) ::
!
      call wf%mem%alloc(g_ij_kl, (wf%n_o)**2, (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_oo_oo(integral_type, g_ij_kl)
!
!     Add two-electron contributions to occupied-occupied block
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            ij = index_two(i, j, wf%n_o)
!
            do k = 1, wf%n_o
!
               kk = index_two(k, k, wf%n_o)
               ik = index_two(i, k, wf%n_o)
               kj = index_two(k, j, wf%n_o)
!
               fock_matrix(i,j) = fock_matrix(i,j) + &
                                    two*g_ij_kl(ij,kk) - &
                                    g_ij_kl(ik,kj)
!
            enddo
!
         enddo
!
      enddo
!
      call wf%mem%dealloc(g_ij_kl, (wf%n_o)**2, (wf%n_o)**2)
!
!     :: Two-electron occupied-virtual blocks ::
!
      call wf%mem%alloc(g_ia_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_oo(integral_type, g_ia_jk)
!
      call wf%mem%alloc(g_ai_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_oo(integral_type, g_ai_jk)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
               ia = index_two(i, a, wf%n_o)
               ja = index_two(j, a, wf%n_o)
!
               ai = index_two(a, i, wf%n_v)
               aj = index_two(a, j, wf%n_v)
!
               jj = index_two(j, j, wf%n_o)
               ji = index_two(j, i, wf%n_o)
               ij = index_two(i, j, wf%n_o)
!
               fock_matrix(i, a + wf%n_o) = fock_matrix(i, a + wf%n_o) + &
                                                 two*g_ia_jk(ia,jj) - g_ia_jk(ja,ij) ! g_ia_jk(ja,ij) = g_jaij = g_ijja
!
               fock_matrix(a + wf%n_o, i) = fock_matrix(a + wf%n_o, i) + &
                                                 two*g_ai_jk(ai,jj) - g_ai_jk(aj,ji)
!
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ia_jk, (wf%n_o)*(wf%n_v), (wf%n_o)**2)
      call wf%mem%dealloc(g_ai_jk, (wf%n_v)*(wf%n_o), (wf%n_o)**2)
!
!     :: Two-electron virtual-virtual block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) ::
!
      call wf%mem%alloc(g_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_vv_oo(integral_type, g_ab_ij)
!
      call wf%mem%alloc(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_ov(integral_type, g_ai_jb)
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            ab = index_two(a, b, wf%n_v)
!
            do i = 1, wf%n_o
!
               ii = index_two(i, i, wf%n_o)
               ai = index_two(a, i, wf%n_v)
               bi = index_two(b, i, wf%n_v)
               ia = index_two(i, a, wf%n_o)
               ib = index_two(i, b, wf%n_o)
!
               fock_matrix(wf%n_o + a, wf%n_o + b) = fock_matrix(wf%n_o + a, wf%n_o + b) &
                                                       + two*g_ab_ij(ab,ii) &
                                                       - g_ai_jb(ai,ib)
!
            enddo
!
         enddo
      enddo
!
     call wf%mem%dealloc(g_ab_ij, (wf%n_v)**2, (wf%n_o)**2)
     call wf%mem%dealloc(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Save the blocks of the Fock matrix in memory (ij, ia, ai, ab)
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            wf%fock_ij(i,j) = fock_matrix(i,j)
!
         enddo
      enddo
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ia(i,a) = fock_matrix(i, wf%n_o + a)
            wf%fock_ai(a,i) = fock_matrix(wf%n_o + a, i)
!
         enddo
      enddo
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            wf%fock_ab(a,b) = fock_matrix(wf%n_o + a, wf%n_o + b)
!
         enddo
      enddo
!
      call wf%mem%dealloc(fock_matrix, wf%n_mo, wf%n_mo)
!
   end subroutine construct_fock_ccs
!
!
   module subroutine one_electron_t1_ccs(wf, h1 ,h1_T1)
!!
!!    One-electron T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    T1-transforms the one-electron MO integrals h_pq
!!
!!       h_p_q_T1 = sum_st x_p_s * y_q_t * h_s_t,
!!
!!    where
!!
!!       x = I - t1,
!!       y = I - t1^T.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo) :: h1
      real(dp), dimension(wf%n_mo, wf%n_mo) :: h1_T1
!
      real(dp), dimension(:,:), allocatable :: x
      real(dp), dimension(:,:), allocatable :: y
      real(dp), dimension(:,:), allocatable :: t1
!
      real(dp), dimension(:,:), allocatable :: Z ! Intermediate for matrix multiplication
!
      integer(i15) :: p = 0, q = 0, a = 0, i = 0
!
!     Allocate the arrays t1, x, and y
!
      call wf%mem%alloc(t1, wf%n_mo, wf%n_mo)
      t1 = zero
!
      call wf%mem%alloc(y, wf%n_mo, wf%n_mo)
      call wf%mem%alloc(x, wf%n_mo, wf%n_mo)
!
      y = zero
      x = zero
!
!     Set t1_p_q = t1am_p_q for p virtual and q occupied, 0 otherwise
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            t1(wf%n_o + a, i) = wf%t1am(a, i)
!
         enddo
      enddo
!
!     Form the x and y arrays
!
      do p = 1, wf%n_mo
         do q = 1, wf%n_mo
!
            if (p .eq. q) then
!
               x(p,q) = 1
               y(p,q) = 1
!
            else
!
               x(p,q) = x(p,q) - t1(p,q)
               y(p,q) = y(p,q) + t1(q,p)
!
            endif
!
         enddo
      enddo
!
!     Deallocate t1 (only x and y are needed below)
!
      call wf%mem%dealloc(t1, wf%n_mo, wf%n_mo)
!
!     Allocate Z intermediate
!
      call wf%mem%alloc(Z, wf%n_mo, wf%n_mo)
!
!     Calculate h1_T1 = x*h1*y^T = x*Z
!
      call dgemm('N','T',  &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  h1,      &
                  wf%n_mo, &
                  y,       &
                  wf%n_mo, &
                  zero,    &
                  Z,       &
                  wf%n_mo)
!
      call dgemm('N','N',  &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  x,       &
                  wf%n_mo, &
                  Z,       &
                  wf%n_mo, &
                  zero,    &
                  h1_T1,   &
                  wf%n_mo)
!
!     Deallocate x and y, and the intermediate Z
!
      call wf%mem%dealloc(Z, wf%n_mo, wf%n_mo)
      call wf%mem%dealloc(y, wf%n_mo, wf%n_mo)
      call wf%mem%dealloc(x, wf%n_mo, wf%n_mo)
!
   end subroutine one_electron_t1_ccs
!
!
end submodule fock
