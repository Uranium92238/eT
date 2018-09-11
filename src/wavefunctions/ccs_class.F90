module ccs_class
!
!!
!!    Coupled cluster singles (ccs) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
   use hf_class
!
   use mo_integral_tool_class
!
   use reordering
   use array_utilities
   use array_analysis
   use interval_class
   use index
!
   implicit none
!
   type, extends(wavefunction):: ccs
!
!     The T1-transformed Fock matrix
!
      real(dp), dimension(:,:), allocatable :: fock_ij
      real(dp), dimension(:,:), allocatable :: fock_ia
      real(dp), dimension(:,:), allocatable :: fock_ai
      real(dp), dimension(:,:), allocatable :: fock_ab
!
      real(dp), dimension(:,:) , allocatable :: fock_diagonal
!
      type(mo_integral_tool) :: integrals
!
      real(dp), dimension(:,:), allocatable :: t1
!
   contains
!
!     Prepare and cleanup wavefunction
!
      procedure :: prepare => prepare_ccs
      procedure :: cleanup   => cleanup_ccs
!
      procedure :: initialize_fock_ij => initialize_fock_ij_ccs
      procedure :: initialize_fock_ia => initialize_fock_ia_ccs
      procedure :: initialize_fock_ai => initialize_fock_ai_ccs
      procedure :: initialize_fock_ab => initialize_fock_ab_ccs
!
      procedure :: initialize_fock_diagonal => initialize_fock_diagonal_ccs
!
      procedure :: initialize_t1 => initialize_t1_ccs
!
      procedure :: destruct_fock_ij => destruct_fock_ij_ccs
      procedure :: destruct_fock_ia => destruct_fock_ia_ccs
      procedure :: destruct_fock_ai => destruct_fock_ai_ccs
      procedure :: destruct_fock_ab => destruct_fock_ab_ccs
!
      procedure :: destruct_fock_diagonal => destruct_fock_diagonal_ccs
!
      procedure :: destruct_t1 => destruct_t1_ccs
!
      procedure :: get_ovov => get_ovov_ccs
!
   end type ccs
!
!
contains
!
!
   subroutine prepare_ccs(wf, ref_wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      class(hf) :: ref_wf
!
      integer(i15) :: p
!
      wf%name = 'ccs'
!
      wf%system = ref_wf%system
!
      wf%n_ao   = ref_wf%n_ao
      wf%n_mo   = ref_wf%n_mo
      wf%n_o    = ref_wf%n_o
      wf%n_v    = ref_wf%n_v
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
      call wf%initialize_fock_diagonal()
!
      wf%fock_ij(:,:) = ref_wf%mo_fock(1 : wf%n_o , 1 : wf%n_o)
      wf%fock_ia(:,:) = ref_wf%mo_fock(1 : wf%n_o , wf%n_o + 1 : wf%n_v)
      wf%fock_ai(:,:) = ref_wf%mo_fock(1 : wf%n_o , wf%n_o + 1 : wf%n_v)
      wf%fock_ab(:,:) = ref_wf%mo_fock( wf%n_o + 1 : wf%n_v , wf%n_o + 1 : wf%n_v)
!
      do p = 1, wf%n_mo
!
         wf%fock_diagonal(p, 1) = ref_wf%mo_fock(p, p)
!
      enddo
!
      call wf%initialize_orbital_coefficients()
      wf%orbital_coefficients = ref_wf%orbital_coefficients
!
   end subroutine prepare_ccs
!
!
   subroutine cleanup_ccs(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
!     Nothing here yet
!
   end subroutine cleanup_ccs
!
!
   subroutine initialize_fock_ij_ccs(wf)
!!
!!    Initialize Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ij)) call mem%alloc(wf%fock_ij, wf%n_o, wf%n_o)
!
   end subroutine initialize_fock_ij_ccs
!
!
   subroutine initialize_fock_ia_ccs(wf)
!!
!!    Initialize Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ia)) call mem%alloc(wf%fock_ia, wf%n_o, wf%n_v)
!
   end subroutine initialize_fock_ia_ccs
!
!
   subroutine initialize_fock_ai_ccs(wf)
!!
!!    Initialize Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ai)) call mem%alloc(wf%fock_ai, wf%n_v, wf%n_o)
!
   end subroutine initialize_fock_ai_ccs
!
!
   subroutine initialize_fock_ab_ccs(wf)
!!
!!    Initialize Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ab)) call mem%alloc(wf%fock_ab, wf%n_v, wf%n_v)
!
   end subroutine initialize_fock_ab_ccs
!
!
   subroutine initialize_fock_diagonal_ccs(wf)
!!
!!    Initialize Fock diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_diagonal)) call mem%alloc(wf%fock_diagonal, wf%n_mo, 1)
!
   end subroutine initialize_fock_diagonal_ccs
!
!
   subroutine initialize_t1_ccs(wf)
!!
!!    Initialize T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%t1)) call mem%alloc(wf%t1, wf%n_v, wf%n_o)
!
   end subroutine initialize_t1_ccs
!
!
   subroutine destruct_fock_ij_ccs(wf)
!!
!!    Destruct Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ij)) call mem%dealloc(wf%fock_ij, wf%n_o, wf%n_o)
!
   end subroutine destruct_fock_ij_ccs
!
!
   subroutine destruct_fock_ia_ccs(wf)
!!
!!    Destruct Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ia)) call mem%dealloc(wf%fock_ia, wf%n_o, wf%n_v)
!
   end subroutine destruct_fock_ia_ccs
!
!
   subroutine destruct_fock_ai_ccs(wf)
!!
!!    Destruct Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ai)) call mem%dealloc(wf%fock_ai, wf%n_v, wf%n_o)
!
   end subroutine destruct_fock_ai_ccs
!
!
   subroutine destruct_fock_ab_ccs(wf)
!!
!!    Destruct Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ab)) call mem%dealloc(wf%fock_ij, wf%n_v, wf%n_v)
!
   end subroutine destruct_fock_ab_ccs
!
!
   subroutine destruct_fock_diagonal_ccs(wf)
!!
!!    Destruct Fock diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_diagonal)) call mem%dealloc(wf%fock_diagonal, wf%n_mo, 1)
!
   end subroutine destruct_fock_diagonal_ccs
!
!
   subroutine destruct_t1_ccs(wf)
!!
!!    Destruct T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1)) call mem%dealloc(wf%t1, wf%n_v, wf%n_o)
!
   end subroutine destruct_t1_ccs
!
!
   subroutine get_ovov_ccs(wf, g_iajb, first_i, last_i, first_a, last_a, &
                                       first_j, last_j, first_b, last_b)
!!
!!    Get ovov 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. If the 
!!    full integral g_pqrs is not on file, the construction is done as 
!!
!!       g_pqrs = sum_J L_pq^J L_rs^J,
!!
!!    where L_pq^J is t1-transformed (either made on-the-fly from the MO L_pq^J or
!!    read from disk, usually before and after the ground state amplitudes have 
!!    converged). On-the-fly t1-transformation can also happen due to disk space 
!!    limitations.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:) :: g_iajb 
!
      integer(i15), optional, intent(in) :: first_i, last_i 
      integer(i15), optional, intent(in) :: first_a, last_a
      integer(i15), optional, intent(in) :: first_j, last_j
      integer(i15), optional, intent(in) :: first_b, last_b
!
!     Special case: the t1-transformed g_iajb integrals are equal to the MO integrals 
!
   !   call wf%integrals%read_ovov(g_iajb, first_i, last_i, first_a, last_a, &
   !                                       first_j, last_j, first_b, last_b)
!
   end subroutine get_ovov_ccs
!
!
end module ccs_class
