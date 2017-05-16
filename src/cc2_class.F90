module cc2_class
!
!!
!!            Coupled cluster perturbative doubles (CC2) class module                                 
!!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
!!                                                                           
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  The ancestor class module (CCS)
!
   use ccs_class
!
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CC2 class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(ccs) :: cc2
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_cc2
      procedure :: drv  => drv_cc2
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_cc2
!
!     Helper routines for construct_omega
!
      procedure :: omega_a1 => omega_a1_cc2 
      procedure :: omega_b1 => omega_b1_cc2 
      procedure :: omega_c1 => omega_c1_cc2 
      procedure :: omega_d1 => omega_d1_cc2      
!
!     Ground state solver helper routines
!
      procedure :: calc_energy => calc_energy_cc2
!
   end type cc2
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
!
      module subroutine construct_omega_cc2(wf)
!!
!!        Construct Omega 
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!        Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!        for the current amplitudes of the object wf 
!!
         implicit none
!
         class(cc2) :: wf 
!
      end subroutine construct_omega_cc2
!
!
      module subroutine omega_a1_cc2(wf, t_kc_di, c_first, c_last, c_length)
!!
!!        Omega A1
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!  
!!        Calculates the A1 term of omega, 
!!  
!!        A1: sum_ckd g_adkc * u_ki^cd,
!! 
!!        and adds it to the projection vector (omega1) of
!!        the wavefunction object wf
!!
!!        u_ki^cd = 2*t_ki^cd - t_ik^cd 
!!
         implicit none
!
         class(cc2) :: wf
!
!        Batching variable for double amplitudes t_kc_di
!  
         integer(i15) :: c_first, c_last, c_length
!
         real(dp), dimension(c_length*(wf%n_o), (wf%n_v)*(wf%n_o)) :: t_kc_di
!
      end subroutine omega_a1_cc2
!
!
      module subroutine omega_b1_cc2(wf, t_lc_ak, c_first, c_last, c_length)
!!
!!        Omega B1
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!        Calculates the B1 term of omega, 
!!
!!        B1: - sum_ckl u_kl^ac * g_kilc,
!! 
!!        and adds it to the projection vector (omeg1) of
!!        the wavefunction object wf
!!
!!        u_kl^ac = 2*t_kl^ac - t_lk^ac 
!!
         implicit none
!
         class(cc2) :: wf 
!
!        Batching variable for double amplitudes t_kc_di
!  
         integer(i15) :: c_first, c_last, c_length
!
         real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_lc_ak
!
      end subroutine omega_b1_cc2
!
!
      module subroutine omega_c1_cc2(wf, t_kc_ai, c_first, c_last, c_length)
!!
!!        Omega C1
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!        Calculates the C1 term of omega,
!! 
!!        C1: sum_ck F_kc*u_ai_ck,
!!
!!        and adds it to the projection vector (omega1) of    
!!        the wavefunction object wf                           
!!
!!        u_ai_ck = 2*t_ck_ai - t_ci_ak
!!
         implicit none
!       
         class(cc2) :: wf 
!
!        Batching variable for double amplitudes t_kc_di
!  
         integer(i15) :: c_first, c_last, c_length
!
         real(dp), dimension(c_length*(wf%n_o),(wf%n_v)*(wf%n_o)) :: t_kc_ai
!
      end subroutine omega_c1_cc2
!
!
      module subroutine omega_d1_cc2(wf)
!!
!!        Omega D1
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!        Calculates the D1 term of omega,
!!
!!        D1: F_ai_T1
!!
!!        and adds it to the projection vector (omega1) of
!!        the wavefunction object wf 
!!
         implicit none
!
         class(cc2) :: wf
!
      end subroutine omega_d1_cc2
!
!
   end interface
!
!
contains
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
   subroutine init_cc2(wf)
!!
!!     Initialize CC2 object
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!     Performs the following tasks
!!
!!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!!        the descendant classes) 
!!     4. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!!     5. Allocate Omega vector
!!
      implicit none
!
      class(cc2)  :: wf
!
      write(unit_output,*)'In init_cc2'
      flush(unit_output)
!
!     Set model name
!
      wf%name = 'CC2'
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
!
!     Allocate Fock matrix and set to zero
!
      call wf%initialize_fock_matrix
!
!     Initialize omega vector
!
      call wf%initialize_omega
!
   end subroutine init_cc2
!
!
   subroutine drv_cc2(wf)
!!
!!     CC2 Driver
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
      implicit none 
!
      class(cc2) :: wf
!
      write(unit_output,*)'In drv_cc2'
      flush(unit_output)
      call wf%ground_state_solver
!
   end subroutine drv_cc2
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine calc_energy_cc2(wf)
!!
!!    Calculate Energy (CC2)
!!
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the CC2 energy, 
!!
!!    E_CC2 = E_HF + sum_aibj L_iajb*(t_ij^ab + t_i^a*t_j^b),
!!   
!!    with t_ij^ab = - g_aibj/(e_a + e_b - e_i - e_j) where 
!!    g_aibj are T1-transformed integrals.
!!    Batching over a.
!!
   implicit none
!
   class(cc2) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_bj_J
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_bj ! = g_aibj
      real(dp), dimension(:,:), allocatable :: g_ia_jb 
!
!     t2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_ia_bj ! = g_aibj/(e_a + e_b - e_i - e_j)
!
!     Batching variables
!  
      integer(i15) :: a_batch, a_first, a_last, a_length
      integer(i15) :: required, available, n_batch, batch_dimension, max_batch_length
!
!     Indices
!
      integer(i15) :: a = 0, b = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ai = 0, bj = 0
      integer(i15) :: ia = 0, ib = 0, jb = 0, ja = 0
!
      integer(i15) :: aibj = 0
!
!     Prepare for batching over index a
!  
      required = (2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                           !    
               + 2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                            ! Needed for g_aibj  
               + 2*((wf%n_v)**2)*(wf%n_J) + ((wf%n_o)**2)*(wf%n_J) &       ! and 't2' amplitudes  
               + 2*(wf%n_v)**2*(wf%n_o)**2)                                !
!      
      required = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
!        :: Calculate cc2 doubles amplitudes ::
!
!        Allocate L_bj_J and L_ia_J (= reordering of L_bj_J constrained to the batch)
!
         call allocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call allocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
         L_bj_J = zero
         L_ia_J = zero
!
         call wf%get_cholesky_ai(L_bj_J)
!
!        Create L_ia_J
!
         do a = 1, a_length
            do i = 1, wf%n_o
               do J = 1, wf%n_J
!
!                 Calculate compound indices
!
                  ia = index_two(i, a, wf%n_o)
                  ai = index_two(a + a_first - 1, i, wf%n_v)
!
                  L_ia_J(ia, J) = L_bj_J(ai, J)
!
               enddo
            enddo
         enddo
!
!        Allocate g_ia_bj
!
         call allocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Construct integral g_ia_bj (= g_aibj for the batch)
!
         call dgemm('N','T',            &
                     a_length*(wf%n_o), &
                     (wf%n_o)*(wf%n_v), &
                     wf%n_J,            &
                     one,               &
                     L_ia_J,            &
                     a_length*(wf%n_o), &
                     L_bj_J,            &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     g_ia_bj,           &
                     a_length*(wf%n_o))
!
!        L_bj_J and L_ia_J
!
         call deallocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call deallocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
!
!        :: Construct the needed integrals for the enegry ::
!
!        Allocate t_ia_bj
!
         call allocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Create t2 amplitudes
!
         do a = 1, a_length
            do i = 1, wf%n_o
!
               ia = index_two(i, a, wf%n_o)
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     t_ia_bj(ia, bj) = - g_ia_bj(ia, bj)/(wf%fock_diagonal(a + wf%n_o, 1) &
                                          + wf%fock_diagonal(b + wf%n_o, 1) &
                                          - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate g_ia_bj
!
         call deallocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
!         Get the Cholesky vector L_ia_J 
!
         call wf%get_cholesky_ia(L_ia_J)
!
!        Allocate g_ia_jb = g_iajb and set it to zero
!
         call allocator(g_ia_jb, (wf%n_o)*a_length, (wf%n_o)*(wf%n_v))
         g_ia_jb = zero
!
!        Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
         call dgemm('N','T',                                   &
                     (wf%n_o)*a_length,                        &
                     (wf%n_o)*(wf%n_v),                        &
                     wf%n_J,                                   &
                     one,                                      &
                     L_ia_J(index_two(1, a_first, wf%n_o), 1), &
                     (wf%n_o)*(wf%n_v),                        &
                     L_ia_J,                                   &
                     (wf%n_o)*(wf%n_v),                        &
                     zero,                                     &
                     g_ia_jb,                                  &
                     (wf%n_o)*a_length)
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Set the initial value of the energy 
!
      wf%energy = wf%scf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
         do i = 1, wf%n_o
            do a = 1, a_length
!
               ai = index_two(a, i, a_length)
               ia = index_two(i, a, wf%n_o)
!
               do j = 1, wf%n_o
!
                  ja = index_two(j, a, wf%n_o)
!
                 do b = 1, wf%n_v
! 
                     bj = index_two(b, j, wf%n_v)
                     jb = index_two(j, b, wf%n_o)
                     ib = index_two(i, b, wf%n_o)
!
                     aibj = index_packed(ai, bj)
!
!                    Add the correlation energy 
!
                     wf%energy = wf%energy + &
                     (two*g_ia_jb(ia,jb) - g_ia_jb(ja, ib))*(t_ia_bj(ia, bj) + (wf%t1am(a,i))*(wf%t1am(b,j)))
                              
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate g_ia_jb
!
         call deallocator(g_ia_jb, (wf%n_o)*a_length, (wf%n_o)*(wf%n_v))
!
!        Deallocate t_ia_bj
!  
         call deallocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
      enddo ! End of batching
!
   end subroutine calc_energy_cc2
!
!
end module cc2_class