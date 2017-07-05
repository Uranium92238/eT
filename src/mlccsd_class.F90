module mlccsd_class
!
!!
!!                 Multi-level CCSD (MLCCSD) class module                                
!!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
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
   use input_reader
!
!  The ancestor class module (MLCC2)
!
   use mlcc2_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the MLCCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::: 
!
   type, extends(mlcc2) :: mlccsd
!
!     ML-variables
!
      integer(i15), dimension(:,:), allocatable :: n_CCSD_o
      integer(i15), dimension(:,:), allocatable :: n_CCSD_v
!
      integer(i15), dimension(:,:), allocatable :: first_CCSD_o
      integer(i15), dimension(:,:), allocatable :: first_CCSD_v
!
      real(dp), dimension(:,:), allocatable :: mo_coef_cc2_ccs ! MO coefficient matrix for mlccsd basis
      real(dp), dimension(:,:), allocatable :: T ! MO transformation matrix for mlccsd basis
!
!     Fock matrix variables
!
      real(dp), dimension(:,:), allocatable :: fock_diagonal_cc2_ccs  ! diagonal vector for mlccsd basis
!
!     Amplitude variables
!
      integer(i15) :: n_t2am = 0                    ! Number of doubles amplitudes
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector
!
!     Schrödinger equation projection vector (the omega vector)
! 
!        < mu | exp(-T) H exp(T) | R >
!
      real(dp), dimension(:,:), allocatable :: omega2 ! Doubles vector
!
   contains
!
      procedure :: init => init_mlccsd
!
      procedure :: initialize_orbital_info => initialize_orbital_info_mlccsd
      procedure :: initialize_amplitudes   => initialize_amplitudes_mlccsd
      procedure :: initialize_omega        => initialize_omega_mlccsd
!
      procedure :: cholesky_localization_drv          => cholesky_localization_drv_mlccsd
      procedure :: cholesky_localization_CCSD_CC2_CCS => cholesky_localization_CCSD_CC2_CCS_mlccsd
      procedure :: cholesky_localization_CCSD_CCS => cholesky_localization_CCSD_CCS_mlccsd
      procedure :: cholesky_localization_CCSD_CC2 => cholesky_localization_CCSD_CC2_mlccsd
!
      procedure :: read_transform_cholesky_for_CC2_amplitude => read_transform_cholesky_for_CC2_amplitude_mlccsd
      procedure :: read_cholesky_ai_for_cc2_amplitudes       => read_cholesky_ai_for_cc2_amplitudes_mlccsd
      procedure :: get_cholesky_ai_for_cc2_amplitudes        => get_cholesky_ai_for_cc2_amplitudes_mlccsd
!
      procedure :: get_CCSD_active_indices => get_CCSD_active_indices_mlccsd
      procedure :: set_n_total_active      => set_n_total_active_mlccsd
!
      procedure :: omega_mlcc2_a1 => omega_mlcc2_a1_mlccsd
      procedure :: omega_mlcc2_b1 => omega_mlcc2_b1_mlccsd
      procedure :: get_mlccsd_s2am => get_mlccsd_s2am_mlccsd
      procedure :: active_space_t2am_offset => active_space_t2am_offset_mlccsd
!
   end type mlccsd
!
   interface
!
!    -::- Orbital partitioning submodule interface -::-
!    :::::::::::::::::::::::::::::::::::::::::::::::::: 
!
!
      module subroutine cholesky_localization_drv_mlccsd(wf)
!!
!!       Cholesky orbital localization. driver,
!!       Written by Sarai D. Folkestad, June 2017
!!
!!       Driver for Cholesky density decomposition.  
!!
!!       - Collects atom and ao-basis information.
!!       - Constructs occupied and vacant densities.
!!       - Constructs AO Fock matrix.  (This is currently an N^5 operation, should be optimized/removed)
!!       - By looping over active spaces, the occupied and virtual densities are Cholesky decomposed
!!         and the cholesky vectors are used to generate new localized MO's.
!!       - New orbitals are tested for orthonormality (Not implemented yet, only need overlap matrix from DALTON)    
!!
!!
         implicit none
!
         class(mlccsd) :: wf
 !
      end subroutine cholesky_localization_drv_mlccsd
!
!
      module subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!       Cholesky orbital localization. driver,
!!       Written by Sarai D. Folkestad, June 2017
!!
!!       Driver for Cholesky density decomposition.  
!!
         implicit none
!
!        Input arguments
!
         class(mlccsd) :: wf
!
         integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
         integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
         real(dp), dimension(:,:) :: ao_fock
!
         integer(i15) :: n_nuclei
         integer(i15) :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd
!
!
      module subroutine cholesky_localization_CCSD_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!       Cholesky orbital localization. driver,
!!       Written by Sarai D. Folkestad, June 2017
!!
!!
         implicit none
!
!        Input arguments
!
         class(mlccsd) :: wf
!
         integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
         integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
         real(dp), dimension(:,:) :: ao_fock
!
         integer(i15) :: n_nuclei
         integer(i15) :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_CCSD_CCS_mlccsd
!
!
      module subroutine cholesky_localization_CCSD_CC2_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!       Cholesky orbital localization driver
!!       Written by Sarai D. Folkestad, June 2017
!!
!!
         implicit none
!
!        Input arguments
!
         class(mlccsd) :: wf
!
         integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
         integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
         real(dp), dimension(:,:) :: ao_fock
!
         integer(i15) :: n_nuclei
         integer(i15) :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_CCSD_CC2_mlccsd
!
   end interface
!
!
   interface
!
!     -::- Cholesky submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::
!
      module subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd(wf)
!
         implicit none
!
         class(mlccsd) :: wf
!
      end subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd
!
!
      module subroutine read_cholesky_ai_for_cc2_amplitudes_mlccsd(wf,L_ai_J, a_first, a_last, i_first, i_last)
!!
!!       Read Cholesky IA 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Reads the MO Cholesky IA (occ-vir) vectors from file and
!!       places them in the incoming L_ia_J matrix
!!
!!
         implicit none
!
         class(mlccsd)            :: wf
         real(dp), dimension(:,:) :: L_ai_J
         integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
!
!
      end subroutine read_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!
      module subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd(wf, L_ai_J, a_first, a_last, i_first, i_last)
!!
!!       Get Cholesky AI
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Read and T1-transform Cholesky AI vectors:
!!    
!!        L_ai_J_T1 = L_ia_J - sum_j  t_aj*L_ji_J 
!!                            + sum_b  t_bi*L_ab_J
!!                            - sum_bj t_aj*t_bi*L_jb_J
!!
!!       Allocations in routine:
!!
!!       (1) n_J*(i_length)*(a_length) + 2*n_J*(a_length)*batch_length  ->  for L_ab_J contribution (batches of b)
!!       (2) n_J*(i_length)*n_v + 2*n_J*n_o*(i_length)                  ->  for L_ij_J contribution
!!       (3) 2*n_J*n_o*n_v                                              ->  for L_jb_J contribution
!!
!!       i_length = i_last - i_first + 1          
!!       a_length = a_last - a_first + 1          
!!
!!       (1) determines memory requirement. 
!!
         implicit none 
!
         class(mlccsd)            :: wf
         real(dp), dimension(:,:) :: L_ai_J
         integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
!
      end subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd
!
   end interface
!
!
   interface
!
!     -::- Omega submodule interface -::-
!     :::::::::::::::::::::::::::::::::::
!
      module subroutine initialize_omega_mlccsd(wf)
!
!        Initialize Omega (MLCCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Allocates the projection vector (omega1, omega2) and sets it
!        to zero.
!
         implicit none 
!
         class(mlccsd) :: wf
!
      end subroutine initialize_omega_mlccsd
!
      module subroutine omega_mlcc2_a1_mlccsd(wf, active_space)
!! 
!!       Omega A1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!       Calculates the A1 term of omega for the active space, 
!!   
!!       A1: sum_bcj g_Abjc * u_ij^bc,
!!  
!!       and adds it to the projection vector (omega1) of
!!       the wavefunction object wf
!! 
!!       u_ij^bc = 2*s_ij^bc - s_ij^cb 
!!
!!       Batching over A
!!
!! 
         implicit none
!
         class(mlccsd)   :: wf
         integer(i15)    :: active_space ! Current active space IS ALWAYS ONE FOR MLCCSD !
!
      end subroutine omega_mlcc2_a1_mlccsd
!
!
      module subroutine omega_mlcc2_b1_mlccsd(wf, active_space)
!! 
!!       Omega B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!       Calculates the B1 term of omega, 
!!   
!!       B1: - sum_bjk u_jk^ab*g_kbjI + sum_bj u_ij^ab F_jb,
!!
!!       with u_ij^ab = 2*s_ij^ab - s_ij^ba. 
!!
         implicit none
!
         class(mlccsd)   :: wf
         integer(i15)   :: active_space 
!
      end subroutine omega_mlcc2_b1_mlccsd
!
!
      module subroutine get_mlccsd_s2am_mlccsd(wf, s_ia_jb)
!!
!!
         implicit none
!
         class(mlccsd) :: wf
!
         real(dp), dimension(:,:) :: s_ia_jb
!
     end subroutine get_mlccsd_s2am_mlccsd
!
   end interface
!
contains
!
!
  subroutine init_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      integer(i15) :: i, j, active_space
      write(unit_output,*)'Init mlccsd'
      flush(unit_output)
!
!     Set model name 
!
      wf%name = 'MLCCSD'
!
!     MLCC sanity check
!
      if(wf%mlcc_settings%CC3) then
         write(unit_output,*)'WARNING: CC3 active spaces not available for MLCCSD'
         stop
      endif
!
!     Set implemented methods
!
!      wf%implemented%ground_state = .true.
!      wf%implemented%excited_state = .true.
!
!     Read Hartree-Fock info
!
      call wf%read_hf_info
!
!     Orbital partitioning - only if we have CCS/CC2 region
!
      if (wf%mlcc_settings%CCS .or. wf%mlcc_settings%CC2) then
!
         call wf%orbital_partitioning
!
      else
!
         write(unit_output,*)'Full CCSD requested, orbital partitioning skipped'
         flush(unit_output)
!
!        Do full space CCSD calculation
!
         wf%n_active_spaces = 1
!
         call wf%initialize_orbital_info
!
         wf%n_CCSD_o = wf%n_o
         wf%n_CCSD_v = wf%n_v
!
         wf%n_CC2_o = 0
         wf%n_CC2_v = 0
!
      endif
!
      call wf%set_n_total_active
!
!     Initialize amplitudes and associated attributes
!
!     Set number of amplitudes for CCSD active space
!
      do active_space = 1, wf%n_active_spaces
!
         wf%n_t2am = wf%n_t2am &
                  + ((wf%n_CCSD_v(active_space,1))*(wf%n_CCSD_o(active_space,1)))&
                   *((wf%n_CCSD_v(active_space,1) )*(wf%n_CCSD_o(active_space,1))+1)/2
! 
      enddo
!
      call wf%initialize_amplitudes
      call wf%initialize_omega
!
!     Set the number of parameters in the wavefunction
!     (that are solved for in the ground and excited state solvers) 
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky_for_CC2_amplitude
!
      call wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call wf%initialize_fock_matrix
!
   end subroutine init_mlccsd
!
!
   subroutine initialize_orbital_info_mlccsd(wf)
!!
!!    Initialize orbital information
!!    Written by Sarai D. Folkestad, June 2017
!!
      implicit none
!
      class(mlccsd) :: wf
!
!     :: CC2 ::
!
      if (wf%mlcc_settings%CCS .and. wf%mlcc_settings%CC2) then ! CCS is inactive method
!
         if (.not. allocated(wf%n_CC2_o)) call allocator_int(wf%n_CC2_o, wf%n_active_spaces, 1)
         if (.not. allocated(wf%n_CC2_v)) call allocator_int(wf%n_CC2_v, wf%n_active_spaces, 1)
!
         if (.not. allocated(wf%first_CC2_o)) call allocator_int(wf%first_CC2_o, wf%n_active_spaces, 1)
         if (.not. allocated(wf%first_CC2_v)) call allocator_int(wf%first_CC2_v, wf%n_active_spaces, 1)
!
      elseif(.not. wf%mlcc_settings%CCS) then ! CC2 is inactive method
!
         if (.not. allocated(wf%n_CC2_o)) call allocator_int(wf%n_CC2_o, 1, 1)
         if (.not. allocated(wf%n_CC2_v)) call allocator_int(wf%n_CC2_v, 1, 1)
!
         if (.not. allocated(wf%first_CC2_o)) call allocator_int(wf%first_CC2_o, 1, 1)
         if (.not. allocated(wf%first_CC2_v)) call allocator_int(wf%first_CC2_v, 1, 1)
!
      endif
!      
      wf%n_CC2_o  = 0
      wf%n_CC2_v  = 0
!
      wf%first_CC2_o  = 0
      wf%first_CC2_v  = 0
!
!     :: CCSD ::
!
      if (.not. allocated(wf%first_CCSD_o)) call allocator_int(wf%first_CCSD_o, wf%n_active_spaces, 1)
      if (.not. allocated(wf%first_CCSD_v)) call allocator_int(wf%first_CCSD_v, wf%n_active_spaces, 1) 
!       
      if (.not. allocated(wf%n_CCSD_o)) call allocator_int(wf%n_CCSD_o, wf%n_active_spaces, 1)
      if (.not. allocated(wf%n_CCSD_v)) call allocator_int(wf%n_CCSD_v, wf%n_active_spaces, 1) 
!
      wf%n_CCSD_o = 0
      wf%n_CCSD_v = 0
!
      wf%first_CCSD_o = 0
      wf%first_CCSD_v = 0
!
   end subroutine initialize_orbital_info_mlccsd
!
!
   subroutine initialize_amplitudes_mlccsd(wf)
!!
!!     Initialize Amplitudes (MLCCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Allocates the amplitudes, sets them to zero, and calculates
!!     the number of amplitudes.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      integer(i15) :: active_space
!
!     Allocate the doubles amplitudes and set to zero
!
      wf%n_t1am = (wf%n_o)*(wf%n_v)
      if (.not. allocated(wf%t1am)) call allocator(wf%t1am, wf%n_t1am, 1)
      wf%t1am = zero
!
      wf%n_t2am = 0
!
      do active_space = 1, wf%n_active_spaces
         wf%n_t2am = wf%n_t2am + (wf%n_CCSD_v(active_space,1))*(wf%n_CCSD_o(active_space,1))&
                  *((wf%n_CCSD_v(active_space,1))*(wf%n_CCSD_o(active_space,1))+1)/2
      enddo
!
      if (.not. allocated(wf%t2am)) call allocator(wf%t2am, wf%n_t2am, 1)
      wf%t2am = zero
!
   end subroutine initialize_amplitudes_mlccsd
!
   subroutine construct_perturbative_doubles_ccsd(wf)
!!
!!    Construct Perturbative Doubles (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Sets the doubles amplitudes (t2am) to its MP2 estimate. This is
!!    the initial guess used in the solver for the ground state amplitude 
!!    equations.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_jb
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0 
!
      integer(i15) :: offset
!
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
      integer(i15) :: active_space ! first active virtual index
!
      offset = 0
      do active_space = 1, wf%n_active_spaces
!
!        Calculate first/last indeces
! 
         call wf%get_CCSD_active_indices(first_active_o, first_active_v, active_space)
!
         n_active_o = wf%n_CCSD_o(active_space, 1)
         n_active_v = wf%n_CCSD_v(active_space, 1)
!
         last_active_o = first_active_o + n_active_o - 1
         last_active_v = first_active_v + n_active_v - 1 
!
!        Allocate L_ia_J and g_ia_jb
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         call allocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         L_ia_J = zero
         g_ia_jb = zero
!
!        Get the Cholesky IA vector 
!
         call wf%get_cholesky_ia(L_ia_J, first_active_o, last_active_o, first_active_v, last_active_v)
!
!        Calculate g_ia_jb = g_iajb
!
         call dgemm('N','T',            &
                     n_active_o*n_active_v, & 
                     n_active_o*n_active_v, &
                     wf%n_J,            &
                     one,               &
                     L_ia_J,            &
                     n_active_o*n_active_v, &
                     L_ia_J,            &
                     n_active_o*n_active_v, &
                     zero,              &
                     g_ia_jb,           &
                     n_active_o*n_active_v)
!
!        Set the doubles amplitudes
!
         do i = 1, n_active_o
            do a = 1, n_active_v
!
               ai = index_two(a, i, n_active_v)
               ia = index_two(i, a, n_active_o)
!
               do j = 1, n_active_o
                  do b = 1, n_active_v
!        
                     jb = index_two(j, b, n_active_o)
                     bj = index_two(b, j, n_active_v)
!
!                    Set the doubles amplitudes
!
                     if (ai .le. bj) then ! To avoid setting the same element twice
!
                        aibj = index_packed(ai,bj)
!
                        wf%t2am(aibj + offset, 1) = - g_ia_jb(ia,jb)/(wf%fock_diagonal(wf%n_o + a + first_active_v - 1, 1) + &
                                                               wf%fock_diagonal(wf%n_o + b + first_active_v - 1, 1) - &
                                                               wf%fock_diagonal(i + first_active_o - 1, 1) - &
                                                               wf%fock_diagonal(j + first_active_o - 1, 1))
!
                     endif
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocations
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), (wf%n_J))
         call deallocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
!
!        Update offset
!
         offset = offset + (n_active_v)*(n_active_o)*((n_active_v)*(n_active_o)+1)/2
      enddo
!
   end subroutine construct_perturbative_doubles_ccsd
!
   subroutine get_CCSD_active_indices_mlccsd(wf, first_o, first_v, active_space)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
   implicit none
!
   class(mlccsd) :: wf
   integer(i15) :: first_o
   integer(i15) :: first_v
   integer(i15) :: active_space
! 
   first_o = wf%first_CCSD_o(active_space, 1)
   first_v = wf%first_CCSD_v(active_space, 1)
!
   end subroutine get_CCSD_active_indices_mlccsd
!
   subroutine set_n_total_active_mlccsd(wf)
!!
!!
      implicit none
!  
      class(mlccsd) :: wf
!
      integer(i15) :: active_space
!
      wf%n_total_active_o = 0
      wf%n_total_active_v = 0
!
      if (wf%mlcc_settings%CCS) then !CCS is inactive method
!
         do active_space = 1, wf%n_active_spaces
            wf%n_total_active_o = wf%n_total_active_o + wf%n_CC2_o(active_space, 1) + wf%n_CCSD_o(active_space, 1)
            wf%n_total_active_v = wf%n_total_active_v + wf%n_CC2_v(active_space, 1) + wf%n_CCSD_v(active_space, 1)
         enddo
!
      else ! CC2 is inactive method
!
          wf%n_total_active_o = wf%n_o
          wf%n_total_active_v = wf%n_v
!
      endif 
!
   end subroutine set_n_total_active_mlccsd
!
!
   function active_space_t2am_offset_mlccsd(wf, active_space)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      integer(i15) :: active_space
!
      integer(i15) :: active_space_t2am_offset_mlccsd
!
      integer(i15) :: i = 0
!
      active_space_t2am_offset_mlccsd = 0
!
      do i = 1, active_space - 1
!
         active_space_t2am_offset_mlccsd = active_space_t2am_offset_mlccsd &
                                       + ((wf%n_CCSD_v(active_space,1))*(wf%n_CCSD_o(active_space,1)))&
                                       *((wf%n_CCSD_v(active_space,1) )*(wf%n_CCSD_o(active_space,1))+1)/2
!
      enddo
!
   end function active_space_t2am_offset_mlccsd
!
!
end module mlccsd_class