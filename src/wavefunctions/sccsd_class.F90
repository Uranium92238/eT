module sccsd_class
!
!!
!!       Similarity constrained CCSD (SCCSD) class module                                 
!!            Written by Eirik F. Kjønstad, June 2017         
!!                                                                           
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- 1. Modules used by the class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  Ancestor class module (CCSD)
!
   use ccsd_class
!
   implicit none 
!
!  ::::::::::::::::::::::::::::::::::::::::::
!  -::- 2. Definition of the SCCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::::::
!
   type, extends(ccsd) :: sccsd
!
!     The triples amplitude that enforces similarity 
!
      real(dp) :: triples = zero ! t = t_IJK^ABC
!
!     The triples excitation indices corresponding to 
!     that amplitude (tau_IJK^ABC)
!
      integer(i15) :: I = 0
      integer(i15) :: J = 0
      integer(i15) :: K = 0
!
      integer(i15) :: A = 0
      integer(i15) :: B = 0
      integer(i15) :: C = 0
!
!     The excited states that are constrained to be similar 
!
      integer(i15) :: state_A = 0
      integer(i15) :: state_B = 0
!
!     The generalized overlap between the two constrained states 
!
      real(dp) :: overlap = zero
      real(dp) :: overlap_threshold = 1.0D-6
!
!     Logical that determines whether to constrain the left or right 
!     eigenstates (note: there is not a rigorous left/right equivalency)
!
      logical :: constrain_right_eigenvectors = .true.
!
   contains
!
!
!     -::- Initialization routines -::-
!     ---------------------------------
!
      procedure :: init => init_sccsd
!
!
!     -::- Omega submodule routine pointers -::-
!     ------------------------------------------
!
      procedure :: construct_omega => construct_omega_sccsd
!
      procedure :: omega_sccsd_a1 => omega_sccsd_a1_sccsd
!
!
!     -::- Jacobian submodule routine pointers -::-
!     ---------------------------------------------
!
      procedure :: jacobian_ccsd_transformation => jacobian_ccsd_transformation_sccsd 
!
      procedure :: jacobian_sccsd_a2 => jacobian_sccsd_a2_sccsd
      procedure :: jacobian_sccsd_b2 => jacobian_sccsd_b2_sccsd 
      procedure :: jacobian_sccsd_c2 => jacobian_sccsd_c2_sccsd
!
!
!     -::- Jacobian transpose submodule routine pointers -::-
!     -------------------------------------------------------
!
      procedure :: jacobian_transpose_ccsd_transformation => jacobian_transpose_ccsd_transformation_sccsd
!
      procedure :: jacobian_transpose_sccsd_a1 => jacobian_transpose_sccsd_a1_sccsd
      procedure :: jacobian_transpose_sccsd_b1 => jacobian_transpose_sccsd_b1_sccsd 
      procedure :: jacobian_transpose_sccsd_c1 => jacobian_transpose_sccsd_c1_sccsd
!
!
!     -::- Metric submodule routine pointers -::-
!     -------------------------------------------
!
      procedure :: metric_transformation => metric_transformation_sccsd
!
!     Helper routines for the metric transformation 
!
      procedure :: construct_q      => construct_q_sccsd       ! q_mu 
      procedure :: Q_transformation => Q_transformation_sccsd  ! Q_mu,nu
      procedure :: S_transformation => S_transformation_sccsd  ! S_mu,nu  
!
!     Routine to calculate the overlap between the similarity constrained 
!     states (zero if equations are satisfied)
!
      procedure :: calc_overlap      => calc_overlap_sccsd      ! Between the right eigenstates
      procedure :: calc_left_overlap => calc_left_overlap_sccsd ! Between the left eigenstates
!
!
!     -::- Excited state submodule routine pointers -::-
!     --------------------------------------------------
!
      procedure :: excited_state_driver => excited_state_driver_sccsd
!
!
!     -::- Other class routine pointers not located in submodules -::-
!     ----------------------------------------------------------------
!
      procedure :: read_triples => read_triples_sccsd
!
   end type sccsd 
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- 3. Interfaces to the submodules of the SCCSD class -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   interface 
!
!
!     -::- Omega submodule interface -::-
!     -----------------------------------
!
      module subroutine construct_omega_sccsd(wf)
!!
!!       Construct Omega (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf
!
      end subroutine construct_omega_sccsd
!
!
      module subroutine omega_sccsd_a1_sccsd(wf)
!!
!!       Omega SCCSD A1 
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine omega_sccsd_a1_sccsd
!
!
   end interface
!
!
   interface 
!
!     -::- Jacobian submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::
!
      module subroutine jacobian_ccsd_transformation_sccsd(wf, c_a_i, c_aibj)
!!
!!       Jacobian transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, May 2017
!!
         implicit none 
!
         class(sccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj ! c_aibj 
!
      end subroutine jacobian_ccsd_transformation_sccsd
!
!
      module subroutine jacobian_sccsd_a2_sccsd(wf, rho, X)
!!
!!       Jacobian SCCSD A2 
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: X                         ! X(l,d) = X_ld 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho ! rho(ai,bj) = rho_aibj 
!
      end subroutine jacobian_sccsd_a2_sccsd
!
!
      module subroutine jacobian_sccsd_b2_sccsd(wf, rho, Y)
!!
!!       Jacobian SCCSD B2
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: Y ! Y(kc,lj) = Y_kclj 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho ! rho(ai,bj) = rho_aibj
!
      end subroutine jacobian_sccsd_b2_sccsd
!
!
      module subroutine jacobian_sccsd_c2_sccsd(wf, rho, Z)
!!
!!       Jacobian SCCSD C2
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)*(wf%n_v)) :: Z ! Z(ac,ld) = Z_acld 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho ! rho(ai,bj) = rho_aibj
!
      end subroutine jacobian_sccsd_c2_sccsd
!
!
   end interface 
!
!
   interface 
!
!     -::- Jacobian transpose submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine jacobian_transpose_ccsd_transformation_sccsd(wf, b_a_i, b_aibj)
!!
!!       Jacobian transpose CCSD transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
         real(dp), dimension(wf%n_t2am, 1)   :: b_aibj 
!
      end subroutine jacobian_transpose_ccsd_transformation_sccsd
!
!
      module subroutine jacobian_transpose_sccsd_a1_sccsd(wf, sigma, b2am, L)
!!
!!       Jacobian transpose SCCSD A1
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      !  sigma(a, i) = sigma_ai 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: L    !    L(ia, jb) = t * L_iajb 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am ! b2am(ai, bj) = b_aibj 
!
      end subroutine jacobian_transpose_sccsd_a1_sccsd
!
!
      module subroutine jacobian_transpose_sccsd_b1_sccsd(wf, sigma, b2am, g)
!!
!!       Jacobian transpose SCCSD B1
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      !  sigma(a, i) = sigma_ai 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: g    !    g(ia, jb) = t * g_iajb 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am ! b2am(ai, bj) = b_aibj 
!
      end subroutine jacobian_transpose_sccsd_b1_sccsd 
!
!
      module subroutine jacobian_transpose_sccsd_c1_sccsd(wf, sigma, b2am, g)
!!
!!       Jacobian transpose SCCSD C1
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      !  sigma(a, i) = sigma_ai 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: g    !    g(ia, jb) = t * g_iajb 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am ! b2am(ai, bj) = b_aibj 
!
      end subroutine jacobian_transpose_sccsd_c1_sccsd  
!
!  
   end interface
!
!
   interface 
!
!     -::- Metric submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::
!
      module subroutine calc_overlap_sccsd(wf)
!!
!!       Calculate Overlap (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine calc_overlap_sccsd
!
!
      module subroutine calc_left_overlap_sccsd(wf)
!!
!!       Calculate Left Overlap (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine calc_left_overlap_sccsd
!
!
      module subroutine metric_transformation_sccsd(wf, r_a_i, r_aibj)
!!
!!       Metric Transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i  ! r_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: r_aibj ! r_aibj 
!
      end subroutine metric_transformation_sccsd
!
!
      module subroutine construct_q_sccsd(wf, q_a_i, q_aibj)
!!
!!       Construct q (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: q_a_i  ! q_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: q_aibj ! q_aibj 
!
      end subroutine construct_q_sccsd
!
!
      module subroutine Q_transformation_sccsd(wf, r_a_i, r_aibj, transpose)
!!
!!       Q transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i  ! r_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: r_aibj ! r_aibj
!
         logical :: transpose ! If true, transform by Q^T instead of Q  
!
      end subroutine Q_transformation_sccsd    
!
!
      module subroutine S_transformation_sccsd(wf, r_a_i, r_aibj)
!!
!!       S Transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         class(sccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i 
         real(dp), dimension(wf%n_t2am, 1)   :: r_aibj
!
      end subroutine S_transformation_sccsd
!
!
   end interface 
!
!
   interface 
!
!     -::- Excited state submodule interface -::-
!     :::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine excited_state_driver_sccsd(wf)
!!
!!       Excited state driver (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
!!       Directs the solution of the excited state problem for SCCSD. 
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine excited_state_driver_sccsd
!
!
   end interface
!
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- 4. Class subroutines and functions -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_sccsd(wf)
!!
!!    Initialize SCCSD object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Performs the following tasks:
!!
!!    - Sets HF orbital and energy information by reading from file (read_hf_info)
!!    - Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!!    - Allocates the Fock matrix and sets it to zero
!!    - Initializes the amplitudes (sets their initial values and associated variables)
!!
      implicit none 
!
      class(sccsd) :: wf
!
!     Set model name 
!
      wf%name = 'SCCSD'
!
!     Set implemented methods
!
      wf%implemented%ground_state  = .true.
      wf%implemented%excited_state = .true.
      wf%implemented%properties    = .false.
!
!     Initialize the triples amplitude 
!
      wf%triples = zero
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky 
!
!     Initialize (singles and doubles) amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
      wf%n_t2am = (wf%n_t1am)*(wf%n_t1am + 1)/2 
!
!     Read the triples amplitude (if it exists)
!
      call wf%read_triples
!
!     The number of parameters solved for in the ground state 
!     and excited state equations is, still, the number of singles
!     and doules - it is a CCSD model, and the triples parameter 
!     is only used to set the non-CCSD part of the model (i.e., the 
!     generalized orthogonality).
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      if (.not. allocated(wf%t1am)) call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
      call wf%initialize_fock_matrix
!
   end subroutine init_sccsd
!
!
   subroutine read_triples_sccsd(wf)
!!
!!    Read Triples (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      logical :: file_exists = .false.
!
      integer(i15) :: unit_triples = -1
!
!     Check to see whether file exists 
!
      inquire(file='triples',exist=file_exists)
!
      if (file_exists) then 
!
         call generate_unit_identifier(unit_triples)
         open(unit_triples, file='triples', status='unknown', form='unformatted')
         rewind(unit_triples)
         read(unit_triples) wf%triples 
!
         write(unit_output,'(/t3,a,f12.10)') &
            'Found a triple amplitude on file. Restarting with amplitude equal to ', wf%triples
!
         close(unit_triples)
!
      endif
!
   end subroutine read_triples_sccsd
!
!
end module sccsd_class