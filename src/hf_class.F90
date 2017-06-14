module hf_class
!
!!
!!                      Hartree-Fock (HF) class module                                 
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
   use calc_procedures_class
   use calc_settings_class
!
   implicit none
!
!  ::::::::::::::::::::::::::::::::::::
!  -::- Definition of the HF class -::-
!  ::::::::::::::::::::::::::::::::::::
!   
   type :: hf
!
!     Model name 
!
      character(len=40) :: name = 'HF'
!
!     Orbital information variables
!
      integer(i15) :: n_o  ! Number of occupied orbitals
      integer(i15) :: n_v  ! Number of virtual orbitals
      integer(i15) :: n_ao ! Number of atomic orbitals (AOs)
      integer(i15) :: n_mo ! Number of molecular orbitals (MOs)
      integer(i15) :: n_J  ! Number of Cholesky vectors
!
      real(dp), dimension(:,:), allocatable :: mo_coef ! MO coefficient matrix
!
!     Fock matrix variables
!
      real(dp), dimension(:,:), allocatable :: fock_diagonal  ! diagonal vector
!
!     Energy variables
!
      real(dp) :: energy            ! Same as scf_energy for HF class, different for descendants
!
      real(dp) :: nuclear_potential ! Nuclear potential energy term
      real(dp) :: scf_energy        ! The Hartree-Fock (HF/SCF) energy
!
!     Calculation settings, tasks, and implemented methods
!
      type(calc_settings)   :: settings
!
      type(calc_procedures) :: tasks
      type(calc_procedures) :: implemented
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_hf
      procedure :: drv  => drv_hf
!
!     Routines to read MO Cholesky vectors from file
!
      procedure, non_overridable :: read_cholesky_ij => read_cholesky_ij_hf ! occ-occ
      procedure, non_overridable :: read_cholesky_ia => read_cholesky_ia_hf ! occ-vir
      procedure, non_overridable :: read_cholesky_ai => read_cholesky_ai_hf ! vir-occ
      procedure, non_overridable :: read_cholesky_ab => read_cholesky_ab_hf ! vir-vir
!
!     Routines needed to initialize HF     
!
!     read_info               : sets variables from file (n_o, n_v, scf_energy,...)
!     read_transform_cholesky : reads AO Cholesky vectors, transforms to MO, and saves to file
!
      procedure, non_overridable :: read_hf_info            => read_hf_info_hf
      procedure, non_overridable :: read_transform_cholesky => read_transform_cholesky_hf 
!
   end type hf
!
!
contains 
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::-
!  ::::::::::::::::::::::::::::::::::::::::::::
! 
   subroutine init_hf(wf)
!!
!!    Initialization of Hartree-Fock object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Performs the following tasks:
!!
!!    1. Sets HF orbital and energy information by reading from file
!!    2. Transforms AO Cholesky vectors to MO basis and saves to file
!!
      implicit none
!
      class(hf) :: wf
!
!     Initialize HF variables
!
      call wf%read_hf_info        
!
!     Initialize Cholesky vectors
!     
      call wf%read_transform_cholesky
!
   end subroutine init_hf
!
!
   subroutine drv_hf(wf)
!!
!!    Driver (HF)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Lets the user know there is no driver for Hartree-Fock and exits
!!    the program if called. The module reads Hartree-Fock information 
!!    from files and contains no independent solver.
!!
      implicit none 
!
      class(hf) :: wf
!
      write(unit_output,*) 'Error: There is no driver for the Hartree-Fock class.'
      stop
!
   end subroutine drv_hf
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::-
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine read_hf_info_hf(wf)
!!
!!    Read HF Info
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the file mlcc_hf_info and sets the following HF variables: 
!!    n_o, n_v, n_mo, orbital_coef, and the fock_diagonal.
!!
!!    The file mlcc_hf_info is written in the mlcc_write_sirifc 
!!    subroutine, which is called from the wr_sirifc subroutine in
!!    the siropt module of the DALTON suite.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: unit_hf = -1 ! Unit identifier for mlcc_hf_info file
!
      integer(i15) :: n_lambda = 0 ! n_ao * n_mo, read but discarded
      integer(i15) :: i = 0, j = 0
!
!     Open the file mlcc_hf_info
!
      call generate_unit_identifier(unit_hf)
      open(unit=unit_hf, file='mlcc_hf_info', status='old', form='formatted')
      rewind(unit_hf)
!
!     Read mlcc_hf_info into HF variables
!  
      read(unit_hf,*) wf%n_mo, wf%n_o, n_lambda, &
                        wf%nuclear_potential, wf%scf_energy
!
!     Set the energy equal to the read SCF energy
!
      wf%energy = wf%scf_energy
!
!     Calculate the number of virtuals
!
      wf%n_v = wf%n_mo - wf%n_o
!      
!     Allocate the Fock diagonal and the MO coefficients
!
      call allocator(wf%fock_diagonal, wf%n_mo, 1)
      wf%fock_diagonal = zero
!
      call allocator(wf%mo_coef, n_lambda, 1)
      wf%mo_coef = zero
!
!     Read in the Fock diagonal and MO coefficients
!
      read(unit_hf,*) (wf%fock_diagonal(i,1), i = 1, wf%n_mo)
      read(unit_hf,*) (wf%mo_coef(i,1), i = 1, n_lambda) 
!
!     Close the mlcc_hf_info file
!    
      close(unit_hf)
!
   end subroutine read_hf_info_hf
!
!
   subroutine read_transform_cholesky_hf(wf)
!!
!!    Read and Transform Cholesky
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 20 Apr 2017
!!
!!    Reads the AO Cholesky vectors from file, transforms the vectors 
!!    to the MO basis, and saves the MO vectors to file
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: unit_chol_ao    = -1 ! Unit identifier for mlcc_cholesky file
      integer(i15) :: unit_chol_mo_ij = -1 ! cholesky_ij file
      integer(i15) :: unit_chol_mo_ia = -1 ! cholesky_ia file
      integer(i15) :: unit_chol_mo_ab = -1 ! cholesky_ab file
      integer(i15) :: unit_chol_mo_ij_direct = -1 ! cholesky_ij direct access file
      integer(i15) :: unit_chol_mo_ia_direct = -1 ! cholesky_ia direct access file
      integer(i15) :: unit_chol_mo_ab_direct = -1 ! cholesky_ab direct access file
      integer(i15) :: ioerror = 0
      integer(i15) :: throw_away_index = 0
      real(dp)     :: throw_away
!
      integer(i15) :: n_ao_sq_packed = 0 ! Packed dimensionality of (n_ao x n_ao) matrix
!
      real(dp), dimension(:,:), allocatable :: chol_ao    ! Packed AO Cholesky vector
      real(dp), dimension(:,:), allocatable :: chol_ao_sq ! Unpacked AO Cholesky vector
      real(dp), dimension(:,:), allocatable :: chol_mo_sq ! Unpacked MO Cholesky vector
!
      real(dp), dimension(:,:), allocatable :: X ! An intermediate matrix
      real(dp), dimension(:,:), allocatable :: L_ij_J
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: L_ab_J
!
      integer(i15) :: i,j,a,b,k,ij,ia, ab
!
!     Batching variables
!
      integer(i15) :: b_batch = 0, b_first = 0, b_last = 0, b_length = 0
      integer(i15) :: required = 0, available = 0, n_batch = 0, batch_dimension = 0
      integer(i15) :: max_batch_length = 0
!
!     Open Dalton file mlcc_cholesky (see mlcc_write_cholesky.F)
! 
      call generate_unit_identifier(unit_chol_ao)
      open(unit=unit_chol_ao, file='mlcc_cholesky', status='old', form='formatted')
      rewind(unit_chol_ao)
!
!     Read the number of Cholesky vectors (n_J) and 
!     the number of atomic orbitals (n_ao)
!
      read(unit_chol_ao,*) wf%n_ao, wf%n_J
!
!     Open files for MO Cholesky vectors 
! 
      call generate_unit_identifier(unit_chol_mo_ij)
      call generate_unit_identifier(unit_chol_mo_ia)
      call generate_unit_identifier(unit_chol_mo_ab)
!
      open(unit_chol_mo_ij, file='cholesky_ij', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ij)
!
      open(unit_chol_mo_ia, file='cholesky_ia', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ia)
!
      open(unit_chol_mo_ab, file='cholesky_ab', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Allocate packed and unpacked Cholesky AO, and 
!     unpacked Cholesky MO vectors
!
      n_ao_sq_packed = packed_size(wf%n_ao)
!
      call allocator(chol_ao, n_ao_sq_packed, 1)
      call allocator(chol_ao_sq, wf%n_ao, wf%n_ao) 
      call allocator(chol_mo_sq, wf%n_mo, wf%n_mo)
!
      chol_ao    = zero
      chol_ao_sq = zero
      chol_mo_sq = zero
!
!     Allocate an intermediate, X
!
      call allocator(X, wf%n_ao, wf%n_mo)
!
      X = zero
!
!     Loop over the number of Cholesky vectors,
!     reading them one by one 
!
      do j = 1, wf%n_J
!
!        Read Cholesky AO vector
!
         read(unit_chol_ao,*) (chol_ao(i,1), i = 1, n_ao_sq_packed)
!
!        Unpack/square up AO vector 
!
         call squareup(chol_ao, chol_ao_sq, wf%n_ao)
!
!        Transform the AO vectors to form the Cholesky MO vectors
!
         call dgemm('N','N',     &
                     wf%n_ao,    &
                     wf%n_mo,    &
                     wf%n_ao,    &
                     one,        &
                     chol_ao_sq, &
                     wf%n_ao,    &
                     wf%mo_coef, &
                     wf%n_ao,    &
                     zero,       &
                     X,          &
                     wf%n_ao)
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
                     chol_mo_sq, &
                     wf%n_mo)
!
!        Write the MO vectors to files in blocks
!
         write(unit_chol_mo_ij) ((chol_mo_sq(i,k), i = 1, wf%n_o), k = 1, i)
         write(unit_chol_mo_ia) ((chol_mo_sq(i,a), i = 1, wf%n_o), a = wf%n_o + 1, wf%n_mo)
         write(unit_chol_mo_ab) ((chol_mo_sq(a,b), a = wf%n_o + 1, wf%n_mo), b = wf%n_o + 1, a)
!
      enddo
!
!     Close files 
!
      close(unit_chol_ao)
      close(unit_chol_mo_ij)
      close(unit_chol_mo_ia)
      close(unit_chol_mo_ab)
! 
!     Rewrite to direct access file, delete sequential file
! 
!     :: L_ij_J ::
!
      call generate_unit_identifier(unit_chol_mo_ij)
      open(unit_chol_mo_ij, file='cholesky_ij', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ij)
!
!     Read L_ij_J
!
      call allocator(L_ij_J, wf%n_o*(wf%n_o+1)/2, wf%n_J)
      do J = 1, wf%n_J
         read(unit_chol_mo_ij) (L_ij_J(ij, J), ij = 1, wf%n_o*(wf%n_o+1)/2)
      enddo
!
!     Close and delete file
!
      close(unit_chol_mo_ij, status='delete')
!
      call generate_unit_identifier(unit_chol_mo_ij_direct)
      open(unit=unit_chol_mo_ij_direct, file='cholesky_ij_direct', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
      do ij = 1, wf%n_o*(wf%n_o+1)/2
         write(unit_chol_mo_ij_direct, rec=ij) (L_ij_J(ij,j), j = 1, wf%n_J)
      enddo
!
      call deallocator(L_ij_J,  wf%n_o*(wf%n_o+1)/2, wf%n_J)
      close(unit_chol_mo_ij_direct)
!
!     :: L_ia_J :: 
!
      call generate_unit_identifier(unit_chol_mo_ia)
      open(unit_chol_mo_ia, file='cholesky_ia', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ia)
!
!     Read L_ia_J
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      do J = 1, wf%n_J
         read(unit_chol_mo_ia) (L_ia_J(ia, J), ia = 1, (wf%n_o)*(wf%n_v))
      enddo
!
!     Close and delete file
!
      close(unit_chol_mo_ia, status='delete')
!
      call generate_unit_identifier(unit_chol_mo_ia_direct)
      open(unit=unit_chol_mo_ia_direct, file='cholesky_ia_direct', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
      do ia = 1, (wf%n_o)*(wf%n_v)
         write(unit_chol_mo_ia_direct, rec=ia) (L_ia_J(ia,j), j = 1, wf%n_J)
      enddo
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      close(unit_chol_mo_ia_direct)
!
!     :: L_ab_J ::
!
      call generate_unit_identifier(unit_chol_mo_ab)
      open(unit_chol_mo_ab, file='cholesky_ab', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ab)
!
      call generate_unit_identifier(unit_chol_mo_ia_direct)
      open(unit=unit_chol_mo_ia_direct, file='cholesky_ia_direct', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
!     Read L_ab_J in batches over b
!
      required = ((wf%n_v)**2)*(wf%n_J)
!
      required = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index b
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do b_batch = 1, n_batch
!
         call batch_limits(b_first, b_last, b_batch, max_batch_length, batch_dimension)
         b_length = b_last - b_first + 1 
!
         call allocator(L_ab_J, (((b_length + 1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length), wf%n_J)
!
         if (b_first .ne. 1) then
!  
!           Calculate index of last element to throw away
!  
            throw_away_index = index_packed(wf%n_v, b_first - 1)
!  
!           Throw away all elements from 1 to throw_away_index, then read from batch start
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab) (throw_away, i = 1, throw_away_index), &
                                    (L_ab_J(a,j), a = 1,(((b_length + 1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length))
!
            enddo
!
         else
!  
!           Read from the start of each entry
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab) (L_ab_J(a,j), a = 1, (((b_length + 1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length))
!
            enddo
!
         endif
!
         do a = 1, wf%n_v
            do b = b_first, b_last
               ab = index_packed(a, b)
               write(unit_chol_mo_ij_direct, rec=ab) (L_ab_J(ab, J), J = 1, wf%n_J)
            enddo
         enddo
         call deallocator(L_ab_J, (((b_length+1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length), wf%n_J)
!
      enddo
      close(unit_chol_mo_ab, status='delete')
      close(unit_chol_mo_ab_direct)
!
   end subroutine read_transform_cholesky_hf
!
   subroutine read_cholesky_ij_hf(wf,L_ij_J)
!!
!!    Read Cholesky IJ 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky IJ (occ-occ) vectors from file and 
!!    places them in the incoming L_ij_J matrix
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_o)**2, wf%n_J) :: L_ij_J ! L_ij^J
!
      integer(i15) :: unit_chol_mo_ij = -1 ! Unit identifier for cholesky_ij file 
      integer(i15) :: i = 0, j = 0
!
!     Prepare for reading: generate unit idientifier, open file, and rewind
!
      call generate_unit_identifier(unit_chol_mo_ij)
      open(unit=unit_chol_mo_ij, file='cholesky_ij', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ij)
!
!     Read the Cholesky vectors into the L_ij_J matrix
!
      do j = 1, wf%n_J
         read(unit_chol_mo_ij) (L_ij_J(i,j), i = 1, (wf%n_o)**2)
      enddo
!
!     Close file
!
      close(unit_chol_mo_ij)    
!   
   end subroutine read_cholesky_ij_hf
!
!
   subroutine read_cholesky_ia_hf(wf,L_ia_J)
!!
!!    Read Cholesky IA 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky IA (occ-vir) vectors from file and
!!    places them in the incoming L_ia_J matrix
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ia_J ! L_ia^J
!
      integer(i15) :: unit_chol_mo_ia = -1 ! Unit identifier for cholesky_ia file
      integer(i15) :: i = 0, j = 0
!
!     Prepare for reading: generate unit idientifier, open, and rewind file
!
      call generate_unit_identifier(unit_chol_mo_ia)
      open(unit=unit_chol_mo_ia, file='cholesky_ia', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ia)
!
!     Read Cholesky vectors into the L_ia_J matrix
!
      do j = 1, wf%n_J
!
         read(unit_chol_mo_ia) (L_ia_J(i,j), i = 1, (wf%n_o)*(wf%n_v))
!
      enddo
!
!     Close file
!
      close(unit_chol_mo_ia)    
!   
   end subroutine read_cholesky_ia_hf
!
!
   subroutine read_cholesky_ai_hf(wf, L_ai_J)
!!
!!    Read Cholesky AI 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky AI (vir-occ) vectors from file and
!!    places them in the incoming L_ai_J matrix
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), wf%n_J) :: L_ai_J ! L_ai^J
!
      real(dp), dimension(:,:), allocatable :: L_ia_J       
!
      integer(i15) :: i = 0, j = 0, a = 0, ia = 0, ai = 0
!
!     Allocation
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ia_J = zero
!     
!     Get Cholesky IA vector 
!
      call wf%read_cholesky_ia(L_ia_J)
!
!     Reorder and save in AI vector 
!      
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
!           Needed indices
!
            ai = index_two(a, i, wf%n_v)
            ia = index_two(i, a, wf%n_o)
!
            do j = 1, wf%n_J
!
               L_ai_J(ai, j) = L_ia_J(ia, j)
!
            enddo
!
         enddo
      enddo
!
!     Deallocate temporary vector 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)   
!
   end subroutine read_cholesky_ai_hf
!   
!
   subroutine read_cholesky_ab_hf(wf, L_ab_J, first, last, ab_dim, reorder)
!!
!!    Read Cholesky AB 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky AB (vir-vir) vectors from file and
!!    places them in the incoming L_ab_J matrix, with batching 
!!    if necessary
!!
!!    If reorder = .true.,  L_ba_J is returned with batching over a
!!    If reorder = .false., L_ab_J is returned with batching over b
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15), intent(in) :: first   ! First index (can differ from 1 when batching)
      integer(i15), intent(in) :: last    ! Last index  (can differ from n_v when batching)
      integer(i15), intent(in) :: ab_dim  ! Dimension of ab index (not n_v^2 when batching) 
!     
      logical, intent(in)      :: reorder ! See description above
!
      real(dp), dimension(ab_dim, wf%n_J) :: L_ab_J ! L_ab^J
!
      integer(i15) :: unit_chol_mo_ab = -1 ! Unit identifier for cholesky_ab file
!
      integer(i15) :: a = 0, b = 0, j = 0, i = 0
      integer(i15) :: batch_length = 0
!
      integer(i15) :: throw_away_index = 0
      real(dp)     :: throw_away = 0
!
!     Prepare for reading: generate unit identifier, open, and rewind file
!  
      call generate_unit_identifier(unit_chol_mo_ab)
      open(unit=unit_chol_mo_ab, file='cholesky_ab', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Calculating batch length
!
      batch_length = last - first + 1
!
      if (.not. reorder) then
!
         if (first .ne. 1) then
!  
!           Calculate index of last element to throw away
!  
            throw_away_index = index_two(wf%n_v, first - 1, wf%n_v)
!  
!           Throw away all elements from 1 to throw_away_index, then read from batch start
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab) (throw_away, i = 1, throw_away_index),&
                                    (L_ab_J(a,j), a = 1, ab_dim)
!
            enddo
!
         else
!  
!           Read from the start of each entry
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab) (L_ab_J(a,j), a = 1, ab_dim)
!
            enddo
!
         endif
!
      else ! Reorder L_ab_J is L_ba_J
!
         if (first .ne. 1) then
            throw_away_index = index_two(wf%n_v, first - 1, wf%n_v)
         else
            throw_away_index = 0
         endif
!
!        Reading vectors
!
         do j = 1, wf%n_J
!
            if (first .eq. 1) then
! 
               read(unit_chol_mo_ab) ((L_ab_J(index_two(b, a, wf%n_v), j), b = 1, wf%n_v), a = 1, batch_length)
!
            else
!
               read(unit_chol_mo_ab) (throw_away, i = 1, throw_away_index),&
                                     ((L_ab_J(index_two(b, a, wf%n_v), j), b = 1, wf%n_v), a = 1, batch_length)
!
            endif
!
         enddo
!     
      endif  ! Reorder
! 
!     Close file
!    
      close(unit_chol_mo_ab)
!
   end subroutine read_cholesky_ab_hf
!
!
end module hf_class
