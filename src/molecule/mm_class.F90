!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module mm_class
!
!!
!!    QM/MM Molecular System class module
!!    Written by Tommaso Giovannini, March 2019
!!
!
   use array_utilities
   use parameters
   use atomic_class
   use io_utilities
   use memory_manager_class
!
   implicit none
!
   type :: mm
!
      integer :: n_atoms     ! number of MM atoms
      integer :: max_mm_mol  ! max number of MM molecules
      integer :: n_variables ! Total number of FQ variables
!      
      character(len=200) :: name    
      character(len=200) :: forcefield  ! MM forcefield
      character(len=200) :: algorithm   ! algorithm for FQ solution
!      
      real(kind=dp), dimension(:,:), allocatable :: coordinates   ! MM coordinates
      real(kind=dp), dimension(:),   allocatable :: charge        ! MM charges 
      real(kind=dp), dimension(:),   allocatable :: chi,eta       ! electronegativity and chemical hardness
      real(kind=dp), dimension(:,:), allocatable :: fq_matrix     ! FQ matrix
!      
      integer, dimension(:), allocatable :: molecule              ! Molecule number
      integer, dimension(:), allocatable :: n_atoms_per_molecule  ! Number of atoms per molecule
!      
      real(dp), dimension(:,:), allocatable :: nopol_h_wx         ! one-electron for non-polarizable QM/MM
      real(dp), dimension(:,:), allocatable :: pol_emb_fock       ! Fock for polarizable QM/MM
      real(dp), dimension(:), allocatable   :: pol_emb_rhs        ! RHS for polarizable QM/MM
      real(dp), dimension(:), allocatable   :: pol_emb_lhs        ! LHS for polarizable QM/MM
!      
      character(len=2), dimension(:), allocatable :: symbol       ! Atomic symbol
!      
   contains
!
      procedure :: prepare                  => prepare_mm
      procedure :: cleanup                  => cleanup_mm
!
      procedure, private :: read_parameters => read_parameters_mm
      procedure, private :: read_geometry   => read_geometry_mm
      procedure, private :: check_molecule  => check_molecule_mm
!
      procedure :: print_system             => print_system_mm
      procedure :: print_geometry           => print_geometry_mm
!
      procedure :: fq_matrix_create         => fq_matrix_create_mm
!
   end type mm
!
contains
!
!
   subroutine prepare_mm(molecule)
!!
!!    Prepare MM system
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none
!
      class(mm) :: molecule
!
!     Read eT.inp and create quantities for QM/MM calculation
!
      molecule%n_atoms = input%get_mm_n_atoms()
      call molecule%read_parameters()
!
      if(.not.allocated(molecule%coordinates)) call mem%alloc(molecule%coordinates, 3,molecule%n_atoms)
!
      if(.not.allocated(molecule%molecule))    allocate(molecule%molecule(molecule%n_atoms))
!
      if(.not.allocated(molecule%symbol))      allocate(molecule%symbol(molecule%n_atoms))
!      
      if(.not.allocated(molecule%charge))      call mem%alloc(molecule%charge, molecule%n_atoms)
!
      call zero_array(molecule%charge,molecule%n_atoms)
!      
      if(trim(molecule%forcefield).eq.'fq') then
!
         if(.not.allocated(molecule%chi)) call mem%alloc(molecule%chi, molecule%n_atoms)
         if(.not.allocated(molecule%eta)) call mem%alloc(molecule%eta, molecule%n_atoms)
!
      endif
!
      call molecule%read_geometry()
!
      molecule%max_mm_mol = maxval(molecule%molecule)
!
      if(.not.allocated(molecule%n_atoms_per_molecule)) allocate(molecule%n_atoms_per_molecule(molecule%max_mm_mol))
!
      call molecule%check_molecule()
!      
      molecule%n_variables = 0
      call molecule%print_system()
!
      if(trim(molecule%forcefield).eq.'fq') then
!
          molecule%n_variables = molecule%n_atoms+molecule%max_mm_mol
!
          if(.not.allocated(molecule%fq_matrix)) call mem%alloc(molecule%fq_matrix,molecule%n_variables, &
                                                                                   molecule%n_variables)
          call zero_array(molecule%fq_matrix,(molecule%n_variables)**2)
!
          call molecule%fq_matrix_create()
!
      endIf
!
   end subroutine prepare_mm
!
!
   subroutine cleanup_mm(molecule)
!!
!!    Clean Up MM Molecular System
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none
!
      class(mm) :: molecule
!
      if(allocated(molecule%coordinates)) call mem%dealloc(molecule%coordinates, 3,molecule%n_atoms)
!
      if(allocated(molecule%molecule))    deallocate(molecule%molecule)
!
      if(allocated(molecule%symbol))      deallocate(molecule%symbol)
!      
      if(allocated(molecule%charge))      call mem%dealloc(molecule%charge, molecule%n_atoms)
!      
      if(allocated(molecule%chi))         call mem%dealloc(molecule%chi, molecule%n_atoms)
!      
      if(allocated(molecule%eta))         call mem%dealloc(molecule%eta, molecule%n_atoms)
!
      if(allocated(molecule%n_atoms_per_molecule)) deallocate(molecule%n_atoms_per_molecule)
!      
      if(allocated(molecule%fq_matrix)) call mem%dealloc(molecule%fq_matrix,molecule%n_variables, &
                                                                            molecule%n_variables)
!                                                                            
      if (allocated(molecule%pol_emb_lhs))  deallocate(molecule%pol_emb_fock)
!
      if (allocated(molecule%nopol_h_wx))   deallocate(molecule%nopol_h_wx)
!      
      if (allocated(molecule%pol_emb_fock)) call mem%dealloc(molecule%pol_emb_rhs, molecule%n_variables)
!      
      if (allocated(molecule%pol_emb_rhs))  call mem%dealloc(molecule%pol_emb_lhs, molecule%n_variables)
!
!
   end subroutine cleanup_mm
!
!
   subroutine read_parameters_mm(molecule)
!!
!!    Read parameters for QMMM calculation
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none
!
      class(mm) :: molecule
!
      molecule%forcefield = 'non-polarizable' ! Standard
      molecule%algorithm  = 'mat_inversion'
!      
      call input%get_keyword_in_section('forcefield','molecular mechanics',molecule%forcefield)
!
      if (trim(molecule%forcefield).eq.'fq') then
!
          call input%get_keyword_in_section('algorithm','molecular mechanics',molecule%algorithm)
!          
          if(trim(molecule%algorithm) .ne. 'mat_inversion') &
             call output%error_msg('Algorithm '//trim(molecule%algorithm)//' not yet implemented')
!
      endIf
!
!
   end subroutine read_parameters_mm
!
!
   subroutine print_system_mm(molecule)
!!
!!    Print MM Details
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none 
!
      class(mm) :: molecule  
!
      write(output%unit, '(/t3,a)')       '- Molecular Mechanics Specifications'
!
      write(output%unit, '(/t6,a14,a)')     'Force Field:  ', trim(molecule%forcefield)
      If(trim(molecule%forcefield).ne.'non-polarizable') &
         write(output%unit, '(t6,a14,a)')      'Algorithm  :  ', trim(molecule%algorithm)
!
      write(output%unit, '(/t6,a34,i4)')     'Number of MM atoms:                ', molecule%n_atoms
      write(output%unit, '(t6,a34,i4)')      'Number of MM molecules:            ', molecule%max_mm_mol

      flush(output%unit)
!
      call molecule%print_geometry()
!
   end subroutine print_system_mm
!
!
   subroutine print_geometry_mm(molecule)
!!
!!    Print MM geometry 
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none 
!
      class(mm) :: molecule  
!
      integer :: I 
!      
      character(len=12) :: frmt0
      character(len=57) :: frmt1
      character(len=61) :: frmt2
!
      write(output%unit,'(a)') ''
!
      frmt0="(t5,68('='))"
      frmt1="(t5,'Atom',4x,'Mol',13x,'X',10x,'Y',10x,'Z',12x,'Charge')"
      frmt2="(t5,'Atom',4x,'Mol',9x,'X',10x,'Y',10x,'Z',9x,'Chi',8x,'Eta')"
!      
      write(output%unit,frmt0 )
      write(output%unit,'(t25,a)') 'MM Geometry (Å) and Parameters'
      write(output%unit,frmt0)
!
      If(trim(molecule%forcefield).eq.'non-polarizable') write(output%unit,frmt1)
      If(trim(molecule%forcefield).eq.'fq') write(output%unit,frmt2)
      write(output%unit,frmt0)
!
      do I = 1, molecule%n_atoms 
!
          If(trim(molecule%forcefield).eq.'non-polarizable') then
!
             write(output%unit, '(t6, a2,4x,i4,6x,3(f11.6),5x,f11.6)') molecule%symbol(I),          &
                                                                       molecule%molecule(I),        &
                                                                       molecule%coordinates(1,I),   &
                                                                       molecule%coordinates(2,I),   &
                                                                       molecule%coordinates(3,I),   &
                                                                       molecule%charge(I)
!
          else if(trim(molecule%forcefield).eq.'fq') then
!
            write(output%unit, '(t6, a2,4x,i4,2x,3(f11.6),2(f11.6))')  &
                  molecule%symbol(I),          &
                  molecule%molecule(I),        &
                  molecule%coordinates(1,I),   &
                  molecule%coordinates(2,I),   &
                  molecule%coordinates(3,I),   &
                  molecule%chi(I),             &
                  molecule%eta(I)
!
          endIf
!  
         flush(output%unit)
!
      enddo 
!
      write(output%unit,frmt0)
!
      if(trim(molecule%forcefield).eq.'non-polarizable') then
!
         write(output%unit,'(/ / t3,a)') '- Electrostatic Embedding'
         write(output%unit,'(/ t6,a)') 'Each atom of the MM portion is endowed with a charge which'
         write(output%unit,'(  t6,a)') 'value is a fixed external parameter.'
         write(output%unit,'(  t6,a)') 'The QM/MM electrostatic interaction energy is defined as:'
         write(output%unit,'(/ t6,a)') '   E^ele_QM/MM = sum_i q_i * V_i(P)'
         write(output%unit,'(/ t6,a)') 'where V_i(P) is the electrostatic potential due to the QM'
         write(output%unit,'(  t6,a)') 'density calculated at the position of the i-th charge q_i.'
         write(output%unit,'(/ t6,a)') 'For further details, see:'
         write(output%unit,'(  t6,a /)') 'Senn & Thiel, Angew. Chem. Int. Ed., 2009, 48, 1198−1229'
         if(input%requested_cc_calculation()) &
            write(output%unit,'(  t6,a/)') 'CC calculation: MM charges only affect MOs and Fock'
!
      else if(trim(molecule%forcefield).eq.'fq') then
!
         write(output%unit,'(/ / t3,a)') '- Polarizable Embedding : Fluctuating Charges (FQ) Force Field'
         write(output%unit,'(/ t6,a)') 'Each atom of the MM portion is endowed with a charge which'
         write(output%unit,'(  t6,a)') 'value can vary in agreement with the Electronegativity'
         write(output%unit,'(  t6,a)') 'Equalization Principle (EEP), which states that at the'
         write(output%unit,'(  t6,a)') 'equilibrium each atom has the same electronegativity.'
         write(output%unit,'(/ t6,a)') 'The force field is defined in terms of electronegativity (Chi)'
         write(output%unit,'(  t6,a)') 'and chemical hardness (Eta), which are specified for each MM'
         write(output%unit,'(  t6,a)') 'atom.'
         write(output%unit,'(  t6,a)') 'The QM/MM electrostatic interaction energy is defined as:'
         write(output%unit,'(/ t6,a)') '   E^ele_QM/MM = sum_i q_i * V_i(P)'
         write(output%unit,'(/ t6,a)') 'where V_i(P) is the electrostatic potential due to the QM'
         write(output%unit,'(  t6,a)') 'density calculated at the position of the i-th charge q_i.'
         write(output%unit,'(  t6,a)') 'The values of the charges are obtained by solving a linear'
         write(output%unit,'(  t6,a)') 'equation:'
         write(output%unit,'(/ t6,a)') '   Dq = -Chi - V(P)'
         write(output%unit,'(/ t6,a)') 'For further details, see:'
         write(output%unit,'(  t6,a /)') 'Lipparini et al. JCTC 2012, 8, 4153–4165.'
         if(input%requested_cc_calculation()) then 
            write(output%unit,'(  t6,a)') 'CC calculation: zero-order approximation'
            write(output%unit,'(  t6,a/)') 'FQ charges only affect MOs and Fock'
         endif
!
      endif
!
      flush(output%unit)
!
!
   end subroutine print_geometry_mm
!
!
   subroutine fq_matrix_create_mm(molecule)
!!
!!    Construction of FQ matrix
!!    Ohno Kernel
!!    Written by Tommaso Giovannini, March 2019
!!
      implicit none
!
      class(mm) :: molecule
!
      integer :: I,J,I2,J2
      real(kind=dp) :: distIJ, chemhard_IJ
!
!
      do I = 1, molecule%n_atoms 
!
         do J = 1, molecule%n_atoms
!
            distIJ = dsqrt((molecule%coordinates(1,I)-molecule%coordinates(1,J))**2 + &
                           (molecule%coordinates(2,I)-molecule%coordinates(2,J))**2 + &
                           (molecule%coordinates(3,I)-molecule%coordinates(3,J))**2 )
!
            distIJ = angstrom_to_bohr*distIJ
!
            chemhard_IJ  = half*(molecule%eta(I)+molecule%eta(J))
!
            molecule%fq_matrix(I,J) = chemhard_IJ/((one+(chemhard_IJ**two)*(distIJ**two))**half)
!
         enddo
!
      enddo
!
!     Add the charge contraints 
!
      J=0
!
      Do I = 1, molecule%max_mm_mol
!
          J = J + 1
!
          J2 = J + molecule%n_atoms
!
          if ( I.eq.1) then
!
             do I2 = 1, molecule%n_atoms_per_molecule(I)
!
               molecule%fq_matrix(I2,J2) = one
               molecule%fq_matrix(J2,I2) = one
!
             enddo
!
          else
!
             do I2 = (molecule%n_atoms_per_molecule(I-1)*(I-1))+1, &
                      molecule%n_atoms_per_molecule(I-1)*(I-1)+molecule%n_atoms_per_molecule(I)
!
               molecule%fq_matrix(I2,J2) = one
               molecule%fq_matrix(J2,I2) = one
!
             enddo
!
          endif
!
      enddo
!
      call output%print_matrix('debug', 'fq_matrix', &
                               molecule%fq_matrix,   &
                               molecule%n_variables, &
                               molecule%n_variables)
!
      if(trim(molecule%algorithm).eq.'mat_inversion') then
!
         call invert(molecule%fq_matrix, &
                     molecule%fq_matrix, &
                     molecule%n_variables) 
!
         call output%print_matrix('debug', 'pol_matrix_inv', &
                                  molecule%fq_matrix,        &
                                  molecule%n_variables,      &
                                  molecule%n_variables)
!
      endif
!
!
   end subroutine fq_matrix_create_mm
!
!
   subroutine check_molecule_mm(molecule)
!!
!!    Check and calculate how many atoms are in each MM fragment/molecule
!!    Written by Tommaso Giovannini, April 2019
!!
      implicit none
!
      class(mm) :: molecule
!
      integer :: i, k
!
      molecule%n_atoms_per_molecule = 1
!
      k = 1
!
      do i = 1, molecule%n_atoms - 1 
!
         if(molecule%molecule(i).eq.molecule%molecule(i+1)) then
!
            molecule%n_atoms_per_molecule(k) = molecule%n_atoms_per_molecule(k) + 1
!
         else 
!
            if(molecule%molecule(i+1).ne.(molecule%molecule(i) + 1)) then
!
               write(output%unit,'(/t6,a)') 'QM/MM Calculation'
               write(output%unit,'(t6,a41,i4)') 'WARNING! Something wrong: missing molecule = ', molecule%molecule(i) + 1
               write(output%unit,'(t6,a41,i4)') '                      but maximum molecule = ', molecule%max_mm_mol
               write(output%unit,'(t6,a)') 'FATAL ERROR'
               call output%error_msg('eT programm will stop.')
!
            endif
!
            k = k + 1
!
         endif
!
      enddo
!
   end subroutine check_molecule_mm
!
!
   subroutine read_geometry_mm(molecule)
!!
!!    Read MM geometry
!!    Written by Tommaso Giovannini, May 2019
!!
!!    Read atoms, molecule to which they belong, their coordinates 
!!    (assumed to be in units of Ångstrøm), and MM parameters
!!
      implicit none
!
      class(mm) :: molecule
!
      if(molecule%forcefield.eq.'non-polarizable') then  
!     
         call input%get_mm_geometry(molecule%n_atoms,      &
                                    molecule%symbol,       &
                                    molecule%molecule,     &
                                    molecule%coordinates,  &
                                    charge=molecule%charge)
!         
      else if(molecule%forcefield.eq.'fq') then
!      
         call input%get_mm_geometry(molecule%n_atoms,     &
                                    molecule%symbol,      &
                                    molecule%molecule,    &
                                    molecule%coordinates, &
                                    chi=molecule%chi,     &
                                    eta=molecule%eta)
!         
      endif
!      
   end subroutine read_geometry_mm
!
!
end module mm_class
