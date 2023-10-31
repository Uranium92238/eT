!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module fq_environment_class
!
!!
!!    FQ fq_environment class module
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2018 and 2020
!!
!!    See C. Cappelli, Int. J. Quantum Chem., 116(21), 1532-1542, (2016).
!!
!
   use parameters
!
   use global_in,             only: input
   use global_out,            only: output
   use memory_manager_class,  only: mem
!
   use environment_class,     only: environment
   use ao_tool_class,         only: ao_tool
   use mm_molecule_class,     only: mm_molecule
   use point_charges_class,   only: point_charges
!
   implicit none
!
   type, extends(environment) :: fq_environment

      integer,                                         private :: n_molecules
      type(mm_molecule), dimension(:),    allocatable, private :: molecules
      real(dp),          dimension(:,:),  allocatable, private :: D ! Ohno kernel
!
   contains
!
!     Public class routines
!
      procedure :: initialize &
                => initialize_fq_environment
!
      procedure :: update &
                => update_fq_environment
!
      procedure :: get_energy &
                => get_energy_fq_environment
!
      procedure :: print_energy &
                => print_energy_fq_environment
!
      procedure :: print_description &
                => print_description_fq_environment
!
!     Private class routines
!
      procedure, private :: set_chi_and_eta &
                         => set_chi_and_eta_fq_environment
!
      procedure, private :: get_q_dot_chi &
                         => get_q_dot_chi_fq_environment
!
      procedure, private :: construct_D &
                         => construct_D_fq_environment
!
      procedure, private :: get_point_charges &
                         => get_point_charges_fq_environment
!
      procedure, private :: print_geometry &
                         => print_geometry_fq_environment
!
      procedure, private :: get_mm_mm_energy &
                         => get_mm_mm_energy_fq_environment
!
      procedure, private :: get_nuclei_mm_energy &
                         => get_nuclei_mm_energy_fq_environment
!
      procedure, private :: get_potential_at_mm_points &
                         => get_potential_at_mm_points_fq_environment
!
      procedure, private :: get_positions &
                         => get_positions_fq_environment
!
      final :: destructor
!
   end type  fq_environment
!
   interface  fq_environment
!
      procedure :: new_fq_environment
!
   end interface  fq_environment
!
contains
!
!
   pure function new_fq_environment() result(embedding)
!!
!!    New fq_environment
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      type(fq_environment) :: embedding
!
      embedding%type_       = 'FQ'
      embedding%n_charges   = 0
      embedding%n_molecules = 0
!
   end function new_fq_environment
!
!
   subroutine initialize_fq_environment(embedding, ao)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    Based on routine written by
!!    Tommaso Giovannini, May 2019
!!
!!    Reads the geometry, chemical hardness (eta) and electronegativity (chi)
!!    Initializes the MM molecules
!!    Prepares the Ohno kernel (J)
!!
      implicit none
!
      class(fq_environment),  intent(inout) :: embedding
      type(ao_tool),          intent(inout) :: ao
!
      real(dp),         dimension(:,:),   allocatable :: coordinates
      real(dp),         dimension(:),     allocatable :: eta, chi
      integer,          dimension(:),     allocatable :: n_charges_in_molecule
      character(len=2), dimension(:),     allocatable :: symbol
!
      integer :: I
      integer :: first_atom
      integer :: last_atom
      integer :: atom_offset
!
!     Read number of atoms and molecules
!
      embedding%n_charges = input%get_n_mm_atoms()
      embedding%n_molecules = input%get_n_mm_molecules()
!
!     Read geometry
!
      call mem%alloc(n_charges_in_molecule, embedding%n_molecules)
      call mem%alloc(coordinates, 3, embedding%n_charges)
      call mem%alloc(chi, embedding%n_charges)
      call mem%alloc(eta, embedding%n_charges)
!
      allocate(symbol(embedding%n_charges))
!
      call input%get_mm_geometry_fq(embedding%n_charges,    &
                                   embedding%n_molecules,   &
                                   n_charges_in_molecule,   &
                                   symbol,                  &
                                   coordinates,             &
                                   chi,                     &
                                   eta)
!
!     Allocate array of molecules
!
      allocate(embedding%molecules(embedding%n_molecules))
!
!     Initialize MM molecules
!
      atom_offset = 0
!
      do I = 1, embedding%n_molecules
!
         embedding%molecules(I) = mm_molecule(n_charges_in_molecule(I))
!
         first_atom = atom_offset + 1
         last_atom  = atom_offset + n_charges_in_molecule(I)
!
         call embedding%molecules(I)%initialize_mm_atoms(coordinates(:, first_atom : last_atom), &
                                                           symbol(first_atom : last_atom))
!
         atom_offset = last_atom
!
      enddo
!
      call mem%dealloc(n_charges_in_molecule, embedding%n_molecules)
      call mem%dealloc(coordinates, 3, embedding%n_charges)
      deallocate(symbol)
!
!     Set charges or electronegativity and hardness (depending on the force field type)
!
      call embedding%set_chi_and_eta(embedding%n_charges, chi, eta)
!
      call mem%dealloc(chi, embedding%n_charges)
      call mem%dealloc(eta, embedding%n_charges)
!
      call embedding%construct_D()
!
      call ao%initialize_oei('electrostatic potential')
!
   end subroutine initialize_fq_environment
!
!
   subroutine update_fq_environment(embedding, ao, density)
!!
!!    Update charges
!!    Written by Tommaso Giovannini
!!
!!    Modified for fq_environment by Sarai D. Folkestad
!!
!!    Solves the equation
!!
!!       J q = - (potential + electronegativity)
!!
!!    where J is the Ohno kernel together with
!!    the charge constraints (see construct_D).
!!
!!    The equation is solved through inversion
!!    J^-1 is stored as the J
!!
!!    See C. Cappelli, Int. J. Quantum Chem., 116(21), 1532-1542, (2016).
!!
      implicit none
!
      class(fq_environment), intent(inout)         :: embedding
      type(ao_tool), intent(inout)                 :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)  :: density
!
      real(dp), dimension(:), allocatable :: rhs ! - (potential + electronegativity)
      real(dp), dimension(:), allocatable :: lhs ! solution (new charges)
!
      integer :: n_variables
      integer :: I, p, atom_offset, m, a, n_atoms
!
      real(dp), dimension(:), allocatable :: potential
!
      real(dp), dimension(:,:), allocatable :: R
      real(dp), dimension(:,:), allocatable :: v
!
      real(dp) :: ddot
!
      call mem%alloc(potential, embedding%n_charges)
!
      call embedding%get_potential_at_mm_points(potential, ao)
!
      call mem%alloc(v, ao%n, ao%n)
!
      I = 0
!
      do m = 1, embedding%n_molecules
         do a = 1, embedding%molecules(m)%n_atoms
!
            I = I + 1
!
            call ao%initialize_external_unit_charges(1, embedding%molecules(m)%get_r_i(a))
            call ao%get_oei('electrostatic potential unit', v)
!
            potential(I) = potential(I) + ddot(ao%n**2, v, 1, density, 1)
!
         enddo
      enddo
!
      call mem%dealloc(v, ao%n, ao%n)
!
      n_variables = embedding%n_charges + embedding%n_molecules
!
      call mem%alloc(rhs, n_variables, set_zero=.true.)
!
!     rhs = - (potential + electronegativity)
!
      atom_offset = 0
!
      do I = 1, embedding%n_molecules
         do p = 1, embedding%molecules(I)%n_atoms
!
            atom_offset = atom_offset + 1
!
            rhs(atom_offset) = - embedding%molecules(I)%get_chi_i(p)
!
         enddo
      enddo
!
      call daxpy(embedding%n_charges, -one, potential, 1, rhs, 1)
      call mem%dealloc(potential, embedding%n_charges)
!
      call mem%alloc(lhs, n_variables) ! solution -> new charges
!
!     q = - J^-1 (potential + electronegativity)
!
      call dgemv('N', &
                  n_variables, &
                  n_variables, &
                  one,         &
                  embedding%D, &
                  n_variables, &
                  rhs,         &
                  1,           &
                  zero,        &
                  lhs,         &
                  1)
!
      call mem%dealloc(rhs, n_variables)
!
      atom_offset = 1
!
      do I = 1, embedding%n_molecules
!
         n_atoms = embedding%molecules(I)%n_atoms
!
         call embedding%molecules(I)%set_q(lhs(atom_offset:atom_offset + n_atoms))
         atom_offset = atom_offset + n_atoms
!
      enddo
!
      call mem%alloc(R, 3, embedding%n_charges)
      call embedding%get_positions(R)
!
      call ao%initialize_external_charges(embedding%n_charges, R, lhs)
!
      call mem%dealloc(R, 3, embedding%n_charges)
      call mem%dealloc(lhs, n_variables)
!
      call ao%construct_stored_oei('electrostatic potential')
!
   end subroutine update_fq_environment
!
!
   function get_energy_fq_environment(embedding, ao, density) result(embedding_energy)
!!
!!    Get energy
!!    Written by Sarai D. Folkestad
!!
!!    Based on QM/MM and PCM routines written by Tommaso Giovannini
!!
!!    1. MM-MM interaction
!!
!!        1/2 sum_I q_I chi_I
!!
!!    2. Nuclei-MM interaction:
!!
!!       1/2 sum_i sum_I q_i Q_I/|r_i - R_I|
!!
!!       Q_I : MM charge
!!       q_i : QM nuclear charge
!!
!!       R_I : MM charge position
!!       r_i : QM nuclei position
!!
!!       Factor 1/2 due to cancelations of MM-MM interaction terms.
!!       See eq. (5), C. Cappelli, Int. J. Quantum Chem., 116(21), 1532-1542, (2016)
!!       and insert the linear equation for the FQs
!!
!!       electrostatic energy = Tr(vD) + sum_i sum_I q_i Q_I/|r_i - R_I|
!!
!!       scf energy = 1/2 Tr(vD) + 1/2 sum_i sum_I q_i Q_I/|r_i - R_I| + sum_i q_i chi_i
!!
      implicit none
!
      class(fq_environment),           intent(in) :: embedding
      type(ao_tool),                   intent(in) :: ao
      real(dp), dimension(ao%n, ao%n), intent(in) :: density
      real(dp) :: embedding_energy
!
      real(dp) :: TrvD, ddot, E_nuc_mm
!
      TrvD = ddot((ao%n)**2, ao%v, 1, density, 1)
!
      E_nuc_mm = embedding%get_nuclei_mm_energy(ao)
!
!     1. 1/2 sum_I q_I chi_I
!
      embedding_energy = embedding%get_mm_mm_energy()
!
!     2. 1/2 sum_i sum_I q_i Q_I/|r_i - R_I| + 1/2 Tr(vD)
!
      embedding_energy = embedding_energy + half*E_nuc_mm + half*TrvD
!
   end function get_energy_fq_environment
!
!
   subroutine print_energy_fq_environment(embedding, ao, density)
!!
!!    Print energy
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(fq_environment)                       :: embedding
      type(ao_tool),                   intent(in) :: ao
      real(dp), dimension(ao%n, ao%n), intent(in) :: density
!
      real(dp) :: scf_energy, electrostatic_energy
      real(dp) :: TrvD, ddot, E_nuc_mm
!
      TrvD = ddot((ao%n)**2, ao%v, 1, density, 1)
      E_nuc_mm = embedding%get_nuclei_mm_energy(ao)
!
!     Electrostatic energy:
!     Tr(vD) + sum_i sum_I q_i Q_I/|r_i - R_I|
!
      electrostatic_energy = E_nuc_mm + TrvD
!
!     SCF energy:
!     1/2 Tr(vD) + 1/2 sum_i sum_I q_i Q_I/|r_i - R_I| + 1/2 sum_i q_i chi_i
!
      scf_energy = half*electrostatic_energy &
                  + embedding%get_mm_mm_energy()
!
      call output%printf('m', '- Summary of QM/MM energetics:', fs='(/t3,a)')
      call output%printf('m', 'a.u.             eV     kcal/mol', fs='(t42,a)')
!
!     SCF energy:
!     1/2 Tr(vD) + 1/2 sum_i sum_I q_i Q_I/|r_i - R_I| + sum_i q_i chi_i
!
      call output%printf('m', 'QM/MM SCF Contribution: (f22.12)(f12.5) (f9.3)', &
                         reals=[scf_energy, &
                                scf_energy*Hartree_to_eV, &
                                scf_energy*Hartree_to_kcalmol], fs='(t6,a)')
!
!     Electrostatic energy:
!     Tr(vD) + sum_i sum_I q_i Q_I/|r_i - R_I|
!
      call output%printf('m', 'QM/MM Electrostatic Energy:(f19.12)(f12.5) (f9.3)', &
                         reals=[electrostatic_energy, &
                                electrostatic_energy*Hartree_to_eV, &
                                electrostatic_energy*Hartree_to_kcalmol], fs='(t6,a)')
!
   end subroutine print_energy_fq_environment
!
!
   subroutine set_chi_and_eta_fq_environment(embedding, n_charges, chi, eta)
!!
!!    Set chi and eta
!!    Written by Sarai D. Folkestad
!!
!!    Sets electronegativity and hardness
!!
      implicit none
!
      class(fq_environment), intent(inout) :: embedding
!
      integer, intent(in) :: n_charges
!
      real(dp), dimension(n_charges), intent(in) :: chi
      real(dp), dimension(n_charges), intent(in) :: eta
!
      integer :: I, first_atom, last_atom, atom_offset
!
      atom_offset = 0
!
      do I = 1, embedding%n_molecules
!
         first_atom = atom_offset + 1
         last_atom  = atom_offset + embedding%molecules(I)%n_atoms
!
         call embedding%molecules(I)%set_chi(chi(first_atom : last_atom))
         call embedding%molecules(I)%set_eta(eta(first_atom : last_atom))
!
         atom_offset = last_atom
!
      enddo
!
   end subroutine set_chi_and_eta_fq_environment
!
!
   subroutine construct_D_fq_environment(embedding)
!!
!!    Construction FQ matrix
!!    Written by Tommaso Giovannini, March 2019
!!
!!    Sets up the Ohno functional/kernel. For more information,
!!    see C. Cappelli, Int. J. Quantum Chem., 116(21), 1532-1542, (2016).
!!
!!       J_pq = eta_pq/(1 + eta_pq^2 R_pq^2)^1/2, p .ne. q
!!       J_pp = 2 eta_p
!!
!!    where
!!
!!       eta_pq = 1/2(eta_p + eta_q)
!!
!!    and
!!
!!       R_pq = |R_p - R_q|
!!
!!
!!    Additionally, the charge constraints for a given molecule
!!    are built into J.
!!
!
      use array_utilities, only: invert_in_place
!
      implicit none
!
      class(fq_environment) :: embedding
!
      integer :: I, J
      integer :: p, q
      integer :: atom_p, atom_q
!
      real(dp) :: R_pq, eta_pq
      real(dp), dimension(3) :: R_p, R_q
!
!     Allocate and construct ohno kernel
!
      call mem%alloc(embedding%D, &
                     embedding%n_charges + embedding%n_molecules, &
                     embedding%n_charges + embedding%n_molecules, &
                     set_zero=.true.)
!
      atom_p = 0
!
      do I = 1, embedding%n_molecules
         do p = 1, embedding%molecules(I)%n_atoms
!
            atom_p = atom_p + 1
            R_p = embedding%molecules(I)%get_r_i(p)
!
            atom_q = 0
!
            do J = 1, embedding%n_molecules
               do q = 1, embedding%molecules(J)%n_atoms
!
                  atom_q = atom_q + 1
!
                  R_q = embedding%molecules(J)%get_r_i(q)
!
                  R_pq = sqrt((R_p(1) - R_q(1))**2 + (R_p(2) - R_q(2))**2 + (R_p(3) - R_q(3))**2)
                  R_pq = R_pq*angstrom_to_bohr
!
                  eta_pq = half*(embedding%molecules(I)%get_eta_i(p) &
                               + embedding%molecules(J)%get_eta_i(q))
!
                  if (atom_p .ne. atom_q) then
!
                     embedding%D(atom_p, atom_q) = eta_pq/((one + (eta_pq**2)*(R_pq**2))**half)
!
                  else
!
                     embedding%D(atom_p, atom_q) = embedding%molecules(I)%get_eta_i(p)
!
                  endif
!
               enddo
!
            enddo
!
         enddo
!
      enddo
!
!     Add the charge contraints
!
      atom_p = 0
!
      do I = 1, embedding%n_molecules
         do p = 1, embedding%molecules(I)%n_atoms
!
            atom_p = atom_p + 1
!
            embedding%D(atom_p, I + embedding%n_charges) = one
            embedding%D(I + embedding%n_charges, atom_p) = one
!
         enddo
      enddo
!
      call invert_in_place(embedding%D,  &
                  embedding%n_charges + embedding%n_molecules)
!
   end subroutine construct_D_fq_environment
!
!
   pure function get_q_dot_chi_fq_environment(embedding) result(q_dot_chi)
!!
!!    Get q dot chi
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    Compute and returns sum_I q_I chi_I
!!
      implicit none
!
      class(fq_environment), intent(in) :: embedding
      real(dp) :: q_dot_chi
!
      integer :: I, p
!
      q_dot_chi = zero
!
      do I = 1, embedding%n_molecules
         do p = 1, embedding%molecules(I)%n_atoms
!
            q_dot_chi = q_dot_chi + embedding%molecules(I)%get_q_i(p) * &
            embedding%molecules(I)%get_chi_i(p)
!
         enddo
      enddo
!
   end function get_q_dot_chi_fq_environment
!
!
   pure function get_mm_mm_energy_fq_environment(embedding) result(E)
!!
!!    Get MM energy contribution
!!    Written by Sarai D. Folkestad
!!
!!    Returns the MM-MM interaction energy
!!
!!       E = 1/2 sum_I q_I chi_I
!!
      implicit none
!
      class(fq_environment), intent(in) :: embedding
      real(dp) :: E
!
      E = half*embedding%get_q_dot_chi()
!
   end function get_mm_mm_energy_fq_environment
!
!
   subroutine print_description_fq_environment(embedding)
!!
!!    Print description
!!    Written by Tommaso Giovannini, April 2019
!!
      implicit none
!
      class(fq_environment), intent(in) :: embedding
!
      call output%printf('m', 'This is a QM/MM calculation', fs='(/t3,a)')
!
      call output%printf('n', 'Polarizable Embedding: (a0)', &
                        chars=[trim(embedding%type_)],  fs='(/t6,a)')
!
      call output%printf('n', 'Each atom of the MM portion is endowed with a charge '        //&
                         'which value can vary in agreement with the Electronegativity' //&
                         'Equalization Principle (EEP), which states that at the '      //&
                         'equilibrium each atom has the same electronegativity.',         &
                          fs='(t6,a)', ffs='(/t6,a)')
!
      call output%printf('n', 'The force field is defined in terms of electronegativity ' //&
                         '(Chi) and chemical hardness (Eta), which are specified '   //&
                         'for each MM atom.', fs='(t6,a)', ffs='(/t6,a)')
!
      call output%printf('n', 'The QM/MM electrostatic interaction energy is defined as:', &
                           fs='(/t6,a)')
      call output%printf('n', 'E^ele_QM/MM = sum_i q_i * V_i(P)',  fs='(/t9,a)')
!
      call output%printf('n', 'where V_i(P) is the electrostatic potential due to the QM'  //&
                         'density calculated at the position of the i-th charge q_i.' //&
                         'The values of the charges are obtained by solving a linear' //&
                         'equation:', fs='(t6,a)', ffs='(/t6,a)')
      call output%printf('n', 'Dq = -Chi - V(P)',  fs='(/t9,a)')
      call output%printf('n', 'For further details, see:',  fs='(/t6,a)')
      call output%printf('n', 'T.Giovannini, F.Egidi, C.Cappelli, Chem. Soc. Rev., 2020, 49, 5664-5677',  fs='(t6,a/)')
!
      if(input%requested_cc_calculation()) then
         call output%printf('n', 'CC calculation: zero-order approximation',  fs='(t6,a)')
         call output%printf('n', 'FQ charges only affect MOs and Fock',  fs='(t6,a/)')
      endif
!
      call embedding%print_geometry('angstrom')
      call embedding%print_geometry('bohr')
!
   end subroutine print_description_fq_environment
!
!
   subroutine print_geometry_fq_environment(embedding, units)
!!
!!    Print geometry
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Prints the MM geometry
!!
      implicit none

      class(fq_environment), intent(in) :: embedding
!
      character(len=*), intent(in) :: units
!
      integer :: m, a, line_length
!
      real(dp) :: conversion_factor
!
      real(dp), dimension(3) :: r
!
!     Determine conversion factor
!
      conversion_factor = get_conversion_factor(from='angstrom', to=units)
!
!     Print header for geometry
!
      line_length = 68
!
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(/t6,a)')
      call output%printf('m', 'MM Geometry ((a0))', &
                        fs='(t33,a)', chars=[get_units_label(units)])
      call output%print_separator('m', line_length,'=', fs='(t6,a)')
!
         call output%printf('m', 'Atom    Mol         X          Y          Z   &
                            &      Chi        Eta', fs='(t6,a)')

      call output%print_separator('m', line_length,'=', fs='(t6,a)')
!
      do m = 1, embedding%n_molecules
         do a = 1, embedding%molecules(m)%n_atoms
!
            r = embedding%molecules(m)%get_r_i(a)
!
            call output%printf('m', '(a2)    (i4)  (f11.6)(f11.6)(f11.6)(f11.6)(f11.6)', &
                           chars=[embedding%molecules(m)%atoms(a)%symbol], &
                           ints=[m], &
                           reals=[r(1)*conversion_factor, &
                           r(2)*conversion_factor, &
                           r(3)*conversion_factor, &
                           embedding%molecules(m)%get_chi_i(a), &
                           embedding%molecules(m)%get_eta_i(a)], fs='(t6,a)')
         enddo
      enddo
!
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
!
   end subroutine print_geometry_fq_environment
!
!
   function get_nuclei_mm_energy_fq_environment(embedding, ao) result(E)
!!
!!    Get MM energy contribution
!!    Written by Sarai D. Folkestad
!!
!!    Returns the nuclei-MM interaction energy
!!
!!       E = sum_i sum_I q_i Q_I/|r_i - R_I|
!!
      implicit none
!
      class(fq_environment),  intent(in)    :: embedding
      type(ao_tool),          intent(in)    :: ao
      real(dp)                              :: E
!
      type(point_charges) :: qm_points, fq_points
!
!     Get QM and MM point charges
!
      call embedding%get_point_charges(fq_points)
      call ao%get_point_charges(qm_points)
!
!     Calculate interaction
!
      E = fq_points%get_coulomb_interaction(qm_points)
!
   end function get_nuclei_mm_energy_fq_environment
!
!
   subroutine get_potential_at_mm_points_fq_environment(embedding, potential, ao)
!!
!!    Get potential at mm points
!!    Written by Sarai D. Folkestad
!!
!!    Returns the potential at the mm charges from the
!!    nuclei
!!
!!      potential(I) = sum_J Q_J/|R_J - r_I|
!!
!!    where Q_J and R_J are the charge and position of the QM nuclei
!!    and r_I is an MM point
!!
      implicit none
!
      class(fq_environment),                    intent(in)  :: embedding
      real(dp), dimension(embedding%n_charges), intent(out) :: potential
      type(ao_tool),                            intent(in)  :: ao
!
      type(point_charges)                                   :: qm_points
      type(point_charges)                                   :: fq_points
!
      call ao%get_point_charges(qm_points)
      call embedding%get_point_charges(fq_points)
!
!     Calculate potential at mm points from qm nuclei
!
      call qm_points%get_potential_at_external_points(fq_points%n_charges,  &
                                                      fq_points%r,          &
                                                      potential)
!
   end subroutine get_potential_at_mm_points_fq_environment
!
!
   subroutine get_point_charges_fq_environment(embedding, pc)
!!
!!    Get point charges
!!    Written by Sarai D.Folkestad, 2020
!!
      implicit none
!
      class(fq_environment), intent(in) :: embedding
!
      type(point_charges), intent(out)  :: pc
!
      integer                           :: m, a, I
!
      pc = point_charges(embedding%n_charges)
      call pc%initialize()
!
      I = 0
!
      do m = 1, embedding%n_molecules
         do a = 1, embedding%molecules(m)%n_atoms
!
            I = I + 1
!
            pc%q(I)   = embedding%molecules(m)%get_q_i(a)
            pc%r(:,I) = embedding%molecules(m)%get_r_i(a)
!
         enddo
      enddo
!
   end subroutine get_point_charges_fq_environment
!
!
   subroutine get_positions_fq_environment(embedding, r)
!!
!!    Get positions
!!    Written by Sarai D.Folkestad, 2020
!!
      implicit none
!
      class(fq_environment),                       intent(in)  :: embedding
      real(dp), dimension(3, embedding%n_charges), intent(out) :: r
!
      integer                                               :: m, a, I
!
      I = 0
!
      do m = 1, embedding%n_molecules
         do a = 1, embedding%molecules(m)%n_atoms
!
            I = I + 1
!
            r(:,I) = embedding%molecules(m)%get_r_i(a)
!
         enddo
      enddo
!
   end subroutine get_positions_fq_environment
!
!
   subroutine destructor(embedding)
!!
!!    Destructor
!!    Written by Sarai D. Folkestad, Sep 2020
!!
      implicit none
!
      type(fq_environment) :: embedding
!
      if (allocated(embedding%D)) then
!
         call mem%dealloc(embedding%D, &
                         (embedding%n_charges + embedding%n_molecules), &
                         (embedding%n_charges + embedding%n_molecules))
!
      endif
!
   end subroutine destructor
!
!
end module fq_environment_class
