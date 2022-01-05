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
module electrostatic_environment_class
!
!!
!!    Electrostatic environment class module
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2018 and 2020
!!
!
   use parameters
!
   use global_in,             only: input 
   use global_out,            only: output
   use point_charges_class,   only: point_charges
   use ao_tool_class,         only: ao_tool
   use memory_manager_class,  only: mem
   use environment_class,     only: environment
!
   implicit none
!
   type, extends(environment) :: electrostatic_environment
!
      type(point_charges), private :: mm_points
!
      character(len=2), dimension(:), allocatable, private  :: symbols
!
   contains
!
!     Public class routines
!
      procedure :: initialize &
                => initialize_electrostatic_environment
!
      procedure :: update &
                => update_electrostatic_environment
!
      procedure :: get_energy &
                => get_energy_electrostatic_environment
!
      procedure :: print_energy &
                => print_energy_electrostatic_environment
!
      procedure :: print_description &
                => print_description_electrostatic_environment
!
!     Private class routines
!
      procedure, private :: get_mm_mm_energy &
                         => get_mm_mm_energy_electrostatic_environment
!
      procedure, private :: get_nuclei_mm_energy &
                         => get_nuclei_mm_energy_electrostatic_environment
!
      procedure, private :: print_geometry &
                         => print_geometry_electrostatic_environment
!
   end type  electrostatic_environment
!
   interface  electrostatic_environment
!
      procedure :: new_electrostatic_environment
!
   end interface  electrostatic_environment
!
contains
!
!
   pure function new_electrostatic_environment() result(embedding)
!!
!!    New environment 
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      type(electrostatic_environment) :: embedding 
!
      embedding%type_      = 'non-polarizable'
      embedding%n_charges  = 0
!
   end function new_electrostatic_environment
!
!
   subroutine initialize_electrostatic_environment(embedding, ao)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    Based on routine written by
!!    Tommaso Giovannini, May 2019
!!
!!    Reads the MM geometry and initializes 
!!    the point charges (R, q) of the
!!    non-polarizable electrostatic embedding
!!
!!    Initializes the integrals
!!
!!       v_wx = <w|q_I/|r - R_I||x>
!!
!!    where r is the electron coordinate.
!!
      implicit none
!
      class(electrostatic_environment), intent(inout) :: embedding
      type(ao_tool),                    intent(inout) :: ao
!
!     Read number of atoms and molecules
!
      embedding%n_charges = input%get_n_mm_atoms()
!
!     Read geometry 
!
      embedding%mm_points = point_charges(embedding%n_charges)
      call embedding%mm_points%initialize()
!
      allocate(embedding%symbols(embedding%n_charges))
!
      call input%get_mm_geometry_non_polarizable(embedding%n_charges,     &
                                                 embedding%symbols,       &
                                                 embedding%mm_points%r,   &
                                                 embedding%mm_points%q)
!
!     Initialize integrals 
!
      call ao%initialize_oei('electrostatic potential')
!
      call ao%initialize_external_charges(embedding%n_charges,    &
                                          embedding%mm_points%r,  &
                                          embedding%mm_points%q)
!
      call ao%construct_stored_oei('electrostatic potential')
!      
   end subroutine initialize_electrostatic_environment
!
!
   subroutine update_electrostatic_environment(embedding, ao, density)
!!
!!    Update
!!    Written by Sarai D. Folkestad
!!
!!    There is no update for electrostatic environment
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(electrostatic_environment),   intent(inout)   :: embedding
      type(ao_tool),                      intent(inout)   :: ao
      real(dp), dimension(ao%n, ao%n),    intent(in)      :: density
!
      call do_nothing(embedding)
      call do_nothing(ao)
      call do_nothing(density)
!
   end subroutine update_electrostatic_environment
!
!
   function get_energy_electrostatic_environment(embedding, ao, density) result(embedding_energy)
!!
!!    Get energy
!!    Written by Sarai D. Folkestad
!!
!!    Based on QM/MM routines written by Tommaso Giovannini 
!!
!!    1. MM-MM interaction:
!!
!!       sum_I<J q_I q_J / R_IJ      
!!
!!    2. Nuclei-MM interaction:
!!
!!       sum_i sum_I q_i Q_I/|r_i - R_I| 
!!
!!       Q_I : QM nuclear charge
!!       q_i : MM charge
!!
!!       R_I : MM charge position
!!       r_i : QM nuclei position
!!
!!    3. 1/2 Tr(Dv) 
!!
!!    Note the difference in contribution 3:
!!    The total energy contribution includes the additional 1/2 Tr(Dv)
!!    that is included in the energy through its contribution to the fock matrix
!!
      implicit none
!
      class(electrostatic_environment),   intent(in) :: embedding
      type(ao_tool),                      intent(in) :: ao
      real(dp), dimension(ao%n, ao%n),    intent(in) :: density
      real(dp)                                       :: embedding_energy
!
      real(dp) :: ddot
!
!     1. sum_I<J q_I q_J / R_IJ 
!
      embedding_energy = embedding%get_mm_mm_energy()
!
!     2. sum_i sum_I q_i Q_I/|r_i - R_I|
!
      embedding_energy = embedding_energy + embedding%get_nuclei_mm_energy(ao)
!
!     3. 1/2 Tr(v D)
!
      embedding_energy = embedding_energy &
                         + half*ddot(ao%n**2, density, 1, ao%v, 1)
!
   end function get_energy_electrostatic_environment
!
!
   subroutine print_energy_electrostatic_environment(embedding, ao, density)
!!
!!    Print energy
!!    Written by Sarai D. Folkestad
!!
!!
!!    SCF energy = Tr(Dv)
!!
!!    electrostatic energy = Tr(vD) + sum_i sum_I q_i Q_I/|r_i - R_I| 
!!                         + sum_I<J Q_I Q_J/|R_I - R_J|
!!
      implicit none
!
      class(electrostatic_environment)             :: embedding
      type(ao_tool),                   intent(in)  :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)  :: density
!
      real(dp) :: scf_energy, electrostatic_energy, ddot
!
      scf_energy           = ddot(ao%n**2, density, 1, ao%v, 1)
      electrostatic_energy = scf_energy &
                           + embedding%get_mm_mm_energy() &
                           + embedding%get_nuclei_mm_energy(ao)
!
      call output%printf('m', '- Summary of QM/MM energetics:', fs='(/t3,a)')
      call output%printf('m', 'a.u.             eV     kcal/mol', fs='(t42,a)')
!
!     SCF energy:  Tr(vD)
!
      call output%printf('m', 'QM/MM SCF Contribution: (f22.12)(f12.5) (f9.3)', &
                         reals=[scf_energy,                                     &
                                scf_energy*Hartree_to_eV,                       &
                                scf_energy*Hartree_to_kcalmol], fs='(t6,a)')
!
!     Electrostatic energy: 
!     Tr(vD) + sum_i sum_I q_i Q_I/|r_i - R_I| + sum_I<J Q_I Q_J/|R_I - R_J|
!
      call output%printf('m', 'QM/MM Electrostatic Energy:(f19.12)(f12.5) (f9.3)',  &
                         reals=[electrostatic_energy,                               &
                                electrostatic_energy*Hartree_to_eV,                 &
                                electrostatic_energy*Hartree_to_kcalmol], fs='(t6,a)')
!
   end subroutine print_energy_electrostatic_environment
!
!
   subroutine print_description_electrostatic_environment(embedding)
!!
!!    Print description
!!    Written by Tommaso Giovannini, April 2019
!!    
      implicit none
!
      class(electrostatic_environment), intent(in) :: embedding
!
      call output%printf('m', 'This is a QM/MM calculation', fs='(/t3,a)')     
!
      call output%printf('n', 'Embedding type: (a0)', chars=[trim(embedding%type_)], fs='(/t6,a)')
      call output%printf('n', 'Each atom of the MM portion is endowed with a charge which &
                         &value is a fixed external parameter.', &
                           ffs='(/t6,a)', fs='(t6,a)')

      call output%printf('n', 'The QM/MM electrostatic interaction energy is defined as:', &
                           fs='(/t6,a)')
                         
      call output%printf('n', 'E^ele_QM/MM = sum_i q_i * V_i(P)',  fs='(/t9,a)')
!
      call output%printf('n', 'where V_i(P) is the electrostatic potential due to the &
                         &QM density calculated at the position of the i-th charge q_i.', &
                           ffs='(/t6,a)', fs='(t6,a)')
      call output%printf('n', 'For further details, see:',  fs='(/t6,a)')
      call output%printf('n', 'Senn & Thiel, Angew. Chem. Int. Ed., 2009, 48, 1198−1229', &
                           fs='(t6,a/)')
!
      if(input%requested_cc_calculation()) &
         call output%printf('n','CC calculation: MM charges only affect MOs and Fock', fs='(t6,a/)')
!
      call embedding%print_geometry('angstrom')
      call embedding%print_geometry('bohr')
!
   end subroutine print_description_electrostatic_environment
!
!
   function get_mm_mm_energy_electrostatic_environment(embedding) result(E) 
!!
!!    Get MM energy contribution
!!    Written by Sarai D. Folkestad
!!
!!    Returns the MM-MM interaction energy 
!!
!!       E = sum_I<J Q_I Q_J / R_IJ 
!!
      implicit none
!
      class(electrostatic_environment), intent(in) :: embedding
!
      real(dp) :: E
!
      E = embedding%mm_points%get_coulomb_interaction()
!
   end function get_mm_mm_energy_electrostatic_environment
!
!
   function get_nuclei_mm_energy_electrostatic_environment(embedding, ao) result(E) 
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
      class(electrostatic_environment),   intent(in) :: embedding
      type(ao_tool),                      intent(in) :: ao
!
      real(dp) :: E
!
      type(point_charges) :: qm_points
!
      call ao%get_point_charges(qm_points)
!
      E = embedding%mm_points%get_coulomb_interaction(qm_points)
!
   end function get_nuclei_mm_energy_electrostatic_environment
!
!
   subroutine print_geometry_electrostatic_environment(embedding, units)
!!
!!    Print geometry
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Prints the MM geometry
!!
      implicit none 
!
      class(electrostatic_environment),                  intent(in)  :: embedding
      character(len=*),                                  intent(in)  :: units
!
      integer  :: k
      integer  :: line_length
      real(dp) :: conversion_factor
!
      conversion_factor = get_conversion_factor(from='angstrom', to=units)
!
!     Print header for geometry
!
      line_length = 68
!
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(/t6,a)')
      call output%printf('m', 'MM Geometry ((a0))', &
                        fs='(t30,a)', chars=[get_units_label(units)])
      call output%print_separator('m', line_length,'=', fs='(t6,a)')
!
      call output%printf('m', 'Atom                    X          Y          &
                            &Z            Charge', fs='(t6,a)')

      call output%print_separator('m', line_length,'=', fs='(t6,a)')
!
      do k = 1, embedding%n_charges
!
         call output%printf('m', '(a2)              (f11.6)(f11.6)(f11.6)    &
                            & (f11.6)', chars=[embedding%symbols(k)], &
                            reals=[embedding%mm_points%r(1, k)*conversion_factor, &
                            embedding%mm_points%r(2, k)*conversion_factor,        &
                            embedding%mm_points%r(3, k)*conversion_factor,        &
                            embedding%mm_points%q(k)], fs='(t6,a)')
      enddo
!     
      call output%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
!
   end subroutine print_geometry_electrostatic_environment
!
!
end module electrostatic_environment_class
