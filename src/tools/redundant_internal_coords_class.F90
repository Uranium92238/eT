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
module redundant_internal_coords_class
!
!!
!!    Redundant internal coordinates class
!!    Written by Eirik F. Kjønstad, 2021
!!
!
   use parameters
   use memory_manager_class, only: mem
   use global_out, only: output
   use atomic_center_class, only: atomic_center
   use array_utilities, only: get_l2_norm
   use math_utilities, only: cross_product_R3, dot_R3, norm_R3
!
   implicit none
!
!  Covalennt radii for elements 1-36
!  From Dalton Trans., 2008, 2832–2838, given in Angstrom
!
   real(dp), dimension(36) :: covalent_radius                                    & 
            = [0.31d0, 0.28d0,                                                   &
               1.28d0, 0.96d0, 0.84d0, 0.73d0, 0.71d0, 0.66d0, 0.57d0, 0.58d0,   &
               1.66d0, 1.41d0, 1.21d0, 1.11d0, 1.07d0, 1.05d0, 1.02d0,           &
               1.06d0, 2.03d0, 1.76d0, 1.70d0, 1.60d0, 1.53d0, 1.39d0, 1.39d0,   &
               1.32d0, 1.26d0, 1.24d0, 1.32d0, 1.22d0, 1.22d0, 1.20d0, 1.19d0,   &
               1.20d0, 1.20d0, 1.16d0]*angstrom_to_bohr
!
   type :: bond 
!
      integer, dimension(2) :: atoms 
      real(dp) :: length
!
   end type bond
!
   type :: angle
!
      real(dp) :: theta
      integer, dimension(3) :: atoms 
!
   end type angle
!
   type :: dihedral 
!
      real(dp) :: theta
      integer, dimension(4) :: atoms 
!
   end type dihedral
!
   type :: redundant_internal_coords
!
      integer, private :: n_atoms, n_bonds, n_angles, n_dihedrals, n_internal, n_cartesian
!
      type(bond), dimension(:), allocatable, private :: bonds  
      integer, dimension(:,:), allocatable, private  :: bond_matrix 
!
      type(angle), dimension(:), allocatable, private    :: angles
      type(dihedral), dimension(:), allocatable, private :: dihedrals
!
      integer, dimension(:), allocatable, private    :: Z ! atomic charges, [n_atoms]
      real(dp), dimension(:,:), allocatable, private :: R ! atomic positions, [3, n_atoms]
!
!     B_ij = dq_i/dx_j (q = internals, x = cartesian)
      real(dp), dimension(:,:), allocatable, private :: B 
!
!     Pseudoinverse of B
      real(dp), dimension(:,:), allocatable, private :: B_inverse 
!
!     Projection matrix: P = G G_inv, where G = B B^T
!     A general change in internals (dq) does not necessarily correspond to some dx 
!     (but dq' = P*dq always corresponds to a dx)
      real(dp), dimension(:,:), allocatable, private :: P 
!
      real(dp), private :: pseudoinverse_threshold
      real(dp), private :: dihedral_skip_threshold
!
      real(dp), private :: internal_to_cartesian_threshold
      integer, private  :: internal_to_cartesian_max_iterations
!
      real(dp), dimension(:), allocatable, private :: q, prev_q 
!
   contains
!
      procedure, public :: initialize
      procedure, public :: set_geometry
!
      procedure, public :: get_n_internal
      procedure, public :: get_approximate_diagonal_Hessian 
!
      procedure, public :: Wilson_project_vector 
      procedure, public :: Wilson_project_matrix 
!
      procedure, public :: calculate_internal_gradient
!
      procedure, public :: compute_next_geometry
!
      procedure, public :: print_geometry 
!
      procedure, private :: identify_internal_coordinates
!
      procedure, private :: identify_bonds 
      procedure, private :: identify_bond_angles 
      procedure, private :: identify_dihedral_angles
!
      procedure, private :: calculate_internal_coordinates
!
      procedure, private :: calculate_bond_lengths
      procedure, private :: calculate_bond_angles
      procedure, private :: calculate_dihedral_angles
!
      procedure, public :: get_internals
      procedure, public :: get_internals_difference
!
      procedure, private :: add_dihedral
      procedure, private :: pseudoinvert_Wilson_B_and_construct_P
!
      procedure, private :: construct_Wilson_B
!
      procedure, private :: construct_Wilson_B_bonds 
      procedure, private :: construct_Wilson_B_bond_angles
      procedure, private :: construct_Wilson_B_dihedral_angles
!
      procedure, private :: perform_cartesian_step
      procedure, private :: remove_multiples_of_360
!
      final :: destructor
!
   end type redundant_internal_coords
!
   interface redundant_internal_coords
!
      procedure :: new_redundant_internal_coords
!
   end interface redundant_internal_coords
!
contains
!
!
   function new_redundant_internal_coords(n_cartesian) result(this)
!!
!!    New redundant internal coords
!!    Written by Eirik F. Kjønstad, 2021  
!!
      implicit none
!
      type(redundant_internal_coords) :: this
!
      integer :: n_cartesian 
!
      this%n_atoms = n_cartesian/3 
      this%n_cartesian = n_cartesian
!
      this%pseudoinverse_threshold = 1.0d-4
!
      this%dihedral_skip_threshold = 1.0d-10
!
      this%internal_to_cartesian_threshold = 1.0d-6
      this%internal_to_cartesian_max_iterations = 25
!
   end function new_redundant_internal_coords
!
!
   subroutine initialize(this, R, Z)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Makes allocations and sets the geometry and atomic charges.
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      integer, dimension(this%n_atoms)     :: Z
      real(dp), dimension(3, this%n_atoms) :: R 
!
      call mem%alloc(this%R, 3, this%n_atoms)
      call mem%alloc(this%Z, this%n_atoms)
!
      this%R = R
      this%Z = Z 
!
      call mem%alloc(this%bond_matrix, this%n_atoms, this%n_atoms)
      this%bond_matrix = 0
!
      call this%identify_internal_coordinates()
!
      call mem%alloc(this%B, this%n_internal, this%n_cartesian)
      call mem%alloc(this%B_inverse, this%n_cartesian, this%n_internal)
      call mem%alloc(this%P, this%n_internal, this%n_internal)
!
      call mem%alloc(this%q, this%n_internal)
      call mem%alloc(this%prev_q, this%n_internal)
!
      this%prev_q = zero 
!
      call this%set_geometry(R)
!
      call this%print_geometry()
!
      call this%get_internals(this%q)
!
   end subroutine initialize
!
!
   subroutine set_geometry(this, R)
!!
!!    Set geometry
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Sets the coordinates R and updates all the quantities that depend on R
!!    (internal coordinates, Wilson, Wilson-inverse, projection, etc.)
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      real(dp), dimension(this%n_cartesian), intent(in) :: R 
!
      call dcopy(this%n_cartesian, R, 1, this%R, 1) 
!
      call this%calculate_internal_coordinates()
!
      call this%construct_Wilson_B()
      call this%pseudoinvert_Wilson_B_and_construct_P()
!
   end subroutine set_geometry
!
!
   subroutine get_internals_difference(this, dq)
!!
!!    Get internals differece
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(redundant_internal_coords) :: this 
!
      real(dp), dimension(this%n_internal) :: dq 
!
      dq = this%q - this%prev_q 
!
      call this%remove_multiples_of_360(dq)
      call this%Wilson_project_vector(dq)
!
   end subroutine get_internals_difference
!
!
   subroutine identify_internal_coordinates(this)
!!
!!    Identify internal coordinates 
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Identifies the atoms involved in bonds (A-B), 
!!    as well as the bond angles (A-B-C) and dihedral angles (A-B-C-D).
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this  
!
!     Count number of bonds and then allocate and set them
!     (same with angles and dihedrals)
!
      call this%identify_bonds(set_bonds=.false.)
      allocate(this%bonds(this%n_bonds))
      call this%identify_bonds(set_bonds=.true.)
!  
      call this%identify_bond_angles(set_angles=.false.)
      allocate(this%angles(this%n_angles))
      call this%identify_bond_angles(set_angles=.true.)
!
      call this%identify_dihedral_angles(set_angles=.false.)
      allocate(this%dihedrals(this%n_dihedrals))
      call this%identify_dihedral_angles(set_angles=.true.)
!
      this%n_internal = this%n_bonds + this%n_angles + this%n_dihedrals
!
      if (this%n_internal .lt. 3*this%n_atoms-6) then
!
         call output%error_msg('Not enough internal coordinates found. Only (i0) found, but expected at least 3N-6 = (i0). &
                                &Try starting from a better geometry.', &
                                ints=[this%n_internal, this%n_atoms*3-6])
!
      endif
!
   end subroutine identify_internal_coordinates
!
!
   subroutine identify_bonds(this, set_bonds)
!!
!!    Identify bonds 
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Counts number of bonds and sets bond matrix array.
!!
!!    If set_bonds = .true., this%bonds is assumed to be allocated 
!!    and is set by the routine
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      logical, intent(in) :: set_bonds 
!
      integer :: bond, A, A_bond, B, Z_A, Z_B 
!
      real(dp), dimension(3) :: R 
!
      logical :: bonded 
!
      integer :: q 
!
      bond = 0 ! Counter for the number of bonds
!
      do A = 1, this%n_atoms 
!
         A_bond = 0
!
         do B = 1, this%n_atoms 
!
            if (A .eq. B) cycle
!
            Z_A = this%Z(A)
            Z_B = this%Z(B)    
!
            do q = 1, 3 
!
               R(q) = this%R(q,A) - this%R(q,B)
!
            end do 
!
            bonded = norm_R3(R) .lt. 1.3d0 * (covalent_radius(Z_A) + covalent_radius(Z_B))
!
            if (bonded) then 
!
               A_bond = A_bond + 1 
               this%bond_matrix(A, A_bond) = B
!
               if (A .gt. B) cycle 
!
               bond = bond + 1
!
               if (set_bonds) &
                  this%bonds(bond)%atoms  = [A, B]
!  
            endif
!
         enddo
      enddo
!
      this%n_bonds = bond 
!
   end subroutine identify_bonds
!
!
   subroutine identify_bond_angles(this, set_angles)
!!
!!    Identify bond angles
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Counts number of bond angles.
!!
!!    If set_angles = .true., this%angles is assumed to be allocated 
!!    and is set by the routine
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      logical, intent(in) :: set_angles 
!
      integer :: angle, B, A, A_red, C_red, C
!
      angle = 0
!
      do B = 1, this%n_atoms
!
         A_red = 1
!
         do while (this%bond_matrix(B, A_red) .gt. 0)
!
            A = this%bond_matrix(B, A_red)
!
            C_red = 1
!
            do while (this%bond_matrix(B, C_red) .gt. 0) 
!
               C = this%bond_matrix(B, C_red)
!
               if (A .lt. C) then ! Avoid  double counting
!
                  angle = angle + 1
!
                  if (set_angles) &
                     this%angles(angle)%atoms = [A, B, C]
!
               endif
!
               C_red = C_red + 1
!
            enddo
!
            A_red = A_red + 1
!
         enddo
!
      enddo
!
      this%n_angles = angle
!
   end subroutine identify_bond_angles
!
!
   subroutine identify_dihedral_angles(this, set_angles)
!!
!!    Identify dihedral angles
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Counts number of dihedral angles.
!!
!!    If set_angles = .true., this%dihedrals is assumed to be allocated 
!!    and is set by the routine
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      logical, intent(in) :: set_angles
!
      integer :: dihedral, angle1, angle2, A, B, C, D 
!
      dihedral = 0
!
      do angle1 = 1, this%n_angles
         do angle2 = angle1 + 1, this%n_angles
!
!           For the angles to define a dihedral, two of the atoms must coincide,
!           and the end-points must be different
!
            if (this%angles(angle1)%atoms(2) .eq. this%angles(angle2)%atoms(3) .and. &
                this%angles(angle1)%atoms(3) .eq. this%angles(angle2)%atoms(2) .and. &
                this%angles(angle1)%atoms(1) .ne. this%angles(angle2)%atoms(1)) then 
!
               A = this%angles(angle1)%atoms(1)
               B = this%angles(angle1)%atoms(2)
               C = this%angles(angle1)%atoms(3)
               D = this%angles(angle2)%atoms(1)
!
               call this%add_dihedral(A, B, C, D, dihedral, set_angles)       
!
            else if (this%angles(angle1)%atoms(2) .eq. this%angles(angle2)%atoms(1) .and. &
                     this%angles(angle1)%atoms(3) .eq. this%angles(angle2)%atoms(2) .and. &
                     this%angles(angle1)%atoms(1) .ne. this%angles(angle2)%atoms(3)) then 
!
               A = this%angles(angle1)%atoms(1)
               B = this%angles(angle1)%atoms(2)
               C = this%angles(angle1)%atoms(3)
               D = this%angles(angle2)%atoms(3)
!
               call this%add_dihedral(A, B, C, D, dihedral, set_angles)    
!
            else if (this%angles(angle1)%atoms(1) .eq. this%angles(angle2)%atoms(2) .and. &
                     this%angles(angle1)%atoms(2) .eq. this%angles(angle2)%atoms(1) .and. &
                     this%angles(angle1)%atoms(3) .ne. this%angles(angle2)%atoms(3)) then  
!
               A = this%angles(angle2)%atoms(3)
               B = this%angles(angle2)%atoms(2)
               C = this%angles(angle2)%atoms(1)
               D = this%angles(angle1)%atoms(3)
!
               call this%add_dihedral(A, B, C, D, dihedral, set_angles)    
!
            else if (this%angles(angle1)%atoms(1) .eq. this%angles(angle2)%atoms(2) .and. &
                     this%angles(angle1)%atoms(2) .eq. this%angles(angle2)%atoms(3) .and. &
                     this%angles(angle1)%atoms(3) .ne. this%angles(angle2)%atoms(1)) then  
!
               A = this%angles(angle2)%atoms(1)
               B = this%angles(angle2)%atoms(2)
               C = this%angles(angle2)%atoms(3)
               D = this%angles(angle1)%atoms(3)
!
               call this%add_dihedral(A, B, C, D, dihedral, set_angles)  
!
            endif
!
         enddo
      enddo
!
      this%n_dihedrals = dihedral
!
   end subroutine identify_dihedral_angles
!
!
   subroutine print_geometry(this)
!!
!!    Print geometry
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(redundant_internal_coords) :: this 
!
      integer :: bond, angle, dihedral
!
      call output%printf('m', 'Molecular geometry (in redundant internal coordinates):', fs='(/t6,a)')
      call output%print_separator('n', 55,'-', fs='(t6,a)')
!
      do bond = 1, this%n_bonds 
!
         call output%printf('m', 'R((i0),(i0))         (f19.12)', &
                                 ints=[this%bonds(bond)%atoms(1), this%bonds(bond)%atoms(2)], &
                                 reals=[this%bonds(bond)%length], fs='(t6,a)')
!
      enddo
!
      do angle = 1, this%n_angles
!
         call output%printf('m', 'A((i0),(i0),(i0))       (f19.12)', &
                                 ints=[this%angles(angle)%atoms(:)], &
                                 reals=[this%angles(angle)%theta], fs='(t6,a)')
!
      enddo
!
      do dihedral = 1, this%n_dihedrals
!
         call output%printf('m', 'D((i0),(i0),(i0),(i0))     (f19.12)', &
                                 ints=[this%dihedrals(dihedral)%atoms(:)], &
                                 reals=[this%dihedrals(dihedral)%theta], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('n', 55,'-', fs='(t6,a)')
!
   end subroutine print_geometry
!
!
   subroutine add_dihedral(this, A, B, C, D, dihedral, set_angles)
!!
!!    Add dihedral
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Adds a dihedral for the quartet A-B-C-D provided it satisfies the 180 degree criterion.
!!
!!    The routine increments dihedral by 1 if the quartet passes the test and is actually added as a dihedral.
!!
      implicit none 
!
      class(redundant_internal_coords) :: this 
!
      integer, intent(in) :: A, B, C, D 
!
      logical, intent(in) :: set_angles
!
      integer, intent(inout) :: dihedral
!
      real(dp), dimension(3) :: AB, BC, CD, ABxBC, BCxCD
!
      real(dp) :: norm_AB, norm_BC, norm_CD, theta_ABC, theta_BCD
!
      integer :: q 
!
      do q = 1, 3 
!
         AB(q) = this%R(q,A) - this%R(q,B)
         BC(q) = this%R(q,B) - this%R(q,C)
         CD(q) = this%R(q,C) - this%R(q,D)
!
      end do 
!
      norm_AB = norm_R3(AB)
      norm_BC = norm_R3(BC)
      norm_CD = norm_R3(CD)
!
      ABxBC = cross_product_R3(AB,BC)
      BCxCD = cross_product_R3(BC,CD)
!
      theta_ABC = acos(dot_R3(AB,BC)/(norm_AB * norm_BC))
      theta_BCD = acos(dot_R3(BC,CD)/(norm_BC * norm_CD))
!
      if ((abs(theta_ABC - pi) .lt. this%dihedral_skip_threshold) .or. &
          (abs(theta_BCD - pi) .lt. this%dihedral_skip_threshold)) then 
!
         call output%printf('v', 'Skipping dihedral between (i0), (i0), (i0), and (i0)', ints=[A,B,C,D])
!
         call output%printf('v', 'theta_ABC: (e9.2), theta_BCD: (e9.2)', reals=[theta_ABC, theta_BCD])
!
      else 
!
         dihedral = dihedral + 1
!
         if (set_angles) &
            this%dihedrals(dihedral)%atoms = [A, B, C, D]      
!
      end if 
!  
   end subroutine add_dihedral
!
!
   pure function get_n_internal(this) result(n_internal)
!!
!!    Get n internal
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      integer :: n_internal
!
      n_internal = this%n_internal
!
   end function get_n_internal
!
!
   subroutine get_internals(this, q)
!!
!!    Get internals
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Returns q = [q_bonds, q_angles, q_dihedrals]
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      real(dp), dimension(this%n_internal), intent(out) :: q 
!
      integer :: bond, angle, dihedral 
!  
      do bond = 1, this%n_bonds
!
         q(bond) = this%bonds(bond)%length
!
      enddo
!
      do angle = 1, this%n_angles
!
         q(this%n_bonds + angle) = this%angles(angle)%theta
!
      enddo
!
      do dihedral = 1, this%n_dihedrals
!
         q(this%n_bonds + this%n_angles + dihedral) = this%dihedrals(dihedral)%theta
!
      enddo
!
   end subroutine get_internals
!
!
   subroutine calculate_internal_coordinates(this)
!!
!!    Calculate internal coordinates
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    q = q(R), internal coordinates corresponding to the cartesian coordinate R    
!!
      implicit none
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      call this%calculate_bond_lengths()
      call this%calculate_bond_angles()
      call this%calculate_dihedral_angles()
!
   end subroutine calculate_internal_coordinates
!
!
   subroutine calculate_bond_lengths(this)
!!
!!    Calculate bond lengths
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      integer :: bond 
!
      integer :: A, B 
!
      real(dp), dimension(3) :: AB 
!
      integer :: q 
!
      do bond = 1, this%n_bonds
!
         A = this%bonds(bond)%atoms(1)
         B = this%bonds(bond)%atoms(2)
!
         do q = 1, 3 
!
            AB(q) = this%R(q,A) - this%R(q,B)
!
         end do 
!
         this%bonds(bond)%length = norm_R3(AB)
!
      enddo
!
   end subroutine calculate_bond_lengths
!
!
   subroutine calculate_bond_angles(this)
!!
!!    Calculate bond angles
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      integer :: angle
!
      integer :: A, B, C 
!
      real(dp), dimension(3) :: AB, CB  
!
      integer :: q 
!
      do angle = 1, this%n_angles
!
         A = this%angles(angle)%atoms(1)
         B = this%angles(angle)%atoms(2)
         C = this%angles(angle)%atoms(3)
!
         do q = 1, 3 
!
            AB(q) = this%R(q,A) - this%R(q,B)
            CB(q) = this%R(q,C) - this%R(q,B)
!
         end do 
!
         this%angles(angle)%theta = acos(dot_R3(AB,CB)/(norm_R3(AB)*norm_R3(CB)))              
!
      enddo
!
   end subroutine calculate_bond_angles
!
!
   subroutine calculate_dihedral_angles(this)
!!
!!    Calculate dihedral angles
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      integer :: dihedral
!
      integer :: A, B, C, D 
!
      real(dp), dimension(3) :: AB, BC, CD, ABxBC, BCxCD
!
      real(dp) :: norm_BC, y, x 
!
      integer :: q 
!
      do dihedral = 1, this%n_dihedrals
!
         A = this%dihedrals(dihedral)%atoms(1)
         B = this%dihedrals(dihedral)%atoms(2)
         C = this%dihedrals(dihedral)%atoms(3)
         D = this%dihedrals(dihedral)%atoms(4)
!
         do q = 1, 3 
!
            AB(q) = this%R(q,A) - this%R(q,B) 
            BC(q) = this%R(q,B) - this%R(q,C) 
            CD(q) = this%R(q,C) - this%R(q,D) 
!
         end do 
!
         norm_BC = norm_R3(BC)
!
         ABxBC = cross_product_R3(AB,BC)
         BCxCD = cross_product_R3(BC,CD)         
!
         y = dot_R3(norm_BC * AB, BCxCD)
         x = dot_R3(ABxBC, BCxCD)       
!
         this%dihedrals(dihedral)%theta = -atan2(y, x)      
!
      enddo
!
   end subroutine calculate_dihedral_angles
!
!
   subroutine construct_Wilson_B(this)
!!
!!    Construct Wilson B
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Constructs the Wilson B matrix, given by
!!
!!       B_ij = dq_i/dx_j 
!!
!!    where q_i are the internal coordinates and x_j are the Cartesian coordinates.
!!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      call this%construct_Wilson_B_bonds()
      call this%construct_Wilson_B_bond_angles()
      call this%construct_Wilson_B_dihedral_angles()
!
      call output%print_matrix('v', 'Wilson B', this%B, this%n_internal, this%n_cartesian)
!
   end subroutine construct_Wilson_B
!
!
   subroutine construct_Wilson_B_bonds(this)
!!
!!    Construct Wilson B (bonds)
!!    Written by Eirik F. Kjønstad, 2021
!!
      use math_utilities, only: delta
!
      implicit none 
!
      class(redundant_internal_coords), target, intent(inout) :: this 
!
      real(dp), dimension(:,:,:), pointer :: B_p 
!
      integer :: m, n, bond, a, i 
!
      real(dp) :: norm_u 
!
      real(dp), dimension(3) :: u 
!
      integer :: q 
!
      B_p(1:this%n_internal, 1:3, 1:this%n_atoms) => this%B
!
      do bond = 1, this%n_bonds
!
         m = this%bonds(bond)%atoms(1)
         n = this%bonds(bond)%atoms(2)
!
         do q = 1, 3 
!
            u(q) = this%R(q,m) - this%R(q,n)
!
         enddo 
!
         norm_u = norm_R3(u)
!
         u = u/norm_u 
!
         do a = 1, this%n_atoms 
            do i = 1, 3
!
               B_p(bond, i, a) = (delta(a, m) - delta(a, n))*u(i)
!
            enddo
         enddo
!
      enddo
!
   end subroutine construct_Wilson_B_bonds
!
!
   subroutine construct_Wilson_B_bond_angles(this)
!!
!!    Construct Wilson B (bond angles)
!!    Written by Eirik F. Kjønstad, 2021
!!
      use math_utilities, only: delta
!
      implicit none 
!
      class(redundant_internal_coords), target, intent(inout) :: this 
!
      real(dp), dimension(:,:,:), pointer :: B_p 
!
      real(dp), dimension(3) :: v1, v2, u, v, w, uxw, wxv
!
      integer :: m, o, n, angle, a, i 
!
      real(dp) :: norm_u, norm_v, norm_w
!
      integer :: q 
!
      B_p(1:this%n_internal, 1:3, 1:this%n_atoms) => this%B
!
      v1 = [one, -one, one]
      v2 = [-one, one, one]  
!
      do angle = 1, this%n_angles
!
         m = this%angles(angle)%atoms(1)
         o = this%angles(angle)%atoms(2)
         n = this%angles(angle)%atoms(3)
!
         do q = 1, 3
!
            u(q) = this%R(q,m) - this%R(q,o)
            v(q) = this%R(q,n) - this%R(q,o)
!
         enddo 
!
         norm_u = norm_R3(u)
         norm_v = norm_R3(v)
!
         u = u/norm_u
         v = v/norm_v
!
         if (dot_R3(u,v) .ne. one) then ! u not parallel with v
!
            w = cross_product_R3(u,v)
!
         elseif (dot_R3(u,v1) .ne. norm_R3(v1)) then ! u not parallel with v1
!
            w = cross_product_R3(u, v1)
!
         else ! u parallel with v1; cross with v2
!
            w = cross_product_R3(u, v2)
!
         endif
!
         norm_w = norm_R3(w)
         w = w/norm_w
!
         uxw = cross_product_R3(u,w)
         wxv = cross_product_R3(w,v)
!
         do a = 1, this%n_atoms
            do i = 1, 3
!
               B_p(this%n_bonds + angle, i, a) = (delta(a,m) - delta(a,o))*uxw(i)/norm_u &
                                               + (delta(a,n) - delta(a,o))*wxv(i)/norm_v
!
            enddo
         enddo
!
      enddo
!
   end subroutine construct_Wilson_B_bond_angles
!
!
   subroutine construct_Wilson_B_dihedral_angles(this)
!!
!!    Construct Wilson B (dihedral angles)
!!    Written by Eirik F. Kjønstad, 2021
!!
      use math_utilities, only: delta
!
      implicit none 
!
      class(redundant_internal_coords), target, intent(inout) :: this 
!
      real(dp), dimension(:,:,:), pointer :: B_p 
!
      integer :: dihedral, m, o, p, n, a, i 
!
      real(dp), dimension(3) :: u, v, w, uxw, vxw
!
      real(dp) :: norm_u, norm_v, norm_w, cos_phiu, cos_phiv, sin_phiu, sin_phiv
!
      integer :: q 
!
      B_p(1:this%n_internal, 1:3, 1:this%n_atoms) => this%B
!
      do dihedral = 1, this%n_dihedrals
!  
         m = this%dihedrals(dihedral)%atoms(1)
         o = this%dihedrals(dihedral)%atoms(2)
         p = this%dihedrals(dihedral)%atoms(3)
         n = this%dihedrals(dihedral)%atoms(4)
!
         do q = 1, 3
!
            u(q) = this%R(q,m) - this%R(q,o)
            v(q) = this%R(q,n) - this%R(q,p)
            w(q) = this%R(q,p) - this%R(q,o)
!
         enddo 
!
         norm_u = norm_R3(u)
         norm_v = norm_R3(v)
         norm_w = norm_R3(w)
!
         u = u/norm_u
         v = v/norm_v
         w = w/norm_w
!
         uxw = cross_product_R3(w,u)
         vxw = cross_product_R3(v,w)
!
         cos_phiu = dot_R3(u,w)
         sin_phiu = sqrt(one - cos_phiu**2)
!
         cos_phiv = -dot_R3(v,w)
         sin_phiv = sqrt(one - cos_phiv**2)
!
         do a = 1, this%n_atoms
            do i = 1, 3
!
               B_p(this%n_bonds + this%n_angles + dihedral, i, a) = -( &
                    + (delta(a,m) - delta(a,o))*uxw(i)/(norm_u*sin_phiu**2) &
                    - (delta(a,p) - delta(a,n))*vxw(i)/(norm_v*sin_phiv**2) &
                    + (delta(a,o) - delta(a,p))*( uxw(i)*cos_phiu/(norm_w*sin_phiu**2) &
                                                - vxw(i)*cos_phiv/(norm_w*sin_phiv**2) ))
!
            enddo
         enddo
!
      enddo
!
   end subroutine construct_Wilson_B_dihedral_angles
!
!
   subroutine pseudoinvert_Wilson_B_and_construct_P(this)
!!
!!    Pseudoinvert Wilson B and construct P 
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Constructs the pseudoinverse of B and P = G G_inv, 
!!    where G = B B^T and G_inv is the pseudoinverse of G
!!
      use array_utilities, only: calculate_pseudoinverse
!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      real(dp), dimension(this%n_internal, this%n_internal) :: G 
      real(dp), dimension(this%n_internal, this%n_internal) :: G_inv 
!
      call dgemm('N', 'T',          &
                  this%n_internal,  &
                  this%n_internal,  &
                  this%n_cartesian, &
                  one,              &
                  this%B,           &
                  this%n_internal,  &
                  this%B,           &
                  this%n_internal,  &
                  zero,             &
                  G,                &
                  this%n_internal)
!
      call calculate_pseudoinverse(G_inv,             &
                                   G,                 &
                                   this%n_internal,   &
                                   this%n_internal,   &
                                   this%pseudoinverse_threshold)
!
      call dgemm('N', 'N',          &
                  this%n_internal,  &
                  this%n_internal,  &
                  this%n_internal,  &
                  one,              &
                  G,                &
                  this%n_internal,  &
                  G_inv,            &
                  this%n_internal,  &
                  zero,             &
                  this%P,           &
                  this%n_internal)
!
      call output%print_matrix('v', 'Wilson P', this%P, this%n_internal, this%n_internal)
!
      call dgemm('T', 'N',             &
                  this%n_cartesian,    &
                  this%n_internal,     &
                  this%n_internal,     &
                  one,                 &
                  this%B,              &
                  this%n_internal,     &
                  G_inv,               &
                  this%n_internal,     &
                  zero,                &
                  this%B_inverse,      &
                  this%n_cartesian)
!
      call output%print_matrix('v', 'Wilson B_inv', &
         this%B_inverse, this%n_cartesian, this%n_internal, fs='(f19.12)')
!
   end subroutine pseudoinvert_Wilson_B_and_construct_P
!
!
   subroutine get_approximate_diagonal_Hessian(this, H_diagonal)
!!
!!    Get approximate diagonal Hessian
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Returns the "simple" Hessian diagonal estimate from Bakken and Helgaker (2002).
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      real(dp), dimension(this%n_internal), intent(out) :: H_diagonal
!
      integer :: k
!
      do k = 1, this%n_bonds 
!
         H_diagonal(k) = 0.5d0
!
      enddo
!
      do k = this%n_bonds + 1, this%n_bonds + this%n_angles
!
         H_diagonal(k) = 0.2d0
!
      enddo
!
      do k = this%n_bonds + this%n_angles + 1, this%n_internal
!
         H_diagonal(k) = 0.1d0
!
      enddo
!
   end subroutine get_approximate_diagonal_Hessian
!
!
   subroutine destructor(this)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      type(redundant_internal_coords) :: this 
!
      if (allocated(this%R)) call mem%dealloc(this%R, 3, this%n_atoms)
      if (allocated(this%Z)) call mem%dealloc(this%Z, this%n_atoms)
!
      if (allocated(this%B)) call mem%dealloc(this%B, this%n_internal, this%n_cartesian)
      if (allocated(this%B_inverse)) call mem%dealloc(this%B_inverse, this%n_cartesian, this%n_internal)
!
      if (allocated(this%P)) call mem%dealloc(this%P, this%n_internal, this%n_internal)
!
      if (allocated(this%q)) call mem%dealloc(this%q, this%n_internal)
      if (allocated(this%prev_q)) call mem%dealloc(this%prev_q, this%n_internal)
!
      if (allocated(this%bond_matrix)) call mem%dealloc(this%bond_matrix, this%n_atoms, this%n_atoms)
!
   end subroutine destructor
!
!
   subroutine Wilson_project_vector(this, q)
!!
!!    Wilson project vector
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    q = P q 
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      real(dp), dimension(this%n_internal), intent(inout) :: q 
!
      real(dp), dimension(this%n_internal) :: q0 
!
      q0 = q 
!
      call dgemm('N', 'N',          &
                  this%n_internal,  &
                  1,                &
                  this%n_internal,  &
                  one,              &
                  this%P,           &
                  this%n_internal,  &
                  q0,               &
                  this%n_internal,  &
                  zero,             &
                  q,                &
                  this%n_internal)
!
   end subroutine Wilson_project_vector
!
!
   subroutine Wilson_project_matrix(this, H)
!!
!!    Wilson project vector
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    H = P H P
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      real(dp), dimension(this%n_internal, this%n_internal), intent(inout) :: H 
!
      real(dp), dimension(this%n_internal, this%n_internal) :: H0, H0P
!
      H0 = H 
!
      call dgemm('N', 'N',          &
                  this%n_internal,  &
                  this%n_internal,  &
                  this%n_internal,  &
                  one,              &
                  H0,               &
                  this%n_internal,  &
                  this%P,           &
                  this%n_internal,  &
                  zero,             &
                  H0P,              &
                  this%n_internal)
!
      call dgemm('N', 'N',          &
                  this%n_internal,  &
                  this%n_internal,  &
                  this%n_internal,  &
                  one,              &
                  this%P,           &
                  this%n_internal,  &
                  H0P,              &
                  this%n_internal,  &
                  zero,             &
                  H,                &
                  this%n_internal)
!
   end subroutine Wilson_project_matrix
!
!
   subroutine calculate_internal_gradient(this, g_x, g_q)
!!
!!    Calculate internal gradient
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    g_q = B_inv g_x
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      real(dp), dimension(this%n_cartesian), intent(in) :: g_x 
      real(dp), dimension(this%n_internal), intent(out) :: g_q 
!
      call dgemm('T', 'N',             &
                  this%n_internal,     &
                  1,                   &
                  this%n_cartesian,    &
                  one,                 &
                  this%B_inverse,      &
                  this%n_cartesian,    &
                  g_x,                 &
                  this%n_cartesian,    &
                  zero,                &
                  g_q,                 &
                  this%n_internal)
!
   end subroutine calculate_internal_gradient
!
!
   subroutine compute_next_geometry(this, x, s_q)
!!
!!    Compute next geometry
!!    Written by Eirik F. Kjønstad and Anna-Kristina Schnack-Petersen, 2021
!!
!!    Updates the cartesian geometry x with the internal coordinate step s_q by
!!    performing a back-transformation through an iterative procedure due to 
!!    Bakken and Helgaker (2002).
!!
      use array_utilities, only: get_root_mean_square
!
      implicit none 
!
      class(redundant_internal_coords), intent(inout) :: this 
!
      real(dp), dimension(this%n_cartesian), intent(inout) :: x 
!
      real(dp), dimension(this%n_internal), intent(in) :: s_q
!
      real(dp), dimension(this%n_cartesian) :: xk_minus_1, xk, dx
      real(dp), dimension(this%n_internal) :: q_original, qk, ddq, dq, s_q_proj
!
      integer :: iteration
      real(dp) :: norm_dx
      logical :: converged 
!
      call this%get_internals(this%prev_q) ! Save current internals as previous internals
!
      call output%printf('m', 'Converting internal step to cartesian step:', fs='(/t6,a/)')
!
      call this%get_internals(q_original)
!
      s_q_proj = s_q 
      call this%Wilson_project_vector(s_q_proj)
      call this%remove_multiples_of_360(s_q_proj)
!
      xk_minus_1 = x
!
      call this%perform_cartesian_step(xk, xk_minus_1, s_q)
!
      call this%set_geometry(xk)
      call this%get_internals(qk)
!
      dx = xk - xk_minus_1
      dq = qk - q_original
!
      call this%remove_multiples_of_360(dq)
!
      ddq = s_q_proj - dq
!
      norm_dx = get_root_mean_square(dx, this%n_cartesian)
!
      call output%printf('m', 'Iteration     RMS error in cartesians', fs='(t6,a)')
      call output%print_separator('n', 37,'-', fs='(t6,a)')
!
      iteration = 0
!
      converged = norm_dx .lt. this%internal_to_cartesian_threshold
!
      do while (.not. converged .and. iteration .lt. this%internal_to_cartesian_max_iterations)
!
         iteration = iteration + 1
!
         call output%printf('m', '(i0)             (e11.5)', &
                        ints=[iteration], reals=[norm_dx], fs='(t6,a)')
!
         xk_minus_1 = xk
!
         call this%Wilson_project_vector(ddq)
!
         call this%perform_cartesian_step(xk, xk_minus_1, ddq)
!
         call this%set_geometry(xk)
!
         dx = xk - xk_minus_1 
!
         norm_dx = get_root_mean_square(dx, this%n_cartesian)
         converged = norm_dx .lt. this%internal_to_cartesian_threshold
!
         if (.not. converged) then 
!
            call this%get_internals(qk)
!
            dq = qk - q_original
            call this%remove_multiples_of_360(dq)
!
            ddq = s_q_proj - dq
!
         endif 
!
      enddo
!
      call output%printf('m', '(i0)             (e11.5)', &
                        ints=[iteration + 1], reals=[norm_dx], fs='(t6,a)')
!
      call output%print_separator('n', 37,'-', fs='(t6,a)')
!
      if (converged) then 
!
         call output%printf('m', 'Converged in (i0) iterations!', ints=[iteration + 1], fs='(t6,a)')
!
      else 
!
         call output%printf('m', 'Could not find cartesian step from internal coordinate step! &
            & Using initial guess instead.')
!
         xk_minus_1 = x 
!
         call this%set_geometry(xk_minus_1)
         call this%perform_cartesian_step(xk, xk_minus_1, s_q)
         call this%set_geometry(xk)
!
      endif 
!
      x = xk
!
      call this%get_internals(this%q) ! Save new current internals
!
   end subroutine compute_next_geometry
!
!
   subroutine perform_cartesian_step(this, x1, x0, q)
!!
!!    Perform cartesian step
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    x1 = x0 + B_inverse q
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this 
!
      real(dp), dimension(this%n_cartesian), intent(out) :: x1 
      real(dp), dimension(this%n_cartesian), intent(in) :: x0 
      real(dp), dimension(this%n_internal), intent(in) :: q 
!
      real(dp), dimension(this%n_cartesian) :: dx 
!
      x1 = x0
!
      call dgemm('N', 'N',             &
                  this%n_cartesian,    &
                  1,                   &
                  this%n_internal,     &
                  one,                 &
                  this%B_inverse,      &
                  this%n_cartesian,    &
                  q,                   &
                  this%n_internal,     &
                  zero,                &
                  dx,                  &
                  this%n_cartesian)
!
      x1 = x0 + dx
!
   end subroutine perform_cartesian_step
!
!
   subroutine remove_multiples_of_360(this, dq)
!!
!!    Remove multiples of 360
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Avoid redundancy in angles by restricting them to [-pi, pi].
!!
      implicit none 
!
      class(redundant_internal_coords), intent(in) :: this
!
      real(dp), dimension(this%n_internal), intent(inout) :: dq 
!
      integer :: k
!
      do k = this%n_bonds + 1, this%n_internal
!  
         if (dq(k) .lt. -pi) then 
!
!           -pi .q_1.. 0 ...q. pi
!
            call output%printf('v', 'Angle less than -pi: (f19.12). Will add 2*pi.', reals=[dq(k)])
!
            dq(k) = dq(k) + two*pi
!
         elseif (dq(k) .gt. pi) then 
!
!           -pi .q.. 0 ...q_1. pi
!
            call output%printf('v', 'Angle greater than pi: (f19.12). Will subtract 2*pi.', reals=[dq(k)])
!
            dq(k) = dq(k) - two*pi
!
         endif
!
      enddo  
!
   end subroutine remove_multiples_of_360

!
!
end module redundant_internal_coords_class
