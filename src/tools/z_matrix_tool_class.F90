!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module z_matrix_tool_class
!
!!
!!    Internal coordinate (Z-matrix) class module
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Always use the following connections:
!!
!!       R_i-1 is the distance connection,
!!       R_i-2 is the angle connection,
!!       R_i-3 is the dihedral connection
!!
!
   use parameters
   use math_utilities, only: cross_product_R3, dot_R3, norm_R3
   use global_out, only: output
   use memory_manager_class, only: mem
!
   implicit none
!
   type :: z_matrix_tool
!
      integer :: n_atoms
!
      character(len=2), dimension(:), allocatable  :: symbols
!
      integer, dimension(:), allocatable           :: distance_connections
      integer, dimension(:), allocatable           :: angle_connections
      integer, dimension(:), allocatable           :: dihedral_connections
!
      real(dp), dimension(:), allocatable          :: distances
      real(dp), dimension(:), allocatable          :: angles
      real(dp), dimension(:), allocatable          :: dihedrals
!
   contains
!
      procedure :: initialize &
                => initialize_z_matrix_tool
!
      procedure :: cleanup_z_matrix &
                => cleanup_z_matrix_z_matrix_tool
!
      procedure :: construct &
                => construct_z_matrix_tool
!
      procedure :: print_ &
                => print_z_matrix_tool
!
      procedure :: convert_to_cartesian &
                => convert_to_cartesian_z_matrix_tool
!
   end type z_matrix_tool
!
   interface z_matrix_tool
!
      procedure :: new_z_matrix_tool
!
   end interface z_matrix_tool
!
contains
!
!
   function new_z_matrix_tool(n_atoms) result(z_matrix)
!!
!!    New internal coordinate converter
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      type(z_matrix_tool) :: z_matrix
!
      integer, intent(in) :: n_atoms
!
      z_matrix%n_atoms = n_atoms
!
   end function new_z_matrix_tool
!
!
   subroutine initialize_z_matrix_tool(z_matrix)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(z_matrix_tool) :: z_matrix
!
      allocate(z_matrix%symbols(z_matrix%n_atoms))
!
      call mem%alloc(z_matrix%distances, z_matrix%n_atoms)
      call mem%alloc(z_matrix%angles, z_matrix%n_atoms)
      call mem%alloc(z_matrix%dihedrals, z_matrix%n_atoms)
!
      z_matrix%distances = zero
      z_matrix%angles    = zero
      z_matrix%dihedrals = zero
!
   end subroutine initialize_z_matrix_tool
!
!
   subroutine cleanup_z_matrix_z_matrix_tool(z_matrix)
!!
!!    Cleanup Z-matrix
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(z_matrix_tool) :: z_matrix
!
      deallocate(z_matrix%symbols)
!
      call mem%dealloc(z_matrix%distances, z_matrix%n_atoms)
      call mem%dealloc(z_matrix%angles, z_matrix%n_atoms)
      call mem%dealloc(z_matrix%dihedrals, z_matrix%n_atoms)
!
   end subroutine cleanup_z_matrix_z_matrix_tool
!
!
   subroutine construct_z_matrix_tool(z_matrix, R, symbols)
!!
!!    Construct Z-matrix
!!    Written by Sarai D. Folkestad, 2020
!!
      use array_utilities, only: get_euclidean_distance
!
      implicit none
!
      class(z_matrix_tool) :: z_matrix
!
      real(dp), dimension(3, z_matrix%n_atoms), intent(in)       :: R
      character(len=2), dimension(z_matrix%n_atoms), intent(in)  :: symbols
!
      real(dp), dimension(3) :: R_kj, R_ij, R_jk, R_lk
      real(dp), dimension(3) :: A, B
!
      real(dp) :: cos_phi, sin_phi
!
      integer :: i
!
      z_matrix%symbols = symbols
!
!     Calculate distances
!
      do i = 2, z_matrix%n_atoms
!
         z_matrix%distances(i) = get_euclidean_distance(R(:, i), R(:, i - 1), 3)
!
      enddo
!
!     Calculate angles
!
!        r sin(theta) = |v1 x v2|
!        r cos(theta) = v1 . v2
!
!        theta = atan2(r sin(theta), r cos(theta))
!
!        v1 = R_ij
!        v2 = R_kj
!
      do i = 3, z_matrix%n_atoms
!
         R_ij = R(:, i - 2) - R(:, i - 1)
         R_kj = R(:, i) - R(:, i - 1)
!
         cos_phi = dot_R3(R_ij, R_kj)
         sin_phi = norm_R3(cross_product_R3(R_ij, R_kj))
!
         z_matrix%angles(i) = atan2(sin_phi, cos_phi)
!
      enddo
!
!     Calculate dihedrals
!
!        See equations (3) and (4)
!        A. Blondel and M. Karplus, J. Comput. Chem., 17(9), 1132-1141, (1996).
!
      do i = 4, z_matrix%n_atoms
!
         R_ij = R(:, i - 3) - R(:, i - 2) ! F, J. Comput. Chem., 17(9), 1132-1141, (1996)
         R_jk = R(:, i - 2) - R(:, i - 1) ! G, J. Comput. Chem., 17(9), 1132-1141, (1996)
         R_lk = R(:, i) - R(:, i - 1)     ! H, J. Comput. Chem., 17(9), 1132-1141, (1996)
!
         A = cross_product_R3(R_ij, R_jk)
         B = cross_product_R3(R_lk, R_jk)

         cos_phi = dot_R3(A, B)/(norm_R3(A)*norm_R3(B))
         sin_phi = dot_R3(cross_product_R3(B, A), R_jk)/(norm_R3(A)*norm_R3(B)*norm_R3(R_jk))
!
         z_matrix%dihedrals(i) = atan2(sin_phi, cos_phi)
!
         if (z_matrix%dihedrals(i) .lt. zero ) &
            z_matrix%dihedrals(i) = z_matrix%dihedrals(i) + two*pi
!
      enddo
!
!
   end subroutine construct_z_matrix_tool
!
!
   subroutine print_z_matrix_tool(z_matrix, file_)
!!
!!    Print Z-matrix
!!    Written by Sarai D. Folkestad
!!
!
      use output_file_class, only: output_file
!
      implicit none
!
      class(z_matrix_tool) :: z_matrix
!
      integer :: i, line_length
      integer :: distance_connect, angle_connect, dihedral_connect
!
      type(output_file), intent(inout) :: file_
!
      line_length = 78
!
      distance_connect = 1
      angle_connect = 1
      dihedral_connect = 1
!
!     Print header for coordinate print
!
      line_length = 78
!
      call file_%print_separator(pl='m', symbol='=', n=line_length, fs='(/t6,a)')
!
      call file_%printf('m', 'Z-matrix (angstrom, degrees)', fs='(t32,a)')
!
      call file_%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
!
      call file_%printf('m',&
            '             Distances             Angles          Dihedrals', &
            fs='(t9,a)')
!
      call file_%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
!
      call file_%printf('m', ' (a2)', chars=[z_matrix%symbols(1)], fs='(t9,a)')

      if (z_matrix%n_atoms .gt. 1) &
         call file_%printf('m', ' (a2) (i4) (f18.12)',   &
                           chars=[z_matrix%symbols(2)],  &
                           reals=[z_matrix%distances(2)],&
                           ints=[distance_connect], fs='(t9,a)')
!
      distance_connect = distance_connect + 1
!
      if (z_matrix%n_atoms .gt. 2) &
         call file_%printf('m', ' (a2) (i4) (f18.12) (i4) (f10.5)',&
                      chars=[z_matrix%symbols(3)],                 &
                      reals=[z_matrix%distances(3),                &
                      radians_to_degrees*z_matrix%angles(3)],      &
                      ints=[distance_connect, angle_connect], fs='(t9,a)')
!
      distance_connect = distance_connect + 1
      angle_connect = angle_connect + 1
!
      do i = 4, z_matrix%n_atoms
!
         call file_%printf('m', ' (a2) (i4) (f18.12) (i4) (f10.5) (i4) (f10.5)', &
            chars=[z_matrix%symbols(i)],                           &
            reals=[z_matrix%distances(i),                          &
                   radians_to_degrees*z_matrix%angles(i),          &
                   radians_to_degrees*z_matrix%dihedrals(i)],      &
            ints=[distance_connect, angle_connect, dihedral_connect], fs='(t9,a)')
!
            distance_connect = distance_connect + 1
            angle_connect = angle_connect + 1
            dihedral_connect = distance_connect + 1
!
      enddo
!
      call file_%print_separator(pl='m', symbol='=', n=line_length, fs='(t6,a)')
!
   end subroutine print_z_matrix_tool
!
!
   subroutine convert_to_cartesian_z_matrix_tool(z_matrix, R)
!!
!!    Convert to cartesian
!!    Written by Sarai D. Folkestad
!!
!!    Converts z matrix to cartesian coordinates
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(z_matrix_tool) :: z_matrix
!
      real(dp), dimension(3, z_matrix%n_atoms), intent(out) :: R
!
      integer :: i
!
      real(dp), dimension(3) :: A, G, F, e_1, e_2, e_3
      real(dp) :: cos_phi, cos_theta, sin_theta, sin_phi, x, y, z
!
      call zero_array(R, z_matrix%n_atoms*3)
!
!     First atom is placed at the origin
!     Place second atom on positive x-axis
!
      R(1, 2) = z_matrix%distances(2)
!
!     Place third atom in xy plane
!
      R(1, 3) = R(1, 2) + z_matrix%distances(3)*cos(pi - z_matrix%angles(3))
      R(2, 3) = R(2, 2) + z_matrix%distances(3)*sin(pi - z_matrix%angles(3))
!
!     Place remaining atoms
!
      do i = 4, z_matrix%n_atoms
!
!        Define coordinate system in R_n-1
!
         F = R(:, i - 3) - R(:, i - 2)
         G = R(:, i - 2) - R(:, i - 1)
         e_1 = G/norm_R3(G) ! G/|G|
!
         A = cross_product_R3(F, G)
         e_2 = A/norm_R3(A) ! (F x G)/|F x G|
!
         e_3 = cross_product_R3(e_1, e_2) ! G x (F x G) / |G| |F x G|
!
         sin_theta = sin(z_matrix%angles(i))
         cos_theta = cos(z_matrix%angles(i))
         sin_phi = sin(z_matrix%dihedrals(i))
         cos_phi = cos(z_matrix%dihedrals(i))
!
!        Projections along the unit vectors e_1, e_2, e_3
!        multiplied by the distance from R_i-1 to R_i
!
         x = z_matrix%distances(i)*cos_theta
         y = z_matrix%distances(i)*sin_phi*sin_theta
         z = z_matrix%distances(i)*cos_phi*sin_theta
!
!        Position relative to R_i-1
!
         R(1, i) = e_1(1)*x + e_2(1)*y + e_3(1)*z
         R(2, i) = e_1(2)*x + e_2(2)*y + e_3(2)*z
         R(3, i) = e_1(3)*x + e_2(3)*y + e_3(3)*z
!
!        Prosition relative to origin (0,0,0)
!
         R(1, i) = R(1, i) + R(1, i - 1)
         R(2, i) = R(2, i) + R(2, i - 1)
         R(3, i) = R(3, i) + R(3, i - 1)
!
      enddo

  end subroutine convert_to_cartesian_z_matrix_tool
!
!
end module z_matrix_tool_class
