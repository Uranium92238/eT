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
module periodic_table
!
!!
!!    Periodic table module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!    Extended from Z=20 to Z=118 by Susi Lehtola, 2020
!!
!
   use kinds
!
   implicit none
!
   integer, parameter :: size_periodic_table = 118
   character(len=2), dimension(size_periodic_table) :: atomic_symbol = &
        [ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
        & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
        & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
        & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
        & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
        & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
        & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
        & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
        & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
        & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
        & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',&
        & 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
!
end module periodic_table
