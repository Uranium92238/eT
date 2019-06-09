//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2019 the authors of eT
//
//  eT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  eT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------
/*
 
   Globals
   Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
 
   Contains the only global variables allowed in the integral calculations:
   the basis set, atomic information, and integral engines.
 
*/
#include "eT_basis.h"

extern eTBasis basis;                                        // The basis set used throughout
extern vector<libint2::Engine> electronic_repulsion_engines; // The ERI engines vector
extern libint2::Engine electronic_repulsion_1der;            // The ERI first derivative engine
extern vector<libint2::Engine> kinetic;                      // The kinetic energy engines vector
extern libint2::Engine kinetic_1der;						       // The kinetic energy first derivative
extern vector<libint2::Engine> nuclear;                      // The nuclear attraction engine vector
extern libint2::Engine nuclear_1der;						       // The nuclear energy first derivative
extern libint2::Engine overlap_1der;                         // The overlap first derivative engine
extern vector<libint2::Engine> overlap;                      // The overlap engine vector
extern vector<libint2::Engine> dipole;                       // The dipole engine
extern vector<libint2::Engine> quadrupole;                   // The quadrupole engine
extern vector<libint2::Atom> atoms;                          // Atoms vector
extern vector<int> shell2atom;                               // Shell center vector
