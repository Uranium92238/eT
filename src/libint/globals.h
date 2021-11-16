//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2021 the authors of eT
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

// Global variables declared in this file

// Basis set, atoms, and AOs
extern eTBasis basis;
extern libint2::BasisSet ri_basis;
extern vector<libint2::Atom> atoms;                          
extern vector<vector<int>> atom_to_shell_list;
extern vector<int> shell_to_first_ao;       

/* Vectors of Libint engines. 
   Length of vectors is equal to number of threads 
*/
extern vector<libint2::Engine> electronic_repulsion;
extern vector<libint2::Engine> electronic_repulsion_2c;
extern vector<libint2::Engine> electronic_repulsion_3c;
extern vector<libint2::Engine> electronic_repulsion_1der;     
extern vector<libint2::Engine> kinetic;                       
extern vector<libint2::Engine> kinetic_1der;
extern vector<libint2::Engine> nuclear;                      
extern vector<libint2::Engine> nuclear_1der;                 
extern vector<libint2::Engine> overlap; 						  
extern vector<libint2::Engine> overlap_1der;                       
extern vector<libint2::Engine> dipole;                        
extern vector<libint2::Engine> quadrupole;                    
extern vector<libint2::Engine> potential;                     
extern vector<libint2::Engine> coulomb_external_charges;      
extern vector<libint2::Engine> coulomb_external_unit_charges; 
