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

 	Atom initialization routines
 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018

*/
#include "atom_init.h"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>

using namespace libint2;
using namespace std; 

void get_shell_numbers(int *atom, int *n_shells){

	auto a2s = basis.atom2shell(atoms); 
	for (std::size_t j = 0; j < a2s[*atom - 1].size(); j++){ 

		auto shell = a2s[*atom - 1][j];

		n_shells[j] = shell + 1;

	}

}

void get_first_ao_in_shells(int *atom, int *first_ao_in_shell){

	auto a2s = basis.atom2shell(atoms); 
	vector<int> stofao = basis.shell2bf(); 
	for (std::size_t j = 0; j < a2s[*atom-1].size(); j++){ 

		auto shell = a2s[*atom - 1][j];
		auto first_ao = stofao[shell];

		first_ao_in_shell[j] = first_ao + 1;

	}

}

void get_n_aos_in_shell(int *atom, int *n_aos_in_shell){

	auto a2s = basis.atom2shell(atoms); 
	for (std::size_t j = 0; j < a2s[*atom - 1].size(); j++){

		auto n = basis[a2s[*atom - 1][j]].size();

		n_aos_in_shell[j] = n;

	}

}

void get_n_shells_on_atom(int *atom, int *n_shells){

	auto a2s = basis.atom2shell(atoms); 
		*n_shells = a2s[*atom - 1].size();

}

void initialize_atom_to_shell_list(){

 	// Don't know why this doesnt work, but we should avoid recalculating these -- see above
	vector<vector<int>> atom_to_shell_list = basis.atom2shell(atoms); 

}

void initialize_shell_to_first_ao(){

 	// Don't know why this doesnt work, but we should avoid recalculating these -- see above
	vector<int> shell_to_first_ao = basis.shell2bf(); 

}
