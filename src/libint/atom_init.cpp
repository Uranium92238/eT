/*
/
/ 	Atom initialization routines
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include "atom_init.h"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace libint2;
using namespace std;

#include "globals.h"

void get_shell_numbers(int *atom, int *sn){

	auto a2s_list = basis.atom2shell(atoms); // Vector of vectors

	for (auto j = 0; j < a2s_list[*atom-1].size(); j++){ // loop over shells on atom

		auto the_shell = a2s_list[*atom-1][j];

		*(sn + j) = the_shell + 1;

	}

}

void get_first_ao_in_shells(int *atom, int *faois){

	auto a2s_list = basis.atom2shell(atoms); // Vector of vectors
	auto shell2bf = basis.shell2bf(); // shell2bf[0] -> first AO index of shell 0
	for (auto j = 0; j < a2s_list[*atom-1].size(); j++){ // loop over shells on atom

		auto the_shell = a2s_list[*atom-1][j];
		auto the_first_ao = shell2bf[the_shell];

		*(faois + j) = the_first_ao + 1;
	}

}

void get_n_basis_in_shells(int *atom, int *nbis){

	auto a2s_list = basis.atom2shell(atoms); // Vector of vectors

	for (auto j = 0; j < a2s_list[*atom-1].size(); j++){

		auto n = basis[a2s_list[*atom-1][j]].size();

		*(nbis + j) = n;

	}

}

void get_n_shells_on_atoms(int *nsoa){

	auto a2s_list = basis.atom2shell(atoms); // Vector of vectors
	int n_shells = 0;

	for (auto j = 0; j < atoms.size(); j++){

		n_shells = a2s_list[j].size();

		*(nsoa + j) = n_shells;

	}

}
