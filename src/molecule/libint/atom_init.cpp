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

void get_shell_numbers(int *atom, int *sn){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet obs("cc-pVDZ", atoms);
	cout.clear();

	auto a2s_list = obs.atom2shell(atoms); // Vector of vectors

	for (auto j = 0; j < a2s_list[*atom-1].size(); j++){ // loop over shells on atom

		auto the_shell = a2s_list[*atom-1][j];

		*(sn + j) = the_shell + 1;
		cout << "The " << j << " th shell on atom " << *atom-1 << " is " << the_shell + 1 << endl;

	}

	finalize();

}

void get_first_ao_in_shells(int *atom, int *faois){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet obs("cc-pVDZ", atoms);
	cout.clear();

	auto a2s_list = obs.atom2shell(atoms); // Vector of vectors
	auto shell2bf = obs.shell2bf(); // shell2bf[0] -> first AO index of shell 0

	for (auto j = 0; j < a2s_list[*atom-1].size(); j++){ // loop over shells on atom

		auto the_shell = a2s_list[*atom-1][j];
		auto the_first_ao = shell2bf[the_shell];

		*(faois + j) = the_first_ao + 1;
	//	cout << "The first AO index in shell " << the_shell << " is " << the_first_ao << "." << endl;

	}

	finalize();

}

void get_n_basis_in_shells(int *atom, int *nbis){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet obs("cc-pVDZ", atoms);
	cout.clear();

	auto a2s_list = obs.atom2shell(atoms); // Vector of vectors

	for (auto j = 0; j < a2s_list[*atom-1].size(); j++){

//		cout << "the " << j << "th shell on atom " << *atom << " is shell nr. " << a2s_list[*atom-1][j] << endl;
		auto n = obs[a2s_list[*atom-1][j]].size();
//		cout << "and the number of basis functions in the shell is " << n << "." << endl;

		*(nbis + j) = n;

	}

	finalize();

}

void get_n_shells_on_atoms(int *nsoa){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet obs("cc-pVDZ", atoms);
	cout.clear();

	auto a2s_list = obs.atom2shell(atoms); // Vector of vectors
	int n_shells = 0;

	for (auto j = 0; j < atoms.size(); j++){

		n_shells = a2s_list[j].size();
		cout << "num shells? " << n_shells << endl;
		*(nsoa + j) = n_shells;

	}

	finalize();

}
