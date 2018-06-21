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

void get_n_shells_on_atom(int *nsoa){

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
