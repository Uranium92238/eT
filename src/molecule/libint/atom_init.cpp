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

int get_n_shells_on_atom(int *i){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet obs("cc-pVDZ", atoms);
	cout.clear();

	int num_aos = 0;
	num_aos = obs.nbf();

 	ofstream outputFile("debug.txt");

	auto a2s_list = obs.atom2shell(atoms); // Vector of vectors

	for (auto j = 0; j < atoms.size(); j++){

		if (j != (*i-1)) {
			// Do nothing; wrong atom
		} else {
			int n_shells = 0;
			n_shells = a2s_list[j].size();
			outputFile << "The number of shells is: " << n_shells << ' or: ' << a2s_list[j].size() << endl;
			return n_shells;
		}

	}

	finalize();

}
