/*
/
/ 	Libint initialization of engines and basis set
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include "libint_initialization.h"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

using namespace libint2;
using namespace std;

BasisSet basis;
extern BasisSet basis;

Engine electronic_repulsion;
extern Engine electronic_repulsion;

void initialize_basis(){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet temporary("cc-pVDZ", atoms);
	cout.clear();

	basis = temporary;

	finalize();

}

void initialize_coulomb(){

	initialize();

	Engine temporary(Operator::coulomb, basis.max_nprim(), basis.max_l());
	electronic_repulsion = temporary;

	finalize();

}
