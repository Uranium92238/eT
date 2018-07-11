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

Engine kinetic;
extern Engine kinetic;

Engine nuclear;
extern Engine nuclear;

Engine overlap;
extern Engine overlap;

vector<Atom> atoms;
extern vector<Atom> atoms;

void initialize_basis(){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> temporary_atoms = read_dotxyz(input_file);
	atoms = temporary_atoms;

	cout.setstate(ios_base::failbit);
	BasisSet temporary("cc-pVQZ", atoms);
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

void initialize_kinetic(){

	initialize();

	Engine temporary(Operator::kinetic, basis.max_nprim(), basis.max_l());
	kinetic = temporary;

	finalize();

}

void initialize_nuclear(){

	initialize();

	Engine temporary(Operator::nuclear, basis.max_nprim(), basis.max_l());
	nuclear = temporary;

	finalize();

}

void initialize_overlap(){

	initialize();

	Engine temporary(Operator::overlap, basis.max_nprim(), basis.max_l());
	overlap = temporary;

	finalize();

}
