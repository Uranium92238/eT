/*
/
/ 	Libint initialization of engines and basis set
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include "libint_initialization.h"
#include "eT_basis.h"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

using namespace libint2;
using namespace std;

eTBasis basis;
extern eTBasis basis;

// BasisSet basis;
// extern BasisSet basis;

vector<Engine> electronic_repulsion_engines(omp_get_max_threads());
extern vector<Engine> electronic_repulsion_engines;

vector<Engine> kinetic(omp_get_max_threads());
extern vector<Engine> kinetic;

vector<Engine> nuclear(omp_get_max_threads());
extern vector<Engine> nuclear;

// Engine kinetic;
// extern Engine kinetic;

// Engine nuclear;
// extern Engine nuclear;

Engine overlap;
extern Engine overlap;

vector<Atom> atoms;
extern vector<Atom> atoms;

vector<long> shell2bf_g;
extern vector<long> shell2bf_g;

void initialize_atoms(char *name){

    string xyzfilename(strcat(name,".xyz"));

    ifstream input_file(xyzfilename);
	vector<Atom> temporary_atoms = read_dotxyz(input_file);
	atoms = temporary_atoms;

}

void initialize_basis(char *basisset, char *filename){
	
	string xyzfilename(strcat(filename,".xyz"));

	//cout << xyzfilename << endl;
	//cout << basisset << endl;

    ifstream input_file(xyzfilename);

	vector<Atom> temporary_atoms = read_dotxyz(input_file);

	cout.setstate(ios_base::failbit);
	BasisSet temporary(basisset, temporary_atoms);
	cout.clear();

	basis.add(temporary);

}


void initialize_libint(){

	initialize();

}

void finalize_libint(){

	finalize();

}


void initialize_coulomb(){


	Engine temporary(Operator::coulomb, basis.max_nprim(), basis.max_l());
	temporary.set_precision(1.0e-25);

	for (int i = 0; i != omp_get_max_threads(); i++){
		electronic_repulsion_engines[i] = temporary; // One engine per thread
	}

	shell2bf_g = basis.shell2bf(); // testing something

}

void set_coulomb_precision(double *prec){

	for (int i = 0; i != omp_get_max_threads(); i++){
		electronic_repulsion_engines[i].set_precision(*prec); // One engine per thread
	}

}

// void initialize_kinetic(){
// 
	// Engine temporary(Operator::kinetic, basis.max_nprim(), basis.max_l());
	// kinetic = temporary;
// 
// }
void initialize_kinetic(){


	Engine temporary(Operator::kinetic, basis.max_nprim(), basis.max_l());

	for (int i = 0; i != omp_get_max_threads(); i++){
		kinetic[i] = temporary; // One engine per thread
	}

}

// void initialize_nuclear(){
// 
	// Engine temporary(Operator::nuclear, basis.max_nprim(), basis.max_l());
	// nuclear = temporary;
// 
// }
void initialize_nuclear(){


	Engine temporary(Operator::nuclear, basis.max_nprim(), basis.max_l());

	for (int i = 0; i != omp_get_max_threads(); i++){
		nuclear[i] = temporary; // One engine per thread
	}

}

void initialize_overlap(){

	Engine temporary(Operator::overlap, basis.max_nprim(), basis.max_l());
	overlap = temporary;

}
