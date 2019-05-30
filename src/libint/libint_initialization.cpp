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
 
    Libint initialization of engines and basis set
    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
 
*/
#include "libint_initialization.h"
#include "eT_basis.h"
#include "omp_control.h"

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

vector<Engine> electronic_repulsion_engines(omp_get_max_threads());
extern vector<Engine> electronic_repulsion_engines;

vector<Engine> kinetic(omp_get_max_threads());
extern vector<Engine> kinetic;

Engine kinetic_1der;
extern Engine kinetic_1der;

vector<Engine> nuclear(omp_get_max_threads());
extern vector<Engine> nuclear;

Engine nuclear_1der;
extern Engine nuclear_1der;

Engine overlap_1der;
extern Engine overlap_1der;

vector<Engine> overlap(omp_get_max_threads());
extern vector<Engine> overlap;

vector<Engine> dipole(omp_get_max_threads());
extern vector<Engine> dipole;

vector<Engine> quadrupole(omp_get_max_threads());
extern vector<Engine> quadrupole;

vector<Atom> atoms;
extern vector<Atom> atoms;

void initialize_atoms(char *name){

    string xyzfilename(strcat(name,".xyz"));

    ifstream input_file(xyzfilename);
    vector<Atom> temporary_atoms = read_dotxyz(input_file);
    atoms = temporary_atoms;

}

void initialize_basis(char *basisset, char *filename){
    
    string xyzfilename(strcat(filename,".xyz"));

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
        electronic_repulsion_engines[i] = temporary;
    }

}

void set_coulomb_precision(double *prec){

    for (int i = 0; i != omp_get_max_threads(); i++){
        electronic_repulsion_engines[i].set_precision(*prec);
    }

}

void initialize_kinetic(){


    Engine temporary(Operator::kinetic, basis.max_nprim(), basis.max_l());

    for (int i = 0; i != omp_get_max_threads(); i++){
        kinetic[i] = temporary; 
    }

    Engine temporary_1der(Operator::kinetic, basis.max_nprim(), basis.max_l(),1);
    kinetic_1der = temporary_1der; 

}

void initialize_nuclear(){


    Engine temporary(Operator::nuclear, basis.max_nprim(), basis.max_l());

    temporary.set_params(make_point_charges(atoms));                           // Tell the engine where the atomic charges are

    for (int i = 0; i != omp_get_max_threads(); i++){
        nuclear[i] = temporary; 
    }

    Engine temporary_1der(Operator::nuclear, basis.max_nprim(), basis.max_l(),1);
    temporary_1der.set_params(make_point_charges(atoms));
    nuclear_1der = temporary_1der; 

}

void initialize_overlap(){

    Engine temporary(Operator::overlap, basis.max_nprim(), basis.max_l());
   
    for (int i = 0; i != omp_get_max_threads(); i++){
        overlap[i] = temporary; 
    }

    Engine temporary_1der(Operator::overlap, basis.max_nprim(), basis.max_l(),1);
    overlap_1der = temporary_1der;

}

void initialize_dipole(){

    Engine temporary(Operator::emultipole1, basis.max_nprim(), basis.max_l()); // Note that expansion point = (0,0,0) by default

    for (int i = 0; i != omp_get_max_threads(); i++){
        dipole[i] = temporary; 
    }

}

void initialize_quadrupole(){

    Engine temporary(Operator::emultipole2, basis.max_nprim(), basis.max_l()); // Note that expansion point = (0,0,0) by default

    for (int i = 0; i != omp_get_max_threads(); i++){
        quadrupole[i] = temporary; 
    }

}
