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

vector<Engine> electronic_repulsion_1der(omp_get_max_threads());
extern vector<Engine> electronic_repulsion_1der;

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

vector<Engine> potential(omp_get_max_threads());
extern vector<Engine> potential;

vector<int> shell2atom;
extern vector<int> shell2atom;

void initialize_atoms(char *name){

    string xyzfilename(strcat(name,".xyz"));

    ifstream input_file(xyzfilename);
    vector<Atom> temporary_atoms = read_dotxyz(input_file);
    atoms = temporary_atoms;

}

void reset_basis(){

    eTBasis temp_basis;
    basis = temp_basis;

}

void initialize_basis(char *basisset, char *filename, bool *cartesian_gaussians){
    
    string xyzfilename(strcat(filename,".xyz"));

    ifstream input_file(xyzfilename);

    vector<Atom> temporary_atoms = read_dotxyz(input_file);

    cout.setstate(ios_base::failbit);
    BasisSet temporary(basisset, temporary_atoms);
    cout.clear();

    temporary.set_pure(!(*cartesian_gaussians));
    basis.add(temporary);

}

void initialize_shell2atom(){

    shell2atom = basis.shell2atom(atoms);

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

    Engine temporary_1der(Operator::coulomb, basis.max_nprim(), basis.max_l(), 1);
    temporary_1der.set_precision(1.0e-25);

    for (int i = 0; i != omp_get_max_threads(); i++){
        electronic_repulsion_1der[i] = temporary_1der;
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

    vector<pair<double, array<double, 3>>> q;
    for (const auto& atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }

    Engine temporary(Operator::nuclear, basis.max_nprim(), basis.max_l());
    temporary.set_params(q);

    for (int i = 0; i != omp_get_max_threads(); i++){
        nuclear[i] = temporary; 
    }

    Engine temporary_1der(Operator::nuclear, basis.max_nprim(), basis.max_l(),1);
    temporary_1der.set_params(q);
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

void initialize_potential(double *charges, double *coordinates, int *n_points){

    Engine temporary(Operator::nuclear, basis.max_nprim(), basis.max_l());
      
    vector<pair<double, array<double, 3>>> points;
    for (int i = 0; i <= *n_points - 1; i++) {
      int offset = (i-1)*3 + 3;
      points.push_back({static_cast<double>(charges[i]),
                   {{*(coordinates+offset),*(coordinates+offset+1),*(coordinates+offset+2)}}});
    }

    temporary.set_params(points); // Tell the engine where the atomic charges are

    for (int i = 0; i != omp_get_max_threads(); i++){
        potential[i] = temporary; 
    }

}
