/*
/
/ 	Globals
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
/ 	Contains the only global variables allowed in the integral calculations,
/ 	namely the basis set, atomic information, and integral engines.
/
*/
#include "eT_basis.h"

extern eTBasis basis; 										 // The basis set used throughout
//extern libint2::Engine electronic_repulsion; 					 // The electronic repulsion engine, deprecated
extern vector<libint2::Engine> electronic_repulsion_engines; // The electronic repulsion engines vector for parallellization
extern libint2::Engine kinetic; 										 // The kinetic energy engine
extern libint2::Engine nuclear; 										 // The nuclear attraction engine
extern libint2::Engine overlap; 										 // The AO overlap engine
extern vector<libint2::Atom> atoms; 								 // Atoms vector
extern vector<long> shell2bf_g;
