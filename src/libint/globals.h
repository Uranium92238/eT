/*
/
/   Globals
/   Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
/   Contains the only global variables allowed in the integral calculations:
/   the basis set, atomic information, and integral engines.
/
*/
#include "eT_basis.h"

extern eTBasis basis;                                        // The basis set used throughout
extern vector<libint2::Engine> electronic_repulsion_engines; // The electronic repulsion engines vector
extern vector<libint2::Engine> kinetic;                      // The kinetic energy engines vector
extern vector<libint2::Engine> nuclear;                      // The nuclear attraction engine vector
extern libint2::Engine overlap;                              // The overlap engine
extern libint2::Engine dipole;                               // The dipole engine
extern vector<libint2::Atom> atoms;                          // Atoms vector
