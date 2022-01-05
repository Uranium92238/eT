//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2022 the authors of eT
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
#include <vector>
#include "globals.h"

#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_shell_numbers(int *atom, int *shells);
void get_first_ao_in_shells(int *atom, int *first_ao_in_shell);
void get_n_shells_on_atom(int *atom, int *n_shells);
void get_n_aos_in_shell(int *atom, int *n_aos_in_shell);
void initialize_shell_to_first_ao();
void initialize_atom_to_shell_list();

#ifdef __cplusplus
}
#endif
