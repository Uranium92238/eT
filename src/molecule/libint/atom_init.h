#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_first_ao_in_shells(int *atom, int *faois);
void get_n_shells_on_atoms(int *nsoa);
void get_n_basis_in_shells(int *atom, int *nbis);

#ifdef __cplusplus
}
#endif
