#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void initialize_coulomb();
void initialize_basis();
void initialize_kinetic();
void initialize_nuclear();
void initialize_overlap();

#ifdef __cplusplus
}
#endif
