#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void initialize_coulomb();
void initialize_basis();
void get_ao_L_xyzw(double *L, int *s1, int *s3);

#ifdef __cplusplus
}
#endif
