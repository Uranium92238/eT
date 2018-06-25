#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void initialize_coulomb();
void initialize_basis();

#ifdef __cplusplus
}
#endif
