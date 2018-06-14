#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_ao_g_xyzw(double *g);

#ifdef __cplusplus
}
#endif
