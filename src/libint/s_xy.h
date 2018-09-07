#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_ao_s_xy(double *s);
void get_n_shells(int *ns);
void construct_ao_s_wx(double *s, long *s1, long *s2);

#ifdef __cplusplus
}
#endif
