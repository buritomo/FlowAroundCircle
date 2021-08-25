#ifndef FDS_KOMEGA_
#define FDS_KOMEGA_

#define EPSILON 1.0
#define KAPPA 1.0 / 3.0

void fds_komega(int dir);
double KWmuscl(double* value, int mg, int kk, int pg, double side);

#endif //FDS_KOMEGA_