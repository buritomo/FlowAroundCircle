#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_

static double* kari;
static double error_max;
static int kkk, lxx, lyy, mmm;

void setAssumedPotential(void);
void releaseAssumedPotential(void);
void rungekutta(void);

void resetError();
void errorCheck(int k, int lx, int ly, int m);
void printErrorMax(void);

#endif //RUNGEKUTTA_H_