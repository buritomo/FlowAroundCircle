#ifndef LUSGS_H_
#define LUSGS_H_

void GaussSeidel(void);
void dQ_initial(void);
void calcPM(int k, int pmflag, double *dQ, double(*flux)[4]);
void calcLX(int k, double(*spe));
void setMatrixs(void);


#ifdef MAIN_C_
double* tmpQ;
double* tmpRHS;
double* dQ;
double* rhsF;
#else
extern double* tmpQ;
extern double* tmpRHS;
extern double* dQ;
extern double* rhsF;
#endif//MAIN_C_

#endif //LUSGS_H_