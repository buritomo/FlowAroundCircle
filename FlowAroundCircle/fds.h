#ifndef FDS_H_
#define FDS_H_

#define EPSILON 1.0
#define KAPPA 1.0 / 3.0

/*
static double b1, b2;
static double Q_L[4];
static double Q_R[4];
static double E_L[4];
static double E_R[4];
static double F_L[4];
static double F_R[4];
static double R[4 * 4];
static double L[4 * 4];
static double lambda[4];
static double w[4 * 4];
static double J_bd;
static double k_x, k_y;
static double theta;
static double Z;
static double rho_L, rho_R;
static double u_L, u_R;
static double v_L, v_R;
static double e_L, e_R;
static double p_L, p_R;
static double H_L, H_R;
static double c_L, c_R;
static double rho_ave;
static double u_ave;
static double v_ave;
static double e_ave;
static double p_ave;
static double H_ave;
static double c_ave;
*/

void fds(int dir);
/*
void RightArray(double* R);
void LeftArray(double* L);
void Lambda(double* lambda);
void musclArray(int k, int dir);
void RoeAverage(void);
void muscl(double* left, double* right, double* value, int k, int dir);
*/
void RightArray(double* R, double k_x, double k_y, double Z, double u_ave, double v_ave, double c_ave, double H_ave);
void LeftArray(double* L, double k_x, double k_y, double Z, double u_ave, double v_ave, double c_ave, double H_ave, double b1, double b2);
void Lambda(double* lambda, double k_x, double k_y, double c_ave, double Z);
void musclArray(int k, int dir, double* Q_L, double* Q_R, double* E_L, double* E_R, double* F_L, double* F_R, double rho_L, double rho_R, double u_L, double u_R, double v_L, double v_R, double p_L, double p_R, double e_L, double e_R, double k_x, double k_y, double J_bd);
void RoeAverage(double* rho_ave, double rho_L, double rho_R, double* u_ave, double u_L, double u_R, double* v_ave, double v_L, double v_R, double* H_ave, double H_L, double H_R, double* c_ave);
void muscl(double* left, double* right, double* value, int k, int dir);
double limiter(double x, double y);
double max(double x, double y);
double min(double x, double y);
double sgn(double x);


#endif //FDS_H_