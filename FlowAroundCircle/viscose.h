#ifndef VISCOSE_H_
#define VISCOSE_H_

/*
static double Jaco_inv;
static double xi_x, xi_y, eta_x, eta_y;    //メトリック
static double ux_dx, vy_dy;                //速度偏微分
static double ux_dy, vy_dx;
static double ux_half, vy_half;            //境界での速度
static double tau_xx, tau_xy, tau_yy;      //せん断応力
static double beta_x, beta_y;
static double mu_ave, kappa_ave;
static double T_dx, T_dy;
*/

void viscose(int dir);
void caocJacobian(int k, int dir, double* Jaco_inv);
void calcMetric(int k, int dir, double* xi_x, double* xi_y, double* eta_x, double* eta_y, double Jaco_inv);
void calcVelocityGra(int k, int dir, double* ux_dx, double* ux_dy, double* vy_dx, double* vy_dy, double xi_x, double xi_y, double eta_x, double eta_y);
void calcVeloBoundary(int k, int dir, double* ux_half, double* vy_half);
void calcViscoseFactor(int k, int dir, double xi_x, double xi_y, double eta_x, double eta_y, double tau_xx, double tau_xy, double tau_yy, double beta_x, double beta_y, double Jaco_inv);
double calcTemp(int k);
void calcCoef(int k, int dir, double* mu_ave, double* kappa_ave);
void calcTGra(int k, int dir, double* T_dx, double* T_dy, double xi_x, double xi_y, double eta_x, double eta_y);

#endif//VISCOSE_H_
