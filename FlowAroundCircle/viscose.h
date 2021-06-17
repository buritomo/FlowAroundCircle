#ifndef VISCOSE_H_
#define VISCOSE_H_

static double Jaco_inv;
static double xi_x, xi_y, eta_x, eta_y;    //���g���b�N
static double ux_dx, vy_dy;                //���x�Δ���
static double ux_half, vy_half;            //���E�ł̑��x
static double tau_xx, tau_xy, tau_yy;      //����f����
static double beta_x, beta_y;
//static double T, T_a, T_b;
static double mu_ave, kappa_ave;
static double T_dx, T_dy;


void viscose(int dir);
void caocJacobian(int k, int dir);
void calcMetric(int k, int dir);
void calcVelocityGra(int k, int dir);
void calcVeloBoundary(int k, int dir);
void calcViscoseFactor(int k, int dir);
double calcTemp(int k);
void calcCoef(int k, int dir);

#endif//VISCOSE_H_
