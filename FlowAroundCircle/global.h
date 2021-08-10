#ifndef GLOBAL_H_
#define GLOBAL_H_

//�����C�v�Z�萔
#define GAMMA 1.4
#define R_CIRCLE 0.1
#define RAIR 287.04
#define MU_0 18.24e-6
#define T_0 293.15
#define C 110.4
#define PR 0.72

//���������萔
#define U_INF 20.0
#define P_INF 101.325e+3
#define RHO_INF 1.184

//�����萔��`
#define II_DIR 1//I������r�����C�y��x����
#define JJ_DIR 2//J������theta, y����

//�X�e�b�v��`
/*	II_STEP, JJ_STEP = STEP+5							*/
/*	Grid	0~STEP+4(II_STEP-1) �� STEP+5��(II_STEP)	*/
/*	Cell	0~STEP+3(II_STEP-2) �� STEP+4��(II_STEP-1)	*/
/*	Cell�̂�������2�_�͉��z�i�q							*/
/*	�v�Z�_�͂䂦��STEP��(II_STEP-5)						*/
//#define DELTA_T 1e-8	//Runge-Kutta
#define DELTA_T 2.5e-7	//LU-SGS
#define STEP 100
#define II_STEP (STEP + 5)
#define JJ_STEP (STEP + 5)
#define DELTA_S 2e-5
#define TIME_MAX 1.00

//�X�e�b�v���E��`
/*�������E����*/
/*���z�i�q�ɏ�������*/
/*�������̓Z�����E��27~77*/
/*����āC�Z�����S�ł�27~76*/
#define J_INLET_MIN 27
#define J_INLET_MAX 77
/*�o�����E����*/
/*���z�i�q�ɏ�������*/
/*�������̓Z�����E��2~27, 77~101*/
/*����āC�Z�����S�ł�2~26, 77~101*/
#define J_OUTLET_START 2
#define J_OUTLET_MAX 27
#define J_OUTLET_MIN 77
#define J_OUTLET_END 101

#define T_NUM 4
#define LAM 1.0

#define PARA

#ifdef MAIN_C_
//������
double* rho;
double* ux;
double* vy;
double* Ur;
double* Vtheta;
double* e;
double* p;
double* H;
double* c;
volatile double time;
volatile double time_flag;
volatile int error_flag;
//�ۑ���
double* Q;
double* E;
double* F;
double* Ehalf;
double* Fhalf;
//�S����
double* Ev;
double* Fv;
//���W
double* x;
double* y;
double* x_cen;//�Z�����S
double* y_cen;//�Z�����S
double* rvec;
double* thetaVec;
//��ʉ����W�֘A�C���g���b�N�ƃ��R�r�A��
double* X_xi_half;
double* Y_xi_half;
double* X_eta_half;
double* Y_eta_half;
double* XI_x;
double* XI_y;
double* Eta_x;
double* Eta_y;
double* J_inv; 
#else
//������
extern double* rho;
extern double* ux;
extern double* vy;
extern double* Ur;
extern double* Vtheta;
extern double* e;
extern double* p;
extern double* H;
extern double* c;
extern volatile double time;
extern volatile double time_flag;
extern volatile int error_flag;
//�ۑ���
extern double* Q;
extern double* E;
extern double* F;
extern double* Ehalf;
extern double* Fhalf;
//�S����
extern double* Ev;
extern double* Fv;
//���W
extern double* x;
extern double* y;
extern double* x_cen;//�Z�����S
extern double* y_cen;//�Z�����S
extern double* rvec;
extern double* thetaVec;
//��ʉ����W�֘A�C���g���b�N�ƃ��R�r�A��
extern double* X_xi_half;
extern double* Y_xi_half;
extern double* X_eta_half;
extern double* Y_eta_half;
extern double* XI_x;
extern double* XI_y;
extern double* Eta_x;
extern double* Eta_y;
extern double* J_inv;
#endif

#endif //GLOBAL_H_