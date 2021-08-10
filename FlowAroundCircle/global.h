#ifndef GLOBAL_H_
#define GLOBAL_H_

//物性，計算定数
#define GAMMA 1.4
#define R_CIRCLE 0.1
#define RAIR 287.04
#define MU_0 18.24e-6
#define T_0 293.15
#define C 110.4
#define PR 0.72

//無限遠方定数
#define U_INF 20.0
#define P_INF 101.325e+3
#define RHO_INF 1.184

//方向定数定義
#define II_DIR 1//I方向はr方向，及びx方向
#define JJ_DIR 2//J方向はtheta, y方向

//ステップ定義
/*	II_STEP, JJ_STEP = STEP+5							*/
/*	Grid	0~STEP+4(II_STEP-1) の STEP+5個(II_STEP)	*/
/*	Cell	0~STEP+3(II_STEP-2) の STEP+4個(II_STEP-1)	*/
/*	Cellのうち両側2点は仮想格子							*/
/*	計算点はゆえにSTEP個(II_STEP-5)						*/
//#define DELTA_T 1e-8	//Runge-Kutta
#define DELTA_T 2.5e-7	//LU-SGS
#define STEP 100
#define II_STEP (STEP + 5)
#define JJ_STEP (STEP + 5)
#define DELTA_S 2e-5
#define TIME_MAX 1.00

//ステップ境界定義
/*入口境界条件*/
/*仮想格子に初期条件*/
/*入口側はセル境界で27~77*/
/*よって，セル中心では27~76*/
#define J_INLET_MIN 27
#define J_INLET_MAX 77
/*出口境界条件*/
/*仮想格子に初期条件*/
/*入口側はセル境界で2~27, 77~101*/
/*よって，セル中心では2~26, 77~101*/
#define J_OUTLET_START 2
#define J_OUTLET_MAX 27
#define J_OUTLET_MIN 77
#define J_OUTLET_END 101

#define T_NUM 4
#define LAM 1.0

#define PARA

#ifdef MAIN_C_
//物理量
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
//保存量
double* Q;
double* E;
double* F;
double* Ehalf;
double* Fhalf;
//粘性項
double* Ev;
double* Fv;
//座標
double* x;
double* y;
double* x_cen;//セル中心
double* y_cen;//セル中心
double* rvec;
double* thetaVec;
//一般化座標関連，メトリックとヤコビアン
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
//物理量
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
//保存量
extern double* Q;
extern double* E;
extern double* F;
extern double* Ehalf;
extern double* Fhalf;
//粘性項
extern double* Ev;
extern double* Fv;
//座標
extern double* x;
extern double* y;
extern double* x_cen;//セル中心
extern double* y_cen;//セル中心
extern double* rvec;
extern double* thetaVec;
//一般化座標関連，メトリックとヤコビアン
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