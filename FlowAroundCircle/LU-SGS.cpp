#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "LU-SGS.h"
#include "fds.h"
#include "boundary.h"
#include "viscose.h"

void GaussSeidel(void) {
	int itr;
	double error_max = 0.0;

	dQ_initial();

	for (itr = 0; itr < T_NUM; itr++) {
		fds(II_DIR);
		fds(JJ_DIR);
		viscose(II_DIR);
		viscose(JJ_DIR);

		for (int ki = 2; ki < II_STEP - 3; ki++) {
			for (int kj = 2; kj < JJ_STEP - 3; kj++) {
				double RHS[4];
				double rhs[4];
				double deltaFlux[2][4];
				double spe_r0[2];
				double LDinv;
				int k = ki + kj * II_STEP;

				for (int m = 0; m < 4; m++) {
					RHS[m] = -DELTA_T * (Ehalf[k + II_STEP * JJ_STEP * m] - Ehalf[k + II_STEP * JJ_STEP * m - 1]);
					RHS[m] += -DELTA_T * (Ev[k + II_STEP * JJ_STEP * m] - Ev[k + II_STEP * JJ_STEP * m - 1]);
					RHS[m] += -DELTA_T * (Fhalf[k + II_STEP * JJ_STEP * m] - Fhalf[k + II_STEP * JJ_STEP * m - II_STEP]);
					RHS[m] += -DELTA_T * (Fv[k + II_STEP * JJ_STEP * m] - Fv[k + II_STEP * JJ_STEP * m - II_STEP]);
				}

				if (itr == 0) {
					for (int m = 0; m < 4; m++) {
						tmpRHS[k + II_STEP * JJ_STEP * m] = RHS[m];
					}
				}

				for (int m = 0; m < 4; m++) {
					rhs[m] = (tmpQ[k + II_STEP * JJ_STEP * m] - Q[k + II_STEP * JJ_STEP * m]) + 0.5 * (RHS[m] + tmpRHS[k + II_STEP * JJ_STEP * m]);
				}

				calcPM(k, -1, dQ, deltaFlux);

				for (int m = 0; m < 4; m++) {
					rhs[m] += LAM * DELTA_T * (deltaFlux[II_DIR - 1][m] + deltaFlux[JJ_DIR - 1][m]);
				}

				calcLX(k, spe_r0);
				LDinv = 1.0 / (1.0 + LAM * DELTA_T * (spe_r0[0] + spe_r0[1]));

				for (int m = 0; m < 4; m++) {
					dQ[k + II_STEP * JJ_STEP * m] = rhs[m] * LDinv;
				}
			}
		}

		for (int ki = 2; ki < II_STEP - 3; ki++) {
			for (int kj = 2; kj < JJ_STEP - 3; kj++) {
				int k = ki + kj * II_STEP;
				double deltaFlux[2][4];
				double rhs[4];
				double spe_r0[2];
				double LDinv;

				for (int m = 0; m < 4; m++) {
					rhsF[k + II_STEP * JJ_STEP * m] = dQ[k + II_STEP * JJ_STEP * m];
				}

				calcPM(k, 1, dQ, deltaFlux);

				for (int m = 0; m < 4; m++) {
					rhs[m] += LAM * DELTA_T * (deltaFlux[II_DIR - 1][m] + deltaFlux[JJ_DIR - 1][m]);
				}
				
				calcLX(k, spe_r0);
				LDinv = 1.0 / (1.0 + LAM * DELTA_T * (spe_r0[0] + spe_r0[1]));

				for (int m = 0; m < 4; m++) {
					rhsF[k + II_STEP * JJ_STEP * m] = dQ[k + II_STEP * JJ_STEP * m] - rhs[m] * LDinv;
				}
			}
		}

		for (int ki = 2; ki < II_STEP - 3; ki++) {
			for (int kj = 2; kj < JJ_STEP - 3; kj++) {
				int k = ki + kj * II_STEP;
				double invQ, Jaco;
				double error4;

				Jaco = 1.0 / J_inv[k];

				for (int m = 0; m < 4; m++) {
					int it = k + II_STEP * JJ_STEP * m;
					Q[it] = Q[it] + rhsF[it];
					dQ[it] = rhsF[it];
				}
				
				invQ = 1.0 / Q[k + II_STEP * JJ_STEP * 0];

				rho[k] = Q[k + II_STEP * JJ_STEP * 0] * Jaco;
				ux[k] = Q[k + II_STEP * JJ_STEP * 1] * invQ;
				vy[k] = Q[k + II_STEP * JJ_STEP * 2] * invQ;
				e[k] = Q[k + II_STEP * JJ_STEP * 3] * invQ;
				p[k] = (e[k] - 0.5 * rho[k] * (ux[k] * ux[k] + vy[k] * vy[k])) * (GAMMA - 1);

				boundaryValue();

				error4 = dQ[k + II_STEP * JJ_STEP * 3] / Q[k + II_STEP * JJ_STEP * 3];
				if (error_max < fabs(error4)) {
					error_max = fabs(error4);
				}

				if (isnan(rho[k])) {
					printf("rho is NaN! at (%d, %d)\n", ki, kj);
					exit(1);
				}
				if (isnan(ux[k])) {
					printf("ux is NaN! at (%d, %d)\n", ki, kj);
					exit(1);
				}
				if (isnan(vy[k])) {
					printf("vy is NaN! at (%d, %d)\n", ki, kj);
					exit(1);
				}
				if (isnan(e[k])) {
					printf("e is NaN! at (%d, %d)\n", ki, kj);
					exit(1);
				}
			}
		}
	}

	printf("Time: %.8f\, error_max %lf", time, error_max);


	return;
}

void dQ_initial(void) {
	for (int ki = 0; ki < II_STEP; ki++) {
		for (int kj = 0; kj < JJ_STEP; kj++) {
			int k = ki + kj * II_STEP;

			dQ[k + II_STEP * JJ_STEP * 0] = 0.0;
			dQ[k + II_STEP * JJ_STEP * 1] = 0.0;
			dQ[k + II_STEP * JJ_STEP * 2] = 0.0;
			dQ[k + II_STEP * JJ_STEP * 3] = 0.0;

			rhsF[k + II_STEP * JJ_STEP * 0] = 0.0;
			rhsF[k + II_STEP * JJ_STEP * 1] = 0.0;
			rhsF[k + II_STEP * JJ_STEP * 2] = 0.0;
			rhsF[k + II_STEP * JJ_STEP * 3] = 0.0;
		}
	}
	return;
}

void calcPM(int k, int pmflag, double *dQ, double(*flux)[4]) {
	double SS[2][2];
	int ddd[2];
	int xi = 0, eta = 1;
	int x = 0, y = 1;
	double alpha = 1.01;
	double Jaco;
	double kx, ky;
	double Sk;
	double invQ;
	double tmprho, tmpux, tmpvy, tmpe, tmpp;
	double ZZ, ZZ0;
	double c, tmp, Mu;
	double P;
	double spe;

	ddd[xi] = 1;
	ddd[eta] = II_STEP;

	for (int mflag = 0; mflag < 2; mflag++) {
		int it = k + pmflag * ddd[mflag];

		SS[xi][x] = 0.5 * (X_eta_half[it] + X_eta_half[it + 1]);
		SS[xi][y] = 0.5 * (X_xi_half[it] + X_xi_half[it + 1]);
		SS[eta][x] = 0.5 * (Y_eta_half[it] + Y_eta_half[it + II_STEP]);
		SS[eta][y] = 0.5 * (Y_xi_half[it] + Y_xi_half[it + II_STEP]);


		if (dQ[it + II_STEP * JJ_STEP * 0] != 0.0) {
			Jaco = 1 / J_inv[k];

			kx = SS[mflag][x];
			ky = SS[mflag][y];
			Sk = sqrt(kx * kx + ky * ky);
			
			invQ = 1.0 / (dQ[it + II_STEP * JJ_STEP * 0] + rho[it] / Jaco);
			tmprho = dQ[it + II_STEP * JJ_STEP * 0] * Jaco + rho[it];
			tmpux = (dQ[it + II_STEP * JJ_STEP * 2] + ux[it] / Jaco * rho[it]) * invQ;
			tmpvy = (dQ[it + II_STEP * JJ_STEP * 3] + vy[it] / Jaco * rho[it]) * invQ;
			tmpe = (dQ[it + II_STEP * JJ_STEP * 4] * Jaco + e[it]);
			tmpp = (GAMMA - 1.0) * (tmpe - 0.5 * tmprho * (tmpux * tmpux + tmpvy * tmpvy));
		
			ZZ = kx * tmpux + ky * tmpvy;
			ZZ0 = kx * ux[it] + ky * vy[it];
			P = (e[it] - 0.5 * rho[it] * (ux[it] * ux[it] + vy[it] * vy[it])) * (GAMMA - 1);
			c = sqrt(GAMMA * P / rho[it]);
			tmp = P / rho[it] * RAIR;
			Mu = MU_0 * pow((tmp / (273.15 + 20.0)), 1.5) * (273.15 + 20.0 + C) / (tmp + C);

			spe = alpha * (fabs(ZZ0) + c * Sk) + 2.0 * (Mu) * Sk * Sk / rho[it];
		
			flux[mflag][0] = 0.5 * ((tmprho * ZZ - rho[it] * ZZ0) / Jaco - pmflag * spe * dQ[it]);
			flux[mflag][1] = 0.5 * (((tmprho * tmpux * ZZ + kx * tmpp) - (rho[it] * ux[it] * ZZ0 + kx * p[it])) / Jaco - pmflag * spe * dQ[it]);
			flux[mflag][2] = 0.5 * (((tmprho * tmpvy * ZZ + ky * tmpp) - (rho[it] * vy[it] * ZZ0 + ky * p[it])) / Jaco - pmflag * spe * dQ[it]);
			flux[mflag][3] = 0.5 * (((tmpe + tmpp) * ZZ - (e[it] + p[it]) * ZZ0) / Jaco - pmflag * spe * dQ[it]);
		
		}
		else {
			flux[mflag][0] = 0.0;
			flux[mflag][1] = 0.0;
			flux[mflag][2] = 0.0;
			flux[mflag][3] = 0.0;
		}
	}

	return;
}

void calcLX(int k, double(*spe)) {
	double Jaco;
	double SS[2][2];
	double kx, ky;
	int xi = 0, eta = 1;
	int x = 0, y = 1;
	double ZZ;
	double Sk;
	double c;
	double tmp, pre;
	double Mu;
	double alpha = 1.01;

	Jaco = 1 / J_inv[k];

	SS[xi][x] = 0.5 * (X_eta_half[k] + X_eta_half[k + 1]);
	SS[xi][y] = 0.5 * (X_xi_half[k] + X_xi_half[k + 1]);
	SS[eta][x] = 0.5 * (Y_eta_half[k] + Y_eta_half[k + II_STEP]);
	SS[eta][y] = 0.5 * (Y_xi_half[k] + Y_xi_half[k + II_STEP]);

	for (int mflag = 0; mflag < 2; mflag++) {
		kx = SS[mflag][x];
		ky = SS[mflag][y];
		pre = (e[k] - 0.5 * rho[k] * (ux[k] * ux[k] + vy[k] * vy[k])) * (GAMMA - 1);
		c = sqrt(GAMMA + pre / rho[k]);
		tmp = pre / rho[k] * RAIR;
		Mu = MU_0 * (pow((tmp) / (273.15 + 20.0), 1.5) * (273.15 + 20.0 + C) / (tmp + C));
		spe[mflag] = alpha * (fabs(ZZ) + c * Sk) + 2.0 * (Mu) * Sk * Sk / rho[k];
	}

	return;
}

void setMatrixs(void) {
	double* tmpQ;
	double* tmpRHS;
	double* dQ;
	double* rhsF;

	tmpQ = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	tmpRHS = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	dQ = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	rhsF = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);

	return;
}