#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#define _USE_MATH_DEFINES//math.h���g����悤�ɂ���R�}���h
#include "math.h"
#include "initial.h"

/*���������͓��G���g���s��������肷��*/
/*�v�Z�̈�̂ݒ�`���Ă��邱�Ƃɒ���*/
void initialValue(void) {
	for (int ki = 2; ki < II_STEP - 3; ki++) {
		for (int kj = 2; kj < JJ_STEP - 3; kj++) {
			int k = ki + kj * II_STEP;
			double r = sqrt(x[k] * x[k] + y[k] * y[k]);
			double theta;

			if (y[k] > 0) {
				theta = acos(x[k] / r);
			}
			else {
				theta = -1 * acos(x[k] / r);
			}

			double Vr = U_INF * (1 - (R_CIRCLE * R_CIRCLE) / (r * r)) * cos(theta);
			double Vtheta = -U_INF * (1 + (R_CIRCLE * R_CIRCLE) / (r * r)) * sin(theta);

			ux[k] = Vr * cos(theta) - Vtheta * sin(theta);
			vy[k] = Vr * sin(theta) + Vtheta * cos(theta);

			rho[k] = RHO_INF;
			p[k] = P_INF + 0.5 * rho[k] * (U_INF * U_INF - (ux[k] * ux[k] + vy[k] * vy[k]));
		}
	}
	return;
}

