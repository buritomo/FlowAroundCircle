#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES//math.hを使えるようにするコマンド
#include <math.h>
#include "global.h"
#include "calc.h"

void calcInternalValues(void) {
    for (int kx = 0; kx < II_STEP - 1; kx++) {
        for (int ky = 0; ky < JJ_STEP - 1; ky++) {
            int k = kx + ky * II_STEP;
            e[k] = rho[k] * (p[k] / (rho[k] * (GAMMA - 1)) + 0.5 * (ux[k] * ux[k] + vy[k] * vy[k]));
            c[k] = sqrt(GAMMA * p[k] / rho[k]);
            H[k] = GAMMA * p[k] / (rho[k] * (GAMMA - 1)) + 0.5 * (ux[k] * ux[k] + vy[k] * vy[k]);
        }
    }
    return;
}

void makePotential(void) {
    for (int ki = 0; ki < II_STEP - 1; ki++) {
        for (int kj = 0; kj < JJ_STEP - 1; kj++) {
            int k = ki + kj * II_STEP;
            Q[k + 0 * II_STEP * JJ_STEP] = J_inv[k] * rho[k];
            Q[k + 1 * II_STEP * JJ_STEP] = J_inv[k] * rho[k] * ux[k];
            Q[k + 2 * II_STEP * JJ_STEP] = J_inv[k] * rho[k] * vy[k];
            Q[k + 3 * II_STEP * JJ_STEP] = J_inv[k] * e[k];
        }
    }
}

void inversePotentialToParams(void) {
    for (int ki = 0; ki < II_STEP - 1; ki++) {
        for (int kj = 0; kj < JJ_STEP - 1; kj++) {
            int k = ki + kj * II_STEP;
            int dim1 = k + 0 * II_STEP * JJ_STEP;
            int dim2 = k + 1 * II_STEP * JJ_STEP;
            int dim3 = k + 2 * II_STEP * JJ_STEP;
            int dim4 = k + 3 * II_STEP * JJ_STEP;
            /*
            if (k == 10713) {
                printf("???");
            }*/

            rho[k] = Q[dim1] / J_inv[k];
            ux[k] = Q[dim2] / rho[k] / J_inv[k];
            vy[k] = Q[dim3] / rho[k] / J_inv[k];
            e[k] = Q[dim4] / J_inv[k];
            p[k] = (e[k] - 0.5 * rho[k] * (ux[k] * ux[k] + vy[k] * vy[k])) * (GAMMA - 1);
            c[k] = sqrt(GAMMA * p[k] / rho[k]);
            H[k] = GAMMA * p[k] / (rho[k] * (GAMMA - 1)) + 0.5 * (ux[k] * ux[k] + vy[k] * vy[k]);
        }
    }
    return;
}

void printTimer(void) {
    if ((time - time_flag) > (TIME_MAX / 1000000)) {
        double rate = time / TIME_MAX * 100.0;
        printf("%.6f %% is finished....\n", rate);
        time_flag = time;
    }
    return;
}