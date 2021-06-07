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