#include <stdio.h>
#include "global.h"
#include "boundary.h"
#define _USE_MATH_DEFINES//math.hを使えるようにするコマンド
#include "math.h"

void boundaryValue(void) {
    /*壁面境界条件*/
    /*仮想格子流束を壁面近傍に揃え，摩擦を無視*/
    for (int kj = 0; kj < JJ_STEP - 1; kj++) {
        int ki1 = 0 + kj * II_STEP;//仮想格子
        int ki2 = 1 + kj * II_STEP;//仮想格子
        int ki = 2 + kj * II_STEP;//計算格子

        double r = sqrt(x[ki] * x[ki] + y[ki] * y[ki]);
        double theta;

        if (y[ki] > 0) {
            theta = acos(x[ki] / r);
        }
        else {
            theta = -1 * acos(x[ki] / r);
        }
        double Vtheta = -ux[ki] * sin(theta) + vy[ki] * cos(theta);

        rho[ki1] = rho[ki];
        p[ki1] = p[ki];
        ux[ki1] = -Vtheta * sin(theta);
        vy[ki1] = Vtheta * cos(theta);

        rho[ki2] = rho[ki];
        p[ki2] = p[ki];
        ux[ki2] = -Vtheta * sin(theta);
        vy[ki2] = Vtheta * cos(theta);
    }


    /*入口境界条件*/
    /*仮想格子に初期条件*/
    /*入口側はセル境界で27~77*/
    /*よって，セル中心では27~76*/
    for (int kj = J_INLET_MIN; kj <= J_INLET_MAX; kj++) {
        int ko1 = (II_STEP - 2) + kj * II_STEP;//仮想格子
        int ko2 = (II_STEP - 3) + kj * II_STEP;//仮想格子

        ux[ko1] = U_INF;
        vy[ko1] = 0;
        p[ko1] = P_INF;
        rho[ko1] = RHO_INF;

        ux[ko2] = U_INF;
        vy[ko2] = 0;
        p[ko2] = P_INF;
        rho[ko2] = RHO_INF;
    }

    /*出口境界条件*/
    /*仮想格子に初期条件*/
    /*入口側はセル境界で2~27, 77~101*/
    /*よって，セル中心では2~26, 77~101*/
    for (int kj = J_OUTLET_START; kj <= J_OUTLET_MAX; kj++) {
        int ko1 = (II_STEP - 2) + kj * II_STEP;//仮想格子
        int ko2 = (II_STEP - 3) + kj * II_STEP;//仮想格子
        int ko3 = (II_STEP - 4) + kj * II_STEP;//計算格子
        int ko4 = (II_STEP - 5) + kj * II_STEP;//計算格子

        ux[ko2] = ux[ko3] + (ux[ko3] - ux[ko4]);
        ux[ko1] = ux[ko2] + (ux[ko2] - ux[ko3]);

        vy[ko2] = vy[ko3] + (vy[ko3] - vy[ko4]);
        vy[ko1] = vy[ko2] + (vy[ko2] - vy[ko3]);

        p[ko2] = p[ko3] + (p[ko3] - p[ko4]);
        p[ko1] = p[ko2] + (p[ko2] - p[ko3]);

        rho[ko2] = rho[ko3] + (rho[ko3] - rho[ko4]);
        rho[ko1] = rho[ko2] + (rho[ko2] - rho[ko3]);
    }

    for (int kj = J_OUTLET_MIN; kj <= J_OUTLET_END; kj++) {
        int ko1 = (II_STEP - 2) + kj * II_STEP;//仮想格子
        int ko2 = (II_STEP - 3) + kj * II_STEP;//仮想格子
        int ko3 = (II_STEP - 4) + kj * II_STEP;//計算格子
        int ko4 = (II_STEP - 5) + kj * II_STEP;//計算格子

        ux[ko2] = ux[ko3] + (ux[ko3] - ux[ko4]);
        ux[ko1] = ux[ko2] + (ux[ko2] - ux[ko3]);

        vy[ko2] = vy[ko3] + (vy[ko3] - vy[ko4]);
        vy[ko1] = vy[ko2] + (vy[ko2] - vy[ko3]);

        p[ko2] = p[ko3] + (p[ko3] - p[ko4]);
        p[ko1] = p[ko2] + (p[ko2] - p[ko3]);

        rho[ko2] = rho[ko3] + (rho[ko3] - rho[ko4]);
        rho[ko1] = rho[ko2] + (rho[ko2] - rho[ko3]);
    }

    /*周期境界条件*/
    for (int ki = 0; ki < II_STEP - 1; ki++) {
        int ks1 = ki + 0 * II_STEP;//スタート側仮想格子
        int ks0 = ki + 1 * II_STEP;//スタート側仮想格子
        int ks = ki + 2 * II_STEP;//スタート側計算格子
        int k1 = ki + 3 * II_STEP;//スタート側計算格子

        int kgb = ki + (JJ_STEP - 5) * II_STEP;//ゴール側計算格子
        int kg = ki + (JJ_STEP - 4) * II_STEP;//ゴール側計算格子
        int kg0 = ki + (JJ_STEP - 3) * II_STEP;//ゴール側仮想格子
        int kg1 = ki + (JJ_STEP - 2) * II_STEP;//ゴール側仮想格子

        ux[ks1] = ux[kgb];
        vy[ks1] = vy[kgb];
        p[ks1] = p[kgb];
        rho[ks1] = rho[kgb];

        ux[ks0] = ux[kg];
        vy[ks0] = vy[kg];
        p[ks0] = p[kg];
        rho[ks0] = rho[kg];

        ux[kg0] = ux[ks];
        vy[kg0] = vy[ks];
        p[kg0] = p[ks];
        rho[kg0] = rho[ks];

        ux[kg1] = ux[k1];
        vy[kg1] = vy[k1];
        p[kg1] = p[k1];
        rho[kg1] = rho[k1];
    }

    return;
}