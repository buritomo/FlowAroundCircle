#include <stdio.h>
#include "global.h"
#include "boundary.h"

void boundaryValue(void) {
    /*壁面境界条件*/
    /*仮想格子流束を壁面近傍に揃え，摩擦を無視*/
    for (int kj = 0; kj < JJ_STEP - 1; kj++) {
        int ki1 = 0 + kj + II_STEP;//仮想格子
        int ki2 = 1 + kj * II_STEP;//仮想格子
        int ki = 2 + kj * II_STEP;//計算格子

        rho[ki1] = rho[ki2] = rho[ki];
        p[ki1] = p[ki2] = p[ki];
        ux[ki1] = ux[ki2] = ux[ki];
        vy[ki1] = vy[ki2] = vy[ki];
    }


    /*入口境界条件*/
    /*仮想格子に初期条件*/
    /*入口側はセル境界で27~77*/
    /*よって，セル中心では27~76*/
    for (int kj = J_INLET_MIN; kj <= J_INLET_MAX; kj++) {
        int ko1 = (II_STEP - 1) + kj * II_STEP;//仮想格子
        int ko2 = (II_STEP - 2) + kj * II_STEP;//仮想格子

        ux[ko1] = ux[ko2] = U_INF;
        vy[ko1] = vy[ko2] = 0;
        p[ko1] = p[ko2] = P_INF;
        rho[ko1] = rho[ko2] = RHO_INF;
    }

    /*出口境界条件*/
    /*仮想格子に初期条件*/
    /*入口側はセル境界で2~27, 77~101*/
    /*よって，セル中心では2~26, 77~101*/
    for (int kj = J_OUTLET_START; kj <= J_OUTLET_MAX; kj++) {
        int ko1 = (II_STEP - 1) + kj * II_STEP;//仮想格子
        int ko2 = (II_STEP - 2) + kj * II_STEP;//仮想格子
        int ko3 = (II_STEP - 3) + kj + II_STEP;//計算格子
        int ko4 = (II_STEP - 4) + kj + II_STEP;//計算格子

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
        rho[ks1] = p[kgb];

        ux[ks0] = ux[kg];
        vy[ks0] = vy[kg];
        p[ks0] = p[kg];
        rho[ks0] = p[kg];

        ux[kg0] = ux[ks];
        vy[kg0] = vy[ks];
        p[kg0] = p[ks];
        rho[kg0] = p[ks];

        ux[kg1] = ux[k1];
        vy[kg1] = vy[k1];
        p[kg1] = p[k1];
        rho[kg1] = p[k1];
    }

    return;
}