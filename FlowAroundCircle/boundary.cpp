#include <stdio.h>
#include "global.h"
#include "boundary.h"

void boundaryValue(void) {
    /*壁面境界条件*/
    /*仮想格子流束を壁面近傍に揃え，摩擦を無視*/
    for (int kj = 0; kj < JJ_STEP - 1; kj++) {
        int ki1 = 0 + kj + JJ_STEP;//仮想格子
        int ki2 = 1 + kj * JJ_STEP;//仮想格子
        int ki = 2 + kj * JJ_STEP;//計算格子

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

    }

    /*出口境界条件*/
    /*仮想格子に初期条件*/
    /*入口側はセル境界で2~27, 77~101*/
    /*よって，セル中心では2~26, 77~101*/


    /*周期境界条件*/
    return;
}