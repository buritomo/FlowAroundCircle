#include <stdio.h>
#include "global.h"
#include "boundary.h"

void boundaryValue(void) {
    /*�ǖʋ��E����*/
    /*���z�i�q������ǖʋߖT�ɑ����C���C�𖳎�*/
    for (int kj = 0; kj < JJ_STEP - 1; kj++) {
        int ki1 = 0 + kj + JJ_STEP;//���z�i�q
        int ki2 = 1 + kj * JJ_STEP;//���z�i�q
        int ki = 2 + kj * JJ_STEP;//�v�Z�i�q

        rho[ki1] = rho[ki2] = rho[ki];
        p[ki1] = p[ki2] = p[ki];
        ux[ki1] = ux[ki2] = ux[ki];
        vy[ki1] = vy[ki2] = vy[ki];
    }


    /*�������E����*/
    /*���z�i�q�ɏ�������*/
    /*�������̓Z�����E��27~77*/
    /*����āC�Z�����S�ł�27~76*/
    for (int kj = J_INLET_MIN; kj <= J_INLET_MAX; kj++) {

    }

    /*�o�����E����*/
    /*���z�i�q�ɏ�������*/
    /*�������̓Z�����E��2~27, 77~101*/
    /*����āC�Z�����S�ł�2~26, 77~101*/


    /*�������E����*/
    return;
}