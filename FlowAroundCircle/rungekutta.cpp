#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "rungekutta.h"

void setAssumedPotential(void) {

    kari = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);

    return;
}

void releaseAssumedPotential(void) {

    free(kari);

    return;
}

void rungekutta(void) {
    double lam[4];
    lam[0] = 0.25;
    lam[1] = 1.0 / 3.0;
    lam[2] = 0.5;
    lam[3] = 1.0;

    double step;

    resetError();

    for (int kx = 0; kx < II_STEP - 1; kx++) {
        for (int ky = 0; ky < JJ_STEP - 1; ky++) {
            int k = kx + ky * II_STEP;
            kari[k + 0 * II_STEP * JJ_STEP] = Q[k + 0 * II_STEP * JJ_STEP];
            kari[k + 1 * II_STEP * JJ_STEP] = Q[k + 1 * II_STEP * JJ_STEP];
            kari[k + 2 * II_STEP * JJ_STEP] = Q[k + 2 * II_STEP * JJ_STEP];
            kari[k + 3 * II_STEP * JJ_STEP] = Q[k + 3 * II_STEP * JJ_STEP];
        }
    }

    for (int k = 0; k < 4; k++) {
        for (int lx = 2; lx < II_STEP - 2; lx++) {
            for (int ly = 0; ly < JJ_STEP - 1; ly++) {
                int l = lx + ly * II_STEP;
                for (int m = 0; m < 4; m++) {
                    kari[l + II_STEP * JJ_STEP * m] = kari[l + II_STEP * JJ_STEP * m] - lam[k] * DELTA_T * (Ehalf[l + II_STEP * JJ_STEP * m] - Ehalf[l + II_STEP * JJ_STEP * m - 1]);
                    //errorCheck(k, lx, ly, m);
                }
            }
        }
    }

    for (int k = 0; k < 4; k++) {
        for (int lx = 0; lx < II_STEP - 1; lx++) {
            for (int ly = 2; ly < JJ_STEP - 2; ly++) {
                int l = lx + ly * II_STEP;
                for (int m = 0; m < 4; m++) {
                    kari[l + II_STEP * JJ_STEP * m] = kari[l + II_STEP * JJ_STEP * m] - lam[k] * DELTA_T * (Fhalf[l + II_STEP * JJ_STEP * m] - Fhalf[l + II_STEP * JJ_STEP * m - II_STEP]);
                    //errorCheck(k, lx, ly, m);
                }

            }
        }
    }

    for (int k = 0; k < 4; k++) {
        for (int lx = 2; lx < II_STEP - 2; lx++) {
            for (int ly = 0; ly < JJ_STEP - 1; ly++) {
                int l = lx + ly * II_STEP;
                for (int m = 0; m < 4; m++) {
                    kari[l + II_STEP * JJ_STEP * m] = kari[l + II_STEP * JJ_STEP * m] + lam[k] * DELTA_T * (Ev[l + II_STEP * JJ_STEP * m] - Ev[l + II_STEP * JJ_STEP * m - 1]);
                    //errorCheck(k, lx, ly, m);
                }
            }
        }
    }

    for (int k = 0; k < 4; k++) {
        for (int lx = 0; lx < II_STEP - 1; lx++) {
            for (int ly = 2; ly < JJ_STEP - 2; ly++) {
                int l = lx + ly * II_STEP;
                for (int m = 0; m < 4; m++) {
                    kari[l + II_STEP * JJ_STEP * m] = kari[l + II_STEP * JJ_STEP * m] + lam[k] * DELTA_T * (Fv[l + II_STEP * JJ_STEP * m] - Fv[l + II_STEP * JJ_STEP * m - II_STEP]);
                    //errorCheck(k, lx, ly, m);
                }

            }
        }
    }

    for (int kx = 0; kx < II_STEP - 1; kx++) {
        for (int ky = 0; ky < JJ_STEP - 1; ky++) {
            int k = kx + ky * II_STEP;
            Q[k + II_STEP * JJ_STEP * 0] = kari[k + II_STEP * JJ_STEP * 0];
            Q[k + II_STEP * JJ_STEP * 1] = kari[k + II_STEP * JJ_STEP * 1];
            Q[k + II_STEP * JJ_STEP * 2] = kari[k + II_STEP * JJ_STEP * 2];
            Q[k + II_STEP * JJ_STEP * 3] = kari[k + II_STEP * JJ_STEP * 3];
        }
    }

    //printErrorMax();

    return;
}

void resetError() {
    error_max = 0;
}

void errorCheck(int k, int lx, int ly, int m) {
    if (lx <= 1 || lx >= II_STEP - 3) {
        if (ly <= 1 || ly >= JJ_STEP - 3) {
            return;
        }
    }
    if (error_max < fabs(kari[lx + ly * II_STEP + II_STEP * JJ_STEP * m] - Q[lx + ly * II_STEP + II_STEP * JJ_STEP * m])) {
        error_max = fabs(kari[lx + ly * II_STEP + II_STEP * JJ_STEP * m] - Q[lx + ly * II_STEP + II_STEP * JJ_STEP * m]);
        kkk = k;
        lxx = lx;
        lyy = ly;
        mmm = m;
    }
    return;
}

void printErrorMax(void) {
    printf("Error Max = %.8f, II = %d, JJ = %d, k = %d, m = %d\n", error_max, lxx, lyy, kkk, mmm);
    return;
}
