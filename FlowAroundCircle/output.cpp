#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "output.h"
#pragma warning(disable : 4996)

void exportf(void) {
    char filename[100], filedata[100];
    int step_cnt = (time / DELTA_T);
    int max_step_cnt = (int)(TIME_MAX / DELTA_T / 10000);

    if (step_cnt % max_step_cnt == 0) {
        FILE* fp;
        sprintf(filename, "Time%.6fsec.txt", time);
        sprintf(filedata, "Time%.6fsec.data", time);

        fp = fopen(filename, "w");

        fprintf(fp, "TITLE = \" %s \" \n", filedata);
        fprintf(fp, "TARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"e\", \"p\", \"H\", \"c\", \"X_xi\", \"Y_xi\", \"X_eta\", \"Y_eta\", \"J_inv\" \n");
        fprintf(fp, "ZONE T=\"th\", I=101,J=101 ,F=POINT\n");

        for (int ki = 2; ki < II_STEP - 2; ki++) {
            for (int kj = 2; kj < JJ_STEP - 2; kj++) {
                int k = ki + kj * II_STEP;
                fprintf(fp, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", x[k], y[k], rho[k], ux[k], vy[k], e[k], p[k], H[k], c[k], X_xi_half[k], Y_xi_half[k], X_eta_half[k], Y_eta_half[k], J_inv[k]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    return;
}