#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "output.h"
#pragma warning(disable : 4996)

void exportf(void) {
    char filename[100], filedata[100];
    int step_cnt = (time / DELTA_T);
    int max_step_cnt = (int)(TIME_MAX / DELTA_T / 100000);

    if (step_cnt % max_step_cnt == 0) {
        FILE* fp;
        sprintf(filename, "Time%.6fsec.dat", time);
        sprintf(filedata, "Time%.6fsec.data", time);

        fp = fopen(filename, "w");

        fprintf(fp, "TITLE = \" %s \" \n", filedata);
        fprintf(fp, "TARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"e\", \"p\", \"H\", \"c\", \"X_xi\", \"Y_xi\", \"X_eta\", \"Y_eta\", \"J_inv\", \"X_cen\", \"Y_cen\", \"r\", \"theta\", \"XI_x\", \"XI_y\", \"ETA_x\", \"ETA_y\" \n");
        fprintf(fp, "ZONE T=\"th\", I=101,J=101 ,F=POINT\n");

        for (int ki = 2; ki < II_STEP - 2; ki++) {
            for (int kj = 2; kj < JJ_STEP - 2; kj++) {
                int k = ki + kj * II_STEP;
                fprintf(fp, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", x[k], y[k], rho[k], ux[k], vy[k], e[k], p[k], H[k], c[k], X_xi_half[k], Y_xi_half[k], X_eta_half[k], Y_eta_half[k], J_inv[k], x_cen[k], y_cen[k], rvec[k], thetaVec[k], XI_x[k], XI_y[k], Eta_x[k], Eta_y[k]);
                //printf("%d, %lf\n", kj, thetaVec[k]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    //return;
}

void ErrorExport(void) {
    char filename[100], filedata[100];
    int step_cnt = (time / DELTA_T);
    int max_step_cnt = (int)(TIME_MAX / DELTA_T / 10000);

    if (error_flag == 1) {
        FILE* fp;
        sprintf(filename, "Time%.6fsec.dat", time);
        sprintf(filedata, "Time%.6fsec.data", time);

        fp = fopen(filename, "w");

        fprintf(fp, "TITLE = \" %s \" \n", filedata);
        fprintf(fp, "TARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"e\", \"p\", \"H\", \"c\", \"X_xi\", \"Y_xi\", \"X_eta\", \"Y_eta\", \"J_inv\", \"X_cen\", \"Y_cen\", \"r\", \"theta\" \n");
        fprintf(fp, "ZONE T=\"th\", I=101,J=101 ,F=POINT\n");

        for (int ki = 2; ki < II_STEP - 2; ki++) {
            for (int kj = 2; kj < JJ_STEP - 2; kj++) {
                int k = ki + kj * II_STEP;
                fprintf(fp, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", x[k], y[k], rho[k], ux[k], vy[k], e[k], p[k], H[k], c[k], X_xi_half[k], Y_xi_half[k], X_eta_half[k], Y_eta_half[k], J_inv[k], x_cen[k], y_cen[k], rvec[k], thetaVec[k]);
                //printf("%d, %lf\n", kj, thetaVec[k]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);

        printf("Error Happen!!\n");
        exit(1);
    }

    return;
}

void setExportBoundary(void) {
    FILE* fp;

    fp = fopen("boundaryChange.csv", "w");

    fprintf(fp, "Time,");
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "ps_%d, ", kj);
    }
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "ux_%d, ", kj);
    }
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "vy_%d, ", kj);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "%.8f, ", time);
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "%.8f, ", p[k]);
    }
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "%.8f, ", ux[k]);
    }
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "%.8f, ", vy[k]);
    }
    fprintf(fp, "\n");

    fclose(fp);
}

void exportBoundary(void) {
    FILE* fp;

    fp = fopen("boundaryChange.csv", "a");

    if (fp == NULL) {
        perror("Cannnot open BoundaryChange!\n");
        exit(1);
    }

    fprintf(fp, "%.8f, ", time);
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "%.8f, ", p[k]);
    }
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "%.8f, ", ux[k]);
    }
    for (int kj = 2; kj < II_STEP - 3; kj++) {
        int k = II_STEP - 4 + kj * II_STEP;
        fprintf(fp, "%.8f, ", vy[k]);
    }
    fprintf(fp, "\n");

    fclose(fp);

}