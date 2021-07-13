#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES//math.hを使えるようにするコマンド
#include <math.h>
#include "global.h"
#include "viscose.h"

void viscose(int dir) {
    //計算範囲
    int kx_max, ky_max;
    int kx_min, ky_min;

    if (dir == II_DIR) {
        kx_min = 1;
        kx_max = II_STEP - 2;
        ky_min = 2;
        ky_max = JJ_STEP - 2;
    }
    else {
        kx_min = 2;
        kx_max = II_STEP - 2;
        ky_min = 1;
        ky_max = JJ_STEP - 2;
    }

    for (int kx = kx_min; kx < kx_max; kx++) {
        for (int ky = ky_min; ky < ky_max; ky++) {
            int k = kx + ky * II_STEP;

            caocJacobian(k, dir);
            calcMetric(k, dir);
            calcVelocityGra(k, dir);
            calcVeloBoundary(k, dir);
            //calcTemp(k, dir);
            calcCoef(k, dir);
            calcTGra(k, dir);

            tau_xx = 2 / 3 * mu_ave * (2 * ux_dx - vy_dy);
            tau_xy = mu_ave * (ux_dy + vy_dx);
            tau_yy = 2 / 3 * mu_ave * (2 * vy_dy - ux_dx);

            beta_x = tau_xx * ux_half + tau_xy * vy_half + kappa_ave * T_dx;
            beta_y = tau_xy * ux_half + tau_yy * vy_half + kappa_ave * T_dy;

            calcViscoseFactor(k, dir);
        }
    }

    return;
}

void caocJacobian(int k, int dir) {
    if (dir == II_DIR) {
        Jaco_inv = ((J_inv[k] + J_inv[k + 1]) / 2);
    }
    else {
        Jaco_inv = ((J_inv[k] + J_inv[k + II_STEP]) / 2);
    }
}

void calcMetric(int k, int dir) {
    //実際に求めたい境界の番号と，使う番号が計算方向に一つずれることに注意
    if (dir == II_DIR) {
        xi_x = Y_eta_half[k + 1] / Jaco_inv;
        xi_y = -X_eta_half[k + 1] / Jaco_inv;
        eta_x = -(Y_xi_half[k + 1] + Y_xi_half[k + 2] + Y_xi_half[k - II_STEP + 1] + Y_xi_half[k - II_STEP + 2]) * 0.25 / Jaco_inv;
        eta_y = (X_xi_half[k + 1] + X_xi_half[k + 2] + X_xi_half[k - II_STEP + 1] + X_xi_half[k - II_STEP + 2]) * 0.25 / Jaco_inv;
    }
    else {
        xi_x = (Y_eta_half[k + II_STEP] + Y_eta_half[k + II_STEP - 1] + Y_eta_half[k + II_STEP * 2] + Y_eta_half[k + II_STEP * 2 - 1]) * 0.25 / Jaco_inv;
        xi_y = -(X_eta_half[k + II_STEP] + X_eta_half[k + II_STEP - 1] + X_eta_half[k + II_STEP * 2] + X_eta_half[k + II_STEP * 2 - 1]) * 0.25 / Jaco_inv;
        eta_x = -Y_xi_half[k + II_STEP] / Jaco_inv;
        eta_y = X_xi_half[k + II_STEP] / Jaco_inv;
    }
    return;
}

void calcVelocityGra(int k, int dir) {
    double ux_xi, ux_eta, vy_xi, vy_eta;
    double k1, k2, k3, k4;
    if (dir == II_DIR) {
        ux_xi = ux[k + 1] - ux[k];

        k1 = ux[k + II_STEP] - ux[k];
        k2 = ux[k] - ux[k - II_STEP];
        k3 = ux[k + II_STEP + 1] - ux[k + 1];
        k4 = ux[k + 1] - ux[k + 1 - II_STEP];
        ux_eta = (k1 + k2 + k3 + k4) / 4;

        vy_xi = vy[k + 1] - vy[k];

        k1 = vy[k + II_STEP] - vy[k];
        k2 = vy[k] - vy[k - II_STEP];
        k3 = vy[k + II_STEP + 1] - vy[k + 1];
        k4 = vy[k + 1] - vy[k + 1 - II_STEP];
        vy_eta = (k1 + k2 + k3 + k4) / 4;
    }
    else {
        k1 = ux[k + 1] - ux[k];
        k2 = ux[k] - ux[k - 1];
        k3 = ux[k + II_STEP + 1] - ux[k + II_STEP];
        k4 = ux[k + II_STEP] - ux[k + II_STEP - 1];
        ux_xi = (k1 + k2 + k3 + k4) / 4;

        ux_eta = ux[k + II_STEP] - ux[k];

        k1 = vy[k + 1] - vy[k];
        k2 = vy[k] - vy[k - 1];
        k3 = vy[k + II_STEP + 1] - vy[k + II_STEP];
        k4 = vy[k + II_STEP] - vy[k + II_STEP - 1];
        vy_xi = (k1 + k2 + k3 + k4) / 4;

        vy_eta = vy[k + II_STEP] - vy[k];
    }

    ux_dx = xi_x * ux_xi + eta_x * ux_eta;
    ux_dy = xi_y * ux_xi + eta_y * ux_eta;
    vy_dy = xi_y * vy_xi + eta_y * vy_eta;
    vy_dx = xi_x * vy_xi + eta_x * vy_eta;

    return;
}

void calcVeloBoundary(int k, int dir) {
    if (dir == II_DIR) {
        ux_half = (ux[k] + ux[k + 1]) / 2;
        vy_half = (vy[k] + vy[k + 1]) / 2;
    }
    else {
        ux_half = (ux[k] + ux[k + II_STEP]) / 2;
        vy_half = (vy[k] + vy[k + II_STEP]) / 2;
    }
    return;
}

void calcViscoseFactor(int k, int dir) {
    if (dir == II_DIR) {
        Ev[k + II_STEP * JJ_STEP * 0] = 0;
        Ev[k + II_STEP * JJ_STEP * 1] = (xi_x * tau_xx + xi_y * tau_xy) * Jaco_inv;
        Ev[k + II_STEP * JJ_STEP * 2] = (xi_x * tau_xy + xi_y * tau_yy) * Jaco_inv;
        Ev[k + II_STEP * JJ_STEP * 3] = (xi_x * beta_x + xi_y * beta_y) * Jaco_inv;
    }
    else {
        Fv[k + II_STEP * JJ_STEP * 0] = 0;
        Fv[k + II_STEP * JJ_STEP * 1] = (eta_x * tau_xx + eta_y * tau_xy) * Jaco_inv;
        Fv[k + II_STEP * JJ_STEP * 2] = (eta_x * tau_xy + eta_y * tau_yy) * Jaco_inv;
        Fv[k + II_STEP * JJ_STEP * 3] = (eta_x * beta_x + eta_y * beta_y) * Jaco_inv;
    }
    return;
}

void calcCoef(int k, int dir) {
    double T_a, T_b;
    double mu_1, mu_2;
    double kappa_1, kappa_2;
    double cp = GAMMA * RAIR / (GAMMA - 1);

    if (dir == II_DIR) {
        T_a = calcTemp(k);
        T_b = calcTemp(k + 1);
    }
    else {
        T_a = calcTemp(k);
        T_b = calcTemp(k + II_STEP);
    }

    mu_1 = MU_0 * pow(T_a / T_0, 1.5) * (T_0 + C) / (T_a + C);
    mu_2 = MU_0 * pow(T_b / T_0, 1.5) * (T_0 + C) / (T_b + C);
    mu_ave = 0.5 * (mu_1 + mu_2);
    kappa_ave = cp * mu_ave / PR;

    return;
}

double calcTemp(int k){
//void calcTemp(int k, int dir) {
    double cp = GAMMA * RAIR / (GAMMA - 1);

    return H[k] / cp;
    /*
    if (dir == II_DIR) {
        T_a = H[k] / cp;
        T_b = H[k + 1] / cp;
    }
    else {
        T_a = H[k] / cp;
        T_b = H[k + II_STEP] / cp;
    }
    T = (T_a + T_b) / 2;
    return;*/
}

void calcTGra(int k, int dir) {/////////
    double T_xi, T_eta;
    double T_1_a, T_1_b;
    double T_2_a, T_2_b, T_2_c, T_2_d;
    double k1, k2, k3, k4;

    if (dir == II_DIR) {
        T_xi = calcTemp(k + 1) - calcTemp(k);

        k1 = calcTemp(k + II_STEP) - calcTemp(k);
        k2 = calcTemp(k) - calcTemp(k - II_STEP);
        k3 = calcTemp(k + II_STEP + 1) - calcTemp(k + 1);
        k4 = calcTemp(k + 1) - calcTemp(k + 1 - II_STEP);
        T_eta = (k1 + k2 + k3 + k4) / 4;
    }
    else {
        k1 = calcTemp(k + 1) - calcTemp(k);
        k2 = calcTemp(k) - calcTemp(k - 1);
        k3 = calcTemp(k + II_STEP + 1) - calcTemp(k + II_STEP);
        k4 = calcTemp(k + II_STEP) - calcTemp(k + II_STEP - 1);
        T_xi = (k1 + k2 + k3 + k4) / 4;

        T_eta = calcTemp(k + II_STEP) - calcTemp(k);
    }

    T_dx = xi_x * T_xi + eta_x * T_eta;
    T_dy = xi_y * T_xi + eta_y * T_eta;

    return;
}