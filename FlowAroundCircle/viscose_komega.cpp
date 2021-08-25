#include <math.h>
#include "global.h"
#include "viscose_komega.h"
#include "SST.h"

void viscouse_komega(int dir) {
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
            int kk = kx + ky * II_STEP;
            int pg;

            double Jaco_inv;
            double xi_x, xi_y, eta_x, eta_y;

            double k_xi, k_eta;
            double omega_xi, omega_eta;
            double k_dx, k_dy;
            double omega_dx, omega_dy;
            double k1, k2, k3, k4;

            double tmpT, Mu;
            double sigma_k, sigma_omega;
            double Mueff_k, Mueff_omega;


            //Jacobian
            if (dir == II_DIR) {
                Jaco_inv = ((J_inv[kk] + J_inv[kk + 1]) / 2);
                pg = kk + 1;
            }
            else {
                Jaco_inv = ((J_inv[kk] + J_inv[kk + II_STEP]) / 2);
                pg = kk + II_STEP;
            }

            //metric
            if (dir == II_DIR) {
                xi_x = Y_eta_half[kk + 1] / Jaco_inv;
                xi_y = -X_eta_half[kk + 1] / Jaco_inv;
                eta_x = -(Y_xi_half[kk + 1] + Y_xi_half[kk + 2] + Y_xi_half[kk - II_STEP + 1] + Y_xi_half[kk - II_STEP + 2]) * 0.25 / Jaco_inv;
                eta_y = (X_xi_half[kk + 1] + X_xi_half[kk + 2] + X_xi_half[kk - II_STEP + 1] + X_xi_half[kk - II_STEP + 2]) * 0.25 / Jaco_inv;
            }
            else {
                xi_x = (Y_eta_half[kk + II_STEP] + Y_eta_half[kk + II_STEP - 1] + Y_eta_half[kk + II_STEP * 2] + Y_eta_half[kk + II_STEP * 2 - 1]) * 0.25 / Jaco_inv;
                xi_y = -(X_eta_half[kk + II_STEP] + X_eta_half[kk + II_STEP - 1] + X_eta_half[kk + II_STEP * 2] + X_eta_half[kk + II_STEP * 2 - 1]) * 0.25 / Jaco_inv;
                eta_x = -Y_xi_half[kk + II_STEP] / Jaco_inv;
                eta_y = X_xi_half[kk + II_STEP] / Jaco_inv;
            }

            //TKE and omega Gradient
            if (dir == II_DIR) {
                k_xi = ux[kk + 1] - ux[kk];

                k1 = ux[kk + II_STEP] - ux[kk];
                k2 = ux[kk] - ux[kk - II_STEP];
                k3 = ux[kk + II_STEP + 1] - ux[kk + 1];
                k4 = ux[kk + 1] - ux[kk + 1 - II_STEP];
                k_eta = (k1 + k2 + k3 + k4) / 4;

                omega_xi = vy[kk + 1] - vy[kk];

                k1 = vy[kk + II_STEP] - vy[kk];
                k2 = vy[kk] - vy[kk - II_STEP];
                k3 = vy[kk + II_STEP + 1] - vy[kk + 1];
                k4 = vy[kk + 1] - vy[kk + 1 - II_STEP];
                omega_eta = (k1 + k2 + k3 + k4) / 4;
            }
            else {
                k1 = ux[kk + 1] - ux[kk];
                k2 = ux[kk] - ux[kk - 1];
                k3 = ux[kk + II_STEP + 1] - ux[kk + II_STEP];
                k4 = ux[kk + II_STEP] - ux[kk + II_STEP - 1];
                k_xi = (k1 + k2 + k3 + k4) / 4;

                k_eta = ux[kk + II_STEP] - ux[kk];

                k1 = vy[kk + 1] - vy[kk];
                k2 = vy[kk] - vy[kk - 1];
                k3 = vy[kk + II_STEP + 1] - vy[kk + II_STEP];
                k4 = vy[kk + II_STEP] - vy[kk + II_STEP - 1];
                omega_xi = (k1 + k2 + k3 + k4) / 4;

                omega_eta = vy[kk + II_STEP] - vy[kk];
            }

            k_dx = xi_x * k_xi + eta_x * k_eta;
            k_dy = xi_y * k_xi + eta_y * k_eta;
            omega_dy = xi_y * omega_xi + eta_y * omega_eta;
            omega_dx = xi_x * omega_xi + eta_x * omega_eta;

            //T, mu
            tmpT = calcTemp(kk);
            Mu = MU_0 * pow(tmpT / T_0, 1.5) * (T_0 + C) / (tmpT + C);;

            sigma_k = mix(kk, SIGMA_K_1, SIGMA_K_2);
            sigma_omega = mix(kk, SIGMA_OMEGA_2, SIGMA_OMEGA_2);

            Mueff_k = Mu + 0.5 * (Mut[kk] + Mut[pg]) * sigma_k;
            Mueff_omega = Mu + 0.5 * (Mut[kk] + Mut[pg]) * sigma_omega;

            if (dir == II_DIR) {
                TurbEv[kk + II_STEP * JJ_STEP * 0] = xi_x * Mueff_k * k_dx + xi_y * Mueff_k * k_dy;
                TurbEv[kk + II_STEP * JJ_STEP * 0] = xi_x * Mueff_omega * omega_dx + xi_y * Mueff_omega * omega_dy;
            }
            else {
                TurbFv[kk + II_STEP * JJ_STEP * 0] = eta_x * Mueff_k * k_dx + eta_y * Mueff_k * k_dy;
                TurbFv[kk + II_STEP * JJ_STEP * 0] = eta_x * Mueff_omega * omega_dx + eta_y * Mueff_omega * omega_dy;
            }
        }
    }
}

double calcTemp(int k) {
    return p[k] / (rho[k] * RAIR);
}