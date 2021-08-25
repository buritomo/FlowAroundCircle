#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
//#include "fds.h"
#include "fds-komega.h"

void fds_komega(int dir) {
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
            double LEFT = 1.0, RIGHT = -1.0;

            double J_bd;
            int mg, pg, ppg;
            double k_x, k_y;

            double rho_L, rho_R;
            double u_L, u_R;
            double v_L, v_R;
            double k_L, k_R;
            double omega_L, omega_R;
            double Vn_L, Vn_R;

            double rho_ave;
            double u_ave;
            double v_ave;
            double k_ave;
            double omega_ave;
            double UU;

            double EE[2];
            double RAL[2][2];


            if (dir == II_DIR) {
                J_bd = 1 / ((J_inv[kk] + J_inv[kk + 1]) / 2);
                mg = kk - 1;
                pg = kk + 1;
                ppg = kk + 2;
                k_x = Y_eta_half[kk + 1] * J_bd;
                k_y = -X_eta_half[kk + 1] * J_bd;
         
            }
            else {
                J_bd = 1 / ((J_inv[kk] + J_inv[kk + II_STEP]) / 2);
                mg = kk - II_STEP;
                pg = kk + II_STEP;
                ppg = kk + 2 * II_STEP;
                k_x = -Y_xi_half[kk + II_STEP] * J_bd;
                k_y = X_xi_half[kk + II_STEP] * J_bd;
            }

            rho_L = KWmuscl(rho, mg, kk, pg, LEFT);
            u_L = KWmuscl(ux, mg, kk, pg, LEFT);
            v_L = KWmuscl(vy, mg, kk, pg, LEFT);
            k_L = KWmuscl(tke, mg, kk, pg, LEFT);
            omega_L = KWmuscl(omega, mg, kk, pg, LEFT);

            rho_R = KWmuscl(rho, kk, pg, ppg, RIGHT);
            u_R = KWmuscl(ux, kk, pg, ppg, RIGHT);
            v_R = KWmuscl(vy, kk, pg, ppg, RIGHT);
            k_R = KWmuscl(tke, kk, pg, ppg, RIGHT);
            omega_R = KWmuscl(omega, kk, pg, ppg, RIGHT);

            rho_ave = sqrt(rho_L * rho_R);
            u_ave = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R) / ((sqrt(rho_L) + sqrt(rho_R)));
            v_ave = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R) / ((sqrt(rho_L) + sqrt(rho_R)));
            k_ave = (sqrt(rho_L) * k_L + sqrt(rho_R) * k_R) / ((sqrt(rho_L) + sqrt(rho_R)));
            omega_ave = (sqrt(rho_L) * omega_L + sqrt(rho_R) * omega_R) / ((sqrt(rho_L) + sqrt(rho_R)));

            Vn_L = (u_L * k_x + v_L * k_y) / J_bd;
            Vn_R = (u_R * k_x + v_R * k_y) / J_bd;

            EE[0] = rho_R * k_R * Vn_R + rho_L * k_L * Vn_L;
            EE[1] = rho_R * omega_R * Vn_R + rho_L * omega_L * Vn_L;

            UU = (k_x * u_ave + k_y * v_ave);

            RAL[0][0] = fabs(UU);
            RAL[0][1] = 0.0;
            RAL[1][0] = 0.0;
            RAL[1][1] = fabs(UU);

            if (dir == II_DIR) {
                for (int l = 0; l < 2; l++) {
                    TurbE[kk + II_STEP * JJ_STEP * l] = 0.5 * EE[l] - RAL[l][0] * (rho_R * k_R - rho_L * k_L) - RAL[l][1] * (rho_R * omega_R - rho_L * omega_L);
                }
            }
            else {
                for (int l = 0; l < 2; l++) {
                    TurbF[kk + II_STEP * JJ_STEP * l] = 0.5 * EE[l] - RAL[l][0] * (rho_R * k_R - rho_L * k_L) - RAL[l][1] * (rho_R * omega_R - rho_L * omega_L);
                }
            }

        }
    }
}

double KWmuscl(double *value, int mg, int kk, int pg, double side) {//fds.h‚Å‚ÌmusclŠÖ”‚Ìo—ˆ‚ªˆ«‚¢‚Ì‚Åì‚è’¼‚µ
    return value[kk] + side * EPSILON / 4 * ((1 - side * KAPPA) * KWlimiter(value[kk] - value[mg], (3 - KAPPA) / (1 - KAPPA) * (value[pg] - value[kk])) + (1 + side * KAPPA) * KWlimiter((value[pg] - value[kk]), (3 - KAPPA) / (1 - KAPPA) * (value[kk] - value[mg])));
}

double KWlimiter(double x, double y) {
    return (sgn(x) * max2f(0, min2f(fabs(x), sgn(x) * y)));
}

double max2f(double x, double y) {
    if (x > y) {
        return x;
    }
    else {
        return y;
    }
}

double min2f(double x, double y) {
    if (x > y) {
        return y;
    }
    else {
        return x;
    }
}

double sgn(double x) {
    if (x > 0) {
        return 1.0;
    }
    else if (x < 0) {
        return -1.0;
    }
    else {
        return 0.0;
    }
}