#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "fds.h"


void fds(int dir) {
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

            //‹«ŠE‚Å‚Ìƒ„ƒRƒrƒAƒ“ŒvŽZ
            if (dir == II_DIR) {
                J_bd = 1 / ((J_inv[k] + J_inv[k + 1]) / 2);
            }
            else {
                J_bd = 1 / ((J_inv[k] + J_inv[k + II_STEP]) / 2);
            }
            /*
            if (k == 316) {
                printf("!!!");
            }*/

            muscl(&rho_L, &rho_R, rho, k, dir);
            muscl(&u_L, &u_R, ux, k, dir);
            muscl(&v_L, &v_R, vy, k, dir);
            muscl(&p_L, &p_R, p, k, dir);
            muscl(&e_L, &e_R, e, k, dir);
            muscl(&H_L, &H_R, H, k, dir);
            muscl(&c_L, &c_R, c, k, dir);

            RoeAverage();
            if (dir == II_DIR) {
                k_x = Y_eta_half[k + 1] * J_bd;
                k_y = -X_eta_half[k + 1] * J_bd;
            }
            else {
                k_x = -Y_xi_half[k + II_STEP] * J_bd;
                k_y = X_xi_half[k + II_STEP] * J_bd;
            }

            Z = k_x * u_ave + k_y * v_ave;
            b1 = (u_ave * u_ave + v_ave * v_ave) / 2 * (GAMMA - 1) / c_ave / c_ave;
            b2 = (GAMMA - 1) / c_ave / c_ave;

            RightArray(R);
            LeftArray(L);
            Lambda(lambda);
            musclArray(k, dir);

            for (int l = 0; l < 4; l++) {
                for (int m = 0; m < 4; m++) {
                    w[l + 4 * m] = 0.0;
                    for (int n = 0; n < 4; n++) {
                        w[l + 4 * m] = w[l + 4 * m] + R[l + 4 * n] * lambda[n] * L[n + 4 * m];
                    }
                }
            }
            if (dir == II_DIR) {
                for (int l = 0; l < 4; l++) {
                    Ehalf[k + II_STEP * JJ_STEP * l] = 0.5 * (E_L[l] + E_R[l]);
                    for (int m = 0; m < 4; m++) {
                        Ehalf[k + II_STEP * JJ_STEP * l] = Ehalf[k + II_STEP * JJ_STEP * l] - 0.5 * w[l + 4 * m] * (Q_R[m] - Q_L[m]);
                    }
                }
            }
            else {
                for (int l = 0; l < 4; l++) {
                    Fhalf[k + II_STEP * JJ_STEP * l] = 0.5 * (F_L[l] + F_R[l]);
                    for (int m = 0; m < 4; m++) {
                        Fhalf[k + II_STEP * JJ_STEP * l] = Fhalf[k + II_STEP * JJ_STEP * l] - 0.5 * w[l + 4 * m] * (Q_R[m] - Q_L[m]);
                    }
                }
            }
        }
    }

    return;
}

void RightArray(double* R) {
    double k_x_var = k_x / sqrt(k_x * k_x + k_y * k_y);
    double k_y_var = k_y / sqrt(k_x * k_x + k_y * k_y);
    double Z_var = Z / sqrt(k_x * k_x + k_y * k_y);
    double theta_var = theta / sqrt(k_x * k_x + k_y * k_y);

    R[0 + 4 * 0] = 1;
    R[0 + 4 * 1] = 1;
    R[0 + 4 * 2] = 1;
    R[0 + 4 * 3] = 0;
    R[1 + 4 * 0] = u_ave - k_x_var * c_ave;
    R[1 + 4 * 1] = u_ave;
    R[1 + 4 * 2] = u_ave + k_x_var * c_ave;
    R[1 + 4 * 3] = -k_y_var;
    R[2 + 4 * 0] = v_ave - k_y_var * c_ave;
    R[2 + 4 * 1] = v_ave;
    R[2 + 4 * 2] = v_ave + k_y_var * c_ave;
    R[2 + 4 * 3] = k_x_var;
    R[3 + 4 * 0] = H_ave - c_ave * Z_var;
    R[3 + 4 * 1] = 0.5 * (u_ave * u_ave + v_ave * v_ave);
    R[3 + 4 * 2] = H_ave + c_ave * Z_var;
    R[3 + 4 * 3] = -(k_y_var * u_ave - k_x_var * v_ave);

    return;
}

void LeftArray(double* L) {
    double k_x_var = k_x / sqrt(k_x * k_x + k_y * k_y);
    double k_y_var = k_y / sqrt(k_x * k_x + k_y * k_y);
    double Z_var = Z / sqrt(k_x * k_x + k_y * k_y);
    double theta_var = theta / sqrt(k_x * k_x + k_y * k_y);

    L[0 + 4 * 0] = 0.5 * (b1 + Z_var / c_ave);
    L[0 + 4 * 1] = -0.5 * (k_x_var / c_ave + b2 * u_ave);
    L[0 + 4 * 2] = -0.5 * (k_y_var / c_ave + b2 * v_ave);
    L[0 + 4 * 3] = 0.5 * b2;
    L[1 + 4 * 0] = 1 - b1;
    L[1 + 4 * 1] = b2 * u_ave;
    L[1 + 4 * 2] = b2 * v_ave;
    L[1 + 4 * 3] = -1 * b2;
    L[2 + 4 * 0] = 0.5 * (b1 - Z_var / c_ave);
    L[2 + 4 * 1] = 0.5 * (k_x_var / c_ave - b2 * u_ave);
    L[2 + 4 * 2] = 0.5 * (k_y_var / c_ave - b2 * v_ave);
    L[2 + 4 * 3] = 0.5 * b2;
    L[3 + 4 * 0] = k_y_var * u_ave - k_x_var * v_ave;
    L[3 + 4 * 1] = -k_y_var;
    L[3 + 4 * 2] = k_x_var;
    L[3 + 4 * 3] = 0;

    return;
}

void Lambda(double* lambda) {
    lambda[0] = fabs(Z - c_ave * sqrt(k_x * k_x + k_y * k_y));
    lambda[1] = fabs(Z);
    lambda[2] = fabs(Z + c_ave * sqrt(k_x * k_x + k_y * k_y));
    lambda[3] = fabs(Z);

    return;
}

void musclArray(int k, int dir) {
    Q_L[0] = rho_L / J_bd;
    Q_L[1] = rho_L * u_L / J_bd;
    Q_L[2] = rho_L * v_L / J_bd;
    Q_L[3] = e_L / J_bd;

    Q_R[0] = rho_R / J_bd;
    Q_R[1] = rho_R * u_R / J_bd;
    Q_R[2] = rho_R * v_R / J_bd;
    Q_R[3] = e_R / J_bd;

    if (dir == JJ_DIR) {
        double U_L = (k_x * u_L + k_y * v_L);
        F_L[0] = rho_L * U_L / J_bd;
        F_L[1] = (rho_L * u_L * U_L + k_x * p_L) / J_bd;
        F_L[2] = (rho_L * v_L * U_L + k_y * p_L) / J_bd;
        F_L[3] = (e_L + p_L) * U_L / J_bd;

        double U_R = (k_x * u_R + k_y * v_R);
        F_R[0] = rho_R * U_R / J_bd;
        F_R[1] = (rho_R * u_R * U_R + k_x * p_R) / J_bd;
        F_R[2] = (rho_R * v_R * U_R + k_y * p_R) / J_bd;
        F_R[3] = (e_R + p_R) * U_R / J_bd;
    }
    else {
        double V_L = (k_x * u_L + k_y * v_L);
        E_L[0] = rho_L * V_L / J_bd;
        E_L[1] = (rho_L * u_L * V_L + k_x * p_L) / J_bd;
        E_L[2] = (rho_L * v_L * V_L + k_y * p_L) / J_bd;
        E_L[3] = (e_L + p_L) * V_L / J_bd;

        double V_R = (k_x * u_R + k_y * v_R);
        E_R[0] = rho_R * V_R / J_bd;
        E_R[1] = (rho_R * u_R * V_R + k_x * p_R) / J_bd;
        E_R[2] = (rho_R * v_R * V_R + k_y * p_R) / J_bd;
        E_R[3] = (e_R + p_R) * V_R / J_bd;
    }

    return;
}

void RoeAverage(void) {
    rho_ave = sqrt(rho_L * rho_R);
    u_ave = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R) / ((sqrt(rho_L) + sqrt(rho_R)));
    v_ave = (sqrt(rho_L) * v_L + sqrt(rho_R) * v_R) / ((sqrt(rho_L) + sqrt(rho_R)));
    H_ave = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R) / ((sqrt(rho_L) + sqrt(rho_R)));
    c_ave = sqrt((GAMMA - 1) * (H_ave - 0.5 * (u_ave * u_ave + v_ave * v_ave)));
    if (isnan(c_ave)) {
        printf("Here!\n");
    }

    return;
}

void muscl(double* left, double* right, double* value, int k, int dir) {
    int delta;
    if (dir == II_DIR) {
        delta = 1;
    }
    else {
        delta = II_STEP;
    }

    double delta_j_plus = value[k + delta] - value[k];
    double delta_j_minus = value[k] - value[k - delta];
    double delta_J1_plus = value[k + 2 * delta] - value[k + delta];
    double b = (3 - KAPPA) / (1 - KAPPA);

    *left = value[k] + EPSILON / 4 * ((1 - KAPPA) * (limiter(delta_j_minus, b * delta_j_plus)) + (1 + KAPPA) * (limiter(delta_j_plus, b * delta_j_minus)));
    *right = value[k + delta] - EPSILON / 4 * ((1 + KAPPA) * (limiter(delta_j_plus, b * delta_J1_plus)) + (1 - KAPPA) * (limiter(delta_J1_plus, b * delta_j_plus)));

    return;
}

double limiter(double x, double y) {
    return (sgn(x) * max(0, min(fabs(x), sgn(x) * y)));
}

double max(double x, double y) {
    if (x > y) {
        return x;
    }
    else {
        return y;
    }
}

double min(double x, double y) {
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