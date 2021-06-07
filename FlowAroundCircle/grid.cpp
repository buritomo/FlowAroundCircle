#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES//math.hを使えるようにするコマンド
#include <math.h>
#include "global.h"

void memorySet(void) {
	rho = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	ux = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	vy = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	e = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	p = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	H = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	c = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);

	Q = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	E = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	F = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	Ehalf = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);
	Fhalf = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP * 4);

	x = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	y = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	x_cen = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	y_cen = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);

	X_xi_half = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	Y_xi_half = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	X_eta_half = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	Y_eta_half = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);
	J_inv = (double*)malloc(sizeof(double) * II_STEP * JJ_STEP);

	return;
}

void cordinateDefine(void) {
	for (int ki = 0; ki < II_STEP; ki++) {
		for (int kj = 0; kj < JJ_STEP; kj++) {
			int k = ki + kj * II_STEP;
			double r = (((double)ki - 2.0) / (double)II_STEP - 5.0) * (((double)ki - 2.0) / (double)II_STEP - 5.0) * DELTA_S;
			double theta = 2 * M_PI * ((double)kj - 2.0) / ((double)JJ_STEP - 5.0);

			x_cen[k] = r * cos(theta);
			y_cen[k] = r * sin(theta);
		}
	}
	return;
}

void metric(void) {
	for (int kx = 0; kx < II_STEP - 1; kx++) {
		for (int kj = 0; kj < JJ_STEP - 1; kj++) {
			int k = kx + kj * II_STEP;
			x[k] = (x_cen[k] + x_cen[k + 1] + x_cen[k + II_STEP] + x_cen[k + II_STEP + 1]) / 4;
			y[k] = (y_cen[k] + y_cen[k + 1] + y_cen[k + II_STEP] + y_cen[k + II_STEP + 1]) / 4;
		}
	}

	for (int kx = 0; kx < II_STEP - 1; kx++) {
		for (int kj = 0; kj < JJ_STEP; kj++) {
			int k = kx + kj * II_STEP;
			X_xi_half[k] = x_cen[k + 1] - x_cen[k];
			Y_xi_half[k] = y_cen[k + 1] - y_cen[k];
		}
	}

	for (int kx = 0; kx < II_STEP; kx++) {
		for (int kj = 0; kj < JJ_STEP - 1; kj++) {
			int k = kx + kj * II_STEP;
			X_eta_half[k] = x_cen[k + II_STEP] - x_cen[k];
			Y_eta_half[k] = y_cen[k + II_STEP] - y_cen[k];
		}
	}

	for (int kx = 0; kx < II_STEP - 1; kx++) {
		for (int kj = 0; kj < JJ_STEP - 1; kj++) {
			int k = kx + kj * II_STEP;
			int kp = (kx + 1) + (kj + 1) * II_STEP;
			J_inv[k] = 0.5 * (x_cen[k + II_STEP + 1] - x_cen[k]) * (y_cen[k + II_STEP] - y_cen[k + 1]);
			J_inv[k] = J_inv[k] - 0.5 * (x_cen[k + II_STEP] - x_cen[k + 1]) * (y_cen[k + II_STEP + 1] - y_cen[k]);
		}
	}
	return;
}