#include <stdio.h>
#include "global.h"
#include "SST.h"

double komega_ini(void) {

}

double mix(int k, double a_1, double a_2) {
	return a_1 * F1(k) + a_2 * (1 - F1(k));
}

double F1(int k) {
	for (int i = 0; i < II_STEP; i++) {
		for(int j = 0; j < JJ_STEP; j++) {

		}
	}
}