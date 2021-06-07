#include <stdio.h>

#define MAIN_C_

#include "global.h"
#include "grid.h"
#include "initial.h"
#include "boundary.h"
#include "calc.h"
#include "fds.h"

int main(void) {
	memorySet();
	printf("%d, %d", II_STEP, JJ_STEP);
	cordinateDefine();
	metric();

	initialValue();
	boundaryValue();
	calcInternalValues();

	while (time <= TIME_MAX) {
		time = time + DELTA_T;
	}

	return 0;
}
