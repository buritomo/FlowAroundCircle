#include <stdio.h>
#include <stdlib.h>

#define MAIN_C_

#include "global.h"
#include "grid.h"
#include "initial.h"
#include "boundary.h"
#include "calc.h"
#include "fds.h"
#include "rungekutta.h"
#include "output.h"

int main(void) {
	error_flag = 0;
	memorySet();
	cordinateDefine();
	metric();

	initialValue();
	boundaryValue();
	calcInternalValues();
	exportf();
	system("pause");
	setAssumedPotential();

	while (time <= TIME_MAX) {
		makePotential();
		fds(II_DIR);
		fds(JJ_DIR);
		rungekutta();
		inversePotentialToParams();
		boundaryValue();
		exportf();
		ErrorExport();
		printTimer();
		time = time + DELTA_T;
		//exportf();
	}
	
	exportf();
	releaseAssumedPotential();
	releaseGrid();
	return 0;
}
