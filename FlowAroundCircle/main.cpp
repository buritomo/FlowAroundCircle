#include <stdio.h>
#include <stdlib.h>

#define MAIN_C_
//#define PARA

#include "global.h"
#include "grid.h"
#include "initial.h"
#include "boundary.h"
#include "calc.h"
#include "fds.h"
#include "rungekutta.h"
#include "output.h"
#include "viscose.h"

int main(void) {
	error_flag = 0;
	memorySet();
	cordinateDefine();
	metric();

	initialValue();
	boundaryValue();
	calcInternalValues();
	exportf();
	setExportBoundary();
	system("pause");
	setAssumedPotential();

	while (time <= TIME_MAX) {
		makePotential();
		fds(II_DIR);
		fds(JJ_DIR);
		viscose(II_DIR);
		viscose(JJ_DIR);
		rungekutta();
		inversePotentialToParams();
		boundaryValue();
		calcInternalValues();
		exportf();
		ErrorExport();
		printTimer();
		time = time + DELTA_T;
		exportBoundary();
		//exportf();
	}
	
	exportf();
	releaseAssumedPotential();
	releaseGrid();
	return 0;
}
