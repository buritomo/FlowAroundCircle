#include <stdio.h>
#include <stdlib.h>

#define MAIN_C_

#include "global.h"
#include "grid.h"
#include "initial.h"
#include "boundary.h"
#include "calc.h"
#include "fds.h"
//#include "rungekutta.h"
#include "LU-SGS.h"
#include "output.h"
#include "viscose.h"

int main(void) {
	error_flag = 0;
	memorySet();
	cordinateDefine();
	metric();
	setMatrixs();

	initialValue();
	boundaryValue();
	calcInternalValues();
	exportf();
	makeErrorFile();
	//setExportBoundary();
	system("pause");
	//setAssumedPotential();

	while (time <= TIME_MAX) {
		makePotential();
		//fds(II_DIR);
		//fds(JJ_DIR);
		//viscose(II_DIR);
		//viscose(JJ_DIR);
		//rungekutta();
		//inversePotentialToParams();
		//boundaryValue();
		GaussSeidel();
		calcInternalValues();
		exportf();
		ErrorExport();
		//printTimer();
		time = time + DELTA_T;
		//exportBoundary();
		//exportf();
	}
	
	exportf();
	//releaseAssumedPotential();
	releaseGrid();
	return 0;
}
