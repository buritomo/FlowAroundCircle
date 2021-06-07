#include <stdio.h>
#include "grid.h"
#include "initial.h"
#include "boundary.h"
#include "calc.h"
#include "fds.h"

int main(void) {
	memorySet();
	cordinateDefine();
	metric();

	initialValue();
	doundaryValue();
	calcInternalValues();

	return 0;
}