#ifndef OUTPUT_H_
#define OUTPUT_H_

#define FNAME "errormax.csv"

void exportf(void);
void ErrorExport(void);
//void ExportTarget(double* data, char* filename, int x_limit, int y_limit);
void setExportBoundary(void);
void exportBoundary(void);
void makeErrorFile(void);
void ExportErrorMax(double errorMax);

#endif //OUTPUT_H_