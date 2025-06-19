#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "time.h"
#include "omp.h"
#define DEBUG 1			/* printf("DEBUG %d\n", DEBUG); */


const char* rename_file(const char *filename, int index);
void decay_law(int *count, double *t2, const char* nome, const char* nome2);
void check_variance(int *count, int *sim, double *t2, const char* nome, const char* nome2, const char* nome3);
