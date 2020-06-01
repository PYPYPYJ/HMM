#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void callocerror(const char errortext[]);
int *ivector(long nl, long nh, const char errortext[]);
double *dvector(long nl, long nh, const char errortext[]);
double **dmatrix(long nrl, long nrh, long ncl, long nch, const char errortext[]);

void freeivector(int* vector, long nl, long nh);
void freedvector(double *vector, long nl, long nh);
void freedmatrix(double **matrix, long nrl, long nrh, long ncl, long nch);