#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>

void callocerror(const char errortext[]);
int* ivector(long nl, long nh, const char callfuncname[]);
double* dvector(long nl, long nh, const char callfuncname[]);
double** dmatrix(long nrl, long nrh, long ncl, long nch, const char callfuncname[]);

void freeivector(int* vector, long nl, long nh);
void freedvector(double* vector, long nl, long nh);
void freedmatrix(double** matrix, long nrl, long nrh, long ncl, long nch);