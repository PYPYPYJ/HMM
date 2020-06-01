#include <malloc.h>
#include "nrutil.h"

void callocerror(const char errortext[])
{
	fprintf(stderr, "malloc function failed : ");
	fprintf(stderr, errortext);
	exit(1);
}

int *ivector(long nl, long nh, const char errortext[])
{
	int *vector = (int *)calloc(nh - nl + 1, sizeof(int));
	if (!vector)
		callocerror(errortext);
	return vector - nl;
}

double *dvector(long nl, long nh, const char errortext[])
{
	double *vector = (double *)calloc(nh-nl+1, sizeof(double));
	if (!vector)
		callocerror(errortext);
	return vector - nl;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch, const char errortext[])
{
	int i;

	double **matrix = (double **)calloc(nrh - nrl + 1, sizeof(double *));
	if (!matrix)
		callocerror(errortext);
	matrix -= nrl;
	for (i = nrl; i <= nrh; ++i) {
		matrix[i] = (double *)calloc(nch - ncl + 1, sizeof(double));
		if (!matrix[i])
			callocerror(errortext);
		matrix[i] -= ncl;
	}
	return matrix;
}

void freeivector(int* vector, long nl, long nh)
{
	vector += nl;
	free(vector);
}

void freedvector(double *vector, long nl, long nh)
{
	vector += nl;
	free(vector);
}

void freedmatrix(double **matrix, long nrl, long nrh, long ncl, long nch)
{
	int i;
	for (i = nrh; i >= nrl; --i)
		free(matrix[i] + ncl);
	free(matrix + nrl);
	//printf("freedmatrix() complete\n");
}