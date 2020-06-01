#include <malloc.h>
#include "nrutil.h"

//���������Ϣ
void callocerror(const char errortext[])
{
	fprintf(stderr, "malloc function failed : ");
	fprintf(stderr, errortext);
	exit(1);
}

//����һ�οռ䣬��СΪ (nh - nl + 1) * sizeof(int)
int* ivector(long nl, long nh, const char callfuncname[])
{
	int* vector = (int*)calloc(nh - nl + 1, sizeof(int));
	if (!vector)
		callocerror(callfuncname);
	return vector - nl;
}

//����һ�οռ䣬��СΪ (nh - nl + 1) * sizeof(double)
double* dvector(long nl, long nh, const char callfuncname[])
{
	double* vector = (double*)calloc(nh - nl + 1, sizeof(double));
	if (!vector)
		callocerror(callfuncname);
	return vector - nl;
}

//����һ�ζ�ά�ռ䣬��һά��СΪ (nrh - nrl + 1) * sizeof(double *)����һάÿ��Ԫ�ش�СΪ (nch - ncl + 1) * sizeof(double)
double** dmatrix(long nrl, long nrh, long ncl, long nch, const char callfuncname[])
{
	int i;

	double** matrix = (double**)calloc(nrh - nrl + 1, sizeof(double*));
	if (!matrix)
		callocerror(callfuncname);
	matrix -= nrl;
	for (i = nrl; i <= nrh; ++i) {
		matrix[i] = (double*)calloc(nch - ncl + 1, sizeof(double));
		if (!matrix[i])
			callocerror(callfuncname);
		matrix[i] -= ncl;
	}
	return matrix;
}

void freeivector(int* vector, long nl, long nh)
{
	vector += nl;
	free(vector);
}

void freedvector(double* vector, long nl, long nh)
{
	vector += nl;
	free(vector);
}

void freedmatrix(double** matrix, long nrl, long nrh, long ncl, long nch)
{
	int i;
	for (i = nrh; i > nrl; i -= 2) {
		matrix[i] += ncl;
		free(matrix[i]);
		matrix[i - 1] += ncl;
		free(matrix[i - 1]);
	}
	for (; i >= nrl; --i) {
		matrix[i] += ncl;
		free(matrix[i]);
	}
	matrix += nrl;
	free(matrix);
	//printf("freedmatrix() complete\n");
}