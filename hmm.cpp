#include "hmm.h"

void ReadHMM(HMM* hmm, FILE* fp)
{
	int i, j, k, M, N;
	double sum, tmp;
	fscanf_s(fp, "M= %d\n", &M);
	fscanf_s(fp, "N= %d\n", &N);
	hmm->M = M;
	hmm->N = N;
	hmm->A = dmatrix(1, N, 1, N, "ReadHMM A alloc");
	fscanf_s(fp, "A:\n");
	for (i = 1; i <= N; ++i) {
		sum = 0.0;
		for (j = 1; j <= N; ++j) {
			fscanf_s(fp, "%lf", &tmp);
			sum += tmp;
			hmm->A[i][j] = tmp;
		}
		if (fabs(sum - 1.0) > 1e-5) {
			fprintf(stderr, "ReadHMM error : A[%d][], sum = %lf", i, sum);
			exit(1);
		}
		fscanf_s(fp, "\n");
	}
	hmm->B = dmatrix(1, N, 1, M, "ReadHMM B alloc");
	fscanf_s(fp, "B:\n");
	for (j = 1; j <= N; ++j) {
		sum = 0.0;
		for (k = 1; k <= M; ++k) {
			fscanf_s(fp, "%lf", &tmp);
			sum += tmp;
			hmm->B[j][k] = tmp;
		}
		if (fabs(sum - 1.0) > 1e-5) {
			fprintf(stderr, "ReadHMM error : B[%d][], sum = %lf", j, sum);
			exit(1);
		}
		fscanf_s(fp, "\n");
	}
	hmm->pi = dvector(1, N, "ReadHMM pi alloc");
	fscanf_s(fp, "pi:\n");
	sum = 0.0;
	for (i = 1; i <= N; ++i) {
		fscanf_s(fp, "%lf", &tmp);
		sum += tmp;
		hmm->pi[i] = tmp;
	}
	if (fabs(sum - 1.0) > 1e-5) {
		fprintf(stderr, "ReadHMM error : pi, sum = %lf", sum);
		exit(1);
	}
	fscanf_s(fp, "\n");
	//printHMM(hmm);
}

void InitHMM(HMM* hmm, int N, int M)
{
	setseed(getseed());		//设置随机数种子
	int i, j, k;
	double sum;		//sum用于归一化处理

	hmm->N = N;
	hmm->M = M;

	sum = 0.0;
	hmm->pi = dvector(1, N, "InitHMM pi alloc");
	for (i = 1; i <= N; ++i) {
		hmm->pi[i] = getrandom();
		sum += hmm->pi[i];
	}
	for (i = 1; i <= N; ++i)
		hmm->pi[i] /= sum;

	//i, j遍历A矩阵
	hmm->A = dmatrix(1, N, 1, N, "InitHMM A alloc");
	for (i = 1; i <= N; ++i) {
		sum = 0.0;
		for (j = 1; j <= N; ++j) {
			hmm->A[i][j] = getrandom();
			sum += hmm->A[i][j];
		}
		for (j = 1; j <= N; ++j)
			hmm->A[i][j] /= sum;
	}

	//j, k遍历B矩阵
	hmm->B = dmatrix(1, N, 1, M, "InitHMM B alloc");
	for (j = 1; j <= N; ++j) {
		sum = 0.0;
		for (k = 1; k <= M; ++k) {
			hmm->B[j][k] = getrandom();
			sum += hmm->B[j][k];
		}
		for (k = 1; k <= M; ++k)
			hmm->B[j][k] /= sum;
	}
}

void freeHMM(HMM* hmm)
{
	freedmatrix(hmm->A, 1, hmm->N, 1, hmm->N);
	freedmatrix(hmm->B, 1, hmm->N, 1, hmm->M);
	free(hmm->pi + 1);
}

void printHMM(HMM* hmm)
{
	int i, j, k;
	printf("number of states : %d\n", hmm->N);
	printf("number of symbols : %d\n", hmm->M);
	printf("pi:\n");
	for (i = 1; i <= hmm->N; ++i) {
		printf(" %f", hmm->pi[i]);
	}
	printf("\n");
	printf("A:\n");
	for (i = 1; i <= hmm->N; ++i) {
		for (j = 1; j <= hmm->N; ++j) {
			printf(" %f", hmm->A[i][j]);
		}
		printf("\n");
	}
	printf("B:\n");
	for (i = 1; i <= hmm->N; ++i) {
		for (k = 1; k <= hmm->M; ++k) {
			printf(" %f", hmm->B[i][k]);
		}
		printf("\n");
	}
}
