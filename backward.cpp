#include "hmm.h"

void Backward(HMM* hmm, long T, int* O, double** beta)
{
	//
	printf("Backward unfinished");
}

//beta是一个T * N矩阵，记录1~T每个时刻中各状态的后向概率
void BackwardWithScale(HMM* hmm, long T, int* O, double** beta, double* scale)
{
	long i, j, t, N;
	double sum, tmp;
	N = hmm->N;

	//初始化第T行beta
	tmp = 1.0 / scale[T];
	for (i = 1; i <= N; ++i)
		beta[T][i] = tmp;

	//其余每行beta[t][i]为状态i转移至下一时刻各状态并输出O[t+1]的概率和
	for (t = T - 1; t >= 1; --t) {
		int next_o = O[t + 1];
		for (i = 1; i <= N; ++i) {
			sum = 0.0;		
			for (j = 1; j <= N; ++j)
				sum += hmm->A[i][j] * hmm->B[j][next_o] * beta[t + 1][j];

			beta[t][i] = sum / scale[t];
		}
	}
}