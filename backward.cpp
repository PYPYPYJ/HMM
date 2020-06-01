#include "hmm.h"

void Backward(HMM *hmm, int T, int *O, double **beta)
{
	//
}

//beta是一个T * N矩阵，记录1~T每个时刻中各状态的后向概率
void BackwardWithScale(HMM *hmm, int T, int *O, double **beta, double *scale)
{
	int i, j, t;
	double sum;

	//初始化第T行beta
	for (i = 1; i <= hmm->N; ++i)
		beta[T][i] = 1.0 / scale[T];

	//其余每行beta[t][i]为状态i转移至下一时刻各状态并输出O[t+1]的概率和
	for (t = T - 1; t >= 1; --t) {
		for (i = 1; i <= hmm->N; ++i) {
			sum = 0.0;
			for (j = 1; j <= hmm->N; ++j)
				sum += hmm->A[i][j] * hmm->B[j][O[t + 1]] * beta[t+1][j];

			beta[t][i] = sum / scale[t];
		}
	}
}