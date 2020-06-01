#include "hmm.h"

void Forward(HMM *hmm, int T, int *O, double **alpha)
{
	//先不写这个，反正不用
}

//alpha是一个T * N矩阵，记录1~T每个时刻中各状态的前向概率
double ForwardWithScale(HMM *hmm, int T, int *O, double **alpha, double *scale)
{
	int i, j, t;
	double sum;

	//初始化第一行alpha
	scale[1] = 0.0;
	for (i = 1; i <= hmm->N; ++i) {
		alpha[1][i] = hmm->pi[i] * (hmm->B[i][O[1]]);
		scale[1] += alpha[1][i];
	}
	for (i = 1; i <= hmm->N; ++i)
		alpha[1][i] /= scale[1];

	//其余每一行alpha[t][j]为上一行各状态转移至状态j并输出O[t]的概率和
	for (t = 2; t <= T; ++t) {
		scale[t] = 0.0;
		for (j = 1; j <= hmm->N; ++j) {
			sum = 0.0;
			for (i = 1; i <= hmm->N; ++i)
				sum += alpha[t - 1][i] * (hmm->A[i][j]);

			alpha[t][j] = sum * (hmm->B[j][O[t]]);
			scale[t] += alpha[t][j];
		}
		for (j = 1; j <= hmm->N; ++j)
			alpha[t][j] /= scale[t];
	}

	//返回值logprobability为当前模型参数下观测序列出现概率的对数
	double logprobability = 0.0;
	for (t = 1; t <= T; ++t)
		logprobability += log(scale[t]);
	return logprobability;
}