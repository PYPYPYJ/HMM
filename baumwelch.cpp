#include "hmm.h"

extern bool first_iteration;
extern double logprobinit;

int BaumWelch(HMM *hmm, int T, int *O, double *logprobprev, double *logprobfinal)
{
	int i, j, k, t;
	double numeratorA, denominatorA, numeratorB, denominatorB;

	int iterationcounter = 0;
	double logprobtmp, logprobcurrent, delta;	//delta记录两次迭代概率之差
	const double limit = 0.001;		//两次迭代概率之差小于此值时迭代停止

	double **alpha, **beta, *scale, **gamma, ***Xi;
	alpha = dmatrix(1, T, 1, hmm->N, "BaumWelch alpha alloc");
	beta = dmatrix(1, T, 1, hmm->N, "BaumWelch beta alloc");
	scale = dvector(1, T, "BaumWelch scale alloc");
	gamma = dmatrix(1, T, 1, hmm->N, "BaumWelch gamma alloc");
	Xi = (double ***)calloc(T, sizeof(double **));
	if (!Xi)
		callocerror("baumwelch() callocXi");
	Xi -= 1;
	for (t = 1; t <= T; ++t)
		Xi[t] = dmatrix(1, hmm->N, 1, hmm->N, "BaumWelch Xi alloc");

	logprobcurrent = *logprobprev = ForwardWithScale(hmm, T, O, alpha, scale);
	if (first_iteration) {
		logprobinit = logprobcurrent;
		first_iteration = !first_iteration;
	}
	BackwardWithScale(hmm, T, O, beta, scale);
	ComputeGamma(hmm, T, alpha, beta, gamma);
	ComputeXi(hmm, T, O, alpha, beta, Xi);

	do
	{
		logprobtmp = logprobcurrent;

		//printHMM(hmm);
		//printABG(hmm->N, T, alpha, beta, gamma);

		//更新模型参数
		//更新pi，初始概率模型更新为t=1时刻各状态的后验概率
		for (i = 1; i <= hmm->N; ++i)
			hmm->pi[i] = 0.001 + 0.999 *  gamma[1][i];
		//更新A和B
		for (i = 1; i <= hmm->N; ++i) {
			//A状态转移模型，A[i][j]更新为
			//1~T-1时刻 （状态i到状态j的联合后验分布和 / 状态i的后验分布和）
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; ++t)
				denominatorA += gamma[t][i];
			for (j = 1; j <= hmm->N; ++j) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; ++t)
					numeratorA += Xi[t][i][j];
				hmm->A[i][j] =0.001 + 0.999 *  numeratorA / denominatorA;
			}
			//B发射概率模型，B[i][k]更新为
			//1~T时刻 （输出为k时状态i的后验分布和 / 状态i的后验分布和）
			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= hmm->M; ++k) {
				numeratorB = 0.0;
				for (t = 1; t <= T; ++t)
					if (O[t] == k)
						numeratorB += gamma[t][i];
				hmm->B[i][k] = 0.001 + 0.999 *  numeratorB / denominatorB;
			}
		}
		//更新完参数后重新计算alpha、beta、gamma、Xi
		logprobcurrent = ForwardWithScale(hmm, T, O, alpha, scale);
		BackwardWithScale(hmm, T, O, beta, scale);
		ComputeGamma(hmm, T, alpha, beta, gamma);
		ComputeXi(hmm, T, O, alpha, beta, Xi);
		//更新迭代次数
		++iterationcounter;
		//计算两次迭代之间的概率差
		delta = logprobcurrent - logprobtmp;
		
		//printABG(hmm->N, T, alpha, beta, gamma);
		//printHMM(hmm);
		//printf("%f\n", logprobprev);
		//printf("%f\n", logprobcurrent);
		//printf("%f\n", delta);

	} while (delta > limit);

	*logprobfinal = logprobcurrent;

	freedmatrix(alpha, 1, T, 1, hmm->N);
	freedmatrix(beta, 1, T, 1, hmm->N);
	freedvector(scale, 1, T);
	freedmatrix(gamma, 1, T, 1, hmm->N);
	freeXi(Xi, T, hmm->N);

	return iterationcounter;
}

//gamma是一个T * N矩阵，gamma[t][i]记录t时刻i状态的后验概率，即alpha[t][i] * beta[t][i]
void ComputeGamma(HMM *hmm, int T, double **alpha, double **beta, double **gamma)
{
	int i, t;
	double sum;

	for (t = 1; t <= T; ++t) {
		sum = 0.0;
		for (i = 1; i <= hmm->N; ++i) {
			gamma[t][i] = alpha[t][i] * beta[t][i];
			sum += gamma[t][i];
		}
		for (i = 1; i <= hmm->N; ++i)
			gamma[t][i] /= sum;
	}
}

//Xi(克西)是一个T * N * N矩阵，Xi[t][i][j]记录t时刻为i状态，t+1时刻为j状态的后验概率
//即alpha[t][i] * hmm->A[i][j] * hmm->B[j][O[t+1]] * beta[t+1][j]
void ComputeXi(HMM *hmm, int T, int *O, double **alpha, double **beta, double ***Xi)
{
	int i, j, t;
	double sum;

	for (t = 1; t <= T - 1; ++t) {
		sum = 0.0;
		for (i = 1; i <= hmm->N; ++i) {
			for (j = 1; j <= hmm->N; ++j) {
				Xi[t][i][j] = alpha[t][i] * hmm->A[i][j] * hmm->B[j][O[t + 1]] * beta[t + 1][j];
				sum += Xi[t][i][j];
			}
		}
		for (i = 1; i <= hmm->N; ++i)
			for (j = 1; j <= hmm->N; ++j)
				Xi[i][i][j] /= sum;
	}
}

void freeXi(double ***Xi, int T, int N)
{
	int t;
	for (t = 1; t <= T; ++t) {
		freedmatrix(Xi[t], 1, N, 1, N);
	}
	++Xi;
	free(Xi);
}

void printABG(int N, int T, double **alpha, double **beta, double **gamma)
{
	int i, t;
	printf("alpha:\n");
	for (t = 1; t <= T; ++t) {
		for (i = 1; i <= N; ++i) {
			printf(" %f", alpha[t][i]);
		}
		printf("\n");
	}
	printf("beta:\n");
	for (t = 1; t <= T; ++t) {
		for (i = 1; i <= N; ++i) {
			printf(" %f", beta[t][i]);
		}
		printf("\n");
	}
	printf("gamma:\n");
	for (t = 1; t <= T; ++t) {
		for (i = 1; i <= N; ++i) {
			printf(" %f", gamma[t][i]);
		}
		printf("\n");
	}
}