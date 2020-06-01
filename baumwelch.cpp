#include "hmm.h"

extern bool first_iteration;
extern double logprobinit;

int BaumWelch(HMM* hmm, long T, int* O, double* logprobprev, double* logprobfinal)
{
	long i, j, k, t, M, N;
	M = hmm->M;
	N = hmm->N;

	double numeratorA, denominatorA, numeratorB, denominatorB;

	int iterationcounter = 0;
	double logprobtmp, logprobcurrent, delta;	//delta记录两次迭代概率之差
	const double limit = 0.001;		//两次迭代概率之差小于此值时迭代停止

	double** alpha, ** beta, * scale, ** gamma, *** Xi;
	alpha = dmatrix(1, T, 1, N, "BaumWelch alpha init");					//前向概率矩阵，T * N
	beta = dmatrix(1, T, 1, N, "BaumWelch beta init");						//后向概率矩阵，T * N
	scale = dvector(1, T, "BaumWelch scale init");							//缩放因子，1 * N
	gamma = dmatrix(1, T, 1, N, "BaumWelch gamma init");					//各状态的后验分布，T * N
	Xi = (double***)calloc(T, sizeof(double**));	//相邻两状态的联合后验分布，T * N * N
	if (!Xi)
		callocerror("BaumWelch Xi init");
	Xi -= 1;
	for (t = 1; t < T; t += 2) {
		Xi[t] = dmatrix(1, N, 1, N, "BaumWelch Xi init 2");
		Xi[t + 1] = dmatrix(1, N, 1, N, "BaumWelch Xi init 2");
	}
	for (; t <= T; ++t)
		Xi[t] = dmatrix(1, N, 1, N, "BaumWelch Xi init 2");

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
		//printABG(N, T, alpha, beta, gamma);

		//更新模型参数
		//更新pi，初始概率模型更新为t=1时刻各状态的后验概率
		for (i = 1; i <= N; ++i)
			hmm->pi[i] = 0.001 + 0.999 * gamma[1][i];
		//更新A和B
		for (i = 1; i <= N; ++i) {
			//A状态转移模型，A[i][j]更新为
			//1~T-1 时刻 （状态i到状态j的联合后验分布和 / 状态i的后验分布和）
			denominatorA = 0.0;
			for (t = 1; t < T - 1; t += 2)		//计算 denominatorA 用到了2 x 1a loop unrolling 技术
				denominatorA += (gamma[t][i] + gamma[t + 1][i]);
			for (; t <= T - 1; ++t)
				denominatorA += gamma[t][i];

			for (j = 1; j <= N; ++j) {
				numeratorA = 0.0;
				for (t = 1; t < T - 1; t += 2)
					numeratorA += (Xi[t][i][j] + Xi[t + 1][i][j]);
				for (; t <= T - 1; ++t)
					numeratorA += Xi[t][i][j];
				hmm->A[i][j] = 0.001 + 0.999 * numeratorA / denominatorA;
			}
			//B发射概率模型，B[i][k]更新为
			//1~T时刻 （输出为k时状态i的后验分布和 / 状态i的后验分布和）
			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= M; ++k) {
				numeratorB = 0.0;
				for (t = 1; t <= T; ++t)
					if (O[t] == k)
						numeratorB += gamma[t][i];
				hmm->B[i][k] = 0.001 + 0.999 * numeratorB / denominatorB;
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

		//printABG(N, T, alpha, beta, gamma);
		//printHMM(hmm);
		//printf("%f\n", logprobprev);
		//printf("%f\n", logprobcurrent);
		//printf("%f\n", delta);

	} while (delta > limit);

	*logprobfinal = logprobcurrent;

	freedmatrix(alpha, 1, T, 1, N);
	freedmatrix(beta, 1, T, 1, N);
	freedvector(scale, 1, T);
	freedmatrix(gamma, 1, T, 1, N);
	freeXi(Xi, T, N);

	return iterationcounter;
}

//gamma是一个T * N矩阵，gamma[t][i]记录t时刻i状态的后验概率，即alpha[t][i] * beta[t][i]
void ComputeGamma(HMM* hmm, long T, double** alpha, double** beta, double** gamma)
{
	long i, t, N;
	double sum;
	N = hmm->N;

	for (t = 1; t <= T; ++t) {
		sum = 0.0;
		for (i = 1; i <= N; ++i) {
			gamma[t][i] = alpha[t][i] * beta[t][i];
			sum += gamma[t][i];
		}

		//归一化
		cblas_dscal(N, 1 / sum, gamma[t] + 1, 1);
	}
}
/*
	Xi(克西)是一个T * N * N矩阵，Xi[t][i][j]记录t时刻为i状态，t+1时刻为j状态的后验概率
	即alpha[t][i] * hmm->A[i][j] * hmm->B[j][O[t+1]] * beta[t+1][j]
*/
void ComputeXi(HMM* hmm, long T, int* O, double** alpha, double** beta, double*** Xi)
{
	long i, j, t, N;
	double sum;
	N = hmm->N;

	for (t = 1; t <= T - 1; ++t) {
		sum = 0.0;
		int next_o = O[t + 1];

		////jout取出发射概率模型的第O[t+1]列，即 各状态输出O[t+1]的概率
		//double* jout = dvector(1, N, "ComputeXi jout alloc");
		//double* tmp = dvector(1, N, "ComputeXi tmp alloc");

		for (i = 1; i <= N; ++i) {

			//for (j = 1; j <= N; ++j)
			//	jout[j] = hmm->B[j][next_o];

			//double* pa, * pb;
			//pa = hmm->A[i];
			//pb = beta[t + 1];
			//for (j = 1; j <= N; ++j) {
			//	tmp[j] = pa[j] * jout[j] * pb[j];
			//}
			//cblas_dscal(N, alpha[t][i], tmp + 1, 1);

			for (j = 1; j <= N; ++j) {
				Xi[t][i][j] = alpha[t][i] * hmm->A[i][j] * hmm->B[j][next_o] * beta[t + 1][j];
				sum += Xi[t][i][j];
				//Xi[t][i][j] = tmp[j];
				//sum += tmp[j];
			}
		}
		//freedvector(jout, 1, N);
		//freedvector(tmp, 1, N);

		for (i = 1; i <= N; ++i)
			cblas_dscal(N, 1 / sum, Xi[t][i], 1);
	}
}

void freeXi(double*** Xi, long T, int N)
{
	long t;
	for (t = 1; t < T; t += 2) {
		freedmatrix(Xi[t], 1, N, 1, N);
		freedmatrix(Xi[t + 1], 1, N, 1, N);
	}
	for (; t <= T; ++t)
		freedmatrix(Xi[t], 1, N, 1, N);
	++Xi;
	free(Xi);
}

/*
	测试用，分别打印alpha、beta、gamma矩阵
*/
void printABG(int N, long T, double** alpha, double** beta, double** gamma)
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