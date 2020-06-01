#include "hmm.h"

extern bool first_iteration;
extern double logprobinit;

int BaumWelch(HMM *hmm, int T, int *O, double *logprobprev, double *logprobfinal)
{
	int i, j, k, t;
	double numeratorA, denominatorA, numeratorB, denominatorB;

	int iterationcounter = 0;
	double logprobtmp, logprobcurrent, delta;	//delta��¼���ε�������֮��
	const double limit = 0.001;		//���ε�������֮��С�ڴ�ֵʱ����ֹͣ

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

		//����ģ�Ͳ���
		//����pi����ʼ����ģ�͸���Ϊt=1ʱ�̸�״̬�ĺ������
		for (i = 1; i <= hmm->N; ++i)
			hmm->pi[i] = 0.001 + 0.999 *  gamma[1][i];
		//����A��B
		for (i = 1; i <= hmm->N; ++i) {
			//A״̬ת��ģ�ͣ�A[i][j]����Ϊ
			//1~T-1ʱ�� ��״̬i��״̬j�����Ϻ���ֲ��� / ״̬i�ĺ���ֲ��ͣ�
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; ++t)
				denominatorA += gamma[t][i];
			for (j = 1; j <= hmm->N; ++j) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; ++t)
					numeratorA += Xi[t][i][j];
				hmm->A[i][j] =0.001 + 0.999 *  numeratorA / denominatorA;
			}
			//B�������ģ�ͣ�B[i][k]����Ϊ
			//1~Tʱ�� �����Ϊkʱ״̬i�ĺ���ֲ��� / ״̬i�ĺ���ֲ��ͣ�
			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= hmm->M; ++k) {
				numeratorB = 0.0;
				for (t = 1; t <= T; ++t)
					if (O[t] == k)
						numeratorB += gamma[t][i];
				hmm->B[i][k] = 0.001 + 0.999 *  numeratorB / denominatorB;
			}
		}
		//��������������¼���alpha��beta��gamma��Xi
		logprobcurrent = ForwardWithScale(hmm, T, O, alpha, scale);
		BackwardWithScale(hmm, T, O, beta, scale);
		ComputeGamma(hmm, T, alpha, beta, gamma);
		ComputeXi(hmm, T, O, alpha, beta, Xi);
		//���µ�������
		++iterationcounter;
		//�������ε���֮��ĸ��ʲ�
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

//gamma��һ��T * N����gamma[t][i]��¼tʱ��i״̬�ĺ�����ʣ���alpha[t][i] * beta[t][i]
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

//Xi(����)��һ��T * N * N����Xi[t][i][j]��¼tʱ��Ϊi״̬��t+1ʱ��Ϊj״̬�ĺ������
//��alpha[t][i] * hmm->A[i][j] * hmm->B[j][O[t+1]] * beta[t+1][j]
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