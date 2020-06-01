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
	double logprobtmp, logprobcurrent, delta;	//delta��¼���ε�������֮��
	const double limit = 0.001;		//���ε�������֮��С�ڴ�ֵʱ����ֹͣ

	double** alpha, ** beta, * scale, ** gamma, *** Xi;
	alpha = dmatrix(1, T, 1, N, "BaumWelch alpha init");					//ǰ����ʾ���T * N
	beta = dmatrix(1, T, 1, N, "BaumWelch beta init");						//������ʾ���T * N
	scale = dvector(1, T, "BaumWelch scale init");							//�������ӣ�1 * N
	gamma = dmatrix(1, T, 1, N, "BaumWelch gamma init");					//��״̬�ĺ���ֲ���T * N
	Xi = (double***)calloc(T, sizeof(double**));	//������״̬�����Ϻ���ֲ���T * N * N
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

		//����ģ�Ͳ���
		//����pi����ʼ����ģ�͸���Ϊt=1ʱ�̸�״̬�ĺ������
		for (i = 1; i <= N; ++i)
			hmm->pi[i] = 0.001 + 0.999 * gamma[1][i];
		//����A��B
		for (i = 1; i <= N; ++i) {
			//A״̬ת��ģ�ͣ�A[i][j]����Ϊ
			//1~T-1 ʱ�� ��״̬i��״̬j�����Ϻ���ֲ��� / ״̬i�ĺ���ֲ��ͣ�
			denominatorA = 0.0;
			for (t = 1; t < T - 1; t += 2)		//���� denominatorA �õ���2 x 1a loop unrolling ����
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
			//B�������ģ�ͣ�B[i][k]����Ϊ
			//1~Tʱ�� �����Ϊkʱ״̬i�ĺ���ֲ��� / ״̬i�ĺ���ֲ��ͣ�
			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= M; ++k) {
				numeratorB = 0.0;
				for (t = 1; t <= T; ++t)
					if (O[t] == k)
						numeratorB += gamma[t][i];
				hmm->B[i][k] = 0.001 + 0.999 * numeratorB / denominatorB;
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

//gamma��һ��T * N����gamma[t][i]��¼tʱ��i״̬�ĺ�����ʣ���alpha[t][i] * beta[t][i]
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

		//��һ��
		cblas_dscal(N, 1 / sum, gamma[t] + 1, 1);
	}
}
/*
	Xi(����)��һ��T * N * N����Xi[t][i][j]��¼tʱ��Ϊi״̬��t+1ʱ��Ϊj״̬�ĺ������
	��alpha[t][i] * hmm->A[i][j] * hmm->B[j][O[t+1]] * beta[t+1][j]
*/
void ComputeXi(HMM* hmm, long T, int* O, double** alpha, double** beta, double*** Xi)
{
	long i, j, t, N;
	double sum;
	N = hmm->N;

	for (t = 1; t <= T - 1; ++t) {
		sum = 0.0;
		int next_o = O[t + 1];

		////joutȡ���������ģ�͵ĵ�O[t+1]�У��� ��״̬���O[t+1]�ĸ���
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
	�����ã��ֱ��ӡalpha��beta��gamma����
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