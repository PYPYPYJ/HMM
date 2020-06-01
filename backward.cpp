#include "hmm.h"

void Backward(HMM* hmm, long T, int* O, double** beta)
{
	//
	printf("Backward unfinished");
}

//beta��һ��T * N���󣬼�¼1~Tÿ��ʱ���и�״̬�ĺ������
void BackwardWithScale(HMM* hmm, long T, int* O, double** beta, double* scale)
{
	long i, j, t, N;
	double sum, tmp;
	N = hmm->N;

	//��ʼ����T��beta
	tmp = 1.0 / scale[T];
	for (i = 1; i <= N; ++i)
		beta[T][i] = tmp;

	//����ÿ��beta[t][i]Ϊ״̬iת������һʱ�̸�״̬�����O[t+1]�ĸ��ʺ�
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