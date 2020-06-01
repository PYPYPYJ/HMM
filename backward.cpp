#include "hmm.h"

void Backward(HMM *hmm, int T, int *O, double **beta)
{
	//
}

//beta��һ��T * N���󣬼�¼1~Tÿ��ʱ���и�״̬�ĺ������
void BackwardWithScale(HMM *hmm, int T, int *O, double **beta, double *scale)
{
	int i, j, t;
	double sum;

	//��ʼ����T��beta
	for (i = 1; i <= hmm->N; ++i)
		beta[T][i] = 1.0 / scale[T];

	//����ÿ��beta[t][i]Ϊ״̬iת������һʱ�̸�״̬�����O[t+1]�ĸ��ʺ�
	for (t = T - 1; t >= 1; --t) {
		for (i = 1; i <= hmm->N; ++i) {
			sum = 0.0;
			for (j = 1; j <= hmm->N; ++j)
				sum += hmm->A[i][j] * hmm->B[j][O[t + 1]] * beta[t+1][j];

			beta[t][i] = sum / scale[t];
		}
	}
}