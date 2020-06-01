#include "hmm.h"

void Forward(HMM *hmm, int T, int *O, double **alpha)
{
	//�Ȳ�д�������������
}

//alpha��һ��T * N���󣬼�¼1~Tÿ��ʱ���и�״̬��ǰ�����
double ForwardWithScale(HMM *hmm, int T, int *O, double **alpha, double *scale)
{
	int i, j, t;
	double sum;

	//��ʼ����һ��alpha
	scale[1] = 0.0;
	for (i = 1; i <= hmm->N; ++i) {
		alpha[1][i] = hmm->pi[i] * (hmm->B[i][O[1]]);
		scale[1] += alpha[1][i];
	}
	for (i = 1; i <= hmm->N; ++i)
		alpha[1][i] /= scale[1];

	//����ÿһ��alpha[t][j]Ϊ��һ�и�״̬ת����״̬j�����O[t]�ĸ��ʺ�
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

	//����ֵlogprobabilityΪ��ǰģ�Ͳ����¹۲����г��ָ��ʵĶ���
	double logprobability = 0.0;
	for (t = 1; t <= T; ++t)
		logprobability += log(scale[t]);
	return logprobability;
}