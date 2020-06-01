#include "hmm.h"

void Forward(HMM* hmm, long T, int* O, double** alpha)
{
	//�Ȳ�д�������������
}

//alpha��һ��T * N���󣬼�¼1~Tÿ��ʱ���и�״̬��ǰ�����
double ForwardWithScale(HMM* hmm, long T, int* O, double** alpha, double* scale)
{
	long i, j, t, N;
	double tmp;
	N = hmm->N;

	//��ʼ����һ��alpha
	scale[1] = tmp = 0.0;
	int o1 = O[1];
	for (i = 1; i <= N; ++i) {
		alpha[1][i] = hmm->pi[i] * (hmm->B[i][o1]);
		tmp += alpha[1][i];
	}
	scale[1] = tmp;

	/*
		ͷ�ļ�<mkl.h>
		����cblas_dscal()����һ���������������Ľ������ ��һ���������з���
		ԭ�ͣ�
		void cblas_dscal (const MKL_INT n, const double a, double *x, const MKL_INT incx);
		n		����Ԫ�ظ���
		a		����
		x		���飬�����С����Ϊ ( 1 + ( n - 1 ) * abs( incx ) )
		incx	x��Ԫ���±������
	*/
	cblas_dscal(N, 1 / tmp, alpha[1] + 1, 1);

	//����ÿһ��alpha[t][j]Ϊ��һ�и�״̬ת����״̬j�����O[t]�ĸ��ʺ�
	for (t = 2; t <= T; ++t) {
		scale[t] = tmp = 0.0;
		int ot = O[t];

		//itoj �� 1 * N ����������¼ ״̬i��״̬j��ת�Ƹ��ʣ��� ״̬ת��ģ�͵ĵ�j��
		double* itoj = dvector(1, N, "ForwardWithScale itoj");

		for (j = 1; j <= N; ++j) {
			for (i = 1; i <= N; ++i)
				itoj[i] = hmm->A[i][j];

			/*
			ͷ�ļ�<mkl.h>
			����cblas_ddot()������������x��y�ĵ�ˣ�ԭ�ͣ�
			double cblas_ddot (const MKL_INT n, const double *x, const MKL_INT incx, const double *y, const MKL_INT incy);
			n		��������Ԫ�ظ���
			x		���飬�����С����Ϊ ( 1 + ( n - 1 ) * abs( incx ) )
			incx	x��Ԫ���±������
			y		���飬�����С����Ϊ ( 1 + ( n - 1 ) * abs( incy ) )
			incy	y��Ԫ���±������

			�� n > 0 ���򷵻ص�˽�������򷵻�0
			*/
			alpha[t][j] = cblas_ddot(N, alpha[t - 1] + 1, 1, itoj + 1, 1) * (hmm->B[j][ot]);
			
			tmp += alpha[t][j];
		}
		freedvector(itoj, 1, N);

		scale[t] = tmp;
		cblas_dscal(N, 1 / tmp, alpha[t] + 1, 1);;
	}

	//����ֵlogprobabilityΪ��ǰģ�Ͳ����¹۲����г��ָ��ʵĶ���
	double logprobability = 0.0;
	for (t = 1; t < T; t += 2)
		logprobability += (log(scale[t]) + log(scale[t + 1]));
	for (; t <= T; ++t)
		logprobability += log(scale[t]);
	return logprobability;
}