#include "hmm.h"

double Prospect(HMM* hmm, long T, int* testO, double **prospectalpha)
{
	setseed(getseed());
	long i, j, t, o, M, N;
	M = hmm->M;
	N = hmm->N;
	double sum;
	
	//prob �� 1 * M ����������¼ÿ��Ԥ��ʱ���۲�״̬���ָ���
	double* prob = dvector(1, M, "Prospect prob alloc");
	//prospectO �� 1 * T �����У���¼Ԥ��Ĺ۲����н��
	int* prospectO = ivector(1, T, "Prospect prospectO alloc");

	FILE* fp;
	errno_t err = 0;
	//err = fopen_s(&fp, "prob.seq", "w");
	//if (err) {
	//	fprintf(stderr, "File prob.seq open failed\n");
	//	exit(1);
	//}
	//fclose(fp);

	for (t = 1; t <= T; ++t) {
		sum = 0.0;
		//itoj �� 1 * N ����������¼ ״̬i��״̬j��ת�Ƹ��ʣ��� ״̬ת��ģ�͵ĵ�j��
		double* itoj = dvector(1, N, "Prospect itoj alloc");

		//������һʱ�̸�״̬����Ԥ����һʱ�̸�״̬����
		for (j = 1; j <= N; ++j) {

			for (i = 1; i <= N; ++i)
				itoj[i] = hmm->A[i][j];
			
			prospectalpha[t][j] = cblas_ddot(N, prospectalpha[t - 1] + 1, 1, itoj + 1, 1);
			sum += prospectalpha[t][j];
		}
		freedvector(itoj, 1, N);
		cblas_dscal(N, 1 / sum, prospectalpha[t] + 1, 1);

		//������۲�״̬���ָ���
		sum = 0.0;
		for (o = 1; o <= M; ++o) {
			prob[o] = cblas_ddot(N, prospectalpha[t] + 1, 1, hmm->B[o] + 1, 1);
			sum += prob[o];
		}
		cblas_dscal(M, 1 / sum, prob + 1, 1);

		double p = getrandom();
		sum = 0.0;
		for (o = 1; o <= M; ++o) {
			sum += prob[o];
			if (p <= sum)
				break;
		}
		prospectO[t] = o;
		/*
			ͷ�ļ�<mkl.h>
			����cblas_idamax()Ѱ�Ҿ���ֵ���Ԫ�ص��±�
			ԭ�ͣ�
			void cblas_dscal (const MKL_INT n, double *x, const MKL_INT incx);
			n		����Ԫ�ظ���
			x		���飬�����С����Ϊ ( 1 + ( n - 1 ) * abs( incx ) )
			incx	x��Ԫ���±������
			���ж�����Ԫ�ؾ���ֵ��ͬ�����ص�һ���±�
			������ NaN Ԫ�أ����� NaN ���±�
		*/
		////Ѱ�������ʳ��ֵĹ۲�ֵ-----------����ˣ������״̬�ֲ�������������������ʼ����һ��״̬����䣬���Բ���ֱ����������������
		//prospectO[t] = 1 + cblas_idamax(M, prob + 1, 1);

		////��ÿ�μ��������ֲ���ӡ�� prob.seq
		//err = fopen_s(&fp, "prob.seq", "a");
		//if (err) {
		//	fprintf(stderr, "File prob.seq open failed\n");
		//	exit(1);
		//}
		//for (o = 1; o <= M; ++o)
		//	fprintf(fp, "\t%lf", prob[o]);
		//fprintf(fp, "\n");
		//fclose(fp);
	}
	freedvector(prob, 1, M);

	//��ԭ������Ԥ�����д�ӡ�� prospect.seq
	err = fopen_s(&fp, "prospect.seq", "w");
	if (err) {
		fprintf(stderr, "File prospect.seq open failed\n");
		exit(1);
	}
	fprintf(fp, "ԭ����\tԤ��ֵ\n");
	for (t = 1; t <= T; ++t) {
		fprintf(fp, "\t%d\t%d", testO[t], prospectO[t]);
		if (testO[t] == prospectO[t])
			fprintf(fp, "\t����\n");
		else
			fprintf(fp, "\n");
	}
	fclose(fp);

	//����Ԥ��ֵ׼ȷ��
	long counter = 0;
	for (t = 1; t <= T; ++t)
		if (prospectO[t] == testO[t])
			++counter;
	
	return (((double)counter) / T);
}