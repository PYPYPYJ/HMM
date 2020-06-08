#include "hmm.h"

double Prospect(HMM* hmm, long T, int* testO, double **prospectalpha)
{
	setseed(getseed());
	long i, j, t, o, M, N;
	M = hmm->M;
	N = hmm->N;
	double sum;
	
	//prob 是 1 * M 的向量，记录每次预测时各观测状态出现概率
	double* prob = dvector(1, M, "Prospect prob alloc");
	//prospectO 是 1 * T 的序列，记录预测的观测序列结果
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
		//itoj 是 1 * N 的向量，记录 状态i到状态j的转移概率，即 状态转移模型的第j列
		double* itoj = dvector(1, N, "Prospect itoj alloc");

		//根据上一时刻各状态概率预测下一时刻各状态概率
		for (j = 1; j <= N; ++j) {

			for (i = 1; i <= N; ++i)
				itoj[i] = hmm->A[i][j];
			
			prospectalpha[t][j] = cblas_ddot(N, prospectalpha[t - 1] + 1, 1, itoj + 1, 1);
			sum += prospectalpha[t][j];
		}
		freedvector(itoj, 1, N);
		cblas_dscal(N, 1 / sum, prospectalpha[t] + 1, 1);

		//计算各观测状态出现概率
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
			头文件<mkl.h>
			函数cblas_idamax()寻找绝对值最大元素的下标
			原型：
			void cblas_dscal (const MKL_INT n, double *x, const MKL_INT incx);
			n		向量元素个数
			x		数组，数组大小至少为 ( 1 + ( n - 1 ) * abs( incx ) )
			incx	x中元素下标的增量
			若有多个最大元素绝对值相同，返回第一个下标
			若包含 NaN 元素，返回 NaN 的下标
		*/
		////寻找最大概率出现的观测值-----------想错了，后面各状态分布会收敛，导致最大概率始终是一个状态不会变，所以不能直接用最大概率来估计
		//prospectO[t] = 1 + cblas_idamax(M, prob + 1, 1);

		////将每次计算的输出分布打印至 prob.seq
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

	//将原序列与预测序列打印至 prospect.seq
	err = fopen_s(&fp, "prospect.seq", "w");
	if (err) {
		fprintf(stderr, "File prospect.seq open failed\n");
		exit(1);
	}
	fprintf(fp, "原序列\t预测值\n");
	for (t = 1; t <= T; ++t) {
		fprintf(fp, "\t%d\t%d", testO[t], prospectO[t]);
		if (testO[t] == prospectO[t])
			fprintf(fp, "\t命中\n");
		else
			fprintf(fp, "\n");
	}
	fclose(fp);

	//计算预测值准确度
	long counter = 0;
	for (t = 1; t <= T; ++t)
		if (prospectO[t] == testO[t])
			++counter;
	
	return (((double)counter) / T);
}