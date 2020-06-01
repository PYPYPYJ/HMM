#include "hmm.h"

void Forward(HMM* hmm, long T, int* O, double** alpha)
{
	//先不写这个，反正不用
}

//alpha是一个T * N矩阵，记录1~T每个时刻中各状态的前向概率
double ForwardWithScale(HMM* hmm, long T, int* O, double** alpha, double* scale)
{
	long i, j, t, N;
	double tmp;
	N = hmm->N;

	//初始化第一行alpha
	scale[1] = tmp = 0.0;
	int o1 = O[1];
	for (i = 1; i <= N; ++i) {
		alpha[1][i] = hmm->pi[i] * (hmm->B[i][o1]);
		tmp += alpha[1][i];
	}
	scale[1] = tmp;

	/*
		头文件<mkl.h>
		函数cblas_dscal()计算一个标量乘以向量的结果，即 对一个向量进行放缩
		原型：
		void cblas_dscal (const MKL_INT n, const double a, double *x, const MKL_INT incx);
		n		向量元素个数
		a		标量
		x		数组，数组大小至少为 ( 1 + ( n - 1 ) * abs( incx ) )
		incx	x中元素下标的增量
	*/
	cblas_dscal(N, 1 / tmp, alpha[1] + 1, 1);

	//其余每一行alpha[t][j]为上一行各状态转移至状态j并输出O[t]的概率和
	for (t = 2; t <= T; ++t) {
		scale[t] = tmp = 0.0;
		int ot = O[t];

		//itoj 是 1 * N 的向量，记录 状态i到状态j的转移概率，即 状态转移模型的第j列
		double* itoj = dvector(1, N, "ForwardWithScale itoj");

		for (j = 1; j <= N; ++j) {
			for (i = 1; i <= N; ++i)
				itoj[i] = hmm->A[i][j];

			/*
			头文件<mkl.h>
			函数cblas_ddot()计算两个向量x、y的点乘，原型：
			double cblas_ddot (const MKL_INT n, const double *x, const MKL_INT incx, const double *y, const MKL_INT incy);
			n		两个向量元素个数
			x		数组，数组大小至少为 ( 1 + ( n - 1 ) * abs( incx ) )
			incx	x中元素下标的增量
			y		数组，数组大小至少为 ( 1 + ( n - 1 ) * abs( incy ) )
			incy	y中元素下标的增量

			若 n > 0 ，则返回点乘结果，否则返回0
			*/
			alpha[t][j] = cblas_ddot(N, alpha[t - 1] + 1, 1, itoj + 1, 1) * (hmm->B[j][ot]);
			
			tmp += alpha[t][j];
		}
		freedvector(itoj, 1, N);

		scale[t] = tmp;
		cblas_dscal(N, 1 / tmp, alpha[t] + 1, 1);;
	}

	//返回值logprobability为当前模型参数下观测序列出现概率的对数
	double logprobability = 0.0;
	for (t = 1; t < T; t += 2)
		logprobability += (log(scale[t]) + log(scale[t + 1]));
	for (; t <= T; ++t)
		logprobability += log(scale[t]);
	return logprobability;
}