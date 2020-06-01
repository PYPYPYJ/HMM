HMM 模型：
typedef struct
{
	int N;			//隐状态数
	int M;			//观测状态数
	double *pi;		//初始概率模型
	double **A;		//转移概率模型
	double **B;		//发射概率模型
} HMM;

MKL版本只改了C版本中以下部分内容：
1.forward.cpp中计算前向概率矩阵可用点乘函数cblas_ddot()，但需要额外时间取出A矩阵的某一列；归一化统一使用cblas_dscal()函数。
2.baumwelch.cpp中基本没变，只用cblas_scal()替换了归一化方法。gamma矩阵(各时刻各状态后验分布)与Xi矩阵(各时刻各状态与下一状态的联合后验分布)因为计算方式原因矩阵计算并不适用，考虑过每次迭代把要计算的向量取出来单独计算再赋值回去，但额外开销似乎更大，所以注释掉了。
