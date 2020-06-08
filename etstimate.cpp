#include "hmm.h"

#define min(x, y) (((x) < (y))?(x):(y))
#define cut_length 262144		//2^18 = 262144
constexpr double train_rate = 0.85;			//可观测序列中用于 训练/测试 的比例

bool first_iteration = true;
double logprobinit;				//初始模型前向概率的对数

int main(int argc, char* argv[])
{
	double starttime, endtime;
	starttime = dsecnd();

	double logprobprev, logprobfinal;		//调用 baumwelch() 前后模型的评价
	HMM hmm;
	long T, totallength;
	int *O, M, N;
	int iterationcounter = 0;

	/*
	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;
	while (c = getopt(argc, argv, ))
	{
		//getopt()是GNU库函数，要include <getopt.h>
		//Linux 处理命令行参数用
	}*/
	FILE* fp;
	errno_t err = 0;
	if (argc != 3 && argc != 6) {
		//fprintf(stderr, "argc = %d\n", argc);
		fprintf(stderr, "Usage error\n");
		fprintf(stderr, "Usage: %s -N <states numbers> -M <observation numbers> <obs.seq>\n", argv[0]);
		fprintf(stderr, "Usage: %s <model.hmm> <obs.seq>\n", argv[0]);
		exit(1);
	}
	else if (argc == 3) {		//从文件读取已设置好的hmm
		err = fopen_s(&fp, argv[1], "r");
		if (err) {
			fprintf(stderr, "File %s open failed\n", argv[1]);
			exit(1);
		}
		ReadHMM(&hmm, fp);
		fclose(fp);

		err = fopen_s(&fp, argv[2], "r");
		if (err) {
			fprintf(stderr, "File %s open failed\n", argv[2]);
			exit(1);
		}
	}
	else {			//随机初始化hmm
		InitHMM(&hmm, atoi(argv[2]), atoi(argv[4]));

		//fprintf(stdout, "number of states : %d\n", atoi(argv[2]));
		//fprintf(stdout, "number of symbols : %d\n", atoi(argv[4]));
		err = fopen_s(&fp, argv[5], "r");
		//fprintf(stdout, "%d\n", err);
		if (err) {
			fprintf(stderr, "File %s open failed\n", argv[5]);
			exit(1);
		}
	}

	//从文件读取观测序列，序列总长度为totallength，其中T用于训练，每次最多读取cut_length个数
	fscanf_s(fp, "T= %d\n", &totallength);
	T = train_rate * totallength;
	M = hmm.M;
	N = hmm.N;
	double** alpha, ** beta, * scale;
	alpha = dmatrix(1, cut_length, 1, N, "BaumWelch alpha alloc");			//前向概率矩阵，T * N
	beta = dmatrix(1, cut_length, 1, N, "BaumWelch beta alloc");			//后向概率矩阵，T * N
	scale = dvector(1, cut_length, "BaumWelch scale alloc");				//缩放因子，1 * T
	long i, tmp;
	long leftlength = T + cut_length;
	do
	{
		leftlength -= cut_length;
		tmp = min(leftlength, cut_length);
		O = ivector(1, tmp, "main O init");

		for (i = 1; i < tmp; i += 2) {
			fscanf_s(fp, "%d", &O[i]);
			fscanf_s(fp, "%d", &O[i + 1]);
		}
		for (; i <= tmp; ++i) {
			fscanf_s(fp, "%d", &O[i]);
		}
		iterationcounter += BaumWelch(&hmm, tmp, O, alpha, beta, scale, &logprobprev, &logprobfinal);
		fprintf(stdout, "iteration: %d\tprobinit:%e | %e\n", iterationcounter, logprobprev, logprobfinal);
		freeivector(O, 1, tmp);
	} while (i <= leftlength);

	fprintf(stdout, "-------------------------\n");
	fprintf(stdout, "Iteration numbers:\t%d\n", iterationcounter);
	fprintf(stdout, "Log probability of initial model:\t%E\n", logprobinit);
	fprintf(stdout, "Log probability of final model:\t%E\n", logprobfinal);
	fprintf(stdout, "-------------------------\n");
	printHMM(&hmm);
	
	endtime = dsecnd() - starttime;
	fprintf(stdout, "------------------------\n");
	fprintf(stdout, "total estimate time : %f s\n", endtime);

	/*
		剩余长度 testlength 用于测试
	*/
	long testlength = totallength - T;
	int* testO = ivector(1, testlength, "main test alloc");
	for (i = 1; i < testlength; i += 2) {
		fscanf_s(fp, "%d", &testO[i]);
		fscanf_s(fp, "%d", &testO[i + 1]);
	}
	for (; i <= testlength; ++i) {
		fscanf_s(fp, "%d", &testO[i]);
	}
	fclose(fp);

	//训练完成后保存 alpha 矩阵最后一行用于预测
	double** prospectalpha = dmatrix(0, testlength, 1, N, "main prospect alloc");
	for (i = 1; i <= N; ++i)
		prospectalpha[0][i] = alpha[tmp][i];
	//fprintf(stdout, "prospect第%ld行：", tmp);
	//for (i = 1; i <= N; ++i)
	//	fprintf(stdout, " %f", prospectalpha[0][i]);
	//fprintf(stdout, "\n");

	freedmatrix(alpha, 1, cut_length, 1, N);
	freedmatrix(beta, 1, cut_length, 1, N);
	freedvector(scale, 1, cut_length);

	double prospectaccuracy = Prospect(&hmm, testlength, testO, prospectalpha);
	prospectaccuracy *= 100;

	fprintf(stdout, "-------------------------\n");
	fprintf(stdout, "prospect accuracy : %f \%\n", prospectaccuracy);

	freeHMM(&hmm);

	
	
	//getchar();
}