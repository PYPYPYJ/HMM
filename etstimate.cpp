#include "hmm.h"

#define min(x, y) (((x) < (y))?(x):(y))
#define cut_length 262144		//2^18 = 262144

bool first_iteration = true;
double logprobinit;				//��ʼģ��ǰ����ʵĶ���

int main(int argc, char* argv[])
{
	double starttime, endtime;
	starttime = dsecnd();

	double logprobprev, logprobfinal;		//���� baumwelch() ǰ��ģ�͵�����
	HMM hmm;
	long T;
	int *O;
	int iterationcounter = 0;

	/*
	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;
	while (c = getopt(argc, argv, ))
	{
		//getopt()��GNU�⺯����Ҫinclude <getopt.h>
		//Linux ���������в�����
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
	else if (argc == 3) {		//���ļ���ȡ�����úõ�hmm
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
	else {			//�����ʼ��hmm
		InitHMM(&hmm, atoi(argv[2]), atoi(argv[4]));

		//printf("number of states : %d\n", atoi(argv[2]));
		//printf("number of symbols : %d\n", atoi(argv[4]));
		err = fopen_s(&fp, argv[5], "r");
		//printf("%d\n", err);
		if (err) {
			fprintf(stderr, "File %s open failed\n", argv[5]);
			exit(1);
		}
	}

	//���ļ���ȡ�۲����У������ܳ���ΪT��ÿ������ȡcut_length����
	fscanf_s(fp, "T= %d\n", &T);
	long i;
	long leftlength = T + cut_length;
	do
	{
		leftlength -= cut_length;
		long tmp = min(leftlength, cut_length);
		O = ivector(1, tmp, "main read O from file");

		for (i = 1; i <= tmp; ++i) {
			fscanf_s(fp, "%d", &O[i]);
		}
		iterationcounter += BaumWelch(&hmm, tmp, O, &logprobprev, &logprobfinal);
		printf("iteration: %d\tprobinit:%e | %e\n", iterationcounter, logprobprev, logprobfinal);
		freeivector(O, 1, tmp);
	} while (i <= leftlength);

	//ReadSequence(fp, &T, &O);
	fclose(fp);

	//printf("sequence length : %d\n", T);
	//printHMM(&hmm);

	printf("-------------------------\n");
	printf("Iteration numbers:\t%d\n", iterationcounter);
	printf("Log probability of initial model:\t%E\n", logprobinit);
	printf("Log probability of final model:\t%E\n", logprobfinal);
	printf("-------------------------\n");
	printHMM(&hmm);
	freeHMM(&hmm);

	endtime = dsecnd() - starttime;
	printf("------------------------\n");
	printf("total execution time : %f s\n", endtime);
	
	//getchar();
}