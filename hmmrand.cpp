#include <time.h>
#include "hmm.h"

int getseed(void)
{
	return ((int)time(0));
}

void setseed(int seed)
{
	srand(seed);
}

/*
rand()返回0~RAND_MAX之间的随机数，除以RAND_MAX则获取0~1
*/
double getrandom(void)
{
	return (double)rand()/RAND_MAX;
}