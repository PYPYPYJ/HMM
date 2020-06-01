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
rand()����0~RAND_MAX֮��������������RAND_MAX���ȡ0~1
*/
double getrandom(void)
{
	return (double)rand()/RAND_MAX;
}